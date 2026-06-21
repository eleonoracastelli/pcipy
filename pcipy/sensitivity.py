#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:41:14 2024

@author: E Castelli, J Baker

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

File: sensitivity.py
Purpose: This module provides functions for estimating the sensitivity of data formatted as TimeData objects.
"""

import numpy as np
import scipy.signal
import healpy as hp

from backgrounds import utils, noise
from backgrounds import StochasticBackgroundResponse
from backgrounds.tdi import convert
from lisaconstants.indexing import LINKS


def compute_welch_matrix_td(y_td, **kwargs):
    """
    Compute the welch estimated PSDs and CSDs of a multivariate time series fro TimeData Class objects.

    Parameters
    ----------
    ydata : ndarray
        array of time series, size n_samples x n_channels
    """

    args=kwargs.copy()
    args['fs']=1/y_td.dt

    ydata=y_td.data.T

    fy, _ = scipy.signal.welch(ydata[:, 0], **kwargs)
    welch_mat = np.zeros((fy.shape[0], ydata.shape[1], ydata.shape[1]), dtype=np.complex128)

    for i in range(ydata.shape[1]):
        for j in range(i, ydata.shape[1]):
            _, welch_mat[:, i, j] = scipy.signal.csd(ydata[:, i], ydata[:, j], **kwargs)
            welch_mat[:, j, i] = np.conjugate(welch_mat[:, i, j])
            welch_mat[:, i, i] = np.abs(welch_mat[:, i, i])

    return fy[1:], welch_mat[1:, : , :]


def estimate_sensitivity(filt, data_n, orbits, e_gw_mat = None, joint=True, welch_kwargs=None, strainPSD=1.0, analytic = False, nsideGWB=12, use_TDI_design_matrix=False,gw_ch_base_name='sci'):
    """
    Estimate the sensitivity of a given dataset against the response to an isotropic SGWB.

    Sensitivity is estimated by explicit computation of the inverse trace of the ratio of the channel simulated GW signal and the channel noise. To scale the sensitivity, we need to know the incident strain PSD value (assumed constant here).

    Parameters
    ----------
    filt : Filter Class.
        If None, data are treated unfilteree.
    data_n : TimeData Class.
        Instrument noise data.
    orbits : str
        Path to orbit file.
    e_gw_mat : GW covariance, optional
        Precalculated GW response from a StochasticBackgroundResponse object from backgrounds. If None, the GW response is evaluated within the class, slowing down calculation time. The default is None.
    joint : bool, optional
        Calculate the joint sensitivity of all channels in the TimeData class. The default is True.
    welch_kwargs : dict, optional
        Dictionary of Welch kwargs for the evaluation of the noise CSD matrix. The default is None.
    strainPSD : float, optional
        Incident strain PSD value (assumed constant here). The default is 1.
    analytic : bool, optional
        Choose between analytic noise covariance matrix, w.r.t. numeric. The default is False.
    nside : float, optional
        Specify the heal-pix number nside value indicating the all-sky GW source density.
    use_TDI_design_matrix : bool, optional
        If true, use backgrounds package code to compute the transfer function from the design matrix instead of computing the transfer function from the filter kernel. Irrelevant if gw_mat is provided.
    gw_ch_base_name: string, optional
        Base name for the input channels with single-link GW response. Default is 'sci'.
    Returns
    -------
    freqs : ndarray
        Array of frequencies.
    sens : ndarray
        Array of sensitivity values. If joint == False, it's an array with as many columns as channels.
    """
    if filt is not None:
        e_n = filt.apply_filter(data_n,method='convolve')
    else:
        e_n = data_n

    if welch_kwargs is None:
        nperseg = e_n.n_samples()
        fs = 1/e_n.dt
        welch_kwargs = {"fs": fs,
                  "window": 'blackman',
                  "nperseg": nperseg,
                  "detrend": 'constant',
                  "return_onesided": True,
                  "scaling": 'density'}

    n_channels=e_n.n_channels()
    #print(n_channels,'channels')

    # Welch spectrum matrix from noise data
    freqs, e_n_mat = compute_welch_matrix_td(e_n, **welch_kwargs)

    # GW response to a SGWB background
    if e_gw_mat is None:
        nside = nsideGWB
        npix = hp.nside2npix(nside)
        # skymap = np.ones(npix) * np.sqrt(4 * np.pi / (2*npix))
        # skymap = np.ones(npix) / np.sqrt(4 * np.pi * npix)
        skymap = np.ones(npix) / np.sqrt(2 * npix)

        #  Instantiate SGWB class
        sgwbcls = StochasticBackgroundResponse(skymap=skymap, orbits=orbits, orbit_interp_order=1)
        sgwbcls_t0 = sgwbcls.orbits.t_interp[0]
        t0 = np.ones(1) * (sgwbcls_t0 + 10)  # need an object with a specific length
        # Now let's do it for all links
        gplus, gcross = sgwbcls.compute_correlations(LINKS, freqs, t0)
        bgresponse = gplus + gcross
        # compute TDI response matrix
        # Time delay interferometry performs a linear combination of delayed single-link measurements to cancel laser frequency noise. It outputs three channels that contain all the GW information.

        # Algebraically, this transformation can be written in the (time-)frequency domain using a matrix operation:
        # If you have already computed the single-link responses like above, you can simply compute the TDI transfer function (for X, Y, Z) and then apply it through the equation above:
        bgresponse_reordered = bgresponse[..., convert, :]
        bgresponse_reordered = bgresponse_reordered[..., convert]

        if use_TDI_design_matrix:
            # tdi transfer function at frequencies freqs and time t0
            # this is an array of 3 x 6 matrices (ncombs x nlinks)        
            designmat = sgwbcls.compute_tdi_design_matrix(freqs, t0, gen = '2.0')
            # compute the response in the frequency domain
            # this is an array of 3 x 3 matrices
            #  Rtdi = Mtdi * Rlinks * Mtdi'
            # (ncombs x ncombs) = (ncombs x nlinks) * (nlinks x nlinks) * (nlinks x ncombs)
            # tdiresponse
            e_gw_mat = utils.transform_covariance(designmat, bgresponse_reordered)
        else:
            #This probably shouldn't be hard coded, but left. We need to ensure that the links are in the right order
            gw_input_names= [gw_ch_base_name +'_'+ link for link in ['12', '23', '31', '13', '21', '32']]
            
            filtGW=filt.restrict(input_names=gw_input_names)
            #Note that we have not addressed the question that the design matrix evolves with time
            #The "tdi_design_matrix" approach uses just a point slightly after t0 as the reference time
            #The piecewise filter approach so far uses the central version of the filter.
            designmat = filtGW.compute_filter_transfer_func(freqs)
            # compute the response in the frequency domain
            # this is an array of 3 x 3 matrices
            #  Rout = Tranf * Rlinks * Transf'
            # (ncombs x ncombs) = (ncombs x nlinks) * (nlinks x nlinks) * (nlinks x ncombs)
            # I expected that I should be using bg_response not bg_response_reordered here, but the latter yields a result which mostly matches theTDI_design_matrix version ... This is probably because of the order of the links in the list above FIXME
            e_gw_mat = utils.transform_covariance(designmat, bgresponse_reordered)

    else:
        ##This looks like an unnecessary hack. Why is it here?
        nside = nsideGWB
        npix = hp.nside2npix(nside)
        # skymap = np.ones(npix) * np.sqrt(4 * np.pi / (2*npix))
        # skymap = np.ones(npix) / np.sqrt(4 * np.pi * npix)
        skymap = np.ones(npix) / np.sqrt(2 * npix)

        #  Instantiate SGWB class
        sgwbcls = StochasticBackgroundResponse(skymap=skymap, orbits=orbits, orbit_interp_order=1)
        t0 = np.ones(1) * (sgwbcls.t0 + 10)  # need an object with a specific length
        

    if analytic==True:
        # OMS amplitudes from LISA Instrument paper
        oms_asds=(6.35e-12, 1.25e-11, 1.42e-12, 3.38e-12, 3.32e-12, 7.90e-12)
        # root sum square of the carrier amplitudes
        Aoms= np.sqrt(oms_asds[0]**2+oms_asds[2]**2+oms_asds[4]**2)
        Atm = 2.24e-15
        noise_classes_analytic = []
        noise_classes_analytic.append(noise.AnalyticOMSNoiseModel(
            freqs, t0, orbits, orbit_interp_order=1, gen="2.0", fs=None, duration=None, oms_isi_carrier_asds=Aoms))
        noise_classes_analytic.append(noise.AnalyticTMNoiseModel(
            freqs, t0, orbits, orbit_interp_order=1, gen="2.0", fs=None, duration=None, tm_isi_carrier_asds=Atm))

        # Compute the TDI spectrum from OMS and TM noise
        cov_tdi_n_comps = [nc.compute_covariances(0.0) for nc in noise_classes_analytic]
        # Compute the full TDI noise spectrum
        e_n_mat = sum(cov_tdi_n_comps)

    # Output sensitivity for each variable, size nfreqs x n_channels
    sens = np.array([np.abs(e_n_mat[:, j, j] / e_gw_mat[:, j, j] * strainPSD) for j in range(n_channels)]).T

    if joint:
        e_n_mat_inv = np.linalg.inv(e_n_mat)
        #print('Calculate trace')
        sens = 1/np.trace(np.linalg.matmul(e_gw_mat,e_n_mat_inv), axis1=1, axis2=2)

    return freqs, sens

def estimate_sensitivity_empirically(filt, data_n, data_gw, gw_input_names=None, joint=True, welch_kwargs=None, strainPSD=1):
    #Sensitivity is estimated by explicit computation of the ratio of the channel noise and channel simulated GW signal
    #To scale the sensitivity, we need to know the incident strain PSD value (assumed constant here)
    #We understand this as the sum of the plus and cross PSDs for each point source in the simulation
    """
    Estimate the sensitivity of a given dataset against the response to a simulated GW signal.

    Sensitivity is estimated by explicit computation of the ratio of the channel simulated GW signal and the channel simulated noise. To scale the sensitivity, we need to know the incident strain PSD value (assumed constant here).

    Parameters
    ----------
    filt : Filter Class.
        If None, data are treated unfilteree.
    data_n : TimeData Class.
        Instrument noise data.
    data_gw : TimeData Class.
        Simulated gw data.
    gw_input_names: list of strings, optional
        In the case that the GW simulation includes fewer channels than the noise simulation, specify the included names here. Must be a subset of the input channel names expected by the filter. Default None indicates all channels are included.
    joint : bool, optional
        Calculate the joint sensitivity of all channels in the TimeData class. The default is True.
    welch_kwargs : dict, optional
        Dictionary of Welch kwargs for the evaluation of the noise CSD matrix. The default is None.
    strainPSD : float, optional
        Incident strain PSD value (assumed constant here). The default is 1.

    Returns
    -------
    freqs : ndarray
        Array of frequencies.
    sens : ndarray
        Array of sensitivity values. If joint == False, it's an array with as many columns as channels.

    """
    
    if filt is not None:
        e_n = filt.apply_filter(data_n,method='convolve')
        if gw_input_names is not None:
            filtGW=filt.restrict(input_names=gw_input_names)
        else:
            filtGW=filt
        e_gw = filtGW.apply_filter(data_gw,method='convolve')
            
    else:
        e_n = data_n
        e_gw = data_gw
        
    n_channels=e_n.n_channels()
    #print(n_channels,'channels')
    
    # Welch spectrum matrix for PCI variables from noise
    freqs, e_n_mat = compute_welch_matrix_td(e_n, **welch_kwargs)

    # Welch spectrum matrix for PCI variables from noise
    freqs, e_gw_mat = compute_welch_matrix_td(e_gw, **welch_kwargs)

    # Orthogonalization
    _, s, vh = np.linalg.svd(e_n_mat)

    # Apply the orthogonal transformation to the GW signal    
    e_gw_mat_ortho = np.einsum("ijk, ikl -> ijl", vh,
                               np.einsum("ijk, ikl -> ijl", e_gw_mat, np.swapaxes(vh, 1, 2).conj()))

    # Apply the orthogonal transformation to the noise covariance
    e_n_mat_ortho = np.einsum("ijk, ikl -> ijl", vh,
                              np.einsum("ijk, ikl -> ijl", e_n_mat, np.swapaxes(vh, 1, 2).conj()))
    
    factor=strainPSD

    # Output sensitivity for each variable, size nfreqs x n_channels
    # pci_sens = np.array([np.abs(s[:, j] / e_pci_gw_mat_ortho[:, j, j]) for j in range(n_channels)]).T
    sens = np.array([np.abs(e_n_mat_ortho[:, j, j] / e_gw_mat_ortho[:, j, j] * factor) for j in range(n_channels)]).T
    
    if joint:
        #print('joining')
        sens = 1 / np.sum(1/np.array(sens), axis=1)
    
    return freqs, sens
