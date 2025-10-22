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


def estimate_sensitivity(filt, data_n, orbits, e_gw_mat = None, joint=True, welch_kwargs=None, strainPSD=1.0, analytic = True):
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
        Choose between analytic noise covariance matrix, w.r.t. numeric. The default is True.

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
    print(n_channels,'channels')

    # Welch spectrum matrix from noise data
    freqs, e_n_mat = compute_welch_matrix_td(e_n, **welch_kwargs)

    # GW response to a SGWB background
    if e_gw_mat is None:
        nside = 12
        npix = hp.nside2npix(nside)
        # skymap = np.ones(npix) * np.sqrt(4 * np.pi / (2*npix))
        # skymap = np.ones(npix) / np.sqrt(4 * np.pi * npix)
        skymap = np.ones(npix) / np.sqrt(2 * npix)

        #  Instantiate SGWB class
        sgwbcls = StochasticBackgroundResponse(skymap=skymap, orbits=orbits, orbit_interp_order=1)
        t0 = np.ones(1) * (sgwbcls.t0 + 10)  # need an object with a specific length
        # Now let's do it for all links
        gplus, gcross = sgwbcls.compute_correlations(sgwbcls.LINKS, freqs, t0)
        bgresponse = gplus + gcross
        # compute TDI response matrix
        # Time delay interferometry performs a linear combination of delayed single-link measurements to cancel laser frequency noise. It outputs three channels that contain all the GW information.

        # Algebraically, this transformation can be written in the (time-)frequency domain using a matrix operation:
        # If you have already computed the single-link responses like above, you can simply compute the TDI transfer function (for X, Y, Z) and then apply it through the equation above:

        bgresponse_ordered = bgresponse[..., convert, :]
        bgresponse_ordered = bgresponse_ordered[..., convert]

        # tdi transfer function at frequencies freqs and time t0
        # this is an array of 3 x 6 matrices (ncombs x nlinks)
        tdimat = sgwbcls.compute_tdi_design_matrix(freqs, t0, gen = '2.0')
        # compute the tdi response in the frequency domain
        # this is an array of 3 x 3 matrices
        #  Rtdi = Mtdi * Rlinks * Mtdi'
        # (ncombs x ncombs) = (ncombs x nlinks) * (nlinks x nlinks) * (nlinks x ncombs)
        # tdiresponse
        e_gw_mat = utils.transform_covariance(tdimat, bgresponse_ordered)
    else:
        nside = 12
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
        print('Calculate trace')
        sens = 1/np.trace(np.linalg.matmul(e_gw_mat,e_n_mat_inv), axis1=1, axis2=2)

    return freqs, sens
