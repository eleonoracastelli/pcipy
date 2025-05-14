#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utils.py

This module provides utility functions for estimating the Power Spectral Density (PSD) of signals.
The functions included in this module facilitate various methods of PSD estimation, including
the periodogram, Welch's method, and multitaper methods. These utilities are designed to assist
researchers and engineers in analyzing the frequency content of signals in various applications.

Classes:
    - PSDStats: class written by Lorenzo Sala in 2021. Used in LPF data analysis and publications: PhysRevLett.120.061101, ...
        Slightly adapted by E Castelli in 2022.    
        Computes the one-sided Cross Power Spectral Density (CPSD) matrix
        at minimally correlated frequencies, averaging over periodograms. Compute cross-coherences.
    
    
Functions:
- estimate_psd: Estimate PSD of a given data stretch, by interpolating the log-periodogram of the data DFT.
- choose_freqs: selection of minimally correlated frequencies according to PhysRevLett.120.061101


@author: E Castelli, 2025
"""

import numpy as np
import scipy.special as sp
from scipy import interpolate

# def estimate_psd(data, dt, Nmax, fmax, axs = None):
#     '''
#     Estimate PSD of a given data stretch, by interpolating the log-periodogram of the data DFT.
#     Uses minimally correlated frequencies according to the CPSD estimation algorithm originally published in PhysRevLett.120.061101
    
#     L Sala, Oct21/Sept22

#     Parameters
#     ----------
#     data : TYPE
#         DESCRIPTION.
#     dt : TYPE
#         DESCRIPTION.
#     Nmax : TYPE
#         DESCRIPTION.
#     fmax : TYPE
#         DESCRIPTION.
#     axs : TYPE, optional
#         DESCRIPTION. The default is None.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     '''

#     # Switching to minimally correlated frequencies seleciton algorithm from LPF
#     # # Number of interior knots
#     # # n_knots = 12 # check number of knots
#     # # # Total number of knots
#     # # n_coeffs = n_knots + 2
#     # # # Frequencies
#     # # logf_knots = np.linspace(np.log(fmin), np.log(fmax), n_coeffs) # check knots location

#     # Switching to minimally correlated frequencies selection algorithm from LPF
#     logf_knots = np.log(choose_freqs(Nmax=Nmax,fmax=fmax,fs=1/dt)[1])
#     n_coeffs = len(logf_knots)
#     # redefine fmin, fmax, ip
#     fmin = np.exp(logf_knots[0])
#     fmax = np.exp(logf_knots[-1])
#     # skip padding
#     skip = 21000   
#     # interp kwargs
#     kwargs = {'kind':'cubic',
#             "axis":-1,
#             "copy":True,
#             "bounds_error": False,
#             "fill_value": "extrapolate",
#             "assume_sorted": True}
#     # basis function
#     basis_func = interpolate.interp1d(logf_knots, np.eye(n_coeffs), **kwargs)
    
#     # set up window mask
#     wd = masking.modified_hann(ub-skip-lb-skip, n_wind=n_wind)
#     # apply windowing
#     k2 = np.sum(wd**2)
#     #create frequency array for the FT of the noise stretch
#     f = np.fft.fftfreq(data.shape[0]) / dt
#     # only select chosen frequencies
#     ip = np.where((f>=fmin) & (f <= fmax))[0]
#     # create array of noise stretch dft for each TDI combination
#     n_dft = np.fft.fft(data * wd, data.shape[0]) * np.sqrt(2 * dt / k2)
    
#     # define design matrix from basis function
#     design_matrix = basis_func(np.log(f[ip])).T
#     # define projector
#     projector = np.linalg.pinv(design_matrix.T.dot(design_matrix)).dot(design_matrix.T)
#     # define Log-periodogram minus mean
#     log_per = np.log(np.asarray([np.abs(n_dft[j][ip])**2 for j in range(3)]).T) + 0.577215
#     # obtain noise coefficients for the stretch
#     noisecoeffs = [np.dot(projector, log_per[:, j]) for j in range(3)]
#     # plot noise stretches and interpolated coefficients
#     if axs is not None:
#         logf_knots = np.log(choose_freqs(Nmax=Nmax,fmax=fmax,fs=1/dt)[1])
#         axs.T[len(noisecoeffs)-1][0].set_title("Noise stretch {n}".format(n = len(noisecoeffs)))
#         for tdi, ax in enumerate(axs.T[len(noisecoeffs)-1]):
#             ax.loglog(f[ip], np.abs(n_dft[tdi][ip])**2, zorder=0, alpha = 0.7)
#             ax.scatter(np.exp(logf_knots), np.exp(noisecoeffs[-1][tdi]), color='tab:orange', zorder=1, s = 10)
#             ax.set_ylabel("PSD")
#         ax.set_xlabel("Frequency [Hz]")
    
#     # average interpolated coefficients
#     coeffs = np.average(np.asarray(noisecoeffs), axis = 0)
#     # interpolate log f knots
#     log_psd_func = [interpolate.interp1d(logf_knots, c, **kwargs) for c in coeffs]
#     # Define noise function
#     def s_function(freqs, **kwargs):
#         s_list = [np.exp(func(np.log(freqs))) for func in log_psd_func]
#         return s_list

#     return noisecoeffs


class PSDstats():
    """
    PSDstats, compute the one-sided Cross Power Spectral Density (CPSD) matrix
        at minimally correlated frequencies, averaging over periodograms.
        Compute cross-coherences.
    L Sala, Oct21/Sept22

    Parameters
    ----------
    datamat : TYPE
        tuple of synchronously sampled time series
    Tmax : TYPE
        max time length of periodograms
    fmax : TYPE
        max frequency to be computed
    fs : TYPE
        sampling frequency
     win:       
        spectral window (function) or None. If None, BH92 is used

    Returns
    -------
    None.

    """
    def __init__(self,datatuple,Tmax,fmax,fs,win=None,units='',olap=50,detrend=False,getperiodograms=False,c=0.68):
        datamat = np.asarray(datatuple)
        if(datamat.ndim==1):
            datamat = np.reshape(datamat,(1,datamat.size))
        if detrend: datamat = datamat - np.mean(datamat,axis=1,keepdims=True)
        win = bh_92 if win is None else win
        initcheck(datamat,Tmax,fmax,fs)
        self.fs = fs
        self.olap = olap
        self.units = units
        self.L, self.freqs = choose_freqs(Nmax=Tmax*fs,fmax=fmax,fs=fs)
        if getperiodograms:
            self.periodograms, self.navs = get_periodograms(datamat=datamat,L=self.L,freqs=self.freqs,fs=fs,win=win,olap=olap)
        self.CPSD, self.navs = get_cpsd(datamat=datamat,L=self.L,freqs=self.freqs,fs=fs,win=win,olap=olap)
        self.cohere = get_cohere(CPSD=self.CPSD)

    def getPSDerrors(self,idx=0,cval=0.68):
        '''
        L Sala, Oct21/Sept22

        Parameters
        ----------
        datamat : TYPE
            DESCRIPTION.
        Tmax : TYPE
            DESCRIPTION.
        fmax : TYPE
            DESCRIPTION.
        fs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        Sxx  = np.real(self.CPSD[:,idx,idx])
        Smin = self.navs*Sxx/sp.gammainccinv(self.navs,(1-cval)/2)
        Smax = self.navs*Sxx/sp.gammainccinv(self.navs,(1+cval)/2)
        return Sxx,Smin,Smax

def bh_lowpass(data, t_win=100,t_sam=5,fs=10):
    '''
    Lowpass data by convolving with a BH92 windowing function.
    L Sala, December 2021
    modified by E Castelli, Sept 2022

    Parameters
    ----------
    data : rec-array
        Synchronously sampled time series with fields t', 'A', 'E', 'T'.
    t_win : TYPE, optional
        Window length, controls cut frequency. The default is 100.
    t_sam : TYPE, optional
        Output sampling time. The default is 5.
    fs : TYPE, optional
        Input sampling frequency. The default is 10.

    Returns
    -------
    datalp : rec-array
        Low-passed time series with fields 't', 'A', 'E', 'T'.
    '''
    names = data.dtype.names[1:]
    dt = 1/fs
    step_win = np.intc(t_win*fs)
    step_sam = np.intc(t_sam*fs)
    assert t_sam>=dt, 'Watch out, do not upsample your data.'
    assert np.isclose(t_sam*fs,int(t_sam*fs),rtol=1e-5), 'Downsampling time must be multiple of sampling time.'
    assert np.isclose(t_win*fs,int(t_win*fs),rtol=1e-5), 'Windowing time must be multiple of sampling time.'
    assert np.isclose(step_win/step_sam,int(step_win/step_sam),rtol=1e-5), 'Watch out, t_win must be multiple of t_sam.'
    
    dtarr = np.diff(data['t'])
    assert np.isclose(dtarr[0],dt,rtol=1e-5), 'Aaargh, sampling frequency is not consistent with data.' #just check fs
    assert np.allclose(dtarr,dt,rtol=1e-5), 'Aaargh, your data are not equally sampled in time.' #just check sampling time

    BHfilt = bh_92(step_win) #build filter
    BHarea = np.sum(BHfilt)
    BHfilt = BHfilt/BHarea
    onearray = np.ones(step_win)/step_win

    #apply filter convolving
    outts = [np.convolve(data['t'],onearray,mode='valid')] #just a simple way to get times, computationally more expensive than linspace, but safer
    for tdi in names:
        outts += [np.convolve(data[tdi], BHfilt,  mode='valid')]
    #downsample it
    for i in range(len(outts)):
        outts[i] = outts[i][::step_sam]
    datalp = np.rec.fromarrays(outts, names = ['t', 'A', 'E', 'T'])
    return datalp

def bh_92(M:int):
    '''
    Mathematical expression of a Blackman-Harris window.
    
    Parameters
    ----------
    M : int
        Length of BH window.

    Returns
    -------
    ndarray
        Array of values for the BH window of desired length.

    '''
    z = np.arange(0,M)*2*np.pi/M
    return 0.35875 - 0.48829 * np.cos(z) + 0.14128 * np.cos(2*z) - 0.01168 * np.cos(3*z)

# def bh_lowpass(data,t_win=100,t_sam=5,fs=10):
#     '''
#     Lowpass data by convolving with a BH92 windowing function.
#     L Sala, December 2021
#     modified by E Castelli, Sept 2022

#     Parameters
#     ----------
#     data : rec-array
#         Synchronously sampled time series with fields t', 'A', 'E', 'T'.
#     t_win : TYPE, optional
#         Window length, controls cut frequency. The default is 100.
#     t_sam : TYPE, optional
#         Output sampling time. The default is 5.
#     fs : TYPE, optional
#         Input sampling frequency. The default is 10.

#     Returns
#     -------
#     datalp : rec-array
#         Low-passed time series with fields 't', 'A', 'E', 'T'.
#     '''

#     names = data.dtype.names[1:]
#     dt = 1/fs
#     step_win = np.intc(t_win*fs)
#     step_sam = np.intc(t_sam*fs)
#     assert t_sam>=dt, 'Watch out, do not upsample your data.'
#     assert np.isclose(t_sam*fs,int(t_sam*fs),rtol=1e-5), 'Downsampling time must be multiple of sampling time.'
#     assert np.isclose(t_win*fs,int(t_win*fs),rtol=1e-5), 'Windowing time must be multiple of sampling time.'
#     assert np.isclose(step_win/step_sam,int(step_win/step_sam),rtol=1e-5), 'Watch out, t_win must be multiple of t_sam.'
    
#     dtarr = np.diff(data['t'])
#     assert np.isclose(dtarr[0],dt,rtol=1e-5), 'Aaargh, sampling frequency is not consistent with data.' #just check fs
#     assert np.allclose(dtarr,dt,rtol=1e-5), 'Aaargh, your data are not equally sampled in time.' #just check sampling time

#     BHfilt = bh92(step_win) #build filter
#     BHarea = np.sum(BHfilt)
#     BHfilt = BHfilt/BHarea
#     onearray = np.ones(step_win)/step_win

#     #apply filter convolving
#     outts = [np.convolve(data['t'],onearray,mode='valid')] #just a simple way to get times, computationally more expensive than linspace, but safer
#     for tdi in names:
#         outts += [np.convolve(data[tdi], BHfilt,  mode='valid')]
#     #downsample it
#     for i in range(len(outts)):
#         outts[i] = outts[i][::step_sam]
#     datalp = np.rec.fromarrays(outts, names = ['t', 'A', 'E', 'T'])
#     return datalp


# def bh92(M:int):
#     '''
#     Mathematical expression of a Blackman-Harris window.
    
#     Parameters
#     ----------
#     M : int
#         Length of BH window.

#     Returns
#     -------
#     ndarray
#         Array of values for the BH window of desired length.

#     '''
#     z = np.arange(0,M)*2*np.pi/M
#     return 0.35875 - 0.48829 * np.cos(z) + 0.14128 * np.cos(2*z) - 0.01168 * np.cos(3*z)

def initcheck(datamat,Tmax,fmax,fs):
    '''
    L Sala, Oct21/Sept22

    Parameters
    ----------
    datamat : TYPE
        DESCRIPTION.
    Tmax : TYPE
        DESCRIPTION.
    fmax : TYPE
        DESCRIPTION.
    fs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    npoints = datamat.shape[1]
    # assert(Tmax*fs<npoints)
    assert(fmax<=fs/2)
    return


def choose_freqs(Nmax:int,fmax:float,fs:float):
    '''
    Minimally correlated frequencies selection, according to the CPSD estimation algorithm originally published in PhysRevLett.120.061101 and Supplemental Material.
    

    Returns 

    Parameters
    ----------
    Nmax : int
        Maximum length of lowest frequency averaging window.
    fmax : float
        Maximum frequency at which to stop evaluation.
    fs : float
        Sampling frequency of the target data.

    Returns
    -------
    ndarray
        Number of averaging periodograms
    ndarray
        Frequency array.
    '''

    M = 4
    r = 3/5
    f = [M/(Nmax/fs)] # []
    L = [Nmax]
    k=2 # 1
    while(1):
        tmpL = r**(k-2) * Nmax
        tmpf = 2*M/(tmpL/fs);
        if (tmpf>fmax):
            break
        L.append(tmpL)
        f.append(tmpf)
        k += 1
    return np.asarray(L,dtype=int),np.asarray(f)

def get_cpsd(datamat,L,freqs,fs,win,olap):
    '''
    L Sala, Oct21/Sept22

    Parameters
    ----------
    datamat : TYPE
        DESCRIPTION.
    Tmax : TYPE
        DESCRIPTION.
    fmax : TYPE
        DESCRIPTION.
    fs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    CPSD = []
    navs = []
    for (tmpL,tmpf) in zip(L,freqs):
        tmpCPSD,tmpnavs,_ = get_cpsd_1freq(datamat=datamat,tmpL=tmpL,tmpf=tmpf,fs=fs,win=win,olap=olap)
        CPSD.append(tmpCPSD)
        navs.append(tmpnavs)
    return np.asarray(CPSD),np.asarray(navs)

def get_periodograms(datamat,L,freqs,fs,win,olap):
    '''
    L Sala, Oct21/Sept22

    Parameters
    ----------
    datamat : TYPE
        DESCRIPTION.
    Tmax : TYPE
        DESCRIPTION.
    fmax : TYPE
        DESCRIPTION.
    fs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    periodograms = []
    navs = []
    for (tmpL,tmpf) in zip(L,freqs):
        _,tmpnavs,tmpPeriodograms = get_cpsd_1freq(datamat=datamat,tmpL=tmpL,tmpf=tmpf,fs=fs,win=win,olap=olap)
        periodograms.append(tmpPeriodograms)
        navs.append(tmpnavs)
    return periodograms,np.asarray(navs)

def get_cpsd_1freq(datamat,tmpL,tmpf,fs,win,olap):
    '''
    L Sala, Oct21/Sept22

    Parameters
    ----------
    datamat : TYPE
        DESCRIPTION.
    Tmax : TYPE
        DESCRIPTION.
    fmax : TYPE
        DESCRIPTION.
    fs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    ndim = datamat.shape[0]
    npoints = datamat.shape[1]

    # DFT coefficients
    winpt = win(tmpL) #spectral window
    p = -2j*np.pi*np.arange(0,tmpL,1)/fs;
    C = np.exp(tmpf*p); #DFT coefficients

    startpoint = 0
    tmpnavs = 0
    Amat = np.zeros((ndim,ndim),dtype=complex)
    periodograms = []
    while(1):
        endpoint = startpoint + tmpL
        if(endpoint > npoints):
            break
        xs0 = datamat[:,startpoint:endpoint] #multivariate stretch
        xs = winpt * xs0 #windowed multivariate stretch
        ax = xs @ C #complex periodogram
        periodograms.append(ax)
        Amat += np.outer(ax,ax.conj()) #fill CPSD matrix
        tmpnavs += 1
        startpoint += int(tmpL*(100-olap)/100)

    wins2 = winpt@winpt
    tmpCPSD  = 2.0*Amat/tmpnavs/fs/wins2; #one-sided CPSD matrix
    periodograms  = np.asarray(periodograms)*np.sqrt(2.0/fs/wins2); #one-sided periodograms
    return tmpCPSD,tmpnavs,periodograms

def get_cohere(CPSD):
    '''
    L Sala, Oct21/Sept22

    Parameters
    ----------
    datamat : TYPE
        DESCRIPTION.
    Tmax : TYPE
        DESCRIPTION.
    fmax : TYPE
        DESCRIPTION.
    fs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    cohere = [ np.diag(np.diag(Sij)**(-1/2)) @ Sij @ np.diag(np.diag(Sij)**(-1/2))
        for Sij in CPSD]
    return np.asarray(cohere)
