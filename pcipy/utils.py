#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:18:38 2024

@author: E Castelli
"""

import numpy as np
from scipy import interpolate

def estimatePSD(f, fmax, dt, axs = None):

    # Switching to minimally correlated frequencies seleciton algorithm from LPF
    # # Number of interior knots
    # # n_knots = 12 # check number of knots
    # # # Total number of knots
    # # n_coeffs = n_knots + 2
    # # # Frequencies
    # # logf_knots = np.linspace(np.log(fmin), np.log(fmax), n_coeffs) # check knots location

    # Switching to minimally correlated frequencies seleciton algorithm from LPF
    min_noise_stretch_len = 7 * 24 * 3600 # 7 days minimum duration for noise estimation
    logf_knots = np.log(choosefreqs(Nmax=min_noise_stretch_len/dt,fmax=fmax,fs=1/dt)[1])
    n_coeffs = len(logf_knots)
    # redefine fmin, fmax, ip
    fmin = np.exp(logf_knots[0])
    fmax = np.exp(logf_knots[-1])
    # Band restriction
    ip = np.where((f>=fmin) & (f <= fmax))[0]    
    # skip padding
    skip = 21000   
    # interp kwargs
    kwargs = {'kind':'cubic',
            "axis":-1,
            "copy":True,
            "bounds_error": False,
            "fill_value": "extrapolate",
            "assume_sorted": True}
    # basis function
    basis_func = interpolate.interp1d(logf_knots, np.eye(n_coeffs), **kwargs)
    
    wd = masking.modified_hann(ub-skip-lb-skip, n_wind=n_wind)
    # apply windowing
    k2 = np.sum(wd**2)
    #create frequency array for the FT of the noise stretch
    fn = np.fft.fftfreq(data[comb].shape[0]) / dt
    # only select chosen frequencies
    ipn = np.where((fn>=fmin) & (fn <= fmax))[0]
    # create array of noise stretch dft for each TDI combination
    n_dft = [np.fft.fft(data[comb] * wd, data[comb].shape[0]) * np.sqrt(2 * dt / k2) for comb in channels]
    
    # define design matrix from basis function
    design_matrix = basis_func(np.log(f[ip])).T
    # define projector
    projector = np.linalg.pinv(design_matrix.T.dot(design_matrix)).dot(design_matrix.T)
    # define Log-periodogram minus mean
    log_per = np.log(np.asarray([np.abs(n_dft[j][ip])**2 for j in range(3)]).T) + 0.577215
    # obtain noise coefficients for the stretch
    noisecoeffs = [np.dot(projector, log_per[:, j]) for j in range(3)]
    # plot noise stretches and interpolated coefficients
    if axs is not None:
        logf_knots = np.log(choosefreqs(Nmax=stretchlen/dt,fmax=fmax,fs=1/dt)[1])
        axs.T[len(noisecoeffs)-1][0].set_title("Noise stretch {n}".format(n = len(noisecoeffs)))
        for tdi, ax in enumerate(axs.T[len(noisecoeffs)-1]):
            ax.loglog(f[ip], np.abs(n_dft[tdi][ip])**2, zorder=0, alpha = 0.7)
            ax.scatter(np.exp(logf_knots), np.exp(noisecoeffs[-1][tdi]), color='tab:orange', zorder=1, s = 10)
            ax.set_ylabel("PSD")
        ax.set_xlabel("Frequency [Hz]")
    
    print("* Noise stretch # {n}: duration {d:.1f} days".format(n=len(noisecoeffs), d=d/86400))
    # average interpolated coefficients
    coeffs = np.average(np.asarray(noisecoeffs), axis = 0)
    # interpolate log f knots
    log_psd_func = [interpolate.interp1d(logf_knots, c, **kwargs) for c in coeffs]
    # Define noise function
    def s_function(freqs, **kwargs):
        s_list = [np.exp(func(np.log(freqs))) for func in log_psd_func]
        return s_list


def interp_noise_stretch(n_dft, f, ip, basis_func, noisecoeffs, axs = None, **kwargs):#logf_knots = logf_knots, channels = channels):
    """
    Interpolate noise stretch periodograms.
    """


    return noisecoeffs

def choosefreqs(Nmax:int,fmax:float,fs:float):
    '''
    Minimally correlated frequencies selection, according to the CPSD estimation algorithm originally published in PhysRevLett.120.061101 and Supplemental Material.

    Returns 

    Parameters
    ----------
    Nmax : int
        DESCRIPTION.
    fmax : float
        DESCRIPTION.
    fs : float
        DESCRIPTION.

    Returns
    -------
    ndarray
        Number of averaging periodograms
    ndarray
        Frequency array.

    '''

    M = 4
    r = 3/5
    f = []
    L = [Nmax]
    k=1
    while(1):
        tmpL = r**(k-2) * Nmax
        tmpf = 2*M/(tmpL/fs);
        if (tmpf>fmax):
            break
        L.append(tmpL)
        f.append(tmpf)
        k += 1
    return np.asarray(L,dtype=int),np.asarray(f)


def bh_lowpass(data,t_win=100,t_sam=5,fs=10):
    '''
    Lowpass data by convolving with a BH92 windowing function.
    L Sala, December 2021
    modified by E Castelli, Sept 2022

    Parameters
    ----------
    data : tuple
        Synchronously sampled time series.
    t_win : TYPE, optional
        Window length, controls cut frequency. The default is 100.
    t_sam : TYPE, optional
        Output sampling time. The default is 5.
    fs : TYPE, optional
        Input sampling frequency. The default is 10.

    Returns
    -------
    datalp : TYPE
        DESCRIPTION.

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

    BHfilt = bh92(step_win) #build filter
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


def bh92(M:int):
    '''
    Parameters
    ----------
    M : int
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    z = np.arange(0,M)*2*np.pi/M
    return 0.35875 - 0.48829 * np.cos(z) + 0.14128 * np.cos(2*z) - 0.01168 * np.cos(3*z)