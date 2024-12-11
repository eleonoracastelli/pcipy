#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 14:18:38 2024

@author: E Castelli
"""

import numpy as np


def interp_noise_stretch(n_dft, f, ip, basis_func, noisecoeffs, axs = None, **kwargs):#logf_knots = logf_knots, channels = channels):
    """
    Interpolate noise stretch periodograms.
    """
    # define design matrix from basis function
    design_matrix = basis_func(np.log(f[ip])).T
    # define projector
    projector = np.linalg.pinv(design_matrix.T.dot(design_matrix)).dot(design_matrix.T)
    # define Log-periodogram minus mean
    log_per = np.log(np.asarray([np.abs(n_dft[j][ip])**2 for j in range(3)]).T) + 0.577215
    # obtain noise coefficients for the stretch
    noisecoeffs += [[np.dot(projector, log_per[:, j]) for j in range(3)]]
    # plot noise stretches and interpolated coefficients
    if axs is not None:
        axs.T[len(noisecoeffs)-1][0].set_title("Noise stretch {n}".format(n = len(noisecoeffs)))
        for tdi, ax in enumerate(axs.T[len(noisecoeffs)-1]):
            ax.loglog(f[ip], np.abs(n_dft[tdi][ip])**2, zorder=0, alpha = 0.7)
            ax.scatter(np.exp(logf_knots), np.exp(noisecoeffs[-1][tdi]), color='tab:orange', zorder=1, s = 10)
            ax.set_ylabel("TDI {t} PSD".format(t=channels[tdi]))
            # for l in range(14
            #     ax.scatter(np.exp(logf_knots[l]), np.exp(noisecoeffs[-1][tdi][l]), zorder=1)
            #     # c = ax.gca().lines[-1].get_color()
            #     ax.loglog((f),): (((design_matrix.T[l]) * np.exp(noisecoeffs[-1][0][l]))), ls = '-')
        ax.set_xlabel("Frequency [Hz]")

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