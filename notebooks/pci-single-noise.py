#!/usr/bin/env python
# coding: utf-8

# # Apply PCI to simulated data 
# ## Investigate the impact of individual noise contributions
# 
# Here we apply PCI to data simulated via the LISA Simulation Suite, following the simulation scripts.
# 
# The simulated datasets are 3 days long. PCI is applied to 12 hours of data, with 4 hours of skipped data at the beginning of the simulation.
# 
# We need the following simulated datasets:
# - full simulation noise dataset (including laser noise and secondary noises), with filename ending in `_measurements_4Hz.h5`
# - secondary noises dataset, with filename ending in `_noise_sec_4Hz.h5`

# # 0. Installations and data generation
# 
# The package dependencies are:
# 
#     pip install numpy scipy sympy h5py matplotlib xarray h5py scikit-learn
#     pip install lisaconstants
#     pip install lisainstrument
#     pip install lisagwresponse
#     pip install pytdi
#     pip install backgrounds
# 
# and after installation, the data generation step is performed by running the simulation scripts: 

#!python ../simulation/noise_simulation.py /Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations --tdi 2 --baseline --individual

#!python ../simulation/signal_simulation.py /Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations --tdi 2

#!python ../simulation/all_sky_signal_simulation.py /Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations --tdi 2

# %% ## 0.1 Settings and imports
# Importing the relevant packages for the notebook. Setting up work directories.

import h5py
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

from datetime import datetime

from pytdi import Data
from pytdi.intervar import ETA_SET

from pcipy import plotting, pci_filter, channel_analysis 

# %% Working directories and paths to files

workdir = "/Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations/"
datadir = '/Users/ecastel2/Documents/research/GSFC/simulation-tests/orbits/'

measpath = '_measurements_4Hz.h5'
secondpath = '_noise_sec_4Hz.h5'

orbits = datadir+"keplerian-orbits.h5"

substring = ['locking_n1_12_laser_tm_oms', 'all_sky_gw', 'point_gw']
sims = ["noise", "all_sky", "point_source"]
datasets = dict(zip(sims, substring))

matchfile={}
dtpath={}

for n,d in zip(substring,datasets):
    timestamp=[]
    matchfile[n] = [f for f in os.listdir(workdir) if n in f]
    for m in matchfile[n]:
        # pick latest date
        timestamp.append(datetime.strptime(m[:17], "%Y-%m-%d_%Hh%M_"))
        #print(n, timestamp[n])
        dtpath[d] = max(timestamp).strftime("%Y-%m-%d_%Hh%M_")

print(dtpath)

# %%


skip_hours=4
pci_hours=12


# %% ## 1. Build data vector of the six LISA single-link channels
# 
# To build the data vector of the six LISA single-link channels $\vec{y} = \left[y_{ij}\right]$, with $i,j=1,2,3$ and $i\neq j$ we resort to the intermediary TDI variables $\eta$, implemented within `pytdi` as `ETA_SET`.
# 
# 
# We build the single link $\vec{y}$ data vector for the full noise simulation and for the secondary noises, ending up with two single link vectors:
# - full simulation single link vector $\vec{y}^{\text{full}}$
# - secondary noises single link vector $\vec{y}^{\text{sec}}$

def build_data_vector(data_noise, skip=300, dtype=np.float64):
    central_freq = 281600000000000.0
    # Conventions 1, 2, 3, 1p, 2p, 3p
    # delays_order = ['23', '31', '12', '32', '13', '21']
    mosas_order = ['12', '23', '31', '13', '21', '32']
    # Form intermediary variables for full data
    ETA_data_set = {key: ETA_SET[key].build(**data_noise.args) for key in ETA_SET.keys()}
    eta_noise_set = {key: ETA_data_set[key](data_noise.measurements) for key in ETA_data_set.keys()}
    # Form the measurement vector for moving arms containing all noises
    y = np.array([eta_noise_set[f'eta_{mosa}'] / central_freq for mosa in mosas_order], dtype=dtype).T
    y_full = y[skip:, :]
    del y

    return y_full

# %% 1.1 Get noise simulation measurements

simpath = workdir + dtpath['noise'] + datasets['noise'] + measpath
print(simpath)
# load hdf5 file to read data attrs
sim = h5py.File(simpath, 'r')
# read attrs
central_freq = sim.attrs['central_freq']
dt = sim.attrs['dt']
n_data = sim.attrs['size']
# load data
data_noise = Data.from_instrument(simpath)
# set sampling frequency
fs = data_noise.fs

# We skip the earliest part of the sim which is not representative
skip = int(skip_hours * 3600 * fs)  
# build data vector of full simulation data
y_full = build_data_vector(data_noise, skip=skip, dtype=np.float64)
# create time vector for plotting
times = np.arange(len(y_full))*dt


# %%

# get some info on Y_full characteristics
print(y_full.shape)
#print(y_full[:ns+2*nhalf, :].shape)
#print(y_full[:ns+2*nhalf, :].T.shape)


# %%

# 1.2 Get secondary noises
secpath = workdir + dtpath['noise']+ datasets['noise'] + secondpath
print(secpath)
# load data
data_sec = Data.from_instrument(secpath)
# build data vector of secondary noises   
y_sec = build_data_vector(data_sec, skip=skip, dtype=np.float64)

# %% ### 1.2 Data quicklook
# 
# Check the stencil size and take a look at the generated $y$s.
# 
# TDI2 fractional delay Lagrange interpolating polynomials are of order `31=1+15*2`
# 
# Overall TDI delays are up to 8x single link delay = `8.34 s * 4 Hz * 8` = about 267 sample. Add 15 on each end: `267 + 2*15 = 297` 
# 
# TDI2 overall stencil is then about 297 samples
# 
# aPCI overall stencil width is `1 + nhalf*2`

# %% Set stencil size and number of samples

nhalf = 45
ns = int(pci_hours * 3600 * fs) 
window = np.ones(ns)

# %% # Time series quicklook

# Plot 4000 data samples:
maxshow=4000
ev=ns//maxshow+1
fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex = True)
for ich in range(len(y_full.T)):
    axes[0].plot(times[skip:skip+ns:ev], y_full[skip:skip+ns:ev, ich])
for ich in range(len(y_sec.T)):
    axes[1].plot(times[skip:skip+ns:ev], y_sec[skip:skip+ns:ev, ich])
axes[1].set_xlabel('Time [s]')
axes[0].set_ylabel(r"$y_\text{full}$ [frac. freq.]")
axes[1].set_ylabel(r"$y_\text{sec}$ [frac. freq.]")
axes[0].grid(), axes[1].grid()



# %%

# Quicklook TDI files and

tdipath2 = workdir +  dtpath['noise'] + datasets['noise'] + '_noise_tdi2_4Hz.h5'

tdi2 = h5py.File(tdipath2, 'r')

x2_noise = tdi2['x'][()] / central_freq
y2_noise = tdi2['y'][()] / central_freq
z2_noise = tdi2['z'][()] / central_freq

nperseg=4e5

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": nperseg,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

# tdi_times = np.arange(ns)*dt

# fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex = False)
# ax[0].plot(tdi_times, x2_noise[skip:skip+ns], label = 'X')
# ax[0].plot(tdi_times, y2_noise[skip:skip+ns], label = 'Y')
# ax[0].plot(tdi_times, z2_noise[skip:skip+ns], label = 'Z')
# ax[0].grid()
# ax[0].legend()
# ax[0].set_xlabel('Time [s]')
# ax[0].set_ylabel('TDI')

# f, xpsd = signal.welch(x2_noise[skip:skip+ns], **kwargs)
# f, ypsd = signal.welch(y2_noise[skip:skip+ns], **kwargs)
# f, zpsd = signal.welch(z2_noise[skip:skip+ns], **kwargs)

# ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
# ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
# ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
# ax[1].grid()
# ax[1].legend()
# ax[1].set_xlabel('Frequency [Hz]')
# ax[1].set_ylabel(r'$S_\text{TDI}(f)$')

# %%## 2. Apply PCI to the data
# 
# We now resort to {class}`PCIFilter` to evaluate PCI from $\vec{y}$.
# 
# An instance of {class}`PCIFilter` has two required inputs:
# - `ydata`: matrix with the single link LISA temporal phase data streams $\vec{y}$ of length `ns`.
# - `fs`: sampling rate of the data streams (Hz).
# 
# The optional parameters are
# - `nhalf`: filter stencil halfwidth in samples. The default is 45.
# - `order`: order of PCI. The default is 1.
# - `maxcompts`: PCA results will be truncated to this length after initial processing. The default is 10.
# - `Tscale`: if dt is None, then dt defaults to Tscale/ns
# 
# The input channels $\vec{y}$ are stretches of data, usually of length `ns+2*nhalf`, but sometimes varying. The variations are:
#   - `:` (full length of matrix)
#   - `0:ns_fixed`
#   - `0:ns+2*nhalf`
#   - `skip:skip+ns+2*nhalf`
#   
# In every case the window is trivial `np.ones([data lenght])`
# 
# Creating an instance of the `PCIFilter` class applies the following methods, in order:
# 
# 1. `pcipy.PCIFilter.build_data_matrix`
# Pre-process $\vec{y}$ data `ydata` to build a matrix of shifted time-series. 
# 
# Output is a matrix $\vec{X}$ of size `[n_c * (2*nhalf+1)] x n_s` , 
# $$
# \vec{X} = \begin{bmatrix}
# y_1(t_0) & y_1(t_1) &y_1(t_2) & \dots & y_1(t_{n_s}) \\
# y_1(t_0-t_\text{st})& y_1(t_1-t_\text{st}) &y_1(t_2-t_\text{st}) & \dots & y_1(t_{n_s}-t_\text{st}) \\
# y_1(t_0-2t_\text{st}) & \dots & \dots & \dots & y_1(t_{n_s}-2t_\text{st}) \\
# \vdots & \vdots & \vdots & \vdots & \vdots \\
# y_1(t-n_\text{half}t_\text{st})&y_1(t_1-n_\text{half}t_\text{st}) &y_1(t_2-n_\text{half}t_\text{st}) & \dots & y_1(t_{n_s}-n_\text{half}t_\text{st}) \\
# y_1(t_0-2t_\text{st}) & \dots & \dots & \dots &\dots \\
# y_1(t_0+t_\text{st})& \dots & \dots & \dots &\dots \\
# y_1(t_0+2t_\text{st}) & \dots & \dots & \dots &\dots \\
# \vdots & \vdots & \vdots &\vdots &\vdots\\
# y_1(t_0+t_\text{nhalf})& \dots & \dots & \dots &\dots \\
# y_2(t_0) & \dots & \dots & \dots & \dots \\
# y_2(t_0-t_\text{st}) & \dots & \dots & \dots &\dots \\
# \vdots & \vdots & \vdots &\vdots &\vdots\\
# y_2(t-n_\text{half}t_\text{st})& \dots & \dots & \dots &\dots \\
# y_2(t_0+t_\text{st})& \dots & \dots & \dots &\dots \\
# \vdots & \vdots & \vdots &\vdots &\vdots\\
# y_2(t_0+n_\text{half}t_\text{st})& \dots & \dots & \dots &\dots \\
# y_3(t_0) & \dots & \dots & \dots &\dots \\
# \vdots & \vdots & \vdots &\vdots &\vdots\\
# y_{n_c}(t_0+n_\text{half}t_\text{st}) & \dots  & \dots & \dots &  y_{n_c}(t_{n_s}+n_\text{half}t_\text{st})
# \end{bmatrix}
# $$
# where $n_c=6$ is the size of the $y$s array, `nhalf=45` is the length of the time-shifting stencil and $n_s$ is the number of data samples in each $y$.
# 
# Option to detrend the data or make them zero-mean.
# 
# 2. `pcipy.PCIFilter.apply_pca`
# Apply Principal Component Analysis (PCA) to matrix $\vec{X}$. PCA is applied using the `PCA` class defined within [scikit learn](https://scikit-learn.org/stable/modules/decomposition.html#pca) `sklearn`, which 
# 
# > is used to decompose a multivariate dataset in a set of successive orthogonal components that explain a maximum amount of the variance. In scikit-learn, PCA is implemented as a transformer object that learns
# components in its fit method, and can be used on new data to project it on these components.
# 
# From `PCA` object, evaluate components and explained variance. Option to sort components by RMS value.
# 
# 4. `pcipy.PCIFilter.set_stencil`

# %% ### 2.1 Evaluate PCI on full data
# #### 2.1.1 Sort components by variance

# import importlib
# importlib.reload(pci_filter)

Tscale=10
pca_list = [pci_filter.PCIFilter(y_full[:ns+2*nhalf, :].T,
                                 fs=fs, 
                                 nhalf=nhalf, 
                                 order=q, 
                                 maxcompts=10, 
                                 Tscale=Tscale,
                                 sort_by_rms=False,
                                verbose=True)            
            for q in range(3)]

#c5demand8 mem peak 77%   ???
#              total        used        free      shared  buff/cache   available
#Mem:       15911316     1540616    13584348        1276      786352    14088684

# get_ipython().system('free')

# %% print shapes of all the involved quantities

p = pca_list[0]
print(p.nsdata, p.nc, p.nhalf, p.nhalf*2+1, p.ns+2*p.nhalf, p.ns, p.maxcompts)

print((p.nhalf*2+1)* p.nc)

p.components.shape, p.explained_variance.shape, p.channels.shape

# %% #### 2.1.2 Sort components by RMS

pca_list_rs = [pci_filter.PCIFilter(y_full[:ns+2*nhalf, :].T,
                                    fs=fs, 
                                    nhalf=nhalf, 
                                    order=q, 
                                    maxcompts=10, 
                                    Tscale=Tscale,
                                    sort_by_rms=True)            
            for q in range(3)]

#c5demand8 mem peak 77%   ???
#              total        used        free      shared  buff/cache   available
#Mem:       15911316     1540616    13584348        1276      786352    14088684

# get_ipython().system('free')

# %% #### 2.1.3 Detrend the components to get zero mean data

pca_list_zm = [pci_filter.PCIFilter(y_full[:ns+2*nhalf, :].T,
                                    fs=fs, 
                                    nhalf=nhalf, 
                                    order=q, 
                                    maxcompts=10, 
                                    Tscale=Tscale,
                                    zero_mean=True)            
            for q in range(3)]


# %% ### 2.2 Plot PCI decomposition
# 
# Compare data sorted by variance with data sorted by RMS.

plotting.plotconfig(lbsize=20, lgsize=16, fsize=18, 
                    ticklabelsize=20, style='publication',
                    fontfamily = 'STIXGeneral')

fig1, ax1 = plt.subplots(nrows=1)

ax1.plot(pca_list[0].explained_variance, 
            #linestyle='dashed',
            label=r'aPCI$_0$',
            linewidth=2,
            rasterized=False)

ax1.plot(pca_list[1].explained_variance, 
            #linestyle='dashed',
            label=r'aPCI$_1$ (1st order)',
            linewidth=2,
            rasterized=False)

ax1.plot(pca_list[2].explained_variance,
            #linestyle='dashed',
            label=r'aPCI$_2$ (2nd order)',
            linewidth=2,
            rasterized=False)

if True:
    ax1.plot(pca_list_rs[0].explained_variance, 
                linestyle='dashed',
                label=r'aPCI$_0$ resorted',
                linewidth=2,
                rasterized=False)

    ax1.plot(pca_list_rs[1].explained_variance, 
                linestyle='dashed',
                label=r'aPCI$_1$ resorted (1st order)',
                linewidth=2,
                rasterized=False)

    ax1.plot(pca_list_rs[2].explained_variance,
                linestyle='dashed',
                label=r'aPCI$_2$ resorted (2nd order)',
                linewidth=2,
                rasterized=False)

if False:
    ax1.plot(pca_list_zm[0].explained_variance, 
                linestyle='dotted',
                label=r'aPCI$_0$ zm',
                linewidth=2,
                rasterized=False)

    ax1.plot(pca_list_zm[1].explained_variance, 
                linestyle='dotted',
                label=r'aPCI$_1$ zm (1st order)',
                linewidth=2,
                rasterized=False)

    ax1.plot(pca_list_zm[2].explained_variance,
                linestyle='dotted',
                label=r'aPCI$_2$ zm (2nd order)',
                linewidth=2,
                rasterized=False)

ax1.set_xscale('linear')
ax1.set_yscale('log')
ax1.set_xlabel(r"Component number", fontsize=20)
ax1.set_ylabel("Variance", fontsize=20)
# ax1.set_ylim([1e-42, 1e-25])
ax1.minorticks_on()
plt.legend(loc='upper right', frameon=False)
plt.title("Duration fixed "+str(int(ns/fs/3600))+ " hours, nh=" + str(nhalf))
plt.show()
# get_ipython().system('free')
# for i in range(3):
#    print(pca_list[i].explained_variance-pca_list_zm[i].explained_variance)


# %% ## 3. Run channel analysis on the PCI output
# 
# The channel analysis we run here does the following:
# - Compare stationarity of the data, for each applied PCI order

order=2
ev=100

# select which order of PCI filter you want to look at
filters=[pca_list_rs[o] for o in range(order+1)]
# select titles for plot
titles = ['rms_sorted order '+str(o) for o in range(order+1)]
# apply_for_channels is a function tasked with applying PCI filters to a specific subset of channels 
# (in this case, the 6 channels with lower variance)
sets=[np.array(xf.apply_for_channels(y_full[:ns+2*nhalf, :].T, 
                                     n_channels=6,
                                     zero_mean=False,
                                     detrend=False)) 
      for xf in filters]
# run channel analysis on selected channels and filter orders
channel_analysis.stationarity_plots(sets,
                                    title=titles)


# %% # - Compare temporal variance of the data, for each applied PCI order

order=2
ev=100

# select which order of PCI filter you want to look at
filters=[pca_list_rs[o] for o in range(order+1)]
# select titles for plot
titles = ['rms_sorted order '+str(o) for o in range(order+1)]

# apply_for_channels is a function tasked with applying PCI filters to a specific subset of channels 
# (in this case, the 6 channels with lower variance)
sets=[np.array(xf.apply_for_channels(y_full[:ns+2*nhalf, :].T,
                                     n_channels=6,
                                     zero_mean=False,detrend=False)) 
      for xf in filters]
# FIX TITLES IN THIS PLOT FUNCTION
channel_analysis.temporal_variance_corr_plots(sets,
                                              nchan=6,
                                              title=titles)


# %% ### 3.1 Frequency Domain analysis

def compute_welch_matrix(ydata, **kwargs):
    """
    Compute the welch estimated PSDs and CSDs of a multivariate time series.

    Parameters
    ----------
    ydata : ndarray
        array of time series, size n_samples x n_channels
    """

    fy, _ = signal.welch(ydata[:, 0], **kwargs)
    welch_mat = np.zeros((fy.shape[0], ydata.shape[1], ydata.shape[1]), dtype=np.complex128)

    for i in range(ydata.shape[1]):
        _, welch_mat[:, i, i] = signal.welch(ydata[:, i], **kwargs)

        for j in range(i+1, ydata.shape[1]):
            _, welch_mat[:, i, j] = signal.csd(ydata[:, i], ydata[:, j], **kwargs)
            welch_mat[:, j, i] = np.conjugate(welch_mat[:, i, j])

    return fy, welch_mat

# %% # Apply Welch to matrix of time-shifted $y$s:

# Welch settings
nperseg = 2**16 # didn't work with some, too short data segment (6hr)?
nperseg = 2**14

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": nperseg,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

# apply Welch
freqs, y_welch_mat = compute_welch_matrix(y_full, **kwargs)
freqs, y_sec_welch_mat = compute_welch_matrix(y_sec, **kwargs)


# %% # Vizualize the single-link measurements in frequency domain

# define ordering of the MOSAs in the constellation
mosas_order = ['12', '23', '31', '13', '21', '32']

# plot config
plotting.plotconfig(lbsize=20, lgsize=16, fsize=18, 
                    ticklabelsize=20, style='publication',
                    fontfamily = 'STIXGeneral')

# plot
fig, axes = plt.subplots(1, 1, figsize=(10, 6))
# For APCI
for i in range(6):
    axes.loglog(freqs, np.sqrt(y_welch_mat[:, i, i].real), 
                linewidth=1, 
                label=r'$y_{\mathrm{'+mosas_order[i]+'}}$',
                rasterized=True)
plt.gca().set_prop_cycle(None)
for i in range(6):
    axes.loglog(freqs, np.sqrt(y_sec_welch_mat[:, i, i].real), 
                linewidth=1,ls='--', 
                label=r'$y_{\mathrm{'+mosas_order[i]+',sec}}$',
                rasterized=True)
axes.legend(loc='upper center', ncol=4, frameon=False)
# axes.grid(linewidth=1, which='both', 
#           color='gray', 
#           linestyle='dotted')
axes.set_xlabel("Frequency [Hz]")
axes.set_ylabel(r"$\mathrm{\sqrt{PSD}}$ [$\mathrm{Hz}^{-1/2}$]")
axes.set_xlim([1e-4, 1])
axes.set_ylim([1e-25, 1e-9])
# axes.set_title("Single-link periodograms")
fig.savefig("single-link-periodogram.pdf", format="pdf", dpi=300)
plt.show()


# %% ## 4. Reconstruct single-link channels

order=2
nchannels=10
ev=100

filters=[pca_list_rs[order-1],pca_list_rs[order]]

[xf.set_stencil(nchannels) for xf in filters]

sets=[np.array(xf.apply_for_channels(y_full[:ns+2*nhalf, :].T, 
                                     n_channels=nchannels,
                                     zero_mean=False,detrend=False)) 
      for xf in filters]

print("Computing single link")
# Compute recovered single-link channel vectors from a set of pri-computed PCI channels using the stencil_compts.
Ysets=[xf.compute_single_links_from_channels(iset) 
       for xf,iset in zip(filters,sets)]

titles=['Y_rs'+str(order-1),'Y_rs'+str(order)]

# the raw single links are the ones obtained from the secondary noises only (no laser noise suppression)
# we're evaluating them for comparison
if True:
    Ysets+=[y_sec[nhalf:ns+nhalf, :].T]
    titles+=['Y_raw']
    
print("Computed single link")
channel_analysis.stationarity_plots(Ysets,
                                    title=titles)


# %% Evaluate temporal variance correlation plots for the reconstructed single links
channel_analysis.temporal_variance_corr_plots(Ysets[:],
                                              nchan=6,
                                              title=titles)
# %% Evaluate Welch CSD matrix

freqs1, y_pci1_welch_mat = compute_welch_matrix(Ysets[0].T, **kwargs)
freqs2, y_pci2_welch_mat = compute_welch_matrix(Ysets[1].T, **kwargs)

# %% Visualize the single-link measurements in frequency domain
plotting.plotconfig(lbsize=20, lgsize=16, fsize=18, 
                    ticklabelsize=20, style='publication',
                    fontfamily = 'STIXGeneral')
fig, axes = plt.subplots(1, 1, figsize=(10, 6))

# For APCI
for i in range(6):
    axes.loglog(freqs, np.sqrt(y_pci1_welch_mat[:, i, i].real), 
                linewidth=1, ls=':',
                label=r'$y_{\mathrm{'+mosas_order[i]+',pci1}}$',
               rasterized=True)
#for i in range(6):
#    axes.loglog(freqs, np.sqrt(y_pci2_welch_mat[:, i, i].real), 
#                linewidth=1, 
#                label=r'$y_{\mathrm{'+mosas_order[i]+',pci2}}$',
#               rasterized=True)
for i in range(6):
    axes.loglog(freqs, np.sqrt(y_sec_welch_mat[:, i, i].real), 
                linewidth=1, ls='--',
                label=r'$y_{\mathrm{'+mosas_order[i]+',sec}}$',
                rasterized=True)
axes.legend(loc='upper center', ncol=4, frameon=False)
# axes.grid(linewidth=1, which='both', 
#           color='gray', 
#           linestyle='dotted')
axes.set_xlabel("Frequency [Hz]")
axes.set_ylabel(r"$\mathrm{\sqrt{PSD}}$ [$\mathrm{Hz}^{-1/2}$]")
axes.set_xlim([1e-4, 1])
axes.set_ylim([1e-25, 1e-9])
# axes.set_title("Single-link periodograms")
fig.savefig("single-link-periodogram.pdf", format="pdf", dpi=300)
plt.show()


# %% Alternative $Y$ reconstruction

def filter_single_link_data(self, ydata,n_channels=None):
    '''
    Compute recovered single-link channel vectors from a set of single-link data. This should yield 
    identical results to appropriately called compute_single_link_from_channels but we first implement
    separately for testing, since the stencil_compts stuff is new. This one is more low-level and 
    doesn't assum the stencil is selected.
    
    Parameters
    ----------
    ydata : ndarray 
        The raw single-link data-set to be transformed.

    Returns
    -------
    reconstructed_ydata : ndarray
        Reconstructed single link data channels
    '''
    if n_channels is None: n_channels=self.maxcompts
    y_transformer=self.components[-n_channels:,:6].T

    chans_data = self.apply_for_channels(ydata, n_channels)
    Z=np.dot(y_transformer,chans_data)
    print('single link data shape', Z.shape)
    return Z


# %%

Zsets=[np.array(xf.filter_single_link_data(y_full[:ns+2*nhalf, :].T, n_channels=nchannels)) for xf in filters]
print("Comparing single link")
diff=[Zsets[i]-Ysets[i] for i in range(order)]
print(diff)


# %% ## 5. Estimate sensitivity
# ### Computation of empirical response with welch periodograms
# Based and empirical sensitivity calculation simulated signal over simulated noise as developed in initial-dev-from-pylisa

# %% # Some options and prep for the welch_matrix calc
nperseg = 1e4
welch_kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": nperseg,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

#and for the orthogonalization...
#multiple_dot = lambda a,b: np.einsum("ijk, ikl -> ijl", a, b)
def multiple_dot(a_mat, b_mat):
    """
    Perform the matrix multiplication of two list of matrices.

    Parameters
    ----------
    a : ndarray
        series of m x n matrices (array of size p x m x n)
    b : ndarray
        series of n x k matrices (array of size p x n x k)

    Returns
    -------
    c : ndarray
        array of size p x m x k containg the dot products of all matrices
        contained in a and b.
    """

    return np.einsum("ijk, ikl -> ijl", a_mat, b_mat)


# %% Estimate sensitivity function

def estimate_sensitivity(pci, data_n, data_gw, n_channels=6, joint=True, single_link=False, welch_kwargs=welch_kwargs):
    
    # PCI transformation vector, size n_channel x p
    #v_pci = pci.v_pci(n_channels)
    # Projection of the data segment t
    #print('compute noise channels')
    #print('for',n_channels,'channels')
    #print('components shape',pci.components.shape)
    if single_link:
        e_pci_n = pci.filter_single_link_data(data_n,n_channels).T
        e_pci_gw = pci.filter_single_link_data(data_gw,n_channels).T
    else:
        e_pci_n = pci.apply_for_channels(data_n,n_channels).T
        e_pci_gw = pci.apply_for_channels(data_gw,n_channels).T
    #print(e_pci_n.shape,e_pci_n.shape)
    #print('compute welch')
    # Welch spectrum matrix for PCI variables from noise
    freqs, e_pci_n_mat = compute_welch_matrix(e_pci_n, **welch_kwargs)
    print(freqs, e_pci_n_mat)
    # Welch spectrum matrix for PCI variables from noise
    freqs, e_pci_gw_mat = compute_welch_matrix(e_pci_gw, **welch_kwargs)
    print(freqs, e_pci_gw_mat)
    # Orthogonalization
    _, s, vh = np.linalg.svd(e_pci_n_mat)
    # Apply the orthogonal transformation to the GW signal    
    e_pci_gw_mat_ortho = multiple_dot(vh, multiple_dot(e_pci_gw_mat, np.swapaxes(vh, 1, 2).conj()))
    # Apply the orthogonal transformation to the noise covariance
    e_pci_n_mat_ortho = multiple_dot(vh, multiple_dot(e_pci_n_mat, np.swapaxes(vh, 1, 2).conj()))
    # Output sensitivity for each variable, size nfreqs x n_channels
    # pci_sens = np.array([np.abs(s[:, j] / e_pci_gw_mat_ortho[:, j, j]) for j in range(n_channels)]).T
    pci_sens = np.array([np.abs(e_pci_n_mat_ortho[:, j, j] / e_pci_gw_mat_ortho[:, j, j]) for j in range(n_channels)]).T
    
    print("Computation completed.")
    if joint:
        pci_sens = 1 / np.sum(1/np.array(pci_sens), axis=1)
    
    return freqs, pci_sens

def process_data(y_noise, y_gw, fs, nhalf, order=1, n_channels=6, pca_y_noise=None, joint=False, pci_kwargs={}, welch_kwargs={}):
    '''
    Apply PCI and compute sensitivity all wrapped up together. If pda_y_noise is provided it is used
    only for computing the PCI and y_noise is used only for the sensitivity application
    '''
    
    if pca_y_noise is None: pca_y_noise=y_noise
        
    print('nhalf:',nhalf, 'sens data_size:',len(y_noise),'pca data_size:',len(pca_y_noise))
    
    # Get the length of the time series
    # data_size = y_noise.shape[0]
    # data_noise = Data.from_instrument(instr_data)

    print('compute PCI')
    pci=pci_filter.PCIFilter(y_noise, fs, maxcompts=10, nhalf=nhalf,order=order,**pci_kwargs)
    
    result=estimate_sensitivity(pci, y_noise, y_gw, n_channels=n_channels, joint=joint, welch_kwargs=welch_kwargs)
    del(pci)
    
    return result


# #### We need to read matching GW simulations results for the empirical sensitivity calculation

# %%


gw_path = workdir+dtpath["point_source"]+datasets["point_source"]+measpath

data_gw = Data.from_gws(gw_path,orbits)

#gw_dataset='tps/y'
#hdf5 = h5py.File(gw_path, 'r')
#dset=hdf5[gw_dataset]

#measurements = {f'isi_{link}': dset[:,ilink] for ilink, link in enumerate(mosas_order)}
#hdf5.close()

y_list = [data_gw.measurements[f'isi_{link}'] for link in mosas_order]
y_gw = np.array(y_list, dtype=np.float64).T[skip:, :]


# #### For comparison, we also need the TDI data

# %%


#tdipath1 = workdir + dtpath + 'noise_tdi1_4Hz.h5'
tdipath2 = workdir +  dtpath['noise'] + datasets['noise'] + '_noise_tdi2_4Hz.h5'
# open hdf5 TDI file
#tdi1 = h5py.File(tdipath1, 'r')
tdi2 = h5py.File(tdipath2, 'r')

x2_noise = tdi2['x'][()] / central_freq
y2_noise = tdi2['y'][()] / central_freq
z2_noise = tdi2['z'][()] / central_freq

tdi2_gw_file = workdir+dtpath["point_source"]+datasets["point_source"]+"_tdi2_4Hz.h5"

hdf5 = h5py.File(tdi2_gw_file, 'r')
x2_gw = hdf5['x'][()]
y2_gw = hdf5['y'][()]
z2_gw = hdf5['z'][()]
hdf5.close()

# Compute welch matrix for the TDI 2.0 noise
e_tdi2_n = np.array([x2_noise[skip:skip+ns],
                     y2_noise[skip:skip+ns],
                     z2_noise[skip:skip+ns]] ).T
freqs, p_tdi2_n_mat = compute_welch_matrix(e_tdi2_n, **welch_kwargs)
# Compute the welch matrix for the TDI 2.0 GW signal
e_tdi2_gw = np.array([x2_gw[skip:skip+ns],
                      y2_gw[skip:skip+ns],
                      z2_gw[skip:skip+ns]]).T
freqs, p_tdi2_gw_mat = compute_welch_matrix(e_tdi2_gw, **welch_kwargs)

# Orthogonalization
u_tdi, s_tdi, vh_tdi = np.linalg.svd(p_tdi2_n_mat)

# Apply the orthogonal transformation to the GW signal
p_tdi2_gw_mat_ortho = multiple_dot(vh_tdi, 

multiple_dot(p_tdi2_gw_mat, np.swapaxes(vh_tdi, 1, 2).conj()))

# Empirical orthogonalization of TDI
mean_tdi2 = 1 / np.sum(
    np.array([np.abs(p_tdi2_gw_mat_ortho[:, j_tdi, j_tdi] / np.abs(s_tdi[:, j_tdi])) for j_tdi in range(3)]), axis=0)


# #### PCI channel sensitivities

# %%


print(ns)
print(y_full[0:ns,:].T)
print(y_gw[0:ns,:].T)
# display(welch_kwargs)


# %%


#Test as it appears in the initial-dev notebook
nh=45
pci_sens_list = [
#    ["PCI-0 std", process_data(y_full[:ns, :].T, y_gw[0:ns, :].T, kwargs['fs'], nh, order=0, n_channels=6, joint=True, pci_kwargs={'sort_by_rms':False},welch_kwargs=kwargs)],
#    ["PCI-1 std", process_data(y_full[:ns, :].T, y_gw[0:ns, :].T, kwargs['fs'], nh, order=1, n_channels=6, joint=True, pci_kwargs={'sort_by_rms':False},welch_kwargs=kwargs)],
#    ["PCI-2 std", process_data(y_full[:ns, :].T, y_gw[0:ns, :].T, kwargs['fs'], nh, order=2, n_channels=6, joint=True, pci_kwargs={'sort_by_rms':False},welch_kwargs=kwargs)],
#    ["PCI-0 rms", process_data(y_full[:ns, :].T, y_gw[0:ns, :].T, kwargs['fs'], nh, order=0, n_channels=6, joint=True, pci_kwargs={'sort_by_rms':True},welch_kwargs=kwargs)],
#    ["PCI-1 rms", process_data(y_full[:ns, :].T, y_gw[0:ns, :].T, kwargs['fs'], nh, order=1, n_channels=6, joint=True, pci_kwargs={'sort_by_rms':True},welch_kwargs=kwargs)],
    ["PCI-2 rms", process_data(y_full[:ns, :].T, y_gw[0:ns, :].T, fs, nh, order=2, n_channels=6, joint=True, pci_kwargs={'sort_by_rms':True},welch_kwargs=welch_kwargs)] 
]


# %%


plotting.plotconfig(lbsize=18, lgsize=16)
_, axes = plt.subplots(1, 1, figsize=(8, 6))
axes.loglog(freqs, np.sqrt(mean_tdi2*ns/fs), 
            linewidth=1, label=r'TDI (empirical)',
            color='black')
for j in range(len(pci_sens_list)):
    #print(j)
    axes.loglog(pci_sens_list[j][1][0], np.sqrt(pci_sens_list[j][1][1]*ns/fs), 
                linewidth=1, label=pci_sens_list[j][0], rasterized=True)
axes.legend(loc='upper left', ncol=2)
axes.set_xlabel("Frequency [Hz]")
axes.set_ylabel(r"Sensitivity $\sqrt{\frac{P_{n}(f)}{P_{\mathrm{GW}}(f)}}$")
axes.set_xlim([1e-3, 1.2])
axes.set_ylim([3e-19, 1e-13])
axes.grid(linewidth=1, which='both', color='gray', linestyle='dotted')
axes.set_title("Sensitivity of combined channels")
plt.show()


# %%


pci_sens_list = []
#pci_sens_list += [
#    ["PCI-"+str(j)+" std", estimate_sensitivity(pca_list[j],y_full[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, welch_kwargs=welch_kwargs)] 
#     for j in range(3)]
pci_sens_list += [
    ["PCI-"+str(j)+" rms", estimate_sensitivity(pca_list_rs[j],y_full[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, welch_kwargs=welch_kwargs)] 
     for j in range(3)]




# %%


plotting.plotconfig(lbsize=18, lgsize=16)
_, axes = plt.subplots(1, 1, figsize=(8, 6))
axes.loglog(freqs, np.sqrt(mean_tdi2*ns/fs), 
            linewidth=1, label=r'TDI (empirical)',
            color='black')
for j in range(len(pci_sens_list)):
    #print(j)
    axes.loglog(pci_sens_list[j][1][0], np.sqrt(pci_sens_list[j][1][1]*ns/fs), 
                linewidth=1, label=pci_sens_list[j][0], rasterized=True)
axes.legend(loc='upper left', ncol=2)
axes.set_xlabel("Frequency [Hz]")
axes.set_ylabel(r"Sensitivity $\sqrt{\frac{P_{n}(f)}{P_{\mathrm{GW}}(f)}}$")
axes.set_xlim([1e-3, 1.2])
axes.set_ylim([3e-19, 1e-13])
axes.grid(linewidth=1, which='both', color='gray', linestyle='dotted')
axes.set_title("Sensitivity of combined channels")
plt.show()


# Reconstructed single link sensitivities

# %%


pci_sens_list = []
#pci_sens_list += [
#    ["PCI-"+str(j)+" std", estimate_sensitivity(pca_list[j],y_full[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, welch_kwargs=welch_kwargs)] 
#     for j in range(3)]
pci_sens_list += [
    ["PCI-"+str(j)+" filtered", estimate_sensitivity(pca_list_rs[j],y_full[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, single_link=True, welch_kwargs=welch_kwargs)] 
     for j in range(3)]


# %%


plotting.plotconfig(lbsize=18, lgsize=16)
_, axes = plt.subplots(1, 1, figsize=(8, 6))
axes.loglog(freqs, np.sqrt(mean_tdi2*ns/fs), 
            linewidth=1, label=r'TDI (raw)',
            color='black')
for j in range(len(pci_sens_list)):
    #print(j)
    axes.loglog(pci_sens_list[j][1][0], np.sqrt(pci_sens_list[j][1][1]*ns/fs), 
                linewidth=1, label=pci_sens_list[j][0], rasterized=True)
axes.legend(loc='upper left', ncol=2)
axes.set_xlabel("Frequency [Hz]")
axes.set_ylabel(r"Sensitivity $\sqrt{\frac{P_{n}(f)}{P_{\mathrm{GW}}(f)}}$")
axes.set_xlim([1e-3, 1.2])
axes.set_ylim([1e-21, 1e-16])
axes.grid(linewidth=1, which='both', color='gray', linestyle='dotted')
axes.set_title("Sensitivity of combined channels")
plt.show()


# ### PCI of single noise components
# 
# #### First we load the single noise data and generate the ys for those

# %%


# Get secondary noises
lockstr = 'locking_n1_12_baseline_noise_'

noises = ['test-mass', 'oms']

y_noise = {}

for nn in noises:
    noisepath = workdir + dtpath['noise']+ lockstr + nn + '_4Hz.h5'
    print(noisepath)
    
    data_noise = Data.from_instrument(noisepath)
    
    y_noise[nn] = build_data_vector(data_noise, skip=skip, dtype=np.float64)


# %%


[y_noise[nn][:ns, :].T, y_full[:ns, :].T]

#[y_noise[nn][:ns, :].T, y_full[nn][:ns, :].T]


# %%


pci_sens_list = []
#pci_sens_list += [
#    ["PCI-"+str(j)+" std", estimate_sensitivity(pca_list[j],y_full[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, welch_kwargs=welch_kwargs)] 
#     for j in range(3)]
labels = ['test-mass', 'oms', 'tm+oms', 'full']

colors = ['']
linestyles = ['']

## add y_full to y_noise in the pci_sens_list to have all of them
## change the labels and the colors

pci_sens_list += [
    ["PCI-"+str(j)+" filtered", estimate_sensitivity(pca_list_rs[j],y_noise[nn][:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, single_link=True, welch_kwargs=welch_kwargs)] 
     for j in range(2,3) for nn in noises];
pci_sens_list += [["PCI-"+str(j)+" filtered", estimate_sensitivity(pca_list_rs[j],y_sec[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, single_link=True, welch_kwargs=welch_kwargs)] 
     for j in range(2,3)];
pci_sens_list += [["PCI-"+str(j)+" filtered", estimate_sensitivity(pca_list_rs[j],y_full[:ns, :].T, y_gw[0:ns, :].T, n_channels=6, joint=True, single_link=True, welch_kwargs=welch_kwargs)] 
     for j in range(2,3)];

plotting.plotconfig(lbsize=18, lgsize=16)
_, axes = plt.subplots(1, 1, figsize=(8, 6))
axes.loglog(freqs, np.sqrt(mean_tdi2*ns/fs), 
            linewidth=2, label=r'TDI (raw)',
            color='grey')
for j in range(len(pci_sens_list)):
    #print(j)
    axes.loglog(pci_sens_list[j][1][0], np.sqrt(pci_sens_list[j][1][1]*ns/fs), 
                linewidth=1.5, label=labels[j], rasterized=True)
axes.legend(loc='upper left', ncol=2)
axes.set_xlabel("Frequency [Hz]")
axes.set_ylabel(r"Sensitivity $\sqrt{\frac{P_{n}(f)}{P_{\mathrm{GW}}(f)}}$")
axes.set_xlim([1e-4, 1.2])
axes.set_ylim([3e-20, 1e-15])
axes.grid(linewidth=1, which='both', color='gray', linestyle='dotted')
axes.set_title("Sensitivity of combined channels")
plt.show()


# %%




