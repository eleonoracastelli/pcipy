#!/usr/bin/env python
# coding: utf-8

"""
Copyright 2025 J Baker, E Castelli

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

File: .py
Purpose: Minimal working exaple for TDI filter with piecewise interpolation.

Next level development of the new PCI infrastructure based on TDI filter.
In the new structure, PCI and TDI filters share a common infrastructure.
We develop that infrastrucutre here using the TDI kernels, as TDI is well
understood, and can be compared with a well-established result.

Key common elements for development and testing are the common infrastructure
for estimating sensitivities and for interpolating results from time-local
kernels for application across a wider temporal stretch of data.
"""

# In[1]: 0. Imports

import os
from datetime import datetime

import h5py
import numpy as np
import scipy
import matplotlib.pyplot as plt

from pytdi import Data

from pcipy.data import TimeData
from pcipy.tdi_filter import DeducedTDIFilter as TDIFilter
from pcipy.filter import PiecewiseFilter as PiecewiseFilter

import pyfftw
from pyfftw.interfaces.numpy_fft import rfft

# %% In[2]: 1. Read in data sets
SKIP = 1000

RANGE_IN_HOURS = 48

DATADIR = "/Users/ecastel2/Documents/research/GSFC/pci-inrep/"
WORKDIR = DATADIR+"/simulations/"

MEASPATH = 'measurements_4Hz.h5'
TDIPATH = 'noise_tdi2_4Hz.h5'
ORBPATH = 'keplerian'

substring = [ORBPATH + '_' +
             s for s in ['locking_N1-12_laser_tm_oms_', 'all_sky_gw', 'gw']]
# substring = ['locking_six_laser_tm_oms', 'all_sky_gw', '2_gw']
sims = ["noise", "all_sky", "point_source"]
datasets = dict(zip(sims, substring))

matchfile = {}
dtpath = {}

for n, d in zip(substring, datasets):
    timestamp = []
    matchfile[n] = [f for f in os.listdir(WORKDIR) if n in f]
    for m in matchfile[n]:
        # pick latest date
        timestamp.append(datetime.strptime(m[:11], "%Y-%m-%d_"))
        # print(n, timestamp[n])
        dtpath[d] = max(timestamp).strftime("%Y-%m-%d_")
print(dtpath)

orbits = WORKDIR + ORBPATH + "-orbits.h5"
noise_file_base = dtpath['noise'] + ORBPATH + "_locking_N1-12_laser_tm_oms_"
# gw_file_base = dtpath['point_source'] + \
#     ORBPATH + '_'  # "2025-04-07_keplerian_"

# %% In[3]: 2. Generate Basic single-link noise data

noise_sim_path = WORKDIR + noise_file_base + MEASPATH

# get the central frequency
with h5py.File(noise_sim_path, 'r') as sim:
    central_freq = sim.attrs['central_freq']
    fs = sim.attrs['fs']

data_noise = Data.from_instrument(noise_sim_path)

in_chans = ["isi"]
in_chans = ["isi", "rfi", "tmi"]
# in_chans=["isi","rfi","tmi","isi_sb","rfi_sb"]

# in_chansGW = ["isi"]

def make_names(in_chans):
    """

    Parameters
    ----------
    in_chans : TYPE
        DESCRIPTION.

    Returns
    -------
    names : TYPE
        DESCRIPTION.

    """
    mosas_order = ['12', '23', '31', '13', '21', '32']
    names = []
    for chan in in_chans:
        for link in mosas_order:
            names.append(f'{chan}_{link}')
    return names

# setup instrumental channels
y_names = make_names(in_chans)
# divide by central_freq because simulation channels are in beat-rate units
y_list = [data_noise.measurements[name]/central_freq for name in y_names]
# create a TimeData object
y_n = TimeData(np.array(y_list, dtype=np.float64)[:, SKIP:],
               dt=1/fs,
               t0=SKIP/fs,
               names=y_names)

# In[4]: 3.2 Matching GW simulations results for the empirical sensitivity calculation

# we need these only for the sensitivity calculation

# gw_path = workdir+gw_file_base+"gw_measurements_4Hz.h5"

# data_gw = Data.from_gws(gw_path,orbits)

# y_namesGW = make_names(in_chansGW)
# #note: the GW mesaurement data are already fractional frequency so we don't divide by central_freq
# y_list = [data_gw.measurements[name] for name in y_namesGW]
# print(y_list)
# np.array(y_list, dtype=np.float64)
# y_gw = TimeData(np.array(y_list, dtype=np.float64)[:,skip:],dt=1/fs,t0=skip/fs,names=y_namesGW)

# In[5]: 3.3 For comparison, we also need the TDI data

# TDI noise from file
tdipath2 = WORKDIR + noise_file_base + TDIPATH

tdi2 = h5py.File(tdipath2, 'r')
x2_noise = tdi2['x'][()] / central_freq
y2_noise = tdi2['y'][()] / central_freq
z2_noise = tdi2['z'][()] / central_freq
tdi2.close()

XYZ_file_noise = TimeData(np.array([x2_noise, y2_noise, z2_noise],
                                   dtype=np.float64)[:, SKIP:],
                          dt=1/fs,
                          t0=SKIP/fs,
                          names='XYZ')

# #GW
# tdi2_gw_file = workdir+gw_file_base+"gw_tdi2_4Hz.h5"

# hdf5 = h5py.File(tdi2_gw_file, 'r')
# x2_gw = hdf5['x'][()] #/ central_freq  (gw simulation already in fractional freq)
# y2_gw = hdf5['y'][()] #/ central_freq
# z2_gw = hdf5['z'][()] #/ central_freq
# hdf5.close()

# XYZ_file_gw = TimeData(np.array([x2_gw,y2_gw,z2_gw],
# dtype=np.float64)[:,skip:],
# dt=1/fs,t0=skip/fs,names='XYZ')

# In[6]: Define data range for this study and set up fourier transforms

ibuff = 250  # edge buffer for kernels

ns = int(RANGE_IN_HOURS * 3600 * fs)  # less for dev

# In[7]: Construct data sets
ytest_n = y_n.get_range(0, ns+2*ibuff)
# ytest_gw = y_gw.get_range(0,ns+2*ibuff)

TDItest_n = XYZ_file_noise.get_range(0, ns+2*ibuff)
# TDItest_gw = XYZ_file_gw.get_range(0,ns+2*ibuff)

# In[8]: Fourier Transform
pyfftw.interfaces.cache.enable()

# Transform PCI variables to Fourier domain
wd = np.blackman(ns)
# Transform simulated TDI variables to Fourier domain
k2 = np.sum(wd**2)

def do_ft(dataset, fs=fs):
    """
    Applies the Fourier Transform to a specific TimeData class dataset.

    Parameters
    ----------
    dataset : TimeData  object.
        Time-series dataset
    fs : float, optional
        Sampling frequency. The default is fs.

    Returns
    -------
    f : ndarray
        Selected frequencies.
    dataset_fft : ndarray
        FT dataset.

    """
    ns = dataset.n_samples()
    wd = np.blackman(ns)
    f = np.fft.rfftfreq(ns) * fs
    sel = f > 0
    # print(sel)
    dataset_fft = [rfft(wd * dataset.data[i]) * np.sqrt(2/(fs*k2))
                   for i in range(dataset.n_channels())]

    return f[sel], np.array(dataset_fft)[:, sel]

print(ytest_n.data.shape)

def do_dec(x, ibuff = ibuff):
    """
    Decimates data under analysis.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """

    return scipy.signal.decimate(x, q=idec)[ibuff:-ibuff]
# epsd=do_dec(eps[:])
# In[12]: Piecewise linear interpolation

tbuff = 0
tstart = ytest_n.t0 - tbuff
tend = ytest_n.t0+ytest_n.dt*ytest_n.n_samples() + tbuff

pwFilter = PiecewiseFilter(tstart, tend, 8, TDIFilter,
                           measurements_data=data_noise,
                           in_chans=in_chans,
                           method='linear')

# In[13]:

XYZpw = pwFilter.apply_filter(ytest_n, check=True)
print(ytest_n.n_samples(), XYZpw.n_samples())
print(ytest_n.dt, XYZpw.dt)
print(ytest_n.t0, XYZpw.t0)

# In[14]: Construct filters

#  test timing
import time

def timed_filter(data_noise, t, in_chans):
    start_time = time.time()  # Start timing
    result = TDIFilter(data_noise, eval_time=t, in_chans=in_chans, method='linear')
    end_time = time.time()      # End timing
    print(f"Time taken for: {end_time - start_time:.4f} seconds")
    return result

nks = [1, 2, 4, 8, 16, 32] # 3.5 s per iteration
tks = [PiecewiseFilter.select_kernel_times(ytest_n, nk) for nk in nks]
filtsets = [[TDIFilter(data_noise, eval_time=t, in_chans=in_chans, method='linear')
             for t in PiecewiseFilter.select_kernel_times(ytest_n, nk)]
            for nk in nks]

# filtsets = [[timed_filter(data_noise, t, in_chans)
#              for t in PiecewiseFilter.select_kernel_times(ytest_n, nk)]
#             for nk in nks]


# In[15]: Apply filters for linear interpolation
XYZset_n = [PiecewiseFilter.piecewise_linear_apply_filter_set(
    ytest_n, filtsets[i], tks[i]) for i in range(len(nks))]
for s in XYZset_n:
    print(s.n_samples())

# In[16]: View data in decimated in the time domain
# Comparison with TDI variables
_, axes = plt.subplots(1, 1, figsize=(8, 6))

idec = 10
iev = 100
ibegin = 0
iend = 1000000

ixyz = 0
# iskip=249
iskip = 282

ioff = len(TDItest_n.data[ixyz])

delta = XYZset_n[0].data[ixyz]-TDItest_n.data[ixyz][iskip:]
print('min', np.min(delta))
print('max', np.max(delta))
print('mean', np.mean(delta))
print('std', np.std(delta))

data = do_dec(TDItest_n.data[ixyz][iskip:])
times = do_dec(TDItest_n.get_times()[iskip:])
refdata = data
reftimes = times
print('refshape:', refdata.shape)

for iset in range(len(XYZset_n)):
    data = do_dec(XYZset_n[iset].data[ixyz])
    times = do_dec(XYZset_n[iset].get_times())
    print('time offest check:', np.mean(times-reftimes))
    print('data shape:', data.shape)
    axes.semilogy(times[ibegin:iend:iev], np.abs(
        data[ibegin:iend:iev]), linewidth=1, linestyle ='-', label='n='+str(nks[iset])+':'+'XYZ'[ixyz])
    # axes.semilogy(times[ibegin:iend:iev], 1e-25+np.abs((data-refdata)[ibegin:iend:iev]),
                  # linewidth=1, label='n='+str(nks[iset])+':'+'XYZ'[ixyz]+' - ref')
    # refdata=data


axes.semilogy(times[ibegin:iend:iev], np.abs(
    data[ibegin:iend:iev]),'k--', linewidth=2, label='TDI2:'+'XYZ'[ixyz])

data = do_dec(XYZpw.data[ixyz])
times = do_dec(XYZpw.get_times())
print('time offest check:', np.mean(times-reftimes))
print('data shape:', data.shape)
axes.semilogy(times[ibegin:iend:iev], np.abs(
    data[ibegin:iend:iev]), 'k', linewidth=1, linestyle ='-', label='piecewise:'+'XYZ'[ixyz])
# axes.semilogy(times[ibegin:iend:iev], 1e-25+np.abs((data-refdata)
#               [ibegin:iend:iev]), linewidth=1, ls='--', label='pw:'+'XYZ'[ixyz]+' - ref')

axes.set_xlabel(r"t")
axes.set_ylabel(r"smoothed chan")
plt.legend(ncols=3)#loc='upper left', fontsize=16)
# axes.set_ylim([1e-24, 1e-17])
# axes.set_xlim([1e-4, 0.1])
# plt.title("PCI 2.0 (" + str(int(ns/fs/3600))+" hours, nh = "+str(nhalf)+")")
plt.grid(linewidth=1.0,
         color='gray',
         linestyle='dotted',
         which='both',
         axis='both')

# %% Comparison with TDI variables
_, axes = plt.subplots(1, 1, figsize=(8, 6))

idec = 10
iev = 100
ibegin = 0
iend = 1000000

ixyz = 0
# iskip=249
iskip = 282

ioff = len(TDItest_n.data[ixyz])

delta = XYZset_n[0].data[ixyz]-TDItest_n.data[ixyz][iskip:]
print('min', np.min(delta))
print('max', np.max(delta))
print('mean', np.mean(delta))
print('std', np.std(delta))

data = do_dec(TDItest_n.data[ixyz][iskip:])
times = do_dec(TDItest_n.get_times()[iskip:])
refdata = data
reftimes = times
print('refshape:', refdata.shape)


for iset in range(len(XYZset_n)):
    data = do_dec(XYZset_n[iset].data[ixyz])
    times = do_dec(XYZset_n[iset].get_times())
    print('time offest check:', np.mean(times-reftimes))
    print('data shape:', data.shape)
    # axes.semilogy(times[ibegin:iend:iev], np.abs(
    #     data[ibegin:iend:iev]), linewidth=1, label='n='+str(nks[iset])+':'+'XYZ'[ixyz])
    axes.semilogy(times[ibegin:iend:iev], 1e-25+np.abs((data-refdata)[ibegin:iend:iev]),
                  linewidth=1, linestyle ='-', label='n='+str(nks[iset])+':'+'XYZ'[ixyz]+' - ref')
    # refdata=data

# axes.semilogy(times[ibegin:iend:iev], np.abs(
#     data[ibegin:iend:iev]), 'k--', linewidth=2, label='TDI2:'+'XYZ'[ixyz])

data = do_dec(XYZpw.data[ixyz])
times = do_dec(XYZpw.get_times())
print('time offest check:', np.mean(times-reftimes))
print('data shape:', data.shape)
# # axes.semilogy(times[ibegin:iend:iev], np.abs(
# #     data[ibegin:iend:iev]), linewidth=1, label='pw:'+'XYZ'[ixyz])
axes.semilogy(times[ibegin:iend:iev], 1e-25+np.abs((data-refdata)
              [ibegin:iend:iev]), 'k--', linewidth=1, ls='-', label='pw:'+'XYZ'[ixyz]+' - ref')

axes.set_xlabel(r"t")
axes.set_ylabel(r"smoothed chan")
plt.legend(ncols=3)#loc='upper left', fontsize=16)
# axes.set_ylim([1e-24, 1e-17])
# axes.set_xlim([1e-4, 0.1])
# plt.title("PCI 2.0 (" + str(int(ns/fs/3600))+" hours, nh = "+str(nhalf)+")")
plt.grid(linewidth=1.0,
         color='gray',
         linestyle='dotted',
         which='both',
         axis='both')

# In[21]: View data in freq domain

# Comparison with TDI variables
_, axes = plt.subplots(1, 1, figsize=(8, 6))

iev = 30
i0 = 0

ixyz = 0

for iset in range(len(XYZset_n)):
    f, data = do_ft(XYZset_n[iset])
    axes.loglog(f[i0::iev], (1+iset*0)*np.abs(data[ixyz, i0::iev]),
                linewidth=1, label='n='+str(nks[iset])+':'+'XYZ'[ixyz])

for iy in range(1):
    f, data = do_ft(y_n)
    axes.loglog(f[i0::iev], (1+iset*0)*np.abs(data[iy, i0::iev]),
                linewidth=1, label="single-link-"+str(iy))


f, data = do_ft(TDItest_n.get_range(iskip, None))
dataref = data
axes.loglog(f[i0::iev], np.abs(dataref[ixyz, i0::iev]),
            'k', linewidth=2, label='TDI2(trim):'+'XYZ'[ixyz])

axes.set_xlabel("f")
axes.set_ylabel(r"fourier transform")
plt.legend()#loc='upper left', fontsize=16)
axes.set_ylim([1e-24, 1e-11])
axes.set_xlim([.3e-4, 1])
# plt.title("PCI 2.0 (" + str(int(ns/fs/3600))+" hours, nh = "+str(nhalf)+")")
plt.grid(linewidth=1.0,
         color='gray',
         linestyle='dotted',
         which='both',
         axis='both')


# %% Comparison with TDI variables
_, axes = plt.subplots(1, 1, figsize=(8, 6))

iev = 30
i0 = 0

ixyz = 0

for iset in range(len(XYZset_n)):
    f, data = do_ft(XYZset_n[iset])
    data = data-dataref
    axes.loglog(f[i0::iev], (1+iset*0)*np.abs(data[ixyz, i0::iev]),
                ls = '--',
                linewidth=1, label='n='+str(nks[iset])+':'+'XYZ'[ixyz]+' err')


# f, data = do_ft(TDItest_n.get_range(iskip, None))
# dataref = data
axes.loglog(f[i0::iev], np.abs(dataref[ixyz, i0::iev]),
            'k', linewidth=2, label='TDI2(trim):'+'XYZ'[ixyz])

axes.set_xlabel("f")
axes.set_ylabel(r"fourier transform")
plt.legend()#loc='upper left', fontsize=16)
axes.set_ylim([1e-24, 1e-11])
axes.set_xlim([.3e-4, 1])
# plt.title("PCI 2.0 (" + str(int(ns/fs/3600))+" hours, nh = "+str(nhalf)+")")
plt.grid(linewidth=1.0,
         color='gray',
         linestyle='dotted',
         which='both',
         axis='both')

# %% Comparison with TDI variables
_, axes = plt.subplots(1, 1, figsize=(8, 6))

iev = 30
i0 = 0

ixyz = 0
for iset in range(len(XYZset_n)):
    f, data = do_ft(XYZset_n[iset])
    data = data-dataref
    axes.loglog(f[i0::iev], (1+iset*0)*np.abs(data[ixyz, i0::iev])/np.abs(dataref[ixyz, i0::iev]),
                ls = '-',
                linewidth=1, label='n='+str(nks[iset])+':'+'XYZ'[ixyz]+' rel err')

axes.loglog(f[i0::iev], np.abs(dataref[ixyz, i0::iev])/np.abs(dataref[ixyz, i0::iev]),
            'k', ls = '-',
            linewidth=2, label='TDI2(trim):'+'XYZ'[ixyz]+' rel err')

axes.set_xlabel("f")
axes.set_ylabel(r"fourier transform")
plt.legend()
# axes.set_ylim([1e-24, 1e-11])
axes.set_xlim([.3e-4, 1])
# plt.title("PCI 2.0 (" + str(int(ns/fs/3600))+" hours, nh = "+str(nhalf)+")")
plt.grid(linewidth=1.0,
         color='gray',
         linestyle='dotted',
         which='both',
         axis='both')
