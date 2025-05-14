#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test PSD estimation script

@author: ecastel2
"""

import os
import numpy as np
from pcipy.utils import PSDstats
from datetime import datetime
import h5py
from pytdi import Data
import scipy.signal as signal
import matplotlib.pyplot as plt

# %% test PSD estimation


workdir = "/Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations/"

measpath = 'measurements_4Hz.h5'
secondpath = 'noise_tdi2_4Hz.h5'

#orbits = workdir+"keplerian-orbits.h5"
orbits = workdir+"equalarm-orbits.h5"

substring = ['equalarm_' + s for s in ['locking_N1-12_laser_tm_oms_', 'all_sky_gw', 'gw']]
#substring = ['locking_six_laser_tm_oms', 'all_sky_gw', '2_gw']
sims = ["noise", "all_sky", "point_source"]
datasets = dict(zip(sims, substring))

matchfile={}
dtpath={}

# %%
for n,d in zip(substring,datasets):
    timestamp=[]
    matchfile[n] = [f for f in os.listdir(workdir) if n in f]
    for m in matchfile[n]:
        # pick latest date
        timestamp.append(datetime.strptime(m[:11], "%Y-%m-%d_"))
        #print(n, timestamp[n])
        dtpath[d] = max(timestamp).strftime("%Y-%m-%d_")
dtpath
# %%
# Get noise simulation measurements
simpath = workdir + dtpath['noise'] + datasets['noise'] + measpath
tdipath = workdir + dtpath['noise'] + datasets['noise'] + secondpath
print(simpath)
# load hdf5 file to read data attrs
sim = h5py.File(simpath, 'r')
# load data
tdi2 = h5py.File(tdipath, 'r')
# %%
central_freq = sim.attrs['central_freq']
dt = sim.attrs['dt']
n_data = sim.attrs['size']
t0 = sim.attrs['t0']

fs = 1/dt

x_noise = tdi2['x'][()] / central_freq
y_noise = tdi2['y'][()] / central_freq
z_noise = tdi2['z'][()] / central_freq



# %%
pytdi_trim=1000

f, xpsd = signal.welch(x_noise[pytdi_trim:], **{"fs": fs,
          "window": 'blackman',
          "nperseg": n_data-pytdi_trim}, )

plt.plot(x_noise[pytdi_trim:])
plt.show()

# %%

dat = tuple(x_noise[pytdi_trim:])

datpsd = PSDstats(dat, Tmax = (n_data-pytdi_trim), fmax = 3e-2, fs=fs)

Sxx, Smin, Smax = datpsd.getPSDerrors(cval=0.68)

# %%
fig, ax = plt.subplots(1,1, dpi=120)
ax.plot(f, xpsd, color='tab:orange', zorder=1)
ax.errorbar(datpsd.freqs, Sxx, 
            yerr=[Sxx-Smin, Smax-Sxx],
            fmt='.', capsize = 3, color="k", zorder=2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
