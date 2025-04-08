#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:39:05 2025

@author: ecastel2
"""

import os
import h5py
from datetime import datetime

# %% LOCKING N1-12 default
# %% equalarm tdi1
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 1 --locking "N1-12" --individual --combined')
# %%
workdir = "/Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations/"

n = 'measurements'

timestamp=[]
matchfile = [f for f in os.listdir(workdir) if n in f]
for m in matchfile:
    # pick latest date
    timestamp.append(datetime.strptime(m[:11], "%Y-%m-%d_"))
    #print(n, timestamp[n])
    dtpath = max(timestamp).strftime("%Y-%m-%d_")

fname = workdir + dtpath + 'equalarm_locking_N1-12_laser_tm_oms_measurements_4Hz.h5'

with h5py.File(fname) as f:
    simseed = int(f.attrs['seed'])
    
simseed

# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --seed {seed} --locking "N1-12" --individual --combined'.format(seed=simseed))

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "N1-12" --baseline'.format(seed=simseed))

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "N1-12" --baseline --individual'.format(seed=simseed))

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "N1-12" --baseline --combined'.format(seed=simseed))

# %% LOCKING 'six'
# %% equalarm tdi1
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 1 --seed {seed} --locking "six" --individual --combined'.format(seed=simseed))

# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --seed {seed} --locking "six" --individual --combined'.format(seed=simseed))

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "six" --baseline'.format(seed=simseed))

# %% keplerian individual
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "six" --baseline --individual'.format(seed=simseed))

# %% keplerian combined
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "six" --baseline --combined'.format(seed=simseed))


# %% Signal simulation
# %% equalarm tdi 1
os.system('python simulation/signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 1')

# %% equalarm tdi 2
os.system('python simulation/signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2')

# %% keplerian
os.system('python simulation/signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2')

# %% All sky simulation

# %% equalarm tdi 1
os.system('python simulation/all_sky_signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 1')

# %% equalarm tdi 2
os.system('python simulation/all_sky_signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2')

# %% keplerian
os.system('python simulation/all_sky_signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2')