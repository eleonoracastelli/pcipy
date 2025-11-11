#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:39:05 2025

Copyright 2025 E Castelli

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

File: simulation_scenatios.py
Purpose: Run all simulation scenarios for benchmarking of TDI/PCI.
"""

import os
import h5py
from datetime import datetime
import json
import lisainstrument

# %% LOCKING N1-12 default
# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --days 14 --locking "N1-12" --individual --combined')

# %% extract seed from the simulation we just ran
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
    if int(lisainstrument.__version__[0])<2:
        simseed = int(f.attrs['seed'])
    else:
        attributes = json.loads(f.attrs['metadata_json'])
        simseed = attributes['seed']

simseed

# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --seed {seed} --days 14 --locking "N1-12" --individual --combined'.format(seed=simseed))

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --days 14 --locking "N1-12" --baseline --individual --combined'.format(seed=simseed))

# %% LOCKING 'six'

# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --seed {seed} --locking "six" --individual --combined'.format(seed=simseed))

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --seed {seed} --locking "six" --baseline --individual --combined'.format(seed=simseed))

# %% Signal simulation

# %% equalarm tdi 2
os.system('python simulation/signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2')

# %% keplerian
os.system('python simulation/signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2')

# %% All sky simulation

# %% equalarm tdi 2
os.system('python simulation/all_sky_signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2')

# %% keplerian
os.system('python simulation/all_sky_signal_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2')
