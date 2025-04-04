#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:39:05 2025

@author: ecastel2
"""

import os

# %% LOCKING N1-12 default
# %% equalarm tdi1
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 1 --locking "N1-12" --individual --combined')

# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --locking "N1-12" --individual --combined')

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --locking "N1-12" --baseline --individual --combined')

# %% LOCKING 'six'
# %% equalarm tdi1
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 1 --locking "six" --individual --combined')

# %% equalarm tdi2
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits equalarm --tdi 2 --locking "six" --individual --combined')

# %% keplerian
os.system('python simulation/noise_simulation.py ../../../research/GSFC/pci-inrep/simulations --orbits keplerian --tdi 2 --locking "six" --baseline --individual --combined')

