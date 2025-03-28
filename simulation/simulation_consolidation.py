#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 16:08:39 2025

Consolidate and test simulation scripts for PCIpy.
@author: ecastel2
"""
# %%
import logging
import h5py
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from scipy import signal
from datetime import datetime
from lisainstrument import Instrument
from lisaorbits import KeplerianOrbits, EqualArmlengthOrbits
from pytdi.michelson import X1, Y1, Z1, X2, Y2, Z2
from pytdi import Data

# %% 

# To print the logs
logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

# create figure objects
figs = []

# define function to export figures to pdf 
# pp = PdfPages('simulation_consolidation_plots.pdf')

# %%
plt.figure(figsize=(10, 8)) 

plt.axis('off')
plt.text(0.5,0.6,"Plots from simulation_consolidation.py.",ha='center',va='top')
plt.text(0.3, 0.45,'''3 days of only noise simulation.
Different noise configurations (laser + tm + oms noise) vs (baseline).
Keplerian orbits, with comparison of old vs new orbit files.
Refer to plot caption for more details.''',
         ha='left',va='bottom')
# pp.savefig()
# plt.close()
    
# if __name__ == "__main__":

# %% 

# To print the logs
logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

# create figure objects
figs = []

# define function to export figures to pdf 
# pp = PdfPages('simulation_consolidation_plots.pdf')

# %%
plt.figure(figsize=(10, 8)) 

plt.axis('off')
plt.text(0.5,0.6,"Plots from simulation_consolidation.py.",ha='center',va='top')
plt.text(0.3, 0.45,'''3 days of only noise simulation.
Different noise configurations (laser + tm + oms noise) vs (baseline).
Keplerian orbits, with comparison of old vs new orbit files.
Refer to plot caption for more details.''',
         ha='left',va='bottom')
# pp.savefig()
# plt.close()
    


# %%
# Sampling time
dt = 0.25
# Sampling frequency
fs = 1 / dt
# Data size: 24 hours or 2 days
tobs = 3 * 24 * 3600
n_data = int(tobs * fs)

t0 = 2173211130.0 # s  datetime.datetime(2038, 11, 12, 16, 45, 30)

print("Data size: " + str(n_data))
print("Data duration: " + str(tobs/3600) + " hours")

# Central frequency
central_freq = 281600000000000.0


X, Y, Z = X2, Y2, Z2

# %% set up proper time grid for simulation

pytdi_trim = 1000
pytdi_t0 = t0 - pytdi_trim * dt
pytdi_size = n_data + pytdi_trim

instrument_t0 = pytdi_t0
instrument_size = pytdi_size

orbits_dt = 100_000
orbits_trim = 100
orbits_t0 = t0 - pytdi_trim * dt - orbits_trim * orbits_dt
orbits_size = np.ceil(3600 * 24 * 365 / orbits_dt) # a year, for illustration purposes

# %% Choose orbit file
# TO DO double check that the keplerian orbits match the orbits simulation 
# parameters in LISA-LCST-SGS-RP-006 
datadir = '/Users/ecastel2/Documents/research/GSFC/simulation-tests/orbits/'


OLD_ORBITS = False
KEPL_ORBITS = True
EQUAL_ORBITS = False



if OLD_ORBITS:
    orbits = datadir+"keplerian-orbits.h5"
elif KEPL_ORBITS:
    orbitsobj = KeplerianOrbits()
    orbits = datadir+"new-orbits.h5"
    orbitsobj.write(orbits, dt=orbits_dt, size=orbits_size, t0=orbits_t0, mode="w")
elif EQUAL_ORBITS:
    orbitsobj = EqualArmlengthOrbits()
    orbits = datadir+"new-equal-orbits.h5"
    orbitsobj.write(orbits, dt=orbits_dt, size=orbits_size, t0=orbits_t0, mode="w")

del orbitsobj


# %% ###################
# Old orbits
########################

# orbits = "/work/SC/lisa/baghiq/orbits/keplerian-orbits.h5"
with h5py.File(orbits) as f:
    orbit_t0 = f.attrs['t0']


# noise parameters to turn selected noises back on
# locking='six'
# oms_asds=(6.35e-12, 1.25e-11, 1.42e-12, 3.38e-12, 3.32e-12, 7.90e-12)        
# tm_asds=2.4E-15
# laser_asds=30


print("*************************************************")
print("Using LISA-LCST-SGS-RP-006 baseline configuration")
print("*************************************************")
locking='N1-12' # default configuration used in LISA-LCST-SGS-RP-006
ranging_asds=3e-9
ranging_b = [ranging_asds * x for x in (2, -1, -1.5, 3, 0.5, 0.75)]
ranging_biases = dict(zip(['12', '23', '31',
        '13', '32', '21'], ranging_b))
# backlink_asds=3e-12
# backlink_fknees=2e-3
# clock_asds=6.32e-14
clock_offsets=(1.5, -0.75, 0.1)
# clock_freqoffsets="default"
# clock_freqlindrifts="default"
# clock_freqquaddrifts="default"
# modulation_asds_left=5.2E-14
# modulation_asds_right=5.2E-13
moc_time_correlation_asds = 0.042

# Instantiate LISA instrument
instr = Instrument(size=instrument_size,
                    dt=dt,
                    t0=orbit_t0+1000,#instrument_t0, 
                    lock=locking, 
                    orbits=orbits)
        
# Disable all noises
instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms'])
instr.simulate()
simseed = instr.seed

# %% evaluate TDI for old orbits

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq
# %% plot for old orbits
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
ax[0].set_title('Laser + tm + oms noise - old orbit file')
=======
ax[0].set_title('Old orbits')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
ax[0].set_title('Laser + tm + oms noise - old orbit file')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([1e-25, 1e-17])
ax[1].set_xlabel('Frequency [Hz]')
<<<<<<< HEAD
<<<<<<< HEAD
ax[1].set_ylabel(r'$S_\text{TDI}(f)$')
# pp.savefig(fig)

=======
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
ax[1].set_ylabel(r'$S_\text{TDI}(f)$')
# pp.savefig(fig)

>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)


# %% ###################
# Old orbits - more settings
########################

# Instantiate LISA instrument
instr = Instrument(size=instrument_size,
                   seed = simseed,
                    dt=dt,
                    t0=orbit_t0+1000, 
                    lock=locking, 
                    orbits=orbits, 
                    aafilter=('kaiser', 240, 1.1, 2.9),
                    clock_offsets={'1':clock_offsets[0],
                                   '2':clock_offsets[1],
                                   '3':clock_offsets[2]},
                    ranging_biases=ranging_biases,
                    moc_time_correlation_asds = moc_time_correlation_asds)
        
# Disable all noises
instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms'])
instr.simulate()

# %% evaluate TDI for old orbits, more settings

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq
# %% plot for old orbits - more settings
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
ax[0].set_title('Laser + tm + oms noise - old orbits - full simulation settings')
=======
ax[0].set_title('Old orbits - more settings')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
ax[0].set_title('Laser + tm + oms noise - old orbits - full simulation settings')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([1e-25, 1e-17])
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)

# %% ###################
# New orbits 
########################

# Instantiate LISA instrument
instr = Instrument(size=instrument_size,
                   seed = simseed,
                    dt=dt,
                    t0=instrument_t0, 
                    lock=locking, 
                    orbits=datadir+"new-orbits.h5", 
                    aafilter=('kaiser', 240, 1.1, 2.9),
                    clock_offsets={'1':clock_offsets[0],
                                   '2':clock_offsets[1],
                                   '3':clock_offsets[2]},
                    ranging_biases=ranging_biases,
                    moc_time_correlation_asds = moc_time_correlation_asds)
        
# Disable all noises
instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms'])
instr.simulate()

# %% evaluate TDI for new orbits

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq
# %% plot for new orbits
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
ax[0].set_title('Laser + tm + oms noise - new orbits')
=======
ax[0].set_title('New orbits')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
ax[0].set_title('Laser + tm + oms noise - new orbits')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([1e-25, 1e-17])
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)


# %% disable laser noise without rerunning the simulation

# Generate secondary noises HDF5 file
# disable laser noise to simulate secondary noises
instr.disable_all_noises(excluding=['test-mass', 'oms'])
instr.simulate()


# %% TDI

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq
# %% plot
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
ax[0].set_title('Laser + tm + oms noise - new orbits - deactivate laser noise from previous sim')
=======
ax[0].set_title('New orbits - deactivate laser no rerunning')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
ax[0].set_title('Laser + tm + oms noise - new orbits - deactivate laser noise from previous sim')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([1e-25, 1e-17])
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)


# %% Individual noise contributions

print("Saving individual noise contribution")

noises = ['oms', 'test-mass']

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
=======

>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
sumpsd = np.ndarray(xpsd.shape)

for nn in noises:
    # Instantiate LISA instrument
    instr = Instrument(size=instrument_size,
                        seed = simseed,
                        dt=dt,
                        t0=instrument_t0, 
                        lock=locking, 
                        orbits=datadir+"new-orbits.h5", 
                        aafilter=('kaiser', 240, 1.1, 2.9),
                        clock_offsets={'1':clock_offsets[0],
                                       '2':clock_offsets[1],
                                       '3':clock_offsets[2]},
                        ranging_biases=ranging_biases,
                        moc_time_correlation_asds = moc_time_correlation_asds
                    )
            
    instr.disable_all_noises(excluding=nn) 
    instr.simulate()
    #

    data_noise = Data.from_instrument(instr)

    # Build other 2.0 Michelson variables
    X_data = X.build(**data_noise.args)
    # Y_data = Y.build(**data_noise.args)
    # Z_data = Z.build(**data_noise.args)
    # Apply TDI 2.0
    x_noise = X_data(data_noise.measurements) / central_freq
    # y_noise = Y_data(data_noise.measurements)  / central_freq
    # z_noise = Z_data(data_noise.measurements)  / central_freq
    
    #
    tdi_times = np.arange(n_data)*dt

# fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)

<<<<<<< HEAD
<<<<<<< HEAD
    ax[0].set_title('Laser + tm + oms noise - new orbits - individual noises')
=======
    ax[0].set_title('Keplerian orbits - individual noises')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
    ax[0].set_title('Laser + tm + oms noise - new orbits - individual noises')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
    ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X {nn}'.format(nn=nn), alpha = 0.5)
    # ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
    # ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
    ax[0].grid()
    ax[0].legend()
    ax[0].set_xlabel('Time [s]')
    ax[0].set_ylabel('TDI')
     
    
    kwargs = {"fs": fs,
              "window": 'blackman',
              "nperseg": instrument_size,
              "detrend": 'constant',
              "return_onesided": True,
              "scaling": 'density'}
    
    f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
    # f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
    # f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)
    
    ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X {nn}'.format(nn=nn))
    # ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
    # ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
    ax[1].grid()
    ax[1].legend()
    # ax[1].set_xlim([1e-5, 1e-1])
    # ax[1].set_ylim([5e-27, 5e-17])
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
    # pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
    # pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
 
    sumpsd += xpsd


# %% disable laser noise by rerunning the simulation

# Instantiate LISA instrument
instr = Instrument(size=instrument_size,
                   seed = simseed,
                    dt=dt,
                    t0=instrument_t0, 
                    lock=locking, 
                    orbits=datadir+"new-orbits.h5", 
                    aafilter=('kaiser', 240, 1.1, 2.9),
                    clock_offsets={'1':clock_offsets[0],
                                   '2':clock_offsets[1],
                                   '3':clock_offsets[2]},
                    ranging_biases=ranging_biases,
                    moc_time_correlation_asds = moc_time_correlation_asds)

instr.disable_all_noises(excluding=['test-mass', 'oms'])
instr.simulate()

# %%

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq
# %%
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
ax[0].set_title('Laser + tm + oms noise - new orbits - only tm and oms')
=======
ax[0].set_title('New orbits - only tm and oms')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
ax[0].set_title('Laser + tm + oms noise - new orbits - only tm and oms')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([5e-27, 5e-12])
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)


# %%

fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
ax.set_title('Laser + tm + oms noise - new orbits')
=======

>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
ax.set_title('Laser + tm + oms noise - new orbits')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax.loglog(f[1:], np.sqrt(xpsd[1:]), label = 'TDI X (tm + oms)')
ax.loglog(f[1:], np.sqrt(sumpsd[1:]), label = 'TDI X tm + TDI X oms')
ax.grid()
ax.legend()
# ax[0].set_xlim([1e-5, 1e-1])
ax.set_ylim([5e-27, 5e-17])
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)


# %%

noises = ['oms', 'test-mass', 'laser']

for nn in noises:
    # Instantiate LISA instrument
    instr = Instrument(size=instrument_size,
                   seed = simseed,
                    dt=dt,
                    t0=instrument_t0, 
                    lock=locking, 
                    orbits=datadir+"new-orbits.h5", 
                    aafilter=('kaiser', 240, 1.1, 2.9),
                    clock_offsets={'1':clock_offsets[0],
                                   '2':clock_offsets[1],
                                   '3':clock_offsets[2]},
                    ranging_biases=ranging_biases,
                    moc_time_correlation_asds = moc_time_correlation_asds
                    )
            
    instr.disable_all_noises(excluding=nn) 
    instr.simulate()
    #

    data_noise = Data.from_instrument(instr)


    # Build other 2.0 Michelson variables
    X_data = X.build(**data_noise.args)
    Y_data = Y.build(**data_noise.args)
    Z_data = Z.build(**data_noise.args)
    # Apply TDI 2.0
    x_noise = X_data(data_noise.measurements) / central_freq
    y_noise = Y_data(data_noise.measurements)  / central_freq
    z_noise = Z_data(data_noise.measurements)  / central_freq
    
    #
    tdi_times = np.arange(n_data)*dt

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
     
    ax[0].set_title('Laser + tm + oms noise - new orbits - {nn} noise'.format(nn=nn))
=======

    ax[0].set_title('Keplerian orbits - {nn} noise'.format(nn=nn))
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
     
    ax[0].set_title('Laser + tm + oms noise - new orbits - {nn} noise'.format(nn=nn))
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
    ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
    ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
    ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
    ax[0].grid()
    ax[0].legend()
    ax[0].set_xlabel('Time [s]')
    ax[0].set_ylabel('TDI')
     
    
    kwargs = {"fs": fs,
              "window": 'blackman',
              "nperseg": instrument_size,
              "detrend": 'constant',
              "return_onesided": True,
              "scaling": 'density'}
    
    f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
    f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
    f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)
    
    ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
    ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
    ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
    ax[1].grid()
    ax[1].legend()
    # ax[1].set_xlim([1e-5, 1e-1])
    ax[1].set_ylim([5e-27, 5e-17])
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
    # pp.savefig(fig)
=======

>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
    # pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)

# %% New orbits - baseline simulation
########################

# Instantiate LISA instrument
instr = Instrument(size=instrument_size,
                   seed = simseed,
                    dt=dt,
                    t0=instrument_t0, 
                    lock=locking, 
                    orbits=datadir+"new-orbits.h5", 
                    aafilter=('kaiser', 240, 1.1, 2.9),
                    clock_offsets={'1':clock_offsets[0],
                                   '2':clock_offsets[1],
                                   '3':clock_offsets[2]},
                    ranging_biases=ranging_biases,
                    moc_time_correlation_asds = moc_time_correlation_asds)
        
# Disable all noises
instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms', 'ranging', 'backlink', 'clock', 'modulation'])
instr.simulate()

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq

#
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 

ax[0].set_title('Baseline noises - new orbits')
=======
ax[0].set_title('New orbits - baseline simulation')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 

ax[0].set_title('Baseline noises - new orbits')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([5e-27, 5e-12])
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)

# %%
# Disable all noises
instr.disable_all_noises(excluding=['test-mass', 'oms', 'ranging', 'backlink', 'clock', 'modulation'])
instr.simulate()

data_noise = Data.from_instrument(instr)

# Build other 2.0 Michelson variables
X_data = X.build(**data_noise.args)
Y_data = Y.build(**data_noise.args)
Z_data = Z.build(**data_noise.args)
# Apply TDI 2.0
x_noise = X_data(data_noise.measurements) / central_freq
y_noise = Y_data(data_noise.measurements) / central_freq
z_noise = Z_data(data_noise.measurements) / central_freq

#
tdi_times = np.arange(n_data)*dt

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 

ax[0].set_title('Baseline noises - new orbits - only secondary noises')
=======
ax[0].set_title('New orbits - baseline simulation only secondary noises')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 

ax[0].set_title('Baseline noises - new orbits - only secondary noises')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X', alpha = 0.5)
ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('TDI')
 

kwargs = {"fs": fs,
          "window": 'blackman',
          "nperseg": instrument_size,
          "detrend": 'constant',
          "return_onesided": True,
          "scaling": 'density'}

f, xpsdsec = signal.welch(x_noise[pytdi_trim:], **kwargs)
f, ypsdsec = signal.welch(y_noise[pytdi_trim:], **kwargs)
f, zpsdsec = signal.welch(z_noise[pytdi_trim:], **kwargs)

ax[1].loglog(f[1:], np.sqrt(xpsdsec[1:]), label = 'X')
ax[1].loglog(f[1:], np.sqrt(ypsdsec[1:]), label = 'Y')
ax[1].loglog(f[1:], np.sqrt(zpsdsec[1:]), label = 'Z')
ax[1].grid()
ax[1].legend()
# ax[1].set_xlim([1e-5, 1e-1])
# ax[1].set_ylim([5e-27, 5e-12])
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
# pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)


# %%


print("Saving individual noise contribution")

noises = ['laser', 'test-mass', 'oms', 'ranging', 'backlink', 'clock', 'modulation']

fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
=======

>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
sumpsd = np.ndarray(xpsd.shape)

for nn in noises:
    # Instantiate LISA instrument
    instr = Instrument(size=instrument_size,
                        seed = simseed,
                        dt=dt,
                        t0=instrument_t0, 
                        lock=locking, 
                        orbits=datadir+"new-orbits.h5", 
                        aafilter=('kaiser', 240, 1.1, 2.9),
                        clock_offsets={'1':clock_offsets[0],
                                       '2':clock_offsets[1],
                                       '3':clock_offsets[2]},
                        ranging_biases=ranging_biases,
                        moc_time_correlation_asds = moc_time_correlation_asds
                    )
            
    instr.disable_all_noises(excluding=nn) 
    instr.simulate()
    #

    data_noise = Data.from_instrument(instr)

    # Build other 2.0 Michelson variables
    X_data = X.build(**data_noise.args)
    # Y_data = Y.build(**data_noise.args)
    # Z_data = Z.build(**data_noise.args)
    # Apply TDI 2.0
    x_noise = X_data(data_noise.measurements) / central_freq
    # y_noise = Y_data(data_noise.measurements)  / central_freq
    # z_noise = Z_data(data_noise.measurements)  / central_freq
    
    #
    tdi_times = np.arange(n_data)*dt

# fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex = False)

<<<<<<< HEAD
<<<<<<< HEAD
    ax[0].set_title('Baseline noises - new orbits - individual noises')
=======
    ax[0].set_title('Keplerian orbits - baseline simulation individual noises')
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
    ax[0].set_title('Baseline noises - new orbits - individual noises')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
    ax[0].plot(tdi_times, x_noise[pytdi_trim:], label = 'X {nn}'.format(nn=nn), alpha = 0.5)
    # ax[0].plot(tdi_times, y_noise[pytdi_trim:], label = 'Y', alpha = 0.5)
    # ax[0].plot(tdi_times, z_noise[pytdi_trim:], label = 'Z', alpha = 0.5)
    ax[0].grid()
    ax[0].legend()
    ax[0].set_xlabel('Time [s]')
    ax[0].set_ylabel('TDI')
     
    
    kwargs = {"fs": fs,
              "window": 'blackman',
              "nperseg": instrument_size,
              "detrend": 'constant',
              "return_onesided": True,
              "scaling": 'density'}
    
    f, xpsd = signal.welch(x_noise[pytdi_trim:], **kwargs)
    # f, ypsd = signal.welch(y_noise[pytdi_trim:], **kwargs)
    # f, zpsd = signal.welch(z_noise[pytdi_trim:], **kwargs)
    
    ax[1].loglog(f[1:], np.sqrt(xpsd[1:]), label = 'X {nn}'.format(nn=nn))
    # ax[1].loglog(f[1:], np.sqrt(ypsd[1:]), label = 'Y')
    # ax[1].loglog(f[1:], np.sqrt(zpsd[1:]), label = 'Z')
    ax[1].grid()
    ax[1].legend()
    # ax[1].set_xlim([1e-5, 1e-1])
    # ax[1].set_ylim([5e-27, 5e-17])
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].set_ylabel(r'$S_\text{TDI}(f)$') 
<<<<<<< HEAD
<<<<<<< HEAD
    # pp.savefig(fig)
=======
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
    # pp.savefig(fig)
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
 
    sumpsd += xpsd
    
    
# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 6), sharex = False)
<<<<<<< HEAD
<<<<<<< HEAD
 
ax.set_title('Baseline noises - new orbits')
=======

>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
 
ax.set_title('Baseline noises - new orbits')
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax.loglog(f[1:], np.sqrt(xpsdsec[1:]), label = 'TDI X baseline secondary noises')
ax.loglog(f[1:], np.sqrt(sumpsd[1:]), label = 'sum of all TDI X noise individual contributions')
ax.grid()
ax.legend()
# ax[0].set_xlim([1e-5, 1e-1])
ax.set_ylim([5e-27, 5e-17])
ax.set_xlabel('Frequency [Hz]')
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
ax.set_ylabel(r'$S_\text{TDI}(f)$') 
# pp.savefig(fig)

# %%
<<<<<<< HEAD
# pp.close()



=======
ax.set_ylabel(r'$S_\text{TDI}(f)$') 
>>>>>>> 8de550e (Add simulation consolidation script with all plots for Issue #23)
=======
# pp.close()
>>>>>>> 2b5f603 (Update simulation consolidation script with all plots for Issue #23)
