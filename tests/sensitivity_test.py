#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 14 16:56:48 2025

1. Instantiate the SGWB class
To compute LISA’s response to a SGWB, we need two things:
    - A sky map defining the background anisotropy (or isotropy) across the sky
    - An orbit file that defines LISA spacecraft trajectories as a function of time.

2. Compute the response for one link

We can compute the single-link responses at frequency
and time . One of these arguments can be a vector. For example, let’s construct a frequency vector and pick a time:

3. Compute all links response matrix

If one wants to compute the full frequency-domain response of all the single-links including cross-correlations, the first argument can be replaced by the list of links:

@author: ecastel2
"""

import numpy as np
import os, sys
import healpy as hp
import h5py
# from backgrounds import priors, sampling, psd
sys.path.append('/Users/ecastel2/Documents/virtual-envs/software-install/backgrounds')
from backgrounds import utils, loadings, plotting, noise, tdi, analysis, signal
from backgrounds import StochasticBackgroundResponse
from lisaorbits import KeplerianOrbits, EqualArmlengthOrbits
import matplotlib.pyplot as plt
from datetime import datetime

import segwo
from lisaconstants import c
from pytdi.michelson import X2_ETA, Y2_ETA, Z2_ETA
from pytdi.core import LISATDICombination


# %% 1. Use Backgrounds
ORB = 'keplerian'
# %% generate isotropic sky map
nside = 12
npix = hp.nside2npix(nside)
# skymap = np.ones(npix) * np.sqrt(4 * np.pi / (2*npix))
skymap = np.ones(npix) / np.sqrt(2*npix)
# %% generate orbit file
workdir = "/Users/ecastel2/Documents/research/GSFC/pci-inrep/simulations/"

orbits = workdir+"sensitivity-test-orbits.h5"
t0 = 2173211130.0 # s  datetime.datetime(2038, 11, 12, 16, 45, 30)
orbits_dt = 100_000
orbits_trim = 100
orbits_t0 = t0 - orbits_trim * orbits_dt
orbits_size = np.ceil(3600 * 24 * 365 / orbits_dt) # a year

if ORB == 'keplerian':
    OrbitsGenerator = KeplerianOrbits
elif ORB == 'equalarm':
    OrbitsGenerator = EqualArmlengthOrbits

print('***************************************************************************')
orbitsobj = OrbitsGenerator()
orbitsobj.write(orbits, dt=orbits_dt, size=orbits_size, t0=orbits_t0, mode="w")
print('***************************************************************************')

# %% Instantiate SGWB class


sgwbcls = StochasticBackgroundResponse(skymap=skymap, orbits=orbits, orbit_interp_order=1)
# %% build frequency vector
fnyq = 0
fmin = -5
nbins = 5000

freqs = np.logspace(fmin, fnyq, nbins)
t0 = np.ones(1) * (sgwbcls.t0 + 10)  #need an object with a specific length
    
# %% To compute the response of the link 12 we can use the following method. It takes a few seconds as integration over the sky is numerial: it sums over all pixels of the sky map.
# gplus, gcross = sgwbcls.compute_correlations([12], freqs, t0)
# # get the full correlation matrix
# response = gplus + gcross

# %% plot

# fig, ax = plt.subplots(1,1)
# ax.loglog(freqs, np.abs(response.flatten()))
# ax.set_xlabel('Frequencies [Hz]')
# ax.set_ylabel(r'$R_{12}(f)$')
# ax.set_title("'12' link SGWB response - sky averaged")
# ax.grid()

# %% Now let's do it for all links
gplus, gcross = sgwbcls.compute_correlations(sgwbcls.LINKS, freqs, t0)

bgresponse = gplus + gcross

# take a look at the shape of the response
# this is an array of 6 x 6 matrices (nlinks x nlinks)
print(bgresponse.shape)
# %% plot
# response of each link
fig, ax = plt.subplots(1,1)
for i, l in enumerate(sgwbcls.LINKS):
    ax.loglog(freqs, np.abs(bgresponse[:, i, i]), label = str(l))
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'$R_{12}(f)$')
ax.set_title('Signal PSD - sky averaged')
ax.grid()
ax.legend()

# correlations among the responses of 12 with all the other links
fig, ax = plt.subplots(1,1)
for i, l in enumerate(sgwbcls.LINKS):
    ax.loglog(freqs, np.abs(bgresponse[:, 0, i]), label = str([sgwbcls.LINKS[0],l]))
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'$R_{12}(f)$')
ax.set_title('Signal CSD - sky averaged')
ax.grid()
ax.legend()

# %% compute TDI response matrix
# Time delay interferometry performs a linear combination of delayed single-link measurements to cancel laser frequency noise. It outputs three channels that contain all the GW information.

# Algebraically, this transformation can be written in the (time-)frequency domain using a matrix operation:
# If you have already computed the single-link responses like above, you can simply compute the TDI transfer function (for X, Y, Z) and then apply it through the equation above:

# tdi transfer function at frequencies freqs and time t0
# this is an array of 3 x 6 matrices (ncombs x nlinks)
tdimat = sgwbcls.compute_tdi_design_matrix(freqs, t0, gen = '2.0')
# compute the tdi response in the frequency domain
# this is an array of 3 x 3 matrices
#  Rtdi = Mtdi * Rlinks * Mtdi'
# (ncombs x ncombs) = (ncombs x nlinks) * (nlinks x nlinks) * (nlinks x ncombs)
tdiresponse = utils.transform_covariance(tdimat, bgresponse)

# look at the shapes
tdimat.shape, tdiresponse.shape

# %% plot tdi response
fig, ax = plt.subplots(1,1)
for i, t in enumerate(['X', 'Y', 'Z']):
    ax.loglog(freqs, np.abs(tdiresponse[:, 0, i]), label = "R_X{t}".format(t = t))
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'$R_{tdi, X}(f)$')
ax.set_ylim([1e-17, 500])
ax.set_xlim([2e-5, 2e0])
ax.grid()
ax.set_title('TDI X, Y, Z response - sky averaged')
ax.legend()

# print(m)
# fr=np.array([0.00010, 0.00100, 0.01000, 0.10000, 1.00000])
# f = np.logspace(faxis[m][0],faxis[m][1], 1000 )
    
# %% Compute SGWB background PSD in TDI

# Now that we have computed the response matrix of any isotropic SGWB, we can derive the PSD of a SGWB given its strain PSD
# . For example, let’s assume that the energy density of the background is a power-law with amplitude and index
# Then the strain PSD at present time is
# You can compute it as follows:

sh = signal.sgwb_psd(freqs, spec_index=0.5, freq0=1e-3, omega_gw=1e-14)

# %% plot strain

fig, ax = plt.subplots(1,1)
ax.loglog(freqs, sh, label = r"$S_h(f)$")
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'PSD Hz$^{-1}$')
ax.set_title('Strain PSD - SGWB')
ax.grid()
ax.legend()

# %% Then we can comp TDI covariance from the SGWB by simply multiplying S_h with the TDI response matrix R_tdi. We obtain a 3x3 matrix for each frequency bin, whose entries are the PSDs and the CSDs of the TDI channels .

tdicovariance = tdiresponse * sh[:, np.newaxis, np.newaxis]

# %%
fig, ax = plt.subplots(1,1)
ax.loglog(freqs, np.abs(tdicovariance[:,0,0]), label=r"$S_{XX}$")
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'TDI PSD Hz$^{-1}$')
ax.set_title('TDI covariance from the SGWB')
ax.grid()
ax.legend()

# %%
dt = 1/4
fs = 4
# Observation time
a_day = 86400
tobs = 7 * a_day
# # Choose frequency segments edges
# f_segments = psd.frequency_grid(100/tobs, 15000/tobs, 1e-5, 2.9e-2)
# # Choose window function
# wd_func = loadings.get_window_function("bh92")
# # Compute full periodogram
# per_xyz_data = psd.periodogram_matrix(xyz_noise.T, fs, wd_func=wd_func)
# # Smooth the periodogram
# fper, per_xyz_smoothed, segment_sizes = psd.smooth(per_xyz_data, fs, f_segments,
#                                                     weights_func=np.ones)
# # We offset the instrument t0 to have orbits information at emission time
# t0 = orbits_t0 + 100.0
# # Restrict the frequency band if not already done
# inds = loadings.restrict_frequencies(orbits,
#                                      fper, t0,
#                                      output_freqs=False,
#                                      fmin=1e-4,
#                                      fmax=2.8e-2,
#                                      df_margin=2e-4)
# # Restrict the frequency band
# finds = fper[inds]

noise_classes_analytic = []
noise_classes_analytic.append(noise.AnalyticOMSNoiseModel(
    freqs, t0, orbits, orbit_interp_order=1, gen="2.0", fs=None, duration=None, oms_isi_carrier_asds=15e-12))
noise_classes_analytic.append(noise.AnalyticTMNoiseModel(
    freqs, t0, orbits, orbit_interp_order=1, gen="2.0", fs=None, duration=None, tm_isi_carrier_asds=3e-15))




# noise_classes_analytic = []
# noise_classes_analytic.append(noise.AnalyticNoiseModel(freqs, t0, orbits, link_psd_func)


# Compute the TDI spectrum from OMS and TM noise
cov_tdi_n_comps = [nc.compute_covariances(0.0) for nc in noise_classes_analytic]
# Compute the full TDI noise spectrum
cov_tdi_n = sum(cov_tdi_n_comps)

# %%
fig, ax = plt.subplots(1,1)
ax.loglog(freqs, (np.abs(cov_tdi_n[:, 0, 0])),
            label="Theoretical PSD X",
            color="tab:orange",
            linestyle='dashed')
# for i, t in enumerate(['X', 'Y', 'Z']):
#     ax.loglog(freqs, np.abs(tdiresponse[:, 0, i])**2, label = "R_X{t}".format(t = t))
ax.set_ylim([1e-45,1e-35])
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'TDI PSD Hz$^{-1}$')
ax.set_title('TDI PSD')
ax.grid()
ax.legend()

# %%

S_hx = np.sqrt((np.abs(cov_tdi_n[:,0,0])/np.abs(tdiresponse[:,0,0])))

fig, ax = plt.subplots(1,1)
ax.loglog(freqs, S_hx, label=r"$S_{h,X}$")
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'TDI X Sensitivity')
ax.grid()
ax.legend()

# %% Use segwo

ltts = np.zeros((1, 6))
positions = np.zeros((1, 3, 3))

with h5py.File(orbits) as f:
    ltts[0,:] = f['tcb']['ltt'][0]
    positions[0,:,:] = f['tcb']['x'][0]

# Compute the signal covariance matrix for eta variables
signal_cov_eta = segwo.response.compute_isotropic_signal_cov(freqs, ltts, positions)

# We get an array of shape (1, 1000, 6, 6) for the single time point (t=0.0),
# our grid of 1000 frequencies, and the 6 eta variables
print(signal_cov_eta.shape)

# %%

fig, ax = plt.subplots(1,1, figsize=(8,6))
for i, l in enumerate(sgwbcls.LINKS):
    ax.loglog(freqs, np.abs(signal_cov_eta[0, :, i, i]), label = str(l)+' segwo')
    ax.loglog(freqs, np.abs(bgresponse[:, i, i]), ls='--', label=str(l)+'backgrounds')
    print(np.abs(bgresponse[:, i, i])/np.abs(signal_cov_eta[0, :, i, i]))
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'Signal PSD for each link')
ax.grid()
ax.legend()


fig, ax = plt.subplots(1,1, figsize=(8,6))
for i, l in enumerate(sgwbcls.LINKS):
    ax.loglog(freqs, np.abs(signal_cov_eta[0, :, 0, i]), label = '12, '+str(l)+' segwo')
    ax.loglog(freqs, np.abs(bgresponse[:, 0, i]), ls='--', label='12, '+str(l)+' backgrounds')
ax.set_xlabel('Frequencies [Hz]')
ax.set_ylabel(r'Signal CSD for 12 with each link')
ax.grid()
ax.legend()

# %% TDI
# We form the ordered list of input variables, i.e., the eta variables
eta_list = [f"eta_{mosa}" for mosa in sgwbcls.LINKS]

# Then we construct our mixing matrix for X2, Y2, and Z2
eta2xyz = segwo.cov.construct_mixing_from_pytdi(
    freqs, eta_list, [X2_ETA, Y2_ETA, Z2_ETA], ltts
)

# It's a 3x6 matrix (transform 6 eta variables into 3 XYZ variables), also given
# for each time point (here, just t=0.0) and frequencies
eta2xyz.shape

# %%
TDI_LABELS = ["X", "Y", "Z"]
for i, t in enumerate(TDI_LABELS):
    fig, axs = plt.subplots(6,1, figsize=(8,6))
    j=0
    ax = axs[j]
    ax.set_title(r'{t} Mixing matrix to each link'.format(t=t))
    for j, l in enumerate(sgwbcls.LINKS):
        ax = axs[j]
        ax.loglog(freqs, np.abs(eta2xyz[0, :, i, j]), label = str(t)+'-'+str(l)+' segwo')
        ax.loglog(freqs, np.abs(tdimat[:, i, j]), ls='--', label=str(t)+'-'+str(l)+' backgrounds')
        # print(np.abs(bgresponse[:, i, i])/np.abs(signal_cov_eta[0, :, i, i]))
        ax.set_xlabel('Frequencies [Hz]')
        ax.grid()
        ax.legend(loc='upper left')

# %%
# We use the mixing matrix to project the eta covariance into the TDI variables
signal_cov_xyz = segwo.cov.project_covariance(signal_cov_eta, eta2xyz)

# It's a 3x3 covariance matrix for the TDI variables, given for each time point
# (here, just t=0.0) and frequencies
signal_cov_xyz.shape

# %%



j=0
fig, ax = plt.subplots(3,1, figsize=(10,8), sharex=True)
# Plotting the signal PSD for each TDI channel
for i, tdi in enumerate(TDI_LABELS[:1]):
    ax[j].loglog(freqs, np.abs(signal_cov_xyz[0, :, i, i]), label=f"PSD TDI {tdi} segwo")
   

# Plotting the signal CSD between TDI channels
for i, tdi_i in enumerate(TDI_LABELS[:1]):
    for j, tdi_j in enumerate(TDI_LABELS):
        if i < j:
            ax[j].loglog(
                freqs,
                np.abs(signal_cov_xyz[0, :, i, j]),
                "-",
                label=f"segwo",
            )
        ax[j].loglog(freqs, np.abs(tdiresponse[:, 0, j]), ls='--', label = "backgrounds".format(t = tdi_j))
        ax[j].set_ylabel("R_X{t}".format(t = tdi_j))
        ax[j].legend()
        ax[j].grid()

ax[j].set_xlabel("Frequency [Hz]")
ax[0].set_title("TDI response for each channel")
fig.tight_layout()

# %%
# Directly compute the signal covariance in the TDI variables
shortcut_signal_cov_xyz = segwo.response.compute_isotropic_signal_cov(
    freqs, ltts, positions, mixings=eta2xyz, nside=nside
)

# We find they give identical results (no exception is raised)
np.testing.assert_allclose(signal_cov_xyz, shortcut_signal_cov_xyz)

# %% Define the noise PSDs, using LISA SciRD noise levels and shapes

# The OMS noise is defined in terms of displacement (meters), which we convert
# to fractional frequency shifts
displ_2_ffd = 2 * np.pi * freqs / c
oms = (15e-12) ** 2 * displ_2_ffd**2 * (1 + ((2e-3) / freqs) ** 4)

# The TM noise is defined in terms of acceleration (m/s^2), which we convert to
# fractional frequency shifts
acc_2_ffd = 1 / (2 * np.pi * freqs * c)
tm = (3e-15) ** 2 * acc_2_ffd**2 * (1 + (0.4e-3 / freqs) ** 2) * (1 + (freqs / 8e-3) ** 4)

# %% 

# Plot the noise PSDs
# As expected, the OMS noise is much larger dominates at high frequencies, while
# the TM noise dominates at low frequencies (below ~ 1 mHz)
plt.figure(figsize=(10, 6))
plt.loglog(freqs, oms, label="OMS")
plt.loglog(freqs, tm, label="TM")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Noise PSD [/Hz]")
plt.title("Noise PSDs")
plt.legend()
plt.show()

# %% Construct the overall noise covariance
noise_cov = segwo.cov.construct_covariance_from_psds([oms] * 6 + [tm] * 6)

# The matrix is diagonal (all noises are uncorrelated) and with shape (1000, 12,
# 12), corresponding to the 1000 frequencies and 12 variables (6 OMS + 6 TM)
noise_cov.shape

# %% Define a dictionary of TDI combinations mapping our noises into each of the 6
# single links (aka, the 6 eta variables)
ETA_COMBS = {}

# We start by defining it for eta 12:
#
# A PyTDI combination is defined as a dictionary, where the keys are the names
# of the input variables (N_tm_12, N_oms_12, etc.), and the values are lists of
# tuples. Each tuple contains a scaling coefficient and a list of delays to be
# applied on the corresponding variable.
ETA_COMBS[12] = LISATDICombination(
    {
        "N_tm_21": [(1, ["D_12"])],
        "N_tm_12": [(1, [])],
        "N_oms_12": [(1, [])],
    }
)

# Use PyTDI symmetry operations to define the other links

# We can do cyclic permuations, rotating indices 1->2->3->1
ETA_COMBS[23] = ETA_COMBS[12].rotated()
ETA_COMBS[31] = ETA_COMBS[23].rotated()

# We can do reflections along the axis going through the respective spacecraft,
# ie., exchanging indices 2<->3, 3<->1, and 1<->2
ETA_COMBS[13] = ETA_COMBS[12].reflected(1)
ETA_COMBS[21] = ETA_COMBS[23].reflected(2)
ETA_COMBS[32] = ETA_COMBS[31].reflected(3)

# %% Define list of noise labels (our input variables)
noise_list = [f"N_oms_{mosa}" for mosa in sgwbcls.LINKS] + [
    f"N_tm_{mosa}" for mosa in sgwbcls.LINKS
]

# Construct the mixing matrix for the noise covariance
noise2eta = segwo.cov.construct_mixing_from_pytdi(
    freqs,
    measurements=noise_list,
    tdi_combinations=[ETA_COMBS[mosa] for mosa in sgwbcls.LINKS],
    ltts=ltts,
)

# %% It's a 6x12 matrix (transforms 12 noise variables into 6 single link
# measurements), with two additional first axes for the time point (t=0.0) and
# frequencies
noise2eta.shape

# We get a 6x6 covariance matrix for the single link measurements, given for
# each time point (here, just t=0.0) and each frequency on our grid
noise_cov_eta = segwo.cov.project_covariance(noise_cov, noise2eta)

noise_cov_eta.shape

# We can then further project the result onto the TDI XYZ observables, producing
# a 3x3 covariance matrix, given for each time point and frequencies
noise_cov_xyz = segwo.cov.project_covariance(noise_cov_eta, eta2xyz)

noise_cov_xyz.shape

# %%
plt.figure(figsize=(10, 6))

# Plot the noise PSD for each TDI channel
for i, tdi in enumerate(TDI_LABELS[:1]):
    plt.loglog(freqs, np.abs(noise_cov_xyz[0, :, i, i]), label=f"segwo {tdi}")


plt.loglog(freqs, (np.abs(cov_tdi_n[:, 0, 0])),
            label="backgrounds X",
            color="tab:orange",
            linestyle='dashed')
# Plot the noise CSD between TDI channels
# for i, tdi_i in enumerate(TDI_LABELS):
#     for j, tdi_j in enumerate(TDI_LABELS):
#         if i < j:
#             plt.loglog(
#                 freqs,
#                 np.abs(noise_cov_xyz[0, :, i, j]),
#                 "--",
#                 label=f"CSD TDI {tdi_i} with {tdi_j}",
#             )

plt.xlabel("Frequency [Hz]")
plt.ylabel("Noise PSD/CSD [/Hz]")
plt.title("Noise PSD/CSD for TDI X")
plt.legend()
plt.grid()
plt.show()

# %%
# First, compose the mixing matrices
noise2xyz_composed = segwo.cov.compose_mixings([noise2eta, eta2xyz])

# And then project only once (from noise to XYZ)
noise_cov_xyz_composed = segwo.cov.project_covariance(noise_cov, noise2xyz_composed)

# Check that the two methods give identical results (no exception is raised)
np.testing.assert_allclose(noise_cov_xyz, noise_cov_xyz_composed)

noise_cov_xyz_composed2 = segwo.cov.project_covariance(noise_cov, [noise2eta, eta2xyz])

# Check that the two methods give identical results (no exception is raised)
np.testing.assert_allclose(noise_cov_xyz_composed, noise_cov_xyz_composed2)

# %%

# Define a dictionary containing the TDI combinations mapping the noises into
# the eta variables, making sure that the keys map the measurement labels used
# in `X2_ETA`, `Y2_ETA`, and `Z2_ETA` (to allow composing them)
ETA_COMBS_SET = {f"eta_{mosa}": ETA_COMBS[mosa] for mosa in sgwbcls.LINKS}

# Define TDI combinations directly mapping noises into XYZ using composition
X_composed = (X2_ETA @ ETA_COMBS_SET).simplified()
Y_composed = (Y2_ETA @ ETA_COMBS_SET).simplified()
Z_composed = (Z2_ETA @ ETA_COMBS_SET).simplified()

# This combination now directly applies delays to the underlying noises,
# as can be seen in the input measurements used by `X_composed` (equivalently
# for `Y_composed` and `Z_composed`)
X_composed.measurements


# %% Construct the mixing matrix for the noise covariance
noise2xyz_pytdi = segwo.cov.construct_mixing_from_pytdi(
    freqs,
    measurements=noise_list,
    tdi_combinations=[X_composed, Y_composed, Z_composed],
    ltts=ltts,
)

# It's a 3x12 matrix (transform 12 noise variables into 3 TDI variables), with
# two additional axis for the time point (t=0.0) and frequencies
noise2xyz_pytdi.shape

# %% Use composed matrix to project the noise covariance into the TDI variables
noise_cov_tdi_pytdi = segwo.cov.project_covariance(noise_cov, noise2xyz_pytdi)

# We can test that all three approaches yield the same result
np.testing.assert_allclose(noise_cov_xyz, noise_cov_xyz_composed)
np.testing.assert_allclose(noise_cov_xyz, noise_cov_tdi_pytdi)

# %% Compute optimal sensitivity curve for TDI XYZ
sensivity_xyz = segwo.sensitivity.compute_sensitivity_from_covariances(
    noise_cov_xyz, signal_cov_xyz
)

# It's of shape (1, 1000), with a single time point (t=0.0) and our grid of 1000
# frequencies
sensivity_xyz.shape

sensitivity_x = segwo.sensitivity.compute_sensitivity_from_covariances(noise_cov_xyz[:, :, :1, 0:1], signal_cov_xyz[:, :, 0:1, 0:1])
    
# %% Plot the optimal sensitivity curve for XYZ
plt.figure(figsize=(10, 6))
# plt.loglog(freqs, np.sqrt(sensivity_xyz[0]), label="segwo sensitivity")
plt.loglog(freqs, np.sqrt(sensitivity_x[0]), label="segwo sensitivity X")
plt.loglog(freqs, S_hx, label="backgrounds sensitivity X")
plt.xlabel("Frequency [Hz]")
plt.ylabel("TDI X sensitivity [1/sqrt(Hz)]")
plt.title("KeplerianOrbits orbit file")
plt.legend()
plt.grid()
plt.show()

