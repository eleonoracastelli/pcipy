#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:56:59 2024

Simulation of SGWB signal using LISA GW Response.
Using the Sept 2024 version of the LISA Simulation Suite: 
    lisagwresponse, pyTDI, lisaorbits

Usage:
    Execute simulation with no TDI computation:
        
        python all_sky_signal_simulation.py path-to-workdir
        
    Execute simulation with TDI computation: use flag --tdi to specify TDI generation 
    
        python all_sky_signal_simulation.py path-to-workdir --tdi 2    

@author: Q Baghi 2021, modified by E Castelli 2024
"""
#
import argparse
import logging
import numpy as np
import healpy as hp
import h5py
import os
from datetime import datetime
from lisaorbits import KeplerianOrbits, EqualArmlengthOrbits
from lisagwresponse import StochasticBackground
from lisagwresponse.psd import white_generator
from pytdi.michelson import X1, Y1, Z1, X2, Y2, Z2
from pytdi import Data

if __name__ == "__main__":

    # To print the logs
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    # Configure program arguments.
    parser = argparse.ArgumentParser(
        description="A wrapper script to run ptemcee's using mpirun"
    )

    parser.add_argument(
        "output_path",
        type=str,
        default=None,
        help="Path of the output file",
    )

    parser.add_argument(
        "-dt",
        "--dt",
        type=float,
        default=1/4,
        help="Sampling time",
    )
    
    parser.add_argument(
        "-orb",
        "--orbits",
        default='keplerian', 
        choices=['keplerian','equalarm'],
        help="Choose orbit type",
    )
    
    parser.add_argument(
        "-tdi",
        "--tdi",
        default=None, 
        choices=[None,'1','2'],
        help="Pass TDI generation: choice between 1 and 2",
    )

    # Parse the input.
    args = parser.parse_args()
    # Sampling time
    dt = args.dt
    # Sampling frequency
    fs = 1 / dt
    # t0
    t0 = 2173211130.0 # s  datetime.datetime(2038, 11, 12, 16, 45, 30)

    # Data size: 24 hours or 2 days
    tobs = 3 * 24 * 3600
    n_data = int(tobs * fs)
    logging.info("Data size: " + str(n_data))
    logging.info("Data duration: " + str(tobs/3600) + " hours")
    # Central frequency
    central_freq = 281600000000000.0

    # set up proper time grid for simulation
    
    pytdi_trim = 1000
    pytdi_t0 = t0 - pytdi_trim * dt
    pytdi_size = n_data + pytdi_trim

    instrument_t0 = pytdi_t0
    instrument_size = pytdi_size

    orbits_dt = 100_000
    orbits_trim = 100
    orbits_t0 = t0 - pytdi_trim * dt - orbits_trim * orbits_dt
    orbits_size = np.ceil(3600 * 24 * 365 / orbits_dt) # a year
    
<<<<<<< HEAD
    if args.orbits == 'keplerian':
        OrbitsGenerator = KeplerianOrbits
    elif args.orbits == 'equalarm':
        OrbitsGenerator = EqualArmlengthOrbits
        
    # Generate new keplerian orbits
    orbits = args.output_path+"/"+args.orbits+"-orbits.h5"
    print('***************************************************************************')
    if not os.path.isfile(orbits):
        print('**** Orbits file not in output path folder. Generating {orb} orbit file.'.format(orb=args.orbits))
        orbitsobj = OrbitsGenerator()
        orbitsobj.write(orbits, dt=orbits_dt, size=orbits_size, t0=orbits_t0, mode="w")
    else:
        print('**** Selecting existing {orb} orbit file.'.format(orb=args.orbits))
    print('***************************************************************************') 
    
=======
    # Generate new keplerian orbits
    orbits = args.output_path+"/keplerian-orbits.h5"
    if not os.path.isfile(orbits):
        print('***************************************************************************')
        print('**** KeplerianOrbits file not in output path folder. Generating orbit file.')
        print('***************************************************************************')
        orbitsobj = KeplerianOrbits()
        orbitsobj.write(orbits, dt=orbits_dt, size=orbits_size, t0=orbits_t0, mode="w")
>>>>>>> 0193744 (Add orbit generation to all sky simulation file)
    # Instantiate GW signal class
    npix = hp.nside2npix(8)
    skymap = np.ones(npix) / np.sqrt(npix)
    generator = white_generator(1)
    
    src_class = StochasticBackground(skymap, 
                                      generator, 
                                      orbits=orbits, 
                                      dt=dt, 
                                      size=n_data, 
                                      t0=instrument_t0, 
                                      optim=True)

    # Choose files' prefixes
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%Y-%m-%d_") + args.orbits + '_'
    
    # Compute and save the GW response
    gw_file = args.output_path + '/' + dt_string + 'all_sky_gw_measurements_'+str(int(fs))+'Hz.h5'
    try:
        os.remove(gw_file)
    except FileNotFoundError:
        pass
    src_class.write(gw_file,
                    dt=dt, 
                    size=n_data, 
                    t0 = instrument_t0)
       
    #  Get data from GW simulation    
    if args.tdi:
        data_signal = Data.from_gws(gw_file, orbits)
        
        if args.tdi == '2':
            X, Y, Z = X2, Y2, Z2
        else:
            X, Y, Z = X1, Y1, Z1
        # Build other 2.0 Michelson variables
        X_data = X.build(**data_signal.args)
        Y_data = Y.build(**data_signal.args)
        Z_data = Z.build(**data_signal.args)
        # Apply TDI 2.0
        x_signal = X_data(data_signal.measurements)
        y_signal = Y_data(data_signal.measurements)
        z_signal = Z_data(data_signal.measurements)
    
        path = args.output_path + '/' + dt_string + 'all_sky_gw_tdi'+args.tdi+'_'+str(int(fs))+'Hz.h5'
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        hdf5 = h5py.File(path, 'w')
        hdf5.create_dataset('x', data=x_signal)
        hdf5.create_dataset('y', data=y_signal)
        hdf5.create_dataset('z', data=z_signal)
        # Closing file
        hdf5.close()