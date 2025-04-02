#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:03:07 2024

Simulation of laser + secondary noises using the LISA Instrument.

Using the Sept 2024 version of the LISA Simulation Suite: 
    lisainstrument, pyTDI, lisaorbits

This script is meant to substitute the older version of noise_simulation.py 
Simulation parameters match the simulation parameters defined in 
LISA-LCST-SGS-RP-006 "End-to-end Demonstration Data Analysis Pipeline"

Usage:
    Execute simulation with no TDI computation:
        
        python noise_simulation.py path-to-workdir
        
    Execute simulation with TDI computation: use flag --tdi to specify TDI generation 
    
        python noise_simulation.py path-to-workdir --tdi 2   
        
    Execute simulation with baseline InRep configuration
    
        python noise_simulation.py path-to-workdir --baseline
        
    Execute simulation with baseline InRep configuration an save all individual noise contributions
    
        python noise_simulation.py path-to-workdir --baseline --individual

@author: Q Baghi 2021, modified by E Castelli 2024
"""

import argparse
import logging
import os
from datetime import datetime
import h5py
import numpy as np
from lisainstrument import Instrument
from lisaorbits import KeplerianOrbits, EqualArmlengthOrbits
from pytdi.michelson import X1, Y1, Z1, X2, Y2, Z2
from pytdi import Data

# %%
if __name__ == "__main__":

    # To print the logs
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)

    # # For parallelization
    # from mpi4py import MPI
    # from schwimmbad import MPIPool

    # Configure program arguments.
    parser = argparse.ArgumentParser(
        description="A wrapper script to run ptemcee's using mpirun"
    )
    # add positional argument
    parser.add_argument(
        "output_path",
        type=str,
        default=None,
        help="Path of the simulation outputs folder",
    )

    # add optional arguments
    parser.add_argument(
        "-dt",
        "--dt",
        type=float,
        default=1/4,
        help="Sampling time",
    )

    parser.add_argument(
        "-f1",
        "--freq1",
        type=float,
        default=1.1,
        help="Kaiser filter first transition frequency (Hz)",
    )

    parser.add_argument(
        "-f2",
        "--freq2",
        type=float,
        default=2.9,
        help="Kaiser filter second transition frequency (Hz)",
    )

    parser.add_argument(
        "-s", 
        "--seed", 
        type=int,
        default = None,
        help = "Specify simulation seed")

    parser.add_argument(
        "-orb",
        "--orbits",
        default='keplerian',
        choices=['keplerian','equalarm'],
        help="Choose orbit type",
    )

    parser.add_argument(
        "-l",
        "--locking",
        default='N1-12',
        choices=['six','N1-12'],
        help="Choose locking configuration",
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

    parser.add_argument(
        "-b",
        "--baseline",
        action="store_true",
        help="Implement baseline simulation configuration from LISA-LCST-SGS-RP-006",
    )

    parser.add_argument(
        "-i",
        "--individual",
        action="store_true",
        help="Save all secondary noises as individual noise sources",
    )

    parser.add_argument(
        "-c",
        "--combined",
        action="store_true",
        help="Save acombinations of laser noise + some individual noise sources",
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
    
    if args.orbits == 'keplerian':
        OrbitsGenerator = KeplerianOrbits
    elif args.orbits == 'equalarm':
        OrbitsGenerator = EqualArmlengthOrbits

    # Generate new orbits
    orbits = args.output_path+"/"+args.orbits+"-orbits.h5"
        
    # Generate new keplerian orbits
    orbits = args.output_path+"/"+args.orbits+"-orbits.h5"
    if not os.path.isfile(orbits):
        print('***************************************************************************')
        print('**** Orbits file not in output path folder. Generating {orb} orbit file.'.format(orb=args.orbits))
        print('***************************************************************************')
        orbitsobj = OrbitsGenerator()
        orbitsobj.write(orbits, dt=orbits_dt, size=orbits_size, t0=orbits_t0, mode="w")
    else:
        print('**** Selecting existing {orb} orbit file.'.format(orb=args.orbits))
    print('***************************************************************************')

    # noise parameters to turn selected noises back on
    locking=args.locking
    print("*************************************************")
    print("Using {locking} locking configuration".format(locking=locking))
    print("*************************************************")
    # default parameters are commented here for reference
    # oms_asds=(6.35e-12, 1.25e-11, 1.42e-12, 3.38e-12, 3.32e-12, 7.90e-12)     
    # tm_asds=2.4E-15
    # laser_asds=30
    clock_offsets=(0,0,0)
    ranging_biases=0
    moc_time_correlation_asds=0.42

    if args.baseline:
        print("*************************************************")
        print("Using LISA-LCST-SGS-RP-006 baseline configuration")
        print("*************************************************")
        # locking='N1-12' # default configuration used in LISA-LCST-SGS-RP-006
        # default parameters are commented here for reference
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
    instr = Instrument(size=n_data,
                        dt=dt,
                        t0=instrument_t0,
                        seed = args.seed,
                        lock=locking,
                        orbits=orbits,
                        aafilter=('kaiser', 240, args.freq1, args.freq2),
                        clock_offsets={'1':clock_offsets[0],
                                       '2':clock_offsets[1],
                                       '3':clock_offsets[2]},
                        ranging_biases=ranging_biases,
                        moc_time_correlation_asds = moc_time_correlation_asds)

    # Disable all noises
    instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms'])
    instr.simulate()

    simseed = instr.seed

    print("*************************************************")
    print("Using seed {seed}".format(seed=simseed))
    print("*************************************************")

    # Choose files' prefixes
    # datetime object containing current date and time
    now = datetime.now()
    # # dd/mm/YY H:M:S
    lockstr = '_locking_'+locking+'_'

    dt_string = now.strftime("%Y-%m-%d_") + args.orbits  +  lockstr + 'laser_tm_oms_'

    writepath = args.output_path + '/' + dt_string + 'measurements_'+str(int(fs))+'Hz.h5'
    # Simulate and save data
    # check if file exists and delete it otherwise
    try:
        os.remove(writepath)
    except FileNotFoundError:
        pass
    instr.write(writepath)

    # Get the single-link outputs and delays

    if args.tdi:
        data_noise = Data.from_instrument(instr)
        print("*************************************************")
        print("Using TDI {tdi}".format(tdi=args.tdi))
        print("*************************************************")
        if args.tdi == '2':
            X, Y, Z = X2, Y2, Z2
        else:
            X, Y, Z = X1, Y1, Z1
        # Build other 2.0 Michelson variables
        X_data = X.build(**data_noise.args)
        Y_data = Y.build(**data_noise.args)
        Z_data = Z.build(**data_noise.args)
        # Apply TDI 2.0
        x_noise = X_data(data_noise.measurements)
        y_noise = Y_data(data_noise.measurements)
        z_noise = Z_data(data_noise.measurements)

        path = args.output_path + '/' + dt_string + 'noise_tdi'+args.tdi+'_'+str(int(fs))+'Hz.h5'
        try:
            os.remove(path)
        except FileNotFoundError:
            pass
        hdf5 = h5py.File(path, 'w')
        hdf5.create_dataset('x', data=x_noise)
        hdf5.create_dataset('y', data=y_noise)
        hdf5.create_dataset('z', data=z_noise)
        # Closing file
        hdf5.close()

    # Generate secondary noises HDF5 file
    # disable laser noise to simulate secondary noises
    instr.disable_all_noises(excluding=['test-mass', 'oms'])
    instr.simulate()
    # Store secondary noises
    writepath = args.output_path + '/' + dt_string + 'noise_sec_'+str(int(fs))+'Hz.h5'
    try:
        os.remove(writepath)
    except FileNotFoundError:
        pass
    instr.write(writepath)

    if args.baseline and (not args.individual and not args.combined):
        print("***** baseline")

        # Instantiate LISA instrument
        instr = Instrument(seed=simseed,
                           size=n_data,
                            dt=dt,
                            t0=instrument_t0,
                            lock=locking,
                            orbits=orbits,
                            aafilter=('kaiser', 240, args.freq1, args.freq2),
                            clock_offsets={'1':clock_offsets[0],
                                           '2':clock_offsets[1],
                                           '3':clock_offsets[2]},
                            ranging_biases=ranging_biases,
                            moc_time_correlation_asds = moc_time_correlation_asds)
        # Disable all noises
        instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms',
                                            'ranging', 'backlink', 'clock', 'modulation'])
        instr.simulate()

        dt_string = now.strftime("%Y-%m-%d_") + args.orbits + lockstr + 'baseline_'

        # Simulate and save data
        writepath = args.output_path + '/' + dt_string + 'measurements_' + str(int(fs)) + 'Hz.h5'
        try:
            os.remove(writepath)
        except FileNotFoundError:
            pass
        instr.write(writepath)

        # Get the single-link outputs and delays
        if args.tdi:
            data_noise = Data.from_instrument(instr)
            print("***** tdi {n}".format(n=args.tdi))

            if args.tdi == '2':
                X, Y, Z = X2, Y2, Z2
            else:
                X, Y, Z = X1, Y1, Z1
            # Build other 2.0 Michelson variables
            X_data = X.build(**data_noise.args)
            Y_data = Y.build(**data_noise.args)
            Z_data = Z.build(**data_noise.args)
            # Apply TDI 2.0
            x_noise = X_data(data_noise.measurements)
            y_noise = Y_data(data_noise.measurements)
            z_noise = Z_data(data_noise.measurements)

            path = args.output_path + '/' + dt_string + 'noise_tdi'+args.tdi+'_'+str(int(fs))+'Hz.h5'
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
            hdf5 = h5py.File(path, 'w')
            hdf5.create_dataset('x', data=x_noise)
            hdf5.create_dataset('y', data=y_noise)
            hdf5.create_dataset('z', data=z_noise)
            # Closing file
            hdf5.close()

        # Generate secondary noises HDF5 file
        # disable laser noise to simulate secondary noises
        instr.disable_all_noises(excluding=['test-mass', 'oms',
                                            'ranging', 'backlink', 'clock', 'modulation'])
        instr.simulate()
        # Store secondary noises
        writepath=args.output_path + '/' + dt_string + 'noise_sec_'+str(int(fs))+'Hz.h5'
        try:
            os.remove(writepath)
        except FileNotFoundError:
            pass
        instr.write(writepath)
    if args.baseline:
        dt_string = now.strftime("%Y-%m-%d_") + args.orbits + lockstr + 'baseline_'

    if args.combined:
        print("***** Saving combined noise contribution")

        noises = ['test-mass', 'oms']

        for n in noises:
            print("***** Simulate laser + {n}".format(n=n))

            # Instantiate LISA instrument
            instr = Instrument(seed=simseed,
                               size=n_data,
                                dt=dt,
                                t0=instrument_t0,
                                lock=locking,
                                orbits=orbits,
                                aafilter=('kaiser', 240, args.freq1, args.freq2))

            instr.disable_all_noises(excluding=["laser", n])
            instr.simulate()
            writepath=args.output_path + '/' + dt_string+'noise_combined_laser_'+n+'_'+str(int(fs))+'Hz.h5'
            try:
                os.remove(writepath)
            except FileNotFoundError:
                pass
            print("***** write laser + {n}".format(n=n))
            instr.write(writepath)

        if args.baseline:
            noises = ['ranging', 'backlink', 'clock', 'modulation']
            for n in noises:
                print("***** Simulate LTO + {n}".format(n=n))
                # Instantiate LISA instrument
                instr = Instrument(seed=simseed,
                                   size=n_data,
                                    dt=dt,
                                    t0=instrument_t0,
                                    lock=locking,
                                    orbits=orbits,
                                    aafilter=('kaiser', 240, args.freq1, args.freq2),
                                    clock_offsets={'1':clock_offsets[0],
                                                   '2':clock_offsets[1],
                                                   '3':clock_offsets[2]},
                                    ranging_biases= ranging_biases,
                                    moc_time_correlation_asds = moc_time_correlation_asds)

                instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms', n])
                instr.simulate()
                writepath=args.output_path + '/' + dt_string+'noise_combined_lto_'+n+'_'+str(int(fs))+'Hz.h5'
                try:
                    os.remove(writepath)
                except FileNotFoundError:
                    pass
                print("***** write lto + {n}".format(n=n))
                instr.write(writepath)

    if args.individual:
        print("***** Saving individual noise contribution")

        noises = ['laser', 'test-mass', 'oms']

        for n in noises:
            print("***** Simulate {n}".format(n=n))

            # Instantiate LISA instrument
            instr = Instrument(seed=simseed,
                               size=n_data,
                                dt=dt,
                                t0=instrument_t0,
                                lock=locking,
                                orbits=orbits,
                                aafilter=('kaiser', 240, args.freq1, args.freq2))

            instr.disable_all_noises(excluding=n)
            instr.simulate()
            writepath=args.output_path + '/' + dt_string+'noise_individual_'+n+'_'+str(int(fs))+'Hz.h5'
            try:
                os.remove(writepath)
            except FileNotFoundError:
                pass
            print("***** write {n}".format(n=n))
            instr.write(writepath)
        if args.baseline:
            noises = ['ranging', 'backlink', 'clock', 'modulation']
            for n in noises:
                print("***** Simulate {n}".format(n=n))

                # Instantiate LISA instrument
                instr = Instrument(seed=simseed,
                                   size=n_data,
                                    dt=dt,
                                    t0=instrument_t0,
                                    lock=locking,
                                    orbits=orbits,
                                    aafilter=('kaiser', 240, args.freq1, args.freq2),
                                    clock_offsets={'1':clock_offsets[0],
                                                   '2':clock_offsets[1],
                                                   '3':clock_offsets[2]},
                                    ranging_biases= ranging_biases,
                                    moc_time_correlation_asds = moc_time_correlation_asds)

                instr.disable_all_noises(excluding=n)
                instr.simulate()
                writepath=args.output_path + '/' + dt_string + 'noise_individual_'+n+'_'+str(int(fs))+'Hz.h5'
                try:
                    os.remove(writepath)
                except FileNotFoundError:
                    pass
                print("***** write {n}".format(n=n))
                instr.write(writepath)
