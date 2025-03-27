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
import h5py
from datetime import datetime
from lisainstrument import Instrument
from pytdi.michelson import X1, Y1, Z1, X2, Y2, Z2
from pytdi import Data

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

    # Parse the input.
    args = parser.parse_args()
    # Sampling time
    dt = args.dt
    # Sampling frequency
    fs = 1 / dt
    # Data size: 24 hours or 2 days
    tobs = 3 * 24 * 3600
    n_data = int(tobs * fs)
    logging.info("Data size: " + str(n_data))
    logging.info("Data duration: " + str(tobs/3600) + " hours")
    # Central frequency
    central_freq = 281600000000000.0

    # Choose orbit file
    # TO DO double check that the keplerian orbits match the orbits simulation 
    # parameters in LISA-LCST-SGS-RP-006 
    datadir = '/Users/ecastel2/Documents/research/GSFC/simulation-tests/orbits/'
    orbits = datadir+"keplerian-orbits.h5"
    # orbits = "/work/SC/lisa/baghiq/orbits/keplerian-orbits.h5"
    with h5py.File(orbits) as f:
        orbit_t0 = f.attrs['t0']
    
    
    # noise parameters to turn selected noises back on
    locking='six'
    # oms_asds=(6.35e-12, 1.25e-11, 1.42e-12, 3.38e-12, 3.32e-12, 7.90e-12)        
    # tm_asds=2.4E-15
    # laser_asds=30
    
    if args.baseline:
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
    instr = Instrument(size=n_data,
                        dt=dt,
                        t0=1000 + orbit_t0, 
                        lock=locking, 
                        orbits=orbits, 
                        aafilter=('kaiser', 240, args.freq1, args.freq2),
                        clock_offsets={'1':clock_offsets[0],'2':clock_offsets[1],'3':clock_offsets[2]},
                        ranging_biases=ranging_biases,
                        moc_time_correlation_asds = moc_time_correlation_asds)
            
    # Disable all noises
    instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms'])
    instr.simulate()
    simseed = instr.seed

    # Choose files' prefixes
    # datetime object containing current date and time
    now = datetime.now()
    # # dd/mm/YY H:M:S
    lockstr = 'locking_n1_12_'
    dt_string = now.strftime("%Y-%m-%d_%Hh%M_") + lockstr + 'laser_tm_oms_'
    
    # Simulate and save data
    instr.write(args.output_path + '/' + dt_string + 'measurements_'+str(int(fs))+'Hz.h5')
    
    # Get the single-link outputs and delays
    
    if args.tdi:
        data_noise = Data.from_instrument(instr)
        
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
        hdf5 = h5py.File(path, 'a')
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
    instr.write(args.output_path + '/' + dt_string + 'noise_sec_'+str(int(fs))+'Hz.h5')
    
    if args.baseline:
        # Instantiate LISA instrument
        instr = Instrument(seed=simseed,
                           size=n_data,
                            dt=dt,
                            t0=1000 + orbit_t0,
                            lock=locking, 
                            orbits=orbits, 
                            aafilter=('kaiser', 240, args.freq1, args.freq2),
                            clock_offsets={'1':clock_offsets[0],'2':clock_offsets[1],'3':clock_offsets[2]},
                            ranging_biases=ranging_biases,
                            moc_time_correlation_asds = moc_time_correlation_asds
                            )
        # Disable all noises
        instr.disable_all_noises(excluding=['laser', 'test-mass', 'oms', 'ranging', 'backlink', 'clock', 'modulation'])
        instr.simulate()
        
        dt_string = now.strftime("%Y-%m-%d_%Hh%M_") + lockstr + 'baseline_'
        
        # Simulate and save data
        instr.write(args.output_path + '/' + dt_string + 'measurements_'+str(int(fs))+'Hz.h5')
        
        # Get the single-link outputs and delays
        if args.tdi:
            data_noise = Data.from_instrument(instr)
            
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
            hdf5 = h5py.File(path, 'a')
            hdf5.create_dataset('x', data=x_noise)
            hdf5.create_dataset('y', data=y_noise)
            hdf5.create_dataset('z', data=z_noise)
            # Closing file
            hdf5.close()
            
        # Generate secondary noises HDF5 file
        # disable laser noise to simulate secondary noises
        instr.disable_all_noises(excluding=['test-mass', 'oms', 'ranging', 'backlink', 'clock', 'modulation'])
        instr.simulate()
        # Store secondary noises
        instr.write(args.output_path + '/' + dt_string + 'noise_sec_'+str(int(fs))+'Hz.h5')
    
    if args.individual:
        print("Saving individual noise contribution")
        
        noises = ['laser', 'test-mass', 'oms']
        
        for n in noises:
            # Instantiate LISA instrument
            instr = Instrument(seed=simseed,
                               size=n_data,
                                dt=dt,
                                t0=1000 + orbit_t0,
                                lock=locking, 
                                orbits=orbits, 
                                aafilter=('kaiser', 240, args.freq1, args.freq2))
                    
            instr.disable_all_noises(excluding=n) 
            instr.simulate()
            instr.write(args.output_path + '/' + dt_string + 'noise_'+n+'_'+str(int(fs))+'Hz.h5')
       
        if args.baseline:
            noises = ['ranging', 'backlink', 'clock', 'modulation']
            for n in noises:
                # Instantiate LISA instrument
                instr = Instrument(seed=simseed,
                                   size=n_data,
                                    dt=dt,
                                    t0=1000 + orbit_t0,
                                    lock=locking, 
                                    orbits=orbits, 
                                    aafilter=('kaiser', 240, args.freq1, args.freq2),
                                    clock_offsets={'1':clock_offsets[0],'2':clock_offsets[1],'3':clock_offsets[2]},
                                    ranging_biases= ranging_biases,
                                    moc_time_correlation_asds = moc_time_correlation_asds)
                        
                instr.disable_all_noises(excluding=n) 
                instr.simulate()
                instr.write(args.output_path + '/' + dt_string + 'noise_'+n+'_'+str(int(fs))+'Hz.h5')
            
            

                