#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:56:59 2024

Simulation of SGWB signal using LISA GW Response.
Using the Sept 2024 version of the LISA Simulation Suite: 
    lisagwresponse, pyTDI, lisaorbits

Usage:
    Execute simulation with no TDI computation:
        
        python signal_simulation.py path-to-workdir
        
    Execute simulation with TDI computation: use flag --tdi to specify TDI generation 
    
        python signal_simulation.py path-to-workdir --tdi 2    

@author: Q Baghi 2021, modified by E Castelli 2024
"""
#
import argparse
import logging
import numpy as np
import h5py
from datetime import datetime
from lisagwresponse import StochasticPointSource
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
    # Data size: 24 hours or 2 days
    tobs = 3 * 24 * 3600
    n_data = int(tobs * fs)
    logging.info("Data size: " + str(n_data))
    logging.info("Data duration: " + str(tobs/3600) + " hours")

    # Choose orbit file
    # orbits = "/data/jgbaker/software/pylisa/data/keplerian-orbits.h5"
    datadir = "/Users/ecastel2/Documents/research/GSFC/simulation-tests/orbits/"
    orbits = datadir+"keplerian-orbits.h5" # datadir + 'new-orbits.h5' #
    
    with h5py.File(orbits) as f:
        orbit_t0 = f.attrs['t0']
    
    
    # Instantiate GW signal class
    src_class = StochasticPointSource(white_generator(1),  #white_noise_generator_at_1,  
                                      orbits=orbits,
                                      gw_beta=np.pi/2, 
                                      gw_lambda=np.pi/2,
                                      dt = dt,
                                      t0 = 1000 + orbit_t0, 
                                      size=n_data)

    # Choose files' prefixes
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%Y-%m-%d_%Hh%M_")
    
    # Compute and save the GW response
    gw_file = args.output_path + '/' + dt_string + 'point_gw_measurements_'+str(int(fs))+'Hz.h5'
    src_class.write(gw_file,   
                    dt=dt, 
                    size=n_data, 
                    t0 = 1000 + orbit_t0)
       
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
    
        path = args.output_path + '/' + dt_string + 'point_gw_tdi'+args.tdi+'_'+str(int(fs))+'Hz.h5'
        hdf5 = h5py.File(path, 'a')
        hdf5.create_dataset('x', data=x_signal)
        hdf5.create_dataset('y', data=y_signal)
        hdf5.create_dataset('z', data=z_signal)
        # Closing file
        hdf5.close()