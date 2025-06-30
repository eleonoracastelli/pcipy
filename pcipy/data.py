#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 16 20:00:05 UTC 2025

@author: J Baker, E Castelli,
"""

import numpy as np

class TimeData:
    '''
    This is a simple base class for temporally sampled data sets.  The base
    base class supports regularly sampled data.

    The idea is to replace the PCI data arguments with these objects.

    Expect to extend this to a version with is segmented to allow for, eg
    gapped data.  

    Probably include support for named channels and for testing against an 
    expected set of names, and possibly reordering the names.  This could 
    include support of various versions of the names to support polymorphism.
    
    Still we will try to keep it simple here!

    It is likely that other special versions of this are useful
    '''

    def __init__(self, data, dt=None, t0=None, names=None):
        '''
        The interface is simple, and the argmuments are optional, but will
        impact the functionality if not provided.

        Parameters
        ----------
        data : ndarray
            The basic data as an (nchannels, nsamples) shaped numpy array.
        
        dt : float, optional
            If provided, this should be the sampling cadence.  It is needed for
            a lot of applications like Fourier transforms.

        t0 : float, optional
            If provided, this is the start-time, the time of the first sample.

        names: list, optional
            If provided, this should be a list of names for the data channels. 
            If it is a simple list (or list-like) then it the names are
            to apply to the channel rows in order. If it is a list of lists, or
            2-d array, then each of those sub-lists should be like the above
            and any of the lists are assumed to be applicable descriptors of the
            data.  The first entry will be taken if names are needed for display
        '''

        self.data = np.array(data)
        self.dt = dt
        self.t0 = t0
        self.names = names
        if names is not None:
            names=np.array(names)
            if len(names.shape)==1:
                names=names.reshape((1,-1))
        self.names=names

    def match_names(self, name_list, verbose=False):
        '''
        Check that a version of this object's names exactly matches the list.
        Returns True or False. If verbose, then print an explanation if False.
        '''
        name_list=np.array(name_list)
        if self.names is None:
            if verbose: print("Failed. No names available")
            return False
        if self.names.shape[1] != len(name_list):
            if verbose: print('Failed. Given',len(name_list),"names to match against object's",names.shape[1],'names.')
            return False
        if np.all(self.names==name_list): return True
        print('Failed. Names did not match.')

    def n_channels(self): return self.data.shape[0]

    def n_samples(self): return self.data.shape[1]

    def get_times(self,t0=None):
        if self.t0 is None:
            assert t0 is not None, "Must have a value for t0 to define times."
        else:
            t0=self.t0
        return np.arange(self.n_samples())*self.dt+t0

    def get_index_at_time(self,time):
        return (time-self.t0)/self.dt

    def check_times(self,other):
        assert self.t0 is not None, "Need a value for t0."
        assert other.t0 is not None, "Need a value for other.t0."
        assert self.t0==other.t0, "t0 values must match"
        assert np.isclose(self.t0,other.t0), "t0 values must agree."
        assert np.isclose(self.dt,other.dt), "dt values must agree."
    def subtract(self,other):
        assert self.data.shape==other.data.shape, "Number of channels and samples must agree to subtract."
        self.check_times(other)
        if self.names is not None or other.names is not None:
            assert self.names==other.names, "Names must match if present."
        return TimeData(self.data-other.data,dt=self.dt,t0=self.t0,names=self.names)
    
    def get_range(self,istart,iend=None):
        #extract a subset of a TD (to be added to TD)
        print(self.data.shape)
        if iend is None:
            iend=self.n_samples()
        newdata = self.data[:,istart:iend]
        #print('-->',newdata.shape)
        newt0=None
        newdt=None
        newnames=None
        if self.t0 is not None:
            newt0 = self.t0 + self.dt*istart
        if self.dt is not None:
            newdt = self.dt
        if self.names is not None:
            newnames=self.names
        #print('old t0,dt:',td.t0,td.dt,'new t0,dt:',newt0,newdt)
        return TimeData(newdata,t0=newt0,dt=newdt,names=newnames)

