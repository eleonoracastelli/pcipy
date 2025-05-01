#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jan 16 20:00:05 UTC 2025

@author: J Baker, E Castelli,
"""

import numpy as np
from scipy.signal import convolve
from .data import TimeData


class LinearFilter:
    '''
    This is a simple base class/interface for PCI as a data filter but it can 
    also apply parallely to TDI or something else.  This base class addresses
    only the functionality for linear filtering itself, not for the filter 
    design and development.  

    The generally support base support is for filtering which takes as input
    a set of several regularly sampled data streams, operates on them with
    some kind of multi-channel quasi-convoluions (which may be time-varying), 
    and generates a set of several output channels.

    The base class encodes a trivial identity-filter operation.

    TBD: Probably the base-class should also support (optional) temporal 
    labeling.  Maybe also support output channel names and some other 
    practical features.  Is support for a time-dependent filter part of the
    interface? Probably so.
    '''

    '''
    Possible stuff to reference from pci_filter (which will inherit)
    def __init__(self, ydata, fs, nhalf=45, order=1, maxcompts=10, t0=None, dt=None,zero_mean=False,detrend=False,sort_by_rms=False,Tscale=1, ref_t0=None, ref_dt=None, verbose=False)

    def i_of_ref_time(self, time)

    def stencil_at_time(self, tref=None)

    def generate_stencil_grid(self, tref_start, dtref, nsamp)

    def apply_stencil_for_channels(self,ydata,tref_start,dtref) ##rename to "apply"
    '''
    def __init__(self, nleft, n_channels, input_names=None, nright=None):
        '''
        Base constructor for LinearFilter.

        This base class addresses only the functionality for linear filtering 
        itself, not any specifics of filter design and development.  It 
        implements a trivial identity filter.

        Parameters
        ----------
        nleft : int
            The left-half width of the applied filter stencil. The output of filtering
            will be shorter than the input series by nleft samples on the left end.
        n_input_channels: int
            The number of channels that the data is expected to be applied to.
        input_names: list
            If provided, the channels of the data provided for filtering can 
            be checked against these names.
        nright : int, OPTIONAL
            The right-half width of the applied filter stencil. The output of filtering
            will be shorter than the input series by nright samples on the right end.
            By default, nright=nelft
        
        Interface data:
            Each of nhalf, n_channels and input_names are stored internally.
            We also generally store the following:
            
            output_names: list or None
                Names of the output channels provided with outptu data.
            dt: float
                In/out data sample rate, verified on filtering and provided
                with filter output data.
            stencil_compts: ndarray
                The components of the convolution stencil realizing the filter.
                To be referenced only internally or supplanted by derived class
                Must be set by derived class if the default apply_filter is
                to work. 
                Shape must be (n_output_channels, 1+nleft+nright, n_input_channels)
            constant_stencil: bool 
                If True then no time info is needed to apply the filter 
                otherwise this but be incorporated in the input_data TimeData 
                object. Must be true for default apply_filter.
        '''
        self.nleft = nleft
        if nright is None: nright=nleft
        self.nright=nright
        self.n_input_channels=n_input_channels
        self.input_names=input_names
        #Default behavior is mathematically trival, but goes through the motions
        self.n_output_channels=n_input_channels
        self.stencil_compts=np.zeros((self.n_output_channels,1+nleft+nright))
        self.stencil_compts[:,nleft]=1
        self.constant_stencil=True
        self.dt=None

    def check_data(self, input_data, verbose=False):
        '''
        Perform basic check on the structure of the filter input data.

        Parameters
        ----------
        input_data : TimeData
            The data to be checked,        
        '''

        if input_data.n_channels()!=self.n_input_channels:
            if verbose: print('Wrong count of channels. Expected '+str(self.n_input_channels)+' but got '+str(input_data.n_channels()))
            return False
        if self.dt is not None:
            if self.dt != input_data.dt:
                if verbose: print('Sampling cadence wrong.')
        if self.constant_stencil:
            if input_data.t0 is None:
                if verbose: print('Unless the stencil is constant the input data must have t0 set to realize the filtering.')
        if self.input_names is not None:
            return self.input_data.match_names(self.input_names,verbose)
        
        return True

    def apply_filter(self, input_data, check=True,method='convolve'):
        '''
        Apply the encoded filter to the input_data.

        The default base version of this computes:
          output_data[ioc,i]=
              sum over j in range(-nleft,nright+1):  
                   stencil[ioc,iic,j] * input_data[iic,i+j]

        Parameters
        ----------
        input_data : TimeData
            The data to be filtered.
        check : bool, optional
            Whether to check the data before applying the filter (def True).
        method : str
            Variants on how to realize the computation, 'dot' for a direct
            approach using np.dot, or 'convolve' using scipy.signal.convolve 
        '''
        assert self.check_data(input_data), 'Data check was:'+str(self.check_data(input_data,verbose=True))
        assert self.constant_stencil, 'Base class apply_filter requires a constant stencil.'

        #Not sure what the fastest implementation of this is
        ns=input_data.n_samples()
        ne=ns-self.nleft-self.nright
        nwid=self.nleft+self.nright+1
        data=np.zeros((self.n_output_channels,ne))
        if method=='dot':
            for ioc in range(self.n_output_channels):
                print(ioc,self.nleft,self.nright)
                for i in range(nwid):
                    data[ioc]+=np.dot(self.stencil_compts[ioc,i],input_data.data[:,i:ne+i])
        elif method=='convolve':
            for ioc in range(self.n_output_channels):
                print(ioc,self.nleft,self.nright,self.stencil_compts.shape,input_data.data.shape)
                for i in range(self.n_input_channels):
                    data[ioc]+=convolve(self.stencil_compts[ioc,:,i][::-1],input_data.data[i,:],mode='valid')

        else: raise ValueError('Invalid value for "method"')
        t0=None
        if input_data.t0 is not None and self.dt is not None: t0=input_data.t0+self.nleft*self.dt
        return TimeData(data, dt=self.dt, t0=t0, names=self.output_names)



            

        
        
            

        
        

        
