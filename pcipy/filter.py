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

    TBD: Possible the base-class should also support (optional) temporal 
    labeling.  Maybe also support output channel names and some other 
    practical features.  Is support for a time-dependent filter part of the
    interface? Probably so.
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
            The left-half width of the applied filter kernel. The output of filtering
            will be shorter than the input series by nleft samples on the left end.
        n_input_channels: int
            The number of channels that the data is expected to be applied to.
        input_names: list
            If provided, the channels of the data provided for filtering can 
            be checked against these names.
        nright : int, OPTIONAL
            The right-half width of the applied filter kernel. The output of filtering
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
            kernel_compts: ndarray
                The components of the convolution kernel realizing the filter.
                To be referenced only internally or supplanted by derived class
                Must be set by derived class if the default apply_filter is
                to work. 
                Shape must be (n_output_channels, 1+nleft+nright, n_input_channels)
            constant_kernel: bool 
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
        self.kernel_compts=np.zeros((self.n_output_channels,1+nleft+nright))
        self.kernel_compts[:,nleft]=1
        self.constant_kernel=True
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
        if self.constant_kernel:
            if input_data.t0 is None:
                if verbose: print('Unless the kernel is constant the input data must have t0 set to realize the filtering.')
        if self.input_names is not None:
            return input_data.match_names(self.input_names,verbose)
        
        return True

    def apply_filter(self, input_data, check=True,method='convolve'):
        '''
        Apply the encoded filter to the input_data.

        The default base version of this computes:
          output_data[ioc,i]=
              sum over j in range(-nleft,nright+1):  
                   kernel[ioc,iic,j] * input_data[iic,i+j]

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
        assert self.constant_kernel, 'Base class apply_filter requires a constant kernel.'

        #Not sure what the fastest implementation of this is
        ns=input_data.n_samples()
        ne=ns-self.nleft-self.nright
        nwid=self.nleft+self.nright+1
        data=np.zeros((self.n_output_channels,ne))
        if method=='dot':
            for ioc in range(self.n_output_channels):
                #print(ioc,self.nleft,self.nright)
                for i in range(nwid):
                    data[ioc]+=np.dot(self.kernel_compts[ioc,i],input_data.data[:,i:ne+i])
        elif method=='convolve':
            for ioc in range(self.n_output_channels):
                #print(ioc,self.nleft,self.nright,self.kernel_compts.shape,input_data.data.shape)
                for i in range(self.n_input_channels):
                    data[ioc]+=convolve(self.kernel_compts[ioc,:,i][::-1],input_data.data[i,:],mode='valid')

        else: raise ValueError('Invalid value for "method"')
        t0=None
        if input_data.t0 is not None and self.dt is not None: t0=input_data.t0+self.nleft*self.dt
        return TimeData(data, dt=self.dt, t0=t0, names=self.output_names)

    
class PiecewiseFilter(LinearFilter):
    '''
    This class builds a time-dependent filter from a set of subfilter kernels 
    which are accurate only locally in time. The current implementation 
    constructs the resultant filter by linear interpolation between the 
    nearest two local subfilters. 

    We expect to develop this further.
    '''

    def __init__(self, tstart, tend, nkern, subfilter_class, **subfilter_kwargs):
        '''
        Parameters:
          tstart  float  
             The start time for the resultant filter's expected range of applications.
          tend    float  
             The end time for the resultant filter's expected range of application
          nkern   int    
             The number of local subfilter kernels.
          subfilter_class  a class inheriting LinearFilter
             This should be the type of class to be applied for the local kernels
          subfilter_kwargs dict
             Additional arguments for the subfilter class constructor. Note that it is expected that there is an eval_time constructor argument which specifies the time at which the kernel is designed to be accurate.
        '''
        
        self.construct_subfilters(subfilter_class, subfilter_kwargs, tstart, tend, nkern)
        
        ## We set some of the internal variables expected with LinearFilter
        
        # nleft is not clearly meaningful
        # nright is not clearly meaningful
        self.n_input_channels=self.subfilters[0].n_input_channels
        self.input_names=self.subfilters[0].input_names
        self.n_output_channels=self.subfilters[0].n_output_channels
        # kernel_compts is not clearly meaningful
        self.constant_kernel=False
        self.dt=self.subfilters[0].dt


                 
    def apply_filter(self, input_data, check=True, method='convolve'):
        '''
        Apply the encoded filter to the input_data.

        In this implementation the filter produces an output channel set, based on the input_data, from the previously 
        generated set of input time-domain filters by piecewise linear interpolation.

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
        return PiecewiseFilter.piecewise_linear_apply_filter_set(input_data, self.subfilters, self.kernel_times, check=check, method=method)
        
    def construct_subfilters(self, subfilter_class, subfilter_kwargs, tstart, tend, nkern):
        '''
        This function assumes that subfilter constructor has an argument "eval_time" which is the only argument which differs as the
        in the construction of the subfilter instances.  This function may be used by a child class, or replaced with something more
        sophisticated.

        Args:
          tstart           double     Start time for filter application range
          tend             double     End time for filter application range
          nkern            int        Number of kernels pieces to use in piece-wise linear combination
          subfilter_class  class      Class for the subfilters
          subfilter_kwargs dict       Additional arguments for the sub_filter constructor
        '''
        #1. First determine the time grid for kernel definition
    
        #first define on unit interval
        ktimes = ( 0.5 + np.arange(nkern) ) / ( nkern )
        #print(ktimes)
        #then scale for the target time-stretch
        ktimes = tstart + (tend - tstart) * ktimes
        self.kernel_times=ktimes

        #2. Next construct the filters. For this example, we use raw LinearFilter objects, but a practical implementation will use a derived class        
        self.subfilters = [ subfilter_class(eval_time=t,**subfilter_kwargs) for t in self.kernel_times]
    
        
    # Internally ulitized functions
    
    def select_kernel_times(td,nkern):
        #define kernel times
    
        #first define on unit interval
        ktimes = ( 0.5 + np.arange(nkern) ) / ( nkern )
        #print(ktimes)
        #then scale for td
        ktimes = td.t0 + td.dt * td.n_samples() * ktimes
        #print(ktimes)
        return ktimes

    def piecewise_linear_apply_filter_set(input_data, filter_set, kernel_times, **filter_kwargs):
        '''
        This function produces an output channel set, based on the input_data, from a set of input time-domain filters by 
        piecewise linear interpolation.
    
        The filters are expected to be optimized at the corresponding kernel times, and time-ordered.
        Then, working across the time domain, the nearest filters are selected pairwise and linear interpolation
        provides the result for the corresponding times. 
        '''
    
        # 1. First some consistency checks. All of these might be obviated if controlled separately
        # Times should be ordered
        assert np.all(kernel_times[:-1] < kernel_times[1:]), "Values in kernel_times must be ordered and increasing"
        # Output channels should match
        for i in range(len(filter_set)-1):
            assert filter_set[i].n_output_channels==filter_set[i+1].n_output_channels and np.all(filter_set[i].output_names==filter_set[i+1].output_names),'Output channels do not make for filters '+str(i)+' and '+str(i+1)
        # For simplicity we demand that the kernel stencils be the same for each filter. This isn't strictly necessary
        # for the interpolation to make sense, but it does simplify the calculation a little.  Could be relaxed later.
        nleft=[filt.nleft for filt in filter_set]
        nright=[filt.nright for filt in filter_set]
        assert np.all(nleft[:-1]==nleft[1:]) and np.all(nright[:-1]==nright[1:]), "Not all filter stencils match."
        nleft=nleft[0]
        nright=nright[0]
        #Require commensurate dt
        dt=[filt.dt for filt in filter_set]
        assert np.all((dt[:-1]==dt[1:])),"Not all filter cadences agree."
        dt=dt[0]
        
        # 2. Identify the extent of the 'valid' part of the output time domain.
        output_t0=input_data.t0+nleft*dt
        output_n=input_data.n_samples()-nleft-nright
        output_times=output_t0+np.arange(output_n)*dt
        output_data=np.zeros((filter_set[0].n_output_channels,output_n))
        
        # 3. Next we need to split the output time domain into a set of subdomains with an assigned pair of nearest filters
        # The subdomains are bounded by the points of nearest kernel times
        kernel_time_indices=np.array([int(x) for x in ((kernel_times-output_t0)/dt+0.5)])
        
        print('kti',kernel_times,kernel_time_indices)
        select_bounds=np.arange(1,len(kernel_times)-1,dtype=int) #No bound associated with outermost filters since we need filter pairs
        #enforce that internal bounds should be strictly internal
        print('sel bounds',select_bounds)
        if len(select_bounds>0): 
            print('ktis',kernel_time_indices[select_bounds])
            select_bounds=select_bounds[kernel_time_indices[select_bounds]>0]
        print('sel bounds',select_bounds)
        if len(select_bounds>0): select_bounds=select_bounds[kernel_time_indices[select_bounds]<output_n-1]
        print('sel bounds',select_bounds)
        istarts=np.concatenate(([0],[kernel_time_indices[i] for i in select_bounds])).astype(int)
        iends=np.concatenate(([kernel_time_indices[i] for i in select_bounds],[output_n])).astype(int)
        k0=0
        if len(select_bounds)>0 and select_bounds[0]>0: k0=select_bounds[0]-1
        left_kernel_inds=np.concatenate(([k0],select_bounds))
        if k0>=len(kernel_times)-1: right_kernel_inds=left_kernel_inds
        else:
            right_kernel_inds=left_kernel_inds+1
        nsegs=len(istarts)
        print('n,starts,ends,kleft,kright',nsegs, istarts, iends, left_kernel_inds, right_kernel_inds)
    
        # 4. Now do apply the filter pairs and linearly interpolate to get the output for each subdomain segment
        #We begin by initiating the left filter version of data
        in_data=input_data.get_range(istarts[0],iends[0]+nleft+nright)
        filt=filter_set[left_kernel_inds[0]]
        left_data=filt.apply_filter(in_data,**filter_kwargs)
        left_istart=istarts[0]
        for j in range(nsegs):
            istart=istarts[j]
            iend=iends[j]
            kleft=left_kernel_inds[j]
            kright=right_kernel_inds[j]
            # A. Compute the right filter version of data 
            jnext=j+1
            if jnext>=nsegs: jnext=nsegs-1
            #we apply the right filter so that it can also be the left filter of next step (if it exists)
            iend_next=iends[jnext] 
            # Note that input data is staggered by nleft from output data
            # and that we also need a buffer by nleft on the left and nright on the right to get data for our output range
            in_data=input_data.get_range(istart,iend_next+nleft+nright)
            filt=filter_set[kright]
            right_data=filt.apply_filter(in_data,**filter_kwargs)
            
            # B. perform linear interpolation on this subdomian
            seg_times = output_times[istart:iend]
            eps = (seg_times-kernel_times[kleft])/(kernel_times[kright]-kernel_times[kleft]+1e-100)
            print(eps)
            print('left',left_data.n_samples(),istart-left_istart,iend-left_istart)
            left_seg_data = left_data.data[:,istart-left_istart:iend-left_istart]
            print('right',right_data.n_samples(),0,iend-istart)
            right_seg_data = right_data.data[:,:iend-istart] 
            print('shapes',left_seg_data.shape,right_seg_data.shape)
            seg_data = left_seg_data + (right_seg_data-left_seg_data)*eps
            output_data[:,istart:iend]=seg_data
    
            # C. Swap right to left
            left_data=right_data
            left_istart=istart
            
        return TimeData(output_data, dt, output_t0, filter_set[0].output_names)

    
class RestrictedFilter(LinearFilter):
    '''
    This class allows to generate a version of a linear filter which is restricted to a apply only to
    a subset of the inputs or outputs.  When applied to a subset of inputs the other input channels are 
    effectively assumed to vanish.  Initial implementation assumes a constant kernel.
    '''

    def __init__(self,other, input_names=None,output_names=None):

        if not other.constant_kernel: raise NotImplementedError("Filter restriction not yet implemented for non-constant filters")
        
        #We copy everything from other except that we restrict the input or output channels
        self.nleft = other.nleft
        self.nright = other.nright
        if input_names is None:
            #Do nothing different if None
            self.input_names = other.input_names
            self.n_input_channels = other.n_input_channels
            input_name_map=None
        else:
            assert other.input_names is not None, "Cannot restrict input channels on a filter which doesn't have input names."
            input_name_map=[]
            for name in input_names:
                locs=np.array(other.input_names)==name
                count=np.sum(locs)==1
                assert count, 'In restricting input channels, found "'+name+'" '+str(count)+' times.'
                input_name_map.append(np.argmax(locs))
            if 1:
                print('Input name mapping:')
                for i in range(len(input_names)):
                    print(str(input_name_map[i])+':'+other.input_names[input_name_map[i]]+' --> '+str(i)+':'+input_names[i])
            self.input_names = input_names
            self.n_input_channels = len(self.input_names)
        if output_names is None:
            #Do nothing different if None
            self.output_names = other.output_names
            self.n_output_channels = other.n_output_channels
            output_name_map=None
        else:
            assert other.output_names is not None, "Cannot restrict output channels on a filter which doesn't have output names."
            output_name_map=[]
            for name in output_names:
                locs=np.array(other.output_names)==name
                count=np.sum(locs)==1
                assert count, 'In restricting output channels, found "'+name+'" '+str(count)+' times.'
                output_name_map.append(np.argmax(locs))
            if 1:
                print('Output name mapping:')
                for i in range(len(out_names)):
                    print(str(output_name_map[i])+':'+other.output_names[output_name_map[i]]+' --> '+str(i)+':'+output_names[i])
            self.output_names = output_names
            self.n_output_channels = len(self.output_names)

        if isinstance(self,PiecewiseFilter):
            self.dt=other.dt
            self.contant_kernel=other.constant_kernel
            self.subfilters=[RestrictedFilter(filt, input_names,output_names) for filt in other.subfilters]
        else:    
            kernel_compts=other.kernel_compts
            if input_name_map is not None: kernel_compts=kernel_compts[:,:,input_name_map]
            if output_name_map is not None: kernel_compts=kernel_compts[output_name_map,:,:]
            self.kernel_compts=kernel_compts
            print('kernel_compts shape changed from',other.kernel_compts.shape,'to',self.kernel_compts)
            self.constant_kernel=other.constant_kernel
            self.dt=other.dt
    

        
        
            

        
        

        
