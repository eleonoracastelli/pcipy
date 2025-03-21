#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:40:27 2024

@author: E Castelli, J Baker
"""

import numpy as np
from sklearn.decomposition import PCA

## Working toward a new class for handling the apci processing
"""
Some plan ideas:
 Want to shift the basic function of PCIFilter class to be for application downstream rather than for the PCI part itself.  This could be via a parent "Filter" class (which could also have a TDI version).  The filter can be used to get a sencil matrix for any regular grid of physical times. For a basic PCI filter this would need to store just the time-referencing info and a set of stencil components.  It would also need the basic generation functions.
 Most of what we have in the current PCIFilter is extra.  We could extract a little from that into a new base class and let the current one become special for experimentation (or drop it).
 For the moment, we keep building in here, but plan to soon split this up into a number of derived class options, each of which will be concrete about the output channels.
"""
class PCIFilter():
    '''
    Class for practical application of automatic pci. 
        
    '''
    def __init__(self, ydata, fs, nhalf=45, order=1, maxcompts=10, t0=None, dt=None,zero_mean=False,detrend=False,sort_by_rms=False,Tscale=1, ref_t0=None, ref_dt=None, verbose=False):
        '''                                   
        Parameters
        ----------
        ydata : ndarray
            Matrix with the single link LISA temporal phase data streams.
        fs : float
            Sampling rate (Hz).
        nhalf : int, optional
            Filter stencil halfwidth. The default is 45.
        order : TYPE, optional
            DESCRIPTION. The default is 1.
        maxcompts : int, optional
            PCA results will be truncated to this length after initial processing.. The default is 10.
        t0 : float, optional
            Zero-time (relative to start of data to use in stencil definition (def: center of data; Note set t0=0 for default of 2023-24). The default is None.
        dt : float, optional
            Timestep size to use for scaling (ie weighting) the higher-order data matrix sectors (def:Tscale/ns). The default is None.
        Tscale : if dt is None, then dt defaults to Tscale/ns
        ref_t0 : external time reference corresponding to the first data sample
        ref_dt : external time reference sampling rate

        Returns
        -------
        None.        
        
        The input channels are stretches of data, usually of length ns+2*nhalf, but sometimes varying. The variations are:
          : (full length of matrix)
          0:ns_fixed
          0:ns+2*nhalf
          skip:skip+ns+2*nhalf
        In every case the window is trivial np.ones([data lenght])  

        '''
        
        # Data sample size
        self.nsdata = ydata[0].shape[0]
        
        # Number of channels
        if type(ydata) == list:
            nc = len(ydata)
        else:
            nc = ydata.shape[0]
            
        self.nc = nc
        self.nhalf = nhalf
        self.ns = self.nsdata-2*nhalf
        self.order = order
        self.maxcompts =  maxcompts
        self.sort_by_rms=sort_by_rms
        
        #external time reference
        #default is zero start and unit cadence
        if ref_t0 is None:
            self.ref_t0=0
        else:
            self.ref_t0=ref_t0
        if ref_dt is None:
            self.ref_dt=1
        else:
            self.ref_dt=ref_dt       
        
        #time weighting params
        self.scaled_dt = dt
        if dt is None: self.scaled_dt = Tscale/self.ns
        self.scaled_t0=t0
        if t0 is None: self.scaled_t0 =  self.ns//2*self.scaled_dt
        
        #define window to avoid zeros        
        self.window=np.ones(self.nsdata,dtype=bool)
        self.window[:nhalf]=0
        self.window[-nhalf:]=0
        
        datamatrix = self.build_data_matrix(ydata, zero_mean=zero_mean,detrend=detrend)
        #print('dm shape',datamatrix.shape)
        
        self.apply_pca(datamatrix, maxcompts)
        #print('pca done')
        
        #stencil default specification
        self.set_stencil(self.maxcompts)
        #print('stencil done')
        
        #print(datamatrix.shape,self.components.shape)
        self.channels = self.components.dot(datamatrix.T).astype(np.float64)
        #print('channels done')
        
        if verbose:
            print("datamatrix mean", np.mean(datamatrix))
            print('channels shape',self.channels.shape)
            print('channel means',np.flip(self.channels.mean(axis=1)))
            print('channel variances:',np.flip(self.explained_variance[-(maxcompts+1):]))
        del(datamatrix)
        
    def i_of_ref_time(self, time):
        '''
        Compute the grid time (as float) from the reference time. This function does not necessarily  generalize to multiple grid versions.
        '''
        time_i = time/self.ref_dt+self.ref_t0
        return time_i
    
    def scaled_t_of_i(self,i):
        return -self.scaled_t0 + self.scaled_dt*i
            
    def build_data_matrix(self, ys, zero_mean=False, detrend=False):
        """
        Pre-process data to build a matrix of features of size n_features x n_data.

        Parameters
        ----------
        y_list : list[ndarrays]
            list of channel data
        nhalf : int
            half-size of shift kernel
        order : int
            order of the Taylor expansion in time

        Returns
        -------
        datamatrix : ndarray
            data matrix of size n_data x (2*nhalf+1)*nc

        Comment: Versions of this func had a "window" option used in periodograms, but this seems not to have been used in any PCI 2.0 notebooks so we leave it out here
        """ 
        detrend_before=True
        detrend_after=False
        
        ###Plan to generalize to work with a sub-stretch of the data???
        ns = len(ys[0])
        nhalf=self.nhalf
        order=self.order

        window=np.ones(ns,dtype=bool)
        window[:nhalf]=0
        window[-nhalf:]=0
        
        # Time vector
        ### Plan to generalize this because the results can depend on how tt is scaled (and offset). 
        ###Control of this is likely crucial for compound-segment treatments
        ###Weighting may also depend on time power
        tt = np.linspace(0, self.scaled_dt*ns, ns) - self.scaled_t0  ##Should be consistent with self.scaled_t_of_i(index)            

        ##Note that we do detrending before timeshifting and windowing and multiplying by T^n.  This doesn't guarantee exactly. Other options can be considered.        
        if detrend_before:
            if detrend:
                from scipy.signal import detrend as scipy_detrend 
                print(ys.shape)
                ys=scipy_detrend(ys)
            elif zero_mean:
                print('ys',ys.shape)
                print('means',np.mean(ys,axis=1))
                ys = (ys.T - np.mean(ys,axis=1).T).T
                print('new means',np.mean(ys,axis=1))
        
        datamatrix_list = [np.hstack([np.array([self.construct_shifted_series(ydata, ishift, 
                                                window=(tt**m))[window]  ###Maybe the optimal window weighting is m-dependent?? 
                                                for ydata in ys]).T 
                                      for ishift in np.arange(-nhalf, nhalf+1)])
                           for m in range(0, order+1)]

        matrix = np.hstack(datamatrix_list)
        #print('shape in build datamatrix',np.shape(matrix))
        if detrend_after:
            if detrend:
                from scipy.signal import detrend as scipy_detrend 
                matrix=scipy_detrend(matrix,axis=0)
            elif zero_mean:
                print('means',np.mean(matrix,axis=0))
                matrix = (matrix - np.mean(matrix,axis=0))
                print('new means',np.mean(matrix,axis=0))
        return matrix
    
    def construct_shifted_series(self, p, ishift, window=1.0):
        """
        Parameters
        ----------
        p : ndarray
            array of size ns
        ishift : int
            delay to apply expressed in number of samples
        window : ndarray or float, optional
            time window of size ns , by default 1.0

        Returns
        -------
        p_shift : ndarray
            Shifted array of size n_data
        """
        ns = p.shape[0]
        pshift = np.zeros_like(p)
        pshift[max([0,ishift]):min(ns,ns+ishift)] = p[max([0,-ishift]):min(ns,ns-ishift)]

        return pshift * window

    def set_stencil(self, n_channels):
        '''
        For use with the stencil functions below.  We specify n_channels to fix the stencil components externally as these functions to support a generic interface where the channels are fully pre-specified.
        The argument is to specify the number of channels. Internally this sets up self.stencil_compts
        which has shape (n_channels, order+1, 6*(1+2*nhalf))
        '''
        self.stencil_n_channels=n_channels
        self.stencil_compts=self.components[-n_channels:].reshape(n_channels,self.order+1,-1)
    
        
    def stencil_at_time(self, tref=None):
        '''
        Evaluate the effective time-dependent PCI[0] stencil for time step 'it'.         
        For higher-order aPCI the effective filter stencil is time-dependent.  If a representative 'constant'
                stencil is needed, one can use the stencil as evaluated at a particular

        Parameters
        ----------
        tref : TYPE, optional
            the external reference time value where the stencil should be evaluated  [default data central time]. 

        Returns
        -------
        v_grid : TYPE
            DESCRIPTION.

        '''
        n_channels=self.stencil_n_channels
        
        if tref is None:
            it = self.ns//2
            tref=self.ref_t0+it*self.ref_dt
        
        t = self.scaled_t_of_i(self.i_of_ref_time(tref)) 
        
        v_compts=self.stencil_compts
        
        v = v_compts[:,0,:]
        tpow=1
        for i in range(order):
            tpow *= t
            v += tpow*v_compts[:,i+1,:]
            
        return v

    def generate_stencil_grid(self, tref_start, dtref, nsamp):
        '''
        Generate the full stencil grid for computing output channelsEvaluate the effective time-dependent PCI[0] stencil for time step 'it'.         
        For higher-order aPCI the effective filter stencil is time-dependent.  If a representative 'constant'
                stencil is needed, one can use the stencil as evaluated at a particular

        Parameters
        ----------
        tref_start : float
            The external reference time value starting the grid where the stencil should be evaluated.
        dtref : float
            The reference time step size for the evaluation grid.
        nsamp : int
            Tthe number of samples in the evaluation grid

        Returns
        -------
        v : ndarray with shape=(nsamp,nchannels,6*(1+2*nhalf))
            Stencil grid
        '''

        times=np.arange(nsamp)*dtref+tref_start
        t0 = self.scaled_t_of_i(self.i_of_ref_time(times)) 
        
        v_compts=self.components[-n_channels:].reshape((nchannels,self.order+1,-1))
        
        v = v_compts[:,0,:]
        tpow=1
        for i in range(self.order):
            tpow *= t0
            v += np.outer(tpow,v_compts[:,i+1,:])
            
        return v                        

    def apply_stencil_for_channels(self,ydata,tref_start,dtref):
        '''
        This function is imcompletely drafted.  When finished, it should provide the same results as
        apply_for_channels, suitably called. This one fits a more generic interface (and may be more efficient)
        '''
        
        nd=len(ydata.T)
        nhalf=self.nhalf
        window=np.ones(nd,dtype=bool)
        window[:nhalf]=0
        window[-nhalf:]=0
        nsamp=sum(window)
        ndelay=self.stencil_compts.shape[2]
        nchan=self.stencil_compts.shape[0]

        # Treating the stincil as an FIR filter, we can avoid constructing the huge data matrix
        # We realize the convolution by hand in the time domain, though for longer kernels
        # it would be faster in Fourier domain
        delays=list(range(-nhalf,nhalf+1))
        times=np.arange(nsamp)*dtref+tref_start
        tt = self.scaled_t_of_i(self.i_of_ref_time(times)) 
        tpow=1
        tpows=[tpow]
        for i in range(self.order):
            tpow*=tt
            tpows+=[tpow]

        X=np.hstack([np.array([
                               self.construct_shifted_series(ydata, ishift, window=window)
                               for ydata in ys]).T 
                                                
                     for ishift in np.arange(-nhalf, nhalf+1)])

        stencil_grid=self.generate_stencil_grid(tref_start+nhalf*dtref, dtref, nsamp)

        return multipledot.something

    
    def apply_for_channels(self,ydata, n_channels=None, zero_mean=False, detrend=False):
        '''
        Apply the precomputed data filters to derive PCI results from distinct data

        Parameters
        ----------
        ydata : TYPE
            The new data.
        n_channels : TYPE, optional
            The number of PCI component channels to compute. The default is None.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''

        datamatrix = self.build_data_matrix(ydata,zero_mean=zero_mean,detrend=detrend)
        if n_channels is None: n_channels=self.maxcompts
        return self.components[-n_channels:].dot(datamatrix.T).astype(np.float64)
    
    def apply_for_channels_split_orders(self,ydata, n_channels=None,zero_mean=False,detrend=False):
        '''
        Apply the precomputed data filters to derive PCI results from distinct data
        but split into separate parts for the different temporal orders.  If the
        channel is defined by
        C_a = H V_a
        C_a = \sum_b C_{ab}
        C_ab = T^b X V_{ab}
        where T is the diagonal time-matrix used in the data matrix construction and 
        V_{ab} is the submatrix of V corresponging to channel a and time-order b
        
        Parameters
        ----------
        ydata : TYPE
            The new data.
        n_channels : TYPE, optional
            The number of PCI component channels to compute. The default is None.

        Returns
        -------
        list
            DESCRIPTION.
            
        
        PCA is always invoked as:
        pca = PCA().fit(datamatrix)
        in some form.  Note that this does not result in some kind of data normalization by default.

        '''
        
        #Undo the X-> H datamatrix stacking:
        
        if n_channels is None: n_channels=self.maxcompt            
        datamatrix = self.build_data_matrix(ydata,zero_mean=zero_mean,detrend=detrend)
        unstacked=datamatrix.reshape(datamatrix.shape[0],self.order+1,-1)
        print('unstacked data mtrix shape',unstacked.shape)
        print('comps shape',self.components[-n_channels:].shape)
        orderedcomps=self.components[-n_channels:].reshape(n_channels,self.order+1,-1)
        print('comps reshaped',orderedcomps.shape)
        
        return [orderedcomps[:,m].dot(unstacked[:,m,:].T).astype(np.float64) for m in range(self.order+1)]
    
    def filter_single_link_data(self, ydata,n_channels=None):
        '''
        Compute recovered single-link channel vectors from a set of single-link data. This should yield 
        identical results to appropriately called compute_single_link_from_channels but we first implement
        separately for testing, since the stencil_compts stuff is new. This one is more low-level and 
        doesn't assum the stencil is selected.
        
        Parameters
        ----------
        ydata : ndarray 
            The raw single-link data-set to be transformed.

        Returns
        -------
        reconstructed_ydata : ndarray
            Reconstructed single link data channels
        '''
        if n_channels is None: n_channels=self.maxcompts
        y_transformer=self.components[-n_channels:,:6].T

        chans_data = self.apply_for_channels(ydata, n_channels)
        Z=np.dot(y_transformer,chans_data)
        print('single link data shape', Z.shape)
        return Z
        
    
    def compute_single_links_from_channels(self,channels_data):
        '''
        Compute recovered single-link channel vectors from a set of pri-computed PCI channels using the stencil_compts.
        
        Parameters
        ----------
        channels_data : ndarray 
            The computed channel data.

        Returns
        -------
        reconstructed_ydata : ndarray
            Reconstructed single link data channels

        The algorithm transform's back to the Y-channel frame by using the transpose (inverse) of the "V" matrix, components of which 
        are used to select the channels. This would identically reconstruct the original data if we used all possible channels. This
        is related to "Zero-phase component analysis" (ZCA) often applied for whitening data.
        '''
        
        assert len(channels_data)==len(self.stencil_compts), "Stencil length ("+str(len(self.stencil_compts))+") doesn't match received channel count."
                                                                                     
        y_transformer = self.stencil_compts[:,0,:6].T  #The first six columns should correspont to the T^0, unshifted Y elems
        Z=np.dot(y_transformer,channels_data)
        print('single link shape', Z.shape)
        return Z


    def apply_pca(self, datamatrix,maxcomponents):
        '''
        Apply principal component analysis to the data matrix and extract relevant results explained_variance and (some) components

        Parameters
        ----------
        datamatrix : TYPE
            DESCRIPTION.
        maxcomponents : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        diagnostics=True
        
        print(datamatrix.shape)
        #We set svd_solver='full' here for consistency, but the default 'auto' may be fine for most purposes
        #A minor issue with auto is that, under some conditions an approximation was used which caused the
        #last component explained variance to come out exactly zero, when the corresponding channel did not vanish
        pca=PCA(svd_solver="full").fit(datamatrix)
        self.explained_variance=pca.explained_variance_
        components=pca.components_
        if self.sort_by_rms:
            ynew=components.dot(datamatrix.T).astype(np.float64)
            rms2=(np.sum(ynew**2,axis=1)/len(ynew.T))
            reordering=np.argsort(-rms2)
            self.explained_variance=self.explained_variance[reordering]
            components=components[reordering]
        self.components=components[-maxcomponents:]
        if diagnostics:
            ynew=components.dot(datamatrix.T).astype(np.float64)
            self.explained_rms=np.sqrt(np.sum(ynew**2,axis=1)/len(ynew.T))
            ynew_zm=(ynew.T-np.mean(ynew,axis=1)).T
            variance_check=np.sum(ynew_zm**2,axis=1)/(len(ynew.T)-1)
            #print('variance check:\n',variance_check, '\n',self.explained_variance)
            print('variance check rms:',np.sqrt(np.sum((variance_check/self.explained_variance-1)**2)))
            chans=ynew_zm[-maxcomponents:]
            nchan=len(chans)
            cov=np.matmul(chans,chans.T)/(len(ynew.T)-1)
            print('shapes of chans, cov, components:',chans.shape,cov.shape,self.components.shape)
            sigmas=np.sqrt([cov[i,i] for i in range(nchan)])
            print(sigmas**2/self.explained_variance[-nchan:])
            corr=np.array([[cov[i,j]/sigmas[i]/sigmas[j] for i in range(nchan)] for j in range(nchan)])
            diag_test=corr-np.identity(nchan)
            print('channel covariance diagonality test on '+str(nchan)+' components:',np.sqrt(np.sum(diag_test**2)/(nchan*(nchan-1))))    
            diag_test=np.matmul(self.components,self.components.T)-np.identity(maxcomponents)
            print('component diagonality test on '+str(maxcomponents)+' components:',np.sqrt(np.sum(diag_test**2)/(maxcomponents*(maxcomponents-1))))
    
    
    def apply_svd(self, datamatrix,maxcomponents):  ##should be equivalent, untested
        '''
        Apply principal component analysis to the data matrix and extract relevant results explained_variance and (some) components
        
         Parameters
         ----------
         datamatrix : TYPE
             DESCRIPTION.
         maxcomponents : TYPE
             DESCRIPTION.

         Returns
         -------
         None.
        '''
        
        #print(datamatrix.shape)
        U, S, Vt = svd(datamatrix, full_matrices=True)
        self.explained_variance=self.explained_variance=S**2/(len(datamatrix)-1)
        self.components=self.components_=Vt[-maxcomponents:]

        
    def v_pci(self, n_channels=None):
        '''
        Get relevant sector of V pci data diagonalizing matrix
        '''
        if n_channels is None: n_channels = self.maxcompts
        #return [self.components_[-n_channels+i, :] 
        #                      for i in range(n_channels)]
        return self.components[-n_channels:] 
    
    def e_pci(self, n_channels=None):
        '''
        Projection of either the original data matrix, or one build from provided y_data. TBD?  See aply_for_channels
        '''
        
        if n_channels is None: n_channels = self.maxcompts
        return self.channels[-n_channels:, :]        

    

class SplitFilter(PCIFilter):
    '''
    Class for practical application of automatic pci. 
    
    
    '''
    def __init__(self, ydata, fs, nhalf=45, order=1, maxcompts=10, t0=None, dt=None,zero_mean=False,detrend=False,sort_by_rms=False,Tscale=1):
        pass
        
#outside class for development; consider moving in     
# def covariance_by_link_psd_model(pci,freqs,s_n,n_components=None):
#     '''
#     Apply a model of the secondary noise PSD to get stantionary FD covariance matrix for the data channels

#     Parameters
#     ----------
#     pci : TYPE
#         DESCRIPTION.
#     freqs : TYPE
#         an freq grid array.
#     s_n : TYPE
#         a corresponding array of PSD values corresponding to link PSDs of secondary noises, all links the same.
#     n_components : TYPE, optional
#         DESCRIPTION. The default is None.

#     Returns
#     -------
#     cov.
#     '''


#     if n_components is None: n_components = pci.maxcompt

#     # Compute 1st order PCI secondary noises PSD from 00 covariance component
#     n_mat = np.array([np.eye(self.nc) for freq in freqs_a])
#     v_pci_const_mat = pci.stencil_at_timestep(pci.ns//2,n_components)

#     # Compute PCI equivalent vector for the noise
#     e_vect = [apci.compute_equivalent_vector(freqs_a, n_mat, v_pci_const_mat, nhalf, i, fs=fs) for i in range(n_components)]

#     # Compute secondary noise PCI covarinace
#     cov = sum([eigen.multiple_dot(e[:, :, np.newaxis], np.conj(e[:, np.newaxis, :]))
#                                      for e in e_vect])

#     cov *= s_n[:, np.newaxis, np.newaxis]  ##Note that we assume equality for link PSDs  ##FIXME
    
#     return cov


