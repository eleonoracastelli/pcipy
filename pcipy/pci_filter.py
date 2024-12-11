#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:40:27 2024

@author: E Castelli, J Baker
"""

## Working toward a new class for handling the apci processing

class apci_data():
    '''
    Class for practical application of automatic pci. 
    
    
    '''
    def __init__(self, ydata, fs, nhalf=45, order=1, maxcompts=10, t0=None, dt=None):
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
            Timestep size to use for scaling (ie weighting) the higher-order data matrix sectors (def:1/ns). The default is None.

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
            
        self.ns = ns
        self.nc = nc
        self.nhalf = nhalf
        self.order = order
        self.maxcompts =  maxcompts
        
        #time weighting params
        self.dt = dt
        if dt is None: self.dt = 1/ns
        self.t0=t0
        if t0 is None: self.t0 =  ns//2*self.dt
        
            
        #define window to avoid zeros        
        self.window=np.ones(ns,dtype=bool)
        self.window[:nhalf]=0
        self.window[-nhalf:]=0
        self.ns=sum(window)
        
        datamatrix = self.build_data_matrix(ydata)
        print(datamatrix.shape)

        self.apply_pca(datamatrix, maxcompts)
        
        print(datamatrix.shape,self.components.shape)
        self.channels = self.components.dot(datamatrix.T).astype(np.float64)
        print('channels shape',self.channels.shape)
        print('channel variances:',np.flip(self.explained_variance[-(maxcompts+1):]))
    
        
    def t_of_i(self,i):
        return -self.t0 + self.dt*i
            
    def build_data_matrix(self, ys):
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
        tt = np.linspace(0, self.dt*ns, ns) - self.t0  ##Should be consistent with self.t_of_i(index)            

        datamatrix_list = [np.hstack([np.array([self.construct_shifted_series(ydata, ishift, 
                                                window=(tt**m))[window]  ###Maybe the optimal window weighting is m-dependent?? 
                                                for ydata in ys]).T 
                                      for ishift in np.arange(-nhalf, nhalf+1)])
                           for m in range(0, order+1)]

        print('shape in build datamatrix',np.shape(datamatrix_list))
        return np.hstack(datamatrix_list)
    
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

    def stencil_at_timestep(self,it=None, n_channels=None):
        '''
        Evaluate the effective time-dependent PCI[0] stencil for time step 'it'.         
        For higher-order aPCI the effective filter stencil is time-dependent.  If a representative 'constant'
                stencil is needed, one can use the stencil as evaluated at a particular

        Parameters
        ----------
        it : TYPE, optional
            the time index (relative to the PCI-generating data)  [default data central time]. The default is None.
        n_channels : TYPE, optional
            number of channel stencils to compute . The default is None.

        Returns
        -------
        v : TYPE
            DESCRIPTION.

        '''
        
        if it is None: it = self.ns//2
        
        t = self.t_of_i(it) 
        
        v_compts=self.components[-n_channels:].reshape((nchannels,self.order+1,-1))
        
        v = v_compts[:,0,:]
        tpow=1
        for i in range(order):
            tpow *= t
            v += tpow*v_compts[:,i+1,:]
            
        return v                        
        
    
    def apply_for_channels(self,ydata, n_channels=None):
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

        datamatrix = self.build_data_matrix(ydata)
        if n_channels is None: n_channels=self.maxcompt
        return self.components[-n_channels:].dot(datamatrix.T).astype(np.float64)
    
    def apply_for_channels_split_orders(self,ydata, n_channels=None):
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
        datamatrix = self.build_data_matrix(ydata)
        unstacked=datamatrix.reshape(datamatrix.shape[0],self.order+1,-1)
        print('unstacked data mtrix shape',unstacked.shape)
        print('comps shape',self.components[-n_channels:].shape)
        orderedcomps=self.components[-n_channels:].reshape(n_channels,self.order+1,-1)
        print('comps reshaped',orderedcomps.shape)
        
        return [orderedcomps[:,m].dot(unstacked[:,m,:].T).astype(np.float64) for m in range(self.order+1)]
        

    

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

        print(datamatrix.shape)
        #We set svd_solver='full' here for consistency, but the default 'auto' may be fine for most purposes
        #A minor issue with auto is that, under some conditions an approximation was used which caused the
        #last component explained variance to come out exactly zero, when the corresponding channel did not vanish
        pca=PCA(svd_solver="full").fit(datamatrix)
        self.explained_variance=self.explained_variance_=pca.explained_variance_
        self.components=self.components_=pca.components_[-maxcomponents:]
        print('diagonality test')
        print(np.matmul(self.components,self.components.T))
    
    
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
        
        print(datamatrix.shape)
        U, S, Vt = svd(datamatrix, full_matrices=True)
        self.explained_variance=self.explained_variance_=S**2/(len(datamatrix)-1)
        self.components=self.components_=Vt[-maxcomponents:]

        
    def v_pci(self, n_channels=None):
        '''
        Get relevant sector of V pci data diagonalizing matrix
        '''
        if n_channels is None: n_channels = self.maxcompts
        #return [self.components_[-n_channels+i, :] 
        #                      for i in range(n_channels)]
        return self.components[-n_channels, :] 
    
    def e_pci(self, n_channels=None):
        '''
        Projection of either the original data matrix, or one build from provided y_data.
        '''
        
        if n_channels is None: n_channels = self.maxcompts
        return self.channels[-n_channels:, :]        
            


#outside class for development; expect to be moved in     
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
