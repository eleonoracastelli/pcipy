# Implementation of a PCI class



```{py:class} #pcipy.PCIFilter
    Class for practical application of automatic pci. 

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

```


