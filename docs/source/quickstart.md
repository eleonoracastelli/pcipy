# Apply PCI to simulated data

Here we apply PCI to data simulated via the LISA Simulation Suite, following the simulation scripts described in <project:#simulation>.

The simulated datasets are 3 days long. PCI is applied to 12 hours of data, with 4 hours of skipped data at the beginning of the simulation.

We need the following simulated datasets:
- full simulation noise dataset (including laser noise and secondary noises), with filename ending in `_measurements_4Hz.h5`
- secondary noises dataset, with filename ending in `_noise_sec_4Hz.h5`

## 1. Build data vector of the six LISA single-link channels

To build the data vector of the six LISA single-link channels $\vec{y} = \left[y_{ij}\right]$, with $i,j=1,2,3$ and $i\neq j$ we resort to the intermediary TDI variables $\eta$, implemented within `pytdi` as `ETA_SET`.

```{eval-rst}
.. note::

    Eleonora: how can we justify being independent from TDI if we use the intermediary TDI variables :math:`\eta` to build the :math:`y` s?
```  

We build the single link $\vec{y}$ data vector for the full noise simulation and for the secondary noises, ending up with two single link vectors:
- full simulation single link vector $\vec{y}^{\text{full}}$
- secondary noises single link vector $\vec{y}^{\text{sec}}$

## 2. Apply PCI to the data

We now resort to {class}`PCIFilter` to evaluate PCI from $\vec{y}$.

An instance of {class}`PCIFilter` has two required inputs:
- `ydata`: matrix with the single link LISA temporal phase data streams $\vec{y}$.
- `fs`: sampling rate of the data streams (Hz).

The optional parameters are
- `nhalf`: filter stencil halfwidth in samples. The default is 45.
- `order`: order of PCI. The default is 1.
- `maxcompts`: PCA results will be truncated to this length after initial processing. The default is 10.
- `Tscale`: if dt is None, then dt defaults to Tscale/ns

The input channels $\vec{y}$ are stretches of data, usually of length `ns+2*nhalf`, but sometimes varying. The variations are:
  : (full length of matrix)
  0:ns_fixed
  0:ns+2*nhalf
  skip:skip+ns+2*nhalf
In every case the window is trivial np.ones([data lenght])  

The PCIFilter class applies the following methods
{method}`self.build_data_matrix(ydata, zero_mean=zero_mean,detrend=detrend)`

{method}`self.apply_pca(datamatrix, maxcompts)`

{method}`self.set_stencil(self.maxcompts)`

{method}`self.channels = self.components.dot(datamatrix.T).astype(np.float64)`

## 3. Run channel analysis on the PCI output

## 4. Reconstruct single-link channels

## 5. Estimate sensitivity