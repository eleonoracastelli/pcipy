# PCIpy documentation

Principal Component Interferometry (PCI) package to process LISA data.

This repo builds on the contents of [`pylisa`](https://github.com/qbaghi/pylisa), but structures the code using Python classes.

It implements the aPCI method outlined in the [first](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.042006) and [second](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.122001) PCI papers.

```{eval-rst}
.. note::

    This project is under active development.
```   
   
## Contents

```{eval-rst}
.. toctree::
    :maxdepth: 1
    :caption: Getting started
   
    Home <self>
    source/PCI
    source/usage
```

```{eval-rst}
.. toctree::
    :maxdepth: 1
    :caption: Code implementation
   
    source/pcifilter
```


## Installation

Install `pcipy` by cloning this repository and unzipping the source code.

Then use this command:

```{eval-rst}
.. code-block::

    python setup.py install
```

### Requirements

Required installations listed in `setup.py`:
- python dependencies: `numpy`, `scipy`,`sympy`, `h5py`, `matplotlib`, `xarray`, `h5py`, `scikit-learn`
- simulation packages: `lisaconstants`, `lisainstrument`, `lisagwresponse`, `pytdi`, `backgrounds`
    
## Authors
This package build on the code available within [`pylisa`](https://github.com/qbaghi/pylisa), authored by Quentin Baghi and John Baker.


- Quentin Baghi 
- John G. Baker
- Eleonora Castelli