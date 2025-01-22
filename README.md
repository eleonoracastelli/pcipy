# PCIpy

Principal Component Interferometry (PCI) package to process LISA data. 

This repo builds on the contents of [`pylisa`](https://github.com/qbaghi/pylisa), but structures the code using Python classes.

It implements the aPCI method outlined in the [first](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.042006) and [second](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.122001) PCI papers.

## Installation 

Install `PCIpy` by cloning this repository and unzipping the source code. Then `cd` into the folder and use this command:

```python
python setup.py install
```

## Simulation

Run simulation using the scripts in the `simulation` folder. Simulations are built using the September 2024 release of the [LISA Simulation suite](https://gitlab.in2p3.fr/lisa-simulation), specifically `lisainstrument, pyTDI, lisaorbits`.

### Simulation of laser + secondary noises using LISA Instrument

The simulation is carried out within `noise_simulation.py`.

Simulation parameters match the simulation parameters defined in 
LISA-LCST-SGS-RP-006 "End-to-end Demonstration Data Analysis Pipeline"

#### Usage:
Execute simulation with no TDI computation:
```
python noise_simulation.py path-to-workdir
```        

Execute simulation with TDI computation: use flag `--tdi` to specify TDI generation 
```    
python noise_simulation.py path-to-workdir --tdi 2   
```        

Execute simulation with baseline InRep configuration
```    
python noise_simulation.py path-to-workdir --baseline
```        

Execute simulation with baseline InRep configuration an save all individual noise contributions
```    
python noise_simulation.py path-to-workdir --baseline --individual
```

### Simulation of SGWB signal using LISA GW Response

The simulation is carried out within `signal_simulation.py` and `all_sky_signal_simulation.py`.

Simulation of SGWB signal using LISA GW Response.
`signal_simulation.py` simulates a stochastic point source located at $\beta = \pi/2$, $\lambda = \pi/2$, while `all_sky_signal_simulation.py` simulates a Stochastic GW background with a white generator.

#### Usage:
Execute simulation with no TDI computation (analogous for `all_sky_signal_simulation.py`):
```
python signal_simulation.py path-to-workdir
```        

Execute simulation with TDI computation (analogous for `all_sky_signal_simulation.py`): use flag `--tdi` to specify TDI generation 
```    
python signal_simulation.py path-to-workdir --tdi 2   
```
## Authors

Q Baghi, J Baker, E Castelli
