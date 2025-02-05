## Simulation

PCI is applied directly to the LISA phasemeter measurements.

The LISA phasemeter measurements are the output of simulations built via the [LISA Simulation suite](https://gitlab.in2p3.fr/lisa-simulation). 

In what follows we use the following packages from the September LISA Simulation Suite:
- `lisainstrument` for the overarching instrument,
- `lisaorbits` for the orbits,
- `lisagwresponse` for the response to GW,
- `pytdi` for the evaluation of TDI variables.

In the package we provide three simulation scripts that can be executed in the terminal:

- `noise_simulation.py` to simulate the laser and secondary noises
- `signal_simulation.py` to simulate a stochastic point source in the sky
- `all_sky_signal_simulation.py` to simualte a stochastic graviational wave background.


### Simulation of laser + secondary noises using LISA Instrument

The simulation is carried out within `noise_simulation.py`.

Simulation parameters match the simulation parameters defined in 
LISA-LCST-SGS-RP-006 "End-to-end Demonstration Data Analysis Pipeline"

#### Usage:
- Execute simulation with no TDI computation:
```
python noise_simulation.py path-to-workdir
```        

- Execute simulation with TDI computation: use flag `--tdi` to specify TDI generation 
```    
python noise_simulation.py path-to-workdir --tdi 2   
```        

- Execute simulation with baseline InRep configuration
```    
python noise_simulation.py path-to-workdir --baseline
```        

= Execute simulation with baseline InRep configuration an save all individual noise contributions
```    
python noise_simulation.py path-to-workdir --baseline --individual
```

### Simulation of SGWB signal using LISA GW Response

The simulation is carried out within `signal_simulation.py` and `all_sky_signal_simulation.py`.

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