## Simulation

PCI is applied directly to the LISA phasemeter measurements.

The LISA phasemeter measurements are the output of simulations built via the [LISA Simulation suite](https://gitlab.in2p3.fr/lisa-simulation). 

In what follows we use the following packages from the September LISA Simulation Suite:
- `lisainstrument >= 1.8` for the overarching instrument simulation,
- `lisaorbits` for the orbits,
- `lisagwresponse` for the response to GW,
- `pytdi` for the evaluation of TDI variables.

In the package we provide three simulation scripts:

- `noise_simulation.py` to simulate the laser and secondary noises for various simulation scenarios;
- `signal_simulation.py` to simulate a stochastic point source in the sky;
- `all_sky_signal_simulation.py` to simulate a stochastic gravitational wave background.


### Simulation scenarios
To allow for benchmarking performance of PCI against TDI, we provide a list of simulation scenarios to test. Various options can be toggled to generate the data.

**Laser locking** Simulation scenarios can be generated with various types of laser locking:`
- `'N1-12'` pre-defined laser locking configuration N1 with 12 as primary laser
or
- `'six'` for 6 lasers locked on cavities,

**Orbits** Choosing between equal armlength and keplerian orbits.

**Secondary noises** Choosing between simpler LTO (laser + TM+ OMS) noise contrbutions, and the baseline noise contributions laid out in the LISA-LCST-SGS-RP-006 "End-to-end Demonstration Data Analysis Pipeline": LTO + ranging + clock + backlink + modulation

Table of simulated scenarios, for each locking configuration:

| `EqualArmlengthOrbits`   | `KeplerianOrbits`     | Comments |
| ------------- | ------------- | ------------- | 
|  | **LTO noise** | |
| only laser noise | only laser noise  | null-hypothesis scenario |
| tm noise | tm noise | individual noise sources |
| oms noise |  oms noise | individual noise sources |
|   laser noise + TM noise |laser noise + TM noise | noise source combined with laser |
|   laser noise + OMS noise | laser noise + OMS noise | noise source combined with laser |
| laser + TM+ OMS noise | laser + TM + OMS noise | full measurement simulation |
| TM + OMS noise | TM + OMS noise | secondary noises |
|  | **LISA-LCST-SGS-RP-006 baseline** |
| ---- |  baseline | full measurement simulation |
| ---- |  baseline secondary noises | secondary noises |
| ---- | ranging  | individual noise sources |
| ---- | clock  | individual noise sources |
| ---- | backlink  | individual noise sources |
| ---- | modulation  | individual noise sources |
| ---- | LTO + ranging  | noise source combined with LTO |
| ---- | LTO + clock | noise source combined with LTO |
| ---- | LTO + backlink  | noise source combined with LTO |
| ---- | LTO + modulation | noise source combined with LTO |


### Simulation of instrumental noise using `LISAInstrument`

The simulation is carried out within `noise_simulation.py`.

The naming convention of the output file is always
 
``` 
yyyy-mm-dd_orbits_locking_noises_4Hz.h5
``` 

#### Usage:
Input arguments for the script:

- `path-to-workdir`

Optional input arguments for the script, with flags:

- `--orbits` specify `equalarm` or `keplerian` (default: `keplerian`)
- `--locking` specify `N1-12` or `six` (default: `N1-12`)
- `--tdi` if desired, specify `1` or `2` (default: `None`): simulate and save TDI combinations for the full measurement simulation
- `--baseline` trigger baseline simulation
- `--individual` simulate and save all individual noise contributions
- `--combined` simulate and save all combined noise contributions


Execute simulation with TDI computation, specifying orbits and locking
```
python noise_simulation.py path-to-workdir --orbits equalarm --locking N1-12 --tdi 2
```             

- Execute simulation with baseline InRep configuration
```    
python noise_simulation.py path-to-workdir --orbits equalarm --locking N1-12 --tdi 2 --baseline
```        

- Execute simulation with baseline InRep configuration an save all individual noise contributions and all the combined noise contributions
```    
python noise_simulation.py path-to-workdir --orbits equalarm --locking N1-12 --tdi 2 --baseline --individual
```

### Simulation of SGWB signal using `LISAGWResponse`

The simulation is carried out within `signal_simulation.py` and `all_sky_signal_simulation.py`.

`signal_simulation.py` simulates a stochastic point source located at $\beta = \pi/2$, $\lambda = \pi/2$, while `all_sky_signal_simulation.py` simulates a Stochastic GW background with a white generator.

#### Usage:
Execute simulation with TDI computation (analogous for `all_sky_signal_simulation.py`):
```    
python signal_simulation.py path-to-workdir --orbits equalarm --tdi 2   
```