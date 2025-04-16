## Simulate LISA data: noise and signal

PCI is applied directly to the LISA phasemeter measurements.

The LISA phasemeter measurements are the output of simulations built via the [LISA Simulation suite](https://gitlab.in2p3.fr/lisa-simulation). 

In what follows we use the following packages from the September LISA Simulation Suite:
- `lisainstrument >= 1.8` for the overarching instrument simulation,
- `lisaorbits` for the orbits,
- `lisagwresponse` for the response to GW,
- `pytdi` for the evaluation of TDI variables.

In the package we provide three simulation scripts that can be executed in the terminal:

- `noise_simulation.py` to simulate the laser and secondary noises
- `signal_simulation.py` to simulate a stochastic point source in the sky
- `all_sky_signal_simulation.py` to simualte a stochastic graviational wave background.


### Simulation scripts
In the package we provide three simulation scripts:

- `noise_simulation.py` to simulate the instrumental noise for various simulation scenarios;
- `signal_simulation.py` to simulate a stochastic point source in the sky;
- `all_sky_signal_simulation.py` to simulate a stochastic gravitational wave background.

A generic simulation scripts `simulation_script.py` can be executed in the terminal and written to path `path-to-workdir` by following the structure
```shell
python simulation_script.py path-to-workdir --flags
```  
where `--flags` are used to pass various options to the script.

The flags implemented in the all the simulation scripts listed here are the following:
- `--dt (float)` setting up sampling time [s] (default: 0.25)
- `--orbits (str)` specify `equalarm` or `keplerian` (default: `keplerian`)
- `--tdi` if desired, specify `1` or `2` (default: `None`): simulate and save TDI combinations for the simulation

#### Simulation scenarios
To allow for benchmarking performance of PCI against TDI, we provide a list of simulation scenarios to test. Various options can be toggled when executing the script to generate the data.

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


#### Simulation of instrumental noise using `LISAInstrument`

The simulation is carried out within `noise_simulation.py`. 

Here we simulate 3 days of data, and provide options to generate the TDI output of the simulated data and specify which simulation scenarios to produce.


The flags implemented in the noise simulation script are the following:
- `--dt (float)` setting up sampling time [s] (default: 0.25)
- `--freq1 (float)` to set up Kaiser filter first transition frequency [Hz] (default: 1.1)
- `--freq2 (float)` to set up Kaiser filter second transition frequency [Hz] (default: 2.9)
- `--orbits (str)` specify `equalarm` or `keplerian` (default: `keplerian`)
- `--locking (str)` specify `N1-12` or `six` (default: `N1-12`)
- `--tdi` if desired, specify `1` or `2` (default: `None`): simulate and save TDI combinations for the full measurement simulation
- `--baseline` trigger baseline simulation to match the simulation parameters defined in LISA-LCST-SGS-RP-006 "End-to-end Demonstration Data Analysis Pipeline", including laser, test-mass, OMS, modulation, ranging, clock and backlink noise 
- `--individual` simulate and save all individual noise contributions
- `--combined` simulate and save all combined noise contributions

The naming convention of the output file is always:
```shell
yyyy-mm-dd_orbits_locking_scenario_data_4Hz.h5
``` 

The outputs of the simulation will be:
- a full simulation dataset (laser noise + all secondary noises), for the selected scenario (`laser_tm_oms` or `baseline`)with filename `'yyyy-mm-dd_orbits_locking_scenario_measurements_4Hz.h5'`
- a secondary noises only simulation dataset, with filename `'yyyy-mm-dd_orbits_locking_scenario_noise_sec_4Hz.h5'`
- if `--individual` is selected, additional simulation datasets for each noise source, with filenames `'yyyy-mm-dd_orbits_locking_scenario_individual_noisename_4Hz.h5'`.
- if `--combined` is selected, additional combined simulation datasets for each noise source, with filenames `'yyyy-mm-dd_orbits_locking_scenario_combined_noisename_4Hz.h5'`.


##### Usage examples:

- Execute simulation in the simple noise configuration containing only laser, test-mass and OMS noise, with equal arm orbits, default laser locking configuration and no TDI computation 
```shell
python noise_simulation.py path-to-workdir --orbits equalarm --locking N1-12
```   

Execute simulation with TDI computation, in the simple noise configuration containing only laser, test-mass and OMS noise, with keplerian orbits, default laser locking configuration and TDI computation 
```shell
python noise_simulation.py path-to-workdir --orbits keplerian --locking N1-12 --tdi 2
```             

- Execute simulation with baseline "End-to-end Demonstration Data Analysis Pipeline" configuration with no TDI computation 
```shell
python noise_simulation.py path-to-workdir --orbits equalarm --locking N1-12 --baseline
```        

- Execute simulation with baseline "End-to-end Demonstration Data Analysis Pipeline" configuration with TDI computation and save all individual noise contributions and all combined noise contributions
```shell  
python noise_simulation.py path-to-workdir --orbits equalarm --locking N1-12 --tdi 2 --baseline --individual
```

### Simulation of SGWB signal using `LISAGWResponse`

The simulation is carried out within `signal_simulation.py` and `all_sky_signal_simulation.py`.

Here we simulate 3 days of data, and provide options to generate the TDI output of the simulated data with different signals depending on the chosen simulation script:
- `signal_simulation.py` simulates a stochastic point source located at $\beta = \pi/2$, $\lambda = \pi/2$
- `all_sky_signal_simulation.py` simulates a stochastic GW background with a white generator.

#### Usage example:

Execute simulation with TDI computation (analogous for `all_sky_signal_simulation.py`):
```    
python signal_simulation.py path-to-workdir --orbits equalarm --tdi 2   
```        
