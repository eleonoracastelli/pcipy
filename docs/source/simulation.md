## Simulate LISA data: noise and signal

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


A generic simulation scripts `simulation_script.py` can be executed in the terminal and written to path `path-to-workdir` by following the structure
```shell
python simulation_script.py path-to-workdir --flags
```  
where `--flags` are used to pass various options to the script.

The flags implemented in the all the simulation scripts listed here are the following:
- `--dt (float)` setting up sampling time (s)
- `--freq1 (float)` to set up Kaiser filter first transition frequency (Hz)
- `--freq2 (float)` to set up Kaiser filter second transition frequency (Hz)
- `--tdi (int)` Pass TDI generation: choice between 1 and 2



### Simulation of noise using LISA Instrument

The simulation is carried out within `noise_simulation.py`. 

Here we simulate 3 days of data, and provide options to generate the TDI output of the simulated data and specify which  noise configurations to use.

The implemented noise configuration can be selected when executing the simulation script:
- simple simulation containing only laser, test-mass and OMS noise (default configuration, no additonal flag needed)
- baseline simulation parameters match the simulation parameters defined in 
LISA-LCST-SGS-RP-006 "End-to-end Demonstration Data Analysis Pipeline", including laser, test-mass, OMS, modulation, ranging, clock and backlink noise (additional flag `--baseline`)

with the additional option of saving each individual noise contribution by adding the flag `--individual`.

The outputs of the simulation will be:
- a full simulation dataset (laser noise + secondary noises), with filename `'yyyy-mm-dd_HHhMM_configuration_measurements_4Hz.h5'`
- a secondary noises only simulation dataset, with filename `'yyyy-mm-dd_HHhMM_configuration_noise_sec_4Hz.h5'`
- if `--individual` is selected, additional simulation datasets for each noise source, with filenames `'yyyy-mm-dd_HHhMM_configuration_individualnoise_4Hz.h5'`.

In the filenames`'configuration'` can either be `'laser_tm_oms'` or `'baseline'` and `'individualnoise'` will be wither of `'laser'`, `'test-mass'`, `'oms'`, `'ranging'`, `'backlink'`, `'clock'`, `'modulation'`

#### Usage examples:
- Execute simulation in the simple noise configuration containing only laser, test-mass and OMS noise with no TDI computation 
```shell
python noise_simulation.py path-to-workdir
```        

- Execute simulation in the simple noise configuration containing only laser, test-mass and OMS noise with TDI computation: use flag `--tdi` to specify TDI generation 
```shell
python noise_simulation.py path-to-workdir --tdi 2   
```        

- Execute simulation with baseline "End-to-end Demonstration Data Analysis Pipeline" configuration with no TDI computation 
```shell   
python noise_simulation.py path-to-workdir --baseline
```        

- Execute simulation with baseline "End-to-end Demonstration Data Analysis Pipeline" configuration with TDI computation and save all individual noise contributions
```shell   
python noise_simulation.py path-to-workdir --tdi 2 --baseline --individual
```

### Simulation of SGWB signal using LISA GW Response

The simulation is carried out within `signal_simulation.py` and `all_sky_signal_simulation.py`.

Here we simulate 3 days of data, and provide options to generate the TDI output of the simulated data with different signals depending on the chosen simulation script:
- `signal_simulation.py` simulates a stochastic point source located at $\beta = \pi/2$, $\lambda = \pi/2$
- `all_sky_signal_simulation.py` simulates a stochastic GW background with a white generator.

#### Usage:
Execute simulation with no TDI computation (analogous for `all_sky_signal_simulation.py`):
```shell
python signal_simulation.py path-to-workdir
```        

Execute simulation with TDI computation (analogous for `all_sky_signal_simulation.py`): use flag `--tdi` to specify TDI generation 
```shell   
python signal_simulation.py path-to-workdir --tdi 2   
```