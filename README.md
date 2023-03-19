# noise-influences-stimulation-effects

Codes to run simulations of the manuscript "Noise improves the association between effects of local stimulation and structural degree of brain networks"

### Language: matlab

## Contents:

Pen201node_par_Euler_HPC1.m: main script used to simulate bain dynamics under different combinations of stimulation sites and noise amplitudes for 30 realizations. Matlab Parallel computing is utilized to increase the speed of computation. The main parameters in this script are G (group-level structural connectivity matrix), D (group-level distance matrix), time (total simulation time), dt (temporal resolution for simulation), c5 (global coupling strength), stim_P (external stimulation), sigma (noise amplitude), res (final temporal resolution), rs (random number stream)

wc_coupled_stochastic1_sd5.m: matlab function, Wilson-Cowan neual mass model to simulate brain activity under a specific condition.

WC_downsampling.m: matlab function, downsampling the resulting neural activity.

parsave_Euler.m: matlab function, saving the data in the parfor loop.



