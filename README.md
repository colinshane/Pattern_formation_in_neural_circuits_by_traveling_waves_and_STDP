Pattern_formation_in_neural_circuits_by_traveling_waves_and_STDP
================================================================

Matlab scripts for journal article: "Refinement and pattern formation in neural circuits by the interaction of traveling waves with spike-timing dependent plasticity".

***********
num_int_w.m
***********
This script is used to numerically integrate w(x) in Eqs. 6 and 8, with non-linearities applied such as hard bounds to the synaptic strengths and alternating wave directions. It can also be used to compute kappa(x) in Eq. 7, the functions that comprise it -- K(x/v), alpha(x/v), epsilon(x/v) -- and the Fourier transforms of these functions, without integrating w(x) over time. 

Example calls to num_int_w:

1) To return the computed kappa(x) and FT{kappa(x)}, as in Eq. 7 and Fig. 2B
w = num_int_w(1,5000,500,0.02,'boxcar',0.1,inf,'asym',1,0.51,0.02,0.04,3,0,50,0.001,0.005,0.1,'uniform',0);

2) To numerically integrate the synaptic weights as in Eq. 8 and Figs. 2D-E: 
w = num_int_w(0,5000,500,0.02,'boxcar',0.1,inf,'asym',1,0.51,0.02,0.04,3,0,50,0.001,0.005,0.1,'uniform',0);

**************
get_spktimes.m
**************
Read spike times from a .nd file into a Matlab cell array. Details of the .nd format are provided at http://www.imodel.org/nd/. Because these .nd files are so large, it is best to read spike times from a limited period of the simulation, as in the third example below.

Example calls to get_spktimes:

1) To read in all spike times from the 64 x 64 inputs:

times = get_spktimes('filename.nd');

2) To read in the first 5000s worth of spike times (note the input layer dimensions, 64-by-64-by-1, bust be included if specifying limits to the spike times):

times = get_spktimes('filename.nd',[64 64 1 5000000]);

2) To read in 500s worth of spike times, starting from 500001ms and ending at 1000000ms into the simulated wave patterns:

times = get_spktimes('filename.nd',[64 64 1 1000000 500001]); 
