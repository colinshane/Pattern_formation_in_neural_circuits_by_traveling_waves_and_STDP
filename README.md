Pattern_formation_in_neural_circuits_by_traveling_waves_and_STDP
================================================================

Matlab scripts for journal article

The script num_int_w.m is included with the journal article entitled "Refinement and pattern formation in neural circuits by the interaction of traveling waves with spike-timing dependent plasticity". It was used to compute kappa(x) and its Fourier transform in Eq. 7 and to solve Eqs. 6 and 8 numerically. The script was last updated using Matlab R2012a.

Example calls to num_int_w:

1) To numerically integrate the synaptic weights as in Eq. 8 and Figs. 2D-E: 
w = num_int_w(0,5000,500,0.02,'boxcar',0.1,inf,'asym',1,0.51,0.02,0.04,3,0,50,0.001,0.005,0.1,'uniform',0);

2) To return the computed kappa(x) and FT{kappa(x)}, as in Eq. 7 and Fig. 2B
w = num_int_w(1,5000,500,0.02,'boxcar',0.1,inf,'asym',1,0.51,0.02,0.04,3,0,50,0.001,0.005,0.1,'uniform',0);
