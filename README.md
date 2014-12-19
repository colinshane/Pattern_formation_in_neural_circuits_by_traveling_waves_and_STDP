Pattern_formation_in_neural_circuits_by_traveling_waves_and_STDP
================================================================

Matlab scripts for journal article

The script compute_kappa.m is included with the journal article entitled "Refinement and pattern formation in neural circuits by the interaction of traveling waves with spike-timing dependent plasticity". It was used to compute kappa(x) and its Fourier transform in Eqs. 6 and 8 and to solve these equations numerically. It was last updated using Matlab R2012a.

An example call to compute_kappa is as follows:

w = compute_kappa(0,5000,500,0.02,'boxcar',0.1,inf,'asym',1,0.51,0.02,0.04,3,0,50,0.001,0.005,0.1,'uniform',0);
