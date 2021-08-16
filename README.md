# Stress-variations-within-fault-damage-zone

This package is to calculate the variations of stress and effective properties in fractured rock mass representative of fault damage zone. Please refer to Zhang and Ma (2021) for more details.


Fault damage zone is idealized as a 2D multilayer model including fractal macroscopic fractures with different densities through layers.


This package mainly includes two steps:

1) Generation of discrete fracture network (DFN), with fracture length following double power law. In the first step, MATLAB scripts "2D_DFN_Kim(2007).m" provided by TAE HYUNG KIM (2007) are used to generate six sets of scale-independent DFNs with six levels of fractue density.

Here are relevant parameters in the generation:
    Fractal dimension D = 2; 
    Power law length exponent a = 3;
    Side width of each representative volume L = 10 m;
    Density-related term alpha = 0.32; 0.8; 1.6; 2.4; 3.2; 4.0


2) Calculation of fracture deformation and thus stress changes according to the specific boundary conditions, given applied far-field stress. In this part, "Main_Func_Stress_Variations.m" is used for simulations. Please go to this file and relevant functions for more details.



References:
Kim, T. H. (2007). Fracture characterization and estimation of fracture porosity of naturally fractured reservoirs with no matrix porosity using stochastic fractal models. Texas A&M University.

Zhang, S., & Ma, X. (2021).How does in situ stress rotate within a fault zone? Insights from explicit modeling of the frictional, fractured rock mass. Journal of Geophysical Research: Solid Earth.
