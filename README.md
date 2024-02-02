This repository provides MATLAB codes for the manuscript "Mathematical modeling of temperature-induced circadian rhythms" by L. Lu, Y. Li, R. Schloss, and I.P. Androulakis.

"CosinorAnalysis_code.m" executes the cosinor analysis for the experimental data (saved in "data for cosinor analysis. xlsx") and generates phase data from the cosine approximation (saved in "data_exp_cosinor.xlsx"), which is used for model calibration and parameter estimation.

"ParamEstim_Temp_HSF1_PCG_code.m" implements parameter estimation and outputs the nominal parameter set for the model. "fminsearchcon.m" is a function file for an optimizer involved in the parameter estimation process that is not a MATLAB built-in routine.

"Temp_HSF1_PCG_code.m" includes the model equations, the generation of virtual cell and individual populations, the calculation of (the dynamics of) their properties under different temperature schedules, whether unperturbed or perturbed (alternating shift schedules), the implementation of linear regression analysis, and the production of figures.

"pop_ALL_samp3%.xlsx", "pop_PCG_samp3%.xlsx", "pop_TS_samp.xlsx", and "pop0_select.txt" record the simulated cells and individuals used to generate the results and figures in this work.
