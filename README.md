# WorkingFromSelfDrivingCar
Collection of GAMS and R-codes and results for the paper Working from Self-Driving Cars (version May 2023) by Georg Hirte, Renee Laes, Regine Gerike (all TU Dresden)

File WFCGenerateParamMC.Rmd is the R-Markdown file that is used to calculate all parameter distributions needed 
	in the Monte Carlo (MC) simulation of the model on working form self-driving cars. 
	The paramters are stored in the files MC_paramDE.gdx and MC_paramUS.gdx.
	
File WFCTRAAllParam202304306.gms is a GAMS file. 
	It includes the simulation code for the impact of single parameter changes on WFC and the Monte Carlo simulations. 
	The GAMS program uses MC_paramDE.gdx and MC_paramUS.gdx as inputs and provides a number of output files OUtMobileWOrk_....gdx. 
	E.g. OutMobileWork_tauw.gdx includes the results form the variation in the wage tax rate tauw,
	OutMobileWork_Benchmark.gdx is the outcome from the benchmark parameter choice, and OutMobileWork_MonteCarlo.gdx is the outcome
	of the Monte Carlo simulation.
	
File WFCTRAResults20230430.RMD is the R-Markdown file to finalize the results form the GAMS simulations, 
	make regressions on the MC outcomes, and produce tables and results used in the paper. Inputs are the outcome *.gdx files.
	
We also added the knitted R-Markdown .pdf files for both R-Markdown files.

Further, GAMS-simulation results are printed as *.gdx files. These are read as inputs in the R-files. *.txt files are results from R calculations.

To reproduce: 1. Run  WFCGenerateParamMC.Rmd, 2. Run WFCTRAAllParam202304306.gms (choose benchmark, next parameters, eventually monte-carlo), 3. Run WFCTRAResults20230430.RMD
