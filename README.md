# Metafoodweb
Simulation code for the publication "Landscape heterogeneity buffers biodiversity of simulated meta-food-webs under global change through rescue and drainage effects"



 	
This repository contains simulation code of a research project. It applies allometric models of animal foraging and network structure extended in a spatial network to simulate the impact of fragmentation and eutrophication on species diversity patterns in a complex landscape.


File List in Code_S.zip

	Code
		FullWeb
			main.cpp
			pdef_dynamics_1.1_space.h
			submit.sh
			summarize.sh
			wrapper.sh
			gcc.sh
		2-patch
			Heatmap.sh
		1-patch
			runlocal.sh	

	read-in
		Nutrients
		Landscapes
		Foodweb
		1/2-patch
	
	R
		FullWeb
			Fullweb.r
		2-patch
			2-patch.r
		1-patch
			nutrient.r
			emigration.r
			Plots.r	


File descriptions

		main.cpp: Main simulation code producing the data for the analysis. The code is an extension of the code used for the article Ryser et al. (2019) ‚The biggest Losers: Habitat fragmentation deconstructs complex food webs from top to bottom, ProcB. A detailed description of parameters and equation used in the model is given in the supplementary file of the manuscript. 

		pdef_dynamics_1.1_space.h: List of libraries needed for simulation run and a list of all functions in the main simulation code. 

		
Note: Working directory has to be adjusted. Input for parallelization of simulation runs on a HPC is adapted for a bash script. Code needs to be compiled beforehand (example for local compilation in gcc.sh).

Libraries:
gsl/1.16-3
sundials/2.7.0-1
gcc/5.2.0-1


