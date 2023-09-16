# InfoTheory-NFkB
The Matlab scripts contained in this repo were utilized to generate the insights presented in the manuscript "Information-theoretic analysis of a model of CAR-4-1BB-mediated NFκB activation" by V. Tseruyan and S. D. Finley (BioRxiv preprint: https://doi.org/10.1101/2023.06.09.544433). The usage of the various script files is as follows:


* SGBN_diagram.nwt/
	A .nwt file containing the SGBN diagram for Supplementary Figure 2.

* ChannelCapacity/noisySimulations/noisySimulatorCodes
 	This folder contains the scripts necessary to reproduce estimates of channel capacity (Fig. 2a). 
	NFkB_model.m — the basic ODE model generated via BioNetGen
	noisySimulator.m — a wrapper function around NFkB_model.m for generating results with stochastically sampled protein concentrations upstream of NFkB, needs the input of extracellular antigen concentration (ant_val) and noise level (sigma)
	noisySimulator_extend.m — a similar wrapper function which also randomizes the cellular NFkB pool.
	last_batch_sims.m — a master script wrapped around noisySimulator.m and noisySimulator_extend.m, which allows to set for different extracellular antigen concentrations and noise levels.


* ./ChannelCapacity/noisySimulations/noisySimulatorresults
 	This folder contains the scripts necessary to reproduce estimates of channel capacity (Fig. 2a). 
	mi_cont_cont.m — a script for computing, an estimate of the mutual information between samples of two continuous distributions
	mi_discrete_cont.m — a script for computing, an estimate of the mutual information between samples of two distributions, one continuous and one discrete
	rangesearch_var_radius.m — help her function for mi_cont_cont.m and mi_discrete_cont.m, computes the number of points in a neighborhood of a given radius
	noisySimulator.m — a wrapper function around NFkB_model.m for generating results with stochastically sampled protein concentrations upstream of NFkB, needs the input of extracellular antigen concentration (ant_val) and noise level (sigma)
	noisySimulator_extend.m — a similar wrapper function which also randomizes the cellular NFkB pool.
	last_batch_sims.m — a master script wrapped around noisySimulator.m and noisySimulator_extend.m, which allows to set for different extracellular antigen concentrations and noise levels.
	mutinfo_analyzer.m and mutinfo_analyzer_extend.m — a wrapper around mi_discrete_cont.m, which computes mutual information given simulation results from ../noisySimulatorCodes


* ./ChannelCapacity_noinact/noisySimulations/noisySimulatorCodes
	This folder contains the scripts necessary to reproduce estimates of channel capacity with disabled IKKbeta deactivation (Fig. 2b). 
	Identical to ChannelCapacity/noisySimulations/noisySimulatorCodes, except the self-inactivation parameter for IKKbeta has been set to zero.


* ./ChannelCapacity_noinact/noisySimulations/noisySimulatorresults
	This folder contains the scripts necessary to reproduce estimates of channel capacity with disabled IKKbeta deactivation (Fig. 2b). 
	Identical to ChannelCapacity/noisySimulations/noisySimulatorCodes, except the self-inactivation parameter for IKKbeta has been set to zero.

* ./Mutinfo_diff_distr
	This folder contains the scripts for evaluating the utilization of channel capacity by different empirically-motivated distributions of antigen concentration (Fig. 4a and 5a). 
	./MutInfo_02_extend — The analysis performed for the distribution with 20% antigen-positive targets and 80% antigen-negative targets
		mi_cont_cont.m — a script for computing, an estimate of the mutual information between samples of two continuous distributions
		mi_discrete_cont.m — a script for computing, an estimate of the mutual information between samples of two distributions, one continuous and one discrete
		rangesearch_var_radius.m — help her function for mi_cont_cont.m and mi_discrete_cont.m, computes the number of points in a neighborhood of a given radius
		NFkB_model.m — the basic ODE model generated via BioNetGen
		mutinfo_analyzer.m — a wrapper around mi_discrete_cont.m, which computes mutual information given simulation results
		wrapper_X.m — a wrapper around NFkB_model.m for generating NFkB activation data for the case where the CAR cell is stimulated by an antigen distribution outlined above and the noise level of intracellular protein concentrations is equal to X (=0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0).
	./MutInfo_05_extend — Similar analysis to above but performed for the distribution with 50% antigen-positive targets and 50% antigen-negative targets
	./MutInfo_08_extend — Similar analysis to above but performed for the distribution with 80% antigen-positive targets and 20% antigen-negative targets
	./MutInfo_10_extend — Similar analysis to above but performed for the distribution with 100% antigen-positive targets


* ./Mutinfo_diff_distr_noinact
	This folder contains the scripts for evaluating the utilization of channel capacity by different empirically-motivated distributions of antigen concentration for the case where self-deactivation of IKKbeta is set to zero  (Fig. 4b and 5b). Contents follow the same structure as ./Mutinfo_diff_distr



		
