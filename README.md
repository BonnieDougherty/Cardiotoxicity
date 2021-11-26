# Identifying metabolic adaptations as a characteristic of cardiotoxicity using paired transcriptomics and metabolomics data integrated with a model of heart metabolism

### Abstract
Improvements in the diagnosis and treatment of cancer has revealed long-term side effects of chemotherapeutics, particularly cardiotoxicity. Here, we present paired transcriptomics and metabolomics data characterizing in vitro cardiotoxicity to three compounds: 5-fluorouracil, acetaminophen, and doxorubicin. Standard gene enrichment and metabolomics approaches identify some commonly affected pathways and metabolites but are not able to readily identify metabolic adaptations in response to cardiotoxicity. The paired data was integrated with a genome-scale metabolic network reconstruction of the heart to identify shifted metabolic functions, unique metabolic reactions, and changes in flux in metabolic reactions in response to these compounds. Using this approach, we confirm known mechanisms of doxorubicin-induced cardiotoxicity and provide hypotheses for metabolic adaptations in cardiotoxicity for 5-fluorouracil, doxorubicin, and acetaminophen.


### Overview

Analyses were run with MATLAB and the COBRA toolbox v3 (accessed 2019-02-18), R and Python. 

	project
	|- README             
  	|
	|- code/              # code used for the presented analysis
 	| |- MATLAB/          # building rat heart-specific model and running TIDEs
	| | |- functions/     # functions necessary for reproducing results
	| |- Python/          # python notebooks for reproducing RIPTIDE integration
 	| |- R/               # scripts for running RNA-seq, metabolomics, dose response, and Seahorse analysis with figure code
 	|
 	|- results/           # figures and supplementary tables
 	| |- figures/
	| |- supplement/
	
	|- data/           # figures and supplementary tables
 	| |- experimental/    # raw Seahorse and dose response data
	| |- metabolomics/    # raw and processed metabolomics data
	| |- RIPTIDE/         # RIPTIDE results for each condition
	| |- RNA-seq/         # DEGs used as input for TIDEs and TPMs for RIPTIDE
	| |- TIDEs/           # output of TIDEs for each condition
