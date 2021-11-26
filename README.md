# Identifying metabolic adaptations as a characteristic of cardiotoxicity using paired transcriptomics and metabolomics data integrated with a model of heart metabolism

### Abstract
Improvements in the diagnosis and treatment of cancer has revealed long-term side effects of chemotherapeutics, particularly cardiotoxicity. Here, we present paired transcriptomics and metabolomics data characterizing in vitro cardiotoxicity to three compounds: 5-fluorouracil, acetaminophen, and doxorubicin. Standard gene enrichment and metabolomics approaches identify some commonly affected pathways and metabolites but are not able to readily identify metabolic adaptations in response to cardiotoxicity. The paired data was integrated with a genome-scale metabolic network reconstruction of the heart to identify shifted metabolic functions, unique metabolic reactions, and changes in flux in metabolic reactions in response to these compounds. Using this approach, we confirm known mechanisms of doxorubicin-induced cardiotoxicity and provide hypotheses for metabolic adaptations in cardiotoxicity for 5-fluorouracil, doxorubicin, and acetaminophen.


### Overview

Analyses were run with MATLAB and the COBRA toolbox v3 (accessed 2019-02-18), R and Python. 

	project
	|- README             
  |
  |- code/              # code used for the presented analysis
  | |- Python/ 
 	| |- MATLAB/          	
	| | |- data/          # necessary models, metabolic tasks, and final results
	| | |- functions/     # functions necessary for reproducing results
	| | |- scripts/       # code for reproducing building draft iCardio models, model curation, and TIDEs analysis
 	| |- R/               # code for microarray DEG analysis
	| | |- data/          # necessary external files to run analyses
	| | |- scripts/       # code for reproducing analyses and figures
 	|
 	|- results/           # figures and supplementary tables
 	| |- figures/
	| |- tables/


### Re-running the TIDEs pipeline for your own data

To run the TIDEs pipeline, you must first run the generateMinRxnList() function to generate a MATLAB structure that contains an entry for each metabolic task, including the reactions necessary for that task. For the analysis presented in the paper, reactions without GPRs were removed (removeNoGPR = 'true'). Next, run the calculateTIDEscores() function with the model, minRxnList, and data for your study. The function will return a structure with a task score for each task in the minRxnList structure with the calculated significance for the task compared to randomly shuffled data. 
