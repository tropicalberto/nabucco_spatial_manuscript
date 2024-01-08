# NABUCCO Spatial
This repository contains the code to reproduce the findings from the 'Spatial relationships in the urothelial and head and neck tumor microenvironment predict response to combination immune checkpoint inhibitors' manuscript.

Title: 
Link to the pre-print: https://www.biorxiv.org/content/10.1101/2023.05.25.542236v1.abstract
Manuscript ID: NCOMMS-23-19326

This repository contains the code necessary to:
* Run the main pipeline

# Dependencies
* R dependencies: `R/install_dependencies.R`

This code has been run on R 3.6.3. The following packages were used in this study: 
●	Spatstat 1.64 28
●	Dplyr 1.0.4
●	Fitdistrplus 1.1.3
●	Patchwork 1.1.1
●	Survival 1.3.24
●	ComplexHeatmap 2.2
●	Circlize 0.4.12
●	Glmnet 4.1.1
●	RColorBrewer 1.1.2
●	Nlme 3.1.144
●	Spastat 1.7.0
●	Ggpubr 0.4.0
●	Ggrepel 0.9.1
●	Plyr 1.8.6
●	Tidyverse 1.3.0
●	Ggplot 2.3.3.3
●	Tibble 3.0.6
●	ggrastr version 1.0.1
●	pROC 1.17.0.1

# Data

The original data used in the manuscript cannot be publicly shared and must be requested as per the `Data availability` section from the manuscript upon reasonable request. For testing purposes, we included a demo dataset in `/data/`. 

The parameters obtained in the original manuscript can be found in `/data/...` (Supplementary Tables from the manuscript).

## Pipeline ##

The generation of the spatial parameters relies on three steps as detailed in the Materials and Methods section: 

1) Compute 1-NN distances and associated 1-NN distance histograms, using: ~/R/1__get1NNdistances.R

2) Estimate initial Weibull parameters, using: ~/R/2__estimateInitialParameters.R

3) Fit a Weibull distribution using a NLME model: ~/R/3__Fitnlmemodel.R

# Reproduce the results

# Contact
Feel free to reach me out (twitter: @_alberto_gil) for further questions
