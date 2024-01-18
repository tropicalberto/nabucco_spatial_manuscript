# Spatial relationship pipeline
This repository contains the code to reproduce the findings from the 'Spatial relationships in the urothelial and head and neck tumor microenvironment predict response to combination immune checkpoint inhibitors' manuscript.

Title: 
Link to the pre-print: https://www.biorxiv.org/content/10.1101/2023.05.25.542236v1.abstract

Manuscript ID: `NCOMMS-23-19326`

This repository contains the code necessary to run the main pipeline

# Dependencies
* R dependencies: `R/install_dependencies.R`

This code has been run on `R 3.6.3`. The following packages were used in this study: 
* Spatstat 1.64 28
* Dplyr 1.0.4
*	Fitdistrplus 1.1.3
*	Patchwork 1.1.1
*	Survival 1.3.24
*	ComplexHeatmap 2.2
*	Circlize 0.4.12
*	Glmnet 4.1.1
*	RColorBrewer 1.1.2
*	Nlme 3.1.144
*	Spastat 1.7.0
*	Ggpubr 0.4.0
*	Ggrepel 0.9.1
*	Plyr 1.8.6
*	Tidyverse 1.3.0
*	Ggplot 2.3.3.3
*	Tibble 3.0.6
*	ggrastr version 1.0.1
*	pROC 1.17.0.1

# Data

The original data used in the manuscript cannot be publicly shared and must be requested as per the `Data availability` section from the manuscript upon reasonable request. For testing purposes, we included a demo dataset in `/data/`. 

The parameters obtained in the original manuscript can be found in `/data/...` (Supplementary Tables from the manuscript).

The example dataset (ovarian dataset) provided in this package has been obtained from the `VectraPolarisData` package from `Bioconductor`
https://bioconductor.org/packages/release/data/experiment/html/VectraPolarisData.html (`DOI: 10.18129/B9.bioc.VectraPolarisData`) by:
`Wrobel J, Ghosh T (2022). VectraPolarisData: Vectra Polaris and Vectra 3 multiplex single-cell imaging data. R package version 1.0.0.`


## Pipeline ##

The generation of the spatial parameters relies on three steps as detailed in the Materials and Methods section: 

### STEP 1) Compute 1-NN distances and associated 1-NN distance histograms 

`Rscript --vanilla ./R/1__get1NNdistances.R <spatial_dataset> <areas_dataset> <output_file>`

The input is: 
* **<spatial_dataset>**: tab-separated file containing the x/y coordinates of the spatial data. Columns should be named as: *sample_id* (sample identifier), *analysisregion* (segmented tissue layer, e.g., tumor/stroma), *Xcenter* (x coordinate), *Ycenter* (y coordinate), *phenotype* (classified cell type, e.g. T-cell, cancer cell)
* **<areas_dataset>**: tab-separated file containing the tissue areas per sample. Columns must be named as: *sample_id* (sample identifier), *Total.Area* (total tissue area per sample)
* **<output_file>**: string for the output file name for step 1. Default is *./results/step1_1nn_output.tsv*

Example run: 

`Rscript --vanilla ./R/1__get1NNdistances.R data/test_spatial_data_ovarian__processed.tsv data/test_spatial_data_areas.tsv ./results/step1_1nn_output.tsv` 


### STEP 2) Estimate initial Weibull parameters

`Rscript --vanilla ./R/2__estimateInitialParameters.R <output_step_1> <output_file>`
The input is: 
* **<output_step_1>**: output file from step 1
* **<output_file>**: string for the output file name for step 2. Default is *./results/step2_weibull_initial_params.tsv*

Example run:

`Rscript --vanilla ./R/2__estimateInitialParameters.R results/step1_1nn_output.tsv ./results/step2_weibull_initial_params.tsv`

### STEP 3) Fit a Weibull distribution using a NLME model

`Rscript --vanilla ./R/3__Fitnlmemodel.R <output_step_1> <output_step_2> <output_file>`
The input is: 
* **<output_step_1>**: output file from step 1
* **<output_step_2>**: output file from step 2
* **<output_file>**: string for the output file name for step 2. Default is *./results/step3_fittednlme_weibull_params.tsv*

Example run:

`Rscript --vanilla ./R/3__Fitnlmemodel.R results/step1_1nn_output.tsv results/step2_weibull_initial_params.tsv  ./results/step3_fittednlme_weibull_params.tsv`

# Notes on the pipeline 

The non-linear mixed effect models trained in **step 3** rely on good initial parameter estimates (initial estimate of shape and scale) to fit a distributtion (Weibull distribution). 
If the pipeline gets stuck at a particular spatial relationship, we advise to inspect how the initial parameters approximate the data. We provide an example on the notebook `rmd/inspect_initial_parameters.Rmd`. 

If required, manually-set initial parameters can be parsed to step 3 so as to force the initial parameters. This ncan be done by parsing an additional argument in step 3
`Rscript --vanilla ./R/3__Fitnlmemodel.R <output_step_1> <output_step_2> <output_file> <forced_initial_params>`
The input is: 
* **<forced_initial_params>**: file with parsed initial parameters forced for the pipeline

# Contact
Feel free to reach me out (twitter: @_alberto_gil) for further questions
