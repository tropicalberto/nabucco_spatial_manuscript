# STEP 3 - Fit a Weibull function using a nonlinear-mixed fixed effect model to the 1-NN distribution curve
library(tidyverse)
library(ggpubr)
library(readxl)
library(reshape2)
library(dbscan)
library(broom)
library(zoo)
library(knitr)
library(magrittr)
library(NMF,quietly = T)
library(MASS)
library(car)
library(doParallel)
library(stats)
library(foreach)
library(spatstat)
library(raster)
library(nlme)
library(dplyr)
library(gridExtra)


spatial_data_directory <- '~/data/' # data directory
spatial_data_file      <- 'xy_data_exc_oct.tsv'
firstNNhistogram_coordinates <- 'AUCScaledSlides_300_tumorandstroma__thresholdStromaMargin_300.tsv' # output from 1__get1NNdistances.R
initial_params_file <- 'initial_params_weib_oct21.tsv'                        # output from 2__estimateInitialParameters.R



# Load initial Weibull parameters (output from 2__estimateInitialParameters.R)
initial_params <- read_delim(paste(spatial_data_directory,initial_params_file, sep=''),delim='\t')


# Calculate number of cells for each cell type / sample (will be used later for filtering purposes)
n_cells <- read_delim(paste(spatial_data_file, number_of_cells_file, sep=''), delim='\t')

n_cells <- n_cells %>%
  mutate(phenotype=gsub('CD3\\+CD8\\+','CD8+CD3+',phenotype))

xy_data_exc_oct <- read_delim(paste(spatial_data_directory, firstNNhistogram_coordinates, sep=''),delim='\t') %>%
  mutate(distance_window = WinMean) %>%
  mutate(x=distance_window, y=`N.per.mm2.scaled`) 


threshold_pixels <- 300
threshold_microns <- threshold_pixels/2

# Define maximum Shape (A_max) and maximum Scale (B_max) parameters for your fits
B_max <- 500 # max scale parameter (0,500]
A_max <- 10  # max shape parameter (0,10]

# Define Weibull and parameterized Weibull functions, and functions to convert between parameterized and
## non parameterized Weibull
## (Parameterization was done to ease the fitting procedure, to allow nlme fitting in a non-restricted manner -- without constraints)
convert_params <- function(a,b){
  A=-log(A_max/a - 1)
  A = as.numeric(A)
  B=-log(B_max/b-1)
  B = as.numeric(B)
  return(c(A=A,B=B))
}

convert_params_reverse <- function(A,B){
  a=A_max/(1+exp(-A))
  b=B_max/(1+exp(-B))
  return(c(a=a,b=b))
}
# Define Weibull distribution function in a parameterized way:
nform2 <- ~(( ( A_max/(1+exp(-A)) )  / ( B_max/(1+exp(-B)) )) * (( x/( B_max/(1+exp(-B))) )^(( A_max/(1+exp(-A)) )-1)) * exp(-1* (x/( B_max/(1+exp(-B))))^( A_max/(1+exp(-A)) )))

weibfunc2 <- function(x,A,B){
  a = A_max/(1+exp(-A))
  b = B_max/(1+exp(-B)) 
  return(((a/b) * ((x/b)^(a-1)) * exp(- (x/b)^a)))
}

nfun2 <- deriv(nform2,namevec=c("A","B"),
               function.arg=c("x","A","B"))



# Define function to fit a Weibull model to x/y coordinates of normalized 1-NN distance
## distribution coordinates, filtering out samples with n_cells < threshold_filter_cells
#' @param xy_coordinates_curves 1-NN distribution coordinate dataframe
#' @param combi Cell combination of interest: 'CellFROM_to_CellTO'
#' @param n_cells_df Dataframe with number of cells per sample (for each cell type)
#' @param threshold_filter_cells Threshold to filter out samples depending on the cell number
#' @param initial_params_df Dataframe with initial Weibull parameters
fit_data <- function(xy_coordinates_curves, combi, n_cells_df, threshold_filter_cells,
                     initial_params_df){
  message(combi)
  xy_data_filt <- xy_coordinates_curves %>% filter(phenotype_combo == combi) %>% as_tibble
  
  xy_data_filt <- xy_data_filt %>% 
    left_join(n_cells_df %>% mutate(n_from=n) %>% dplyr::select(-n), by=c('tnumber','phenotype_from'='phenotype')) %>%
    left_join(n_cells_df %>% mutate(n_to=n) %>% dplyr::select(-n), by=c('tnumber','phenotype_to'='phenotype'))
  
  threshold_cells <- threshold_filter_cells
  message(xy_data_filt %>% filter(n_from<=threshold_cells & n_to <= threshold_cells) %>%
            mutate(tnumber=paste(tnumber, 'nFROM', as.character(n_from), 'nTO',as.character(n_to))) %>% 
            pull(tnumber) %>% unique %>% paste(collapse=' / '))
  
  xy_data_filt <- xy_data_filt %>% filter(n_from>threshold_cells & n_to > threshold_cells)
  
  # Add missing observations in data and set to 0
  ## i.e. if no distances were observed beyond 100 microns, the normalized scaled counts 
  ## from 100 to 300 microns will be set to 0 
  # include 0,0 or only 1,0??
  xy_data_filt <- expand.grid(x=c(1,xy_data_filt %>%pull(distance_window) %>% unique), 
                              tnumber=xy_data_filt %>%pull(tnumber) %>% unique,
                              phenotype_combo=combi) %>%
    left_join(xy_data_filt, by=c('x','tnumber', 'phenotype_combo')) %>%
    mutate(y=ifelse(is.na(y),0,y))
  
  mean_params <- initial_params_df %>% filter(combo == combi) %>%
    filter(tnumber %in% unique(xy_data_filt$tnumber)) %>% 
    group_by(term) %>% dplyr::summarise(mean_estimate = mean(estimate))
  
  start_params <- c(a=mean_params %>% filter(term == 'shape') %>% pull(mean_estimate),
                    b=mean_params %>% filter(term == 'scale') %>% pull(mean_estimate))
  
  if(combi == 'PanCK+_to_PanCK+'){
    start_params <- c(a=3, b=10) # rewrite exception imsicion for pck-pck
  }
  
  fit52_oct <- tryCatch({
    nlme(y~nfun2(x,A,B),
         fixed=A+B~1,
         random=list(tnumber=pdDiag(A+B~1)), # or pdDiag, pdLogChol
         data=xy_data_filt %>% mutate(tnumber = factor(tnumber)),
         start=convert_params(start_params['a'], start_params['b']),
         method='REML', control = nlmeControl(maxIter = 1000, 
                                              returnObject = F,
                                              msMaxIter=200,
                                              pnlsTol=1e-1, 
                                              tolerance=1e-6, opt="nlm"))
  },
  error=function(cond){
    message(paste('failed',combi))
    return(combi)
  })
  # }
  return(list(model=fit52_oct, xy_data=xy_data_filt))
  
}



# Implement pipeline:
# * First try to fit models for ALL data
# * If for a particular phenotype combo it doesn't work, fit:
# ** Fit all data with n_cells > 20
# ** If not, increase the threshold from 20 to 100
# 
# This aided the fit of the following phenotype combinations: 
# * failed CD3+FoxP3+_to_CD3+FoxP3+
# * failed CD20+_to_CD20+
# * failed CD20+_to_CD3+FoxP3+
# * failed CD3+FoxP3+_to_CD20+

# Define phenotype combinations to study in data
combinations_to_study <- unique(xy_data_exc_oct %>% pull(phenotype_combo))

success_models <- list() # Fitted models

# Dataframe to store fitted Weibull coefficients
weib_coefs <- data.frame(matrix(nrow=0,ncol=4))
colnames(weib_coefs) <- c('a','b','sample','phenotype_combo')
weib_coefs$a %<>% as.numeric()  # shape
weib_coefs$b %<>% as.numeric()  # scale
weib_coefs$sample %<>% as.character()
weib_coefs$phenotype_combo %<>% as.character()

# Function to store x/y coordinates used for modelling
xy_data_weib_imcision <- xy_data_exc_oct %>% slice(0)
xy_data_weib_imcision <- xy_data_weib_imcision %>% add_column(yhat=as.numeric())

# Function to iterate through pheno combos
subset_combinations <- combinations_to_study
# Iterate through different number of threshold minimum number of cells sample
## (in some cases, all samples word, in other strenghter filter need to be 
## applied to samples with low number of cells)
for(threshold_cells in c(0,20,50,70,100)){
  model_objects <- list()
  if(length(subset_combinations) > 0){
    for(combi in subset_combinations){
      message(threshold_cells)
      # Fit Weibull distribution
      ## If fit succeds, a nlme model is returned, otherwise a character
      ## with the name of the failed phenotype combination is returned
      pipeline_run <- fit_data(xy_data_exc_oct, combi, n_cells, 
                               threshold_cells, initial_params)
      model_objects <- append(model_objects,
                              list(combi=pipeline_run[['model']]))
      fitted_model <- tail(model_objects, 1)[[1]] # retrieve fitted model in the for loop
      # If the model successfully fitted, append estimated parameters of the model                        
      if(class(fitted_model)  != 'character'){
        # Weibull parameters:
        new_data_weib_param <- as_tibble(ranef(fitted_model)) %>%
          mutate(sample=rownames(ranef(fitted_model))) %>%
          mutate(A=A + fixef(fitted_model)['A'],
                 B=B+ fixef(fitted_model)['B']) %>%
          mutate(phenotype_combo=combi) %>%
          mutate(a=10/(1+exp(-A)),b=500/(1+exp(-B)))
        weib_coefs <- bind_rows(weib_coefs, new_data_weib_param)
        
        # x/y coordinates of input / output model: 
        xy_data_weib_imcision <- bind_rows(xy_data_weib_imcision,
                                           pipeline_run[['xy_data']] %>% mutate(yhat=fitted(fitted_model)))
      }
    }
  }
  names(model_objects) <- subset_combinations
  # Update subset combinations. Include only phenotype combinations whose modelling failed
  subset_combinations <- names(model_objects)[lapply(model_objects, class) == 'character']
  
  # Update list of succesfully fitted models
  success_models <- append(success_models,
                                    model_objects[lapply(model_objects, class) != 'character'])
  
}


# Process parameters (extract weibull paramters from fixed/random effects)
final_weibull_parameters <- data.frame()
for(combi in names(success_models)){
  fixefect <- fixef(success_models[[combi]])
  ranef <- ranef(success_models[[combi]])
  
  
  
  final_weibull_parameters <- final_weibull_parameters %>% 
    bind_rows(  ranef %>%
                  mutate(tnumber= rownames(ranef),
                         phenotype_combo=combi) %>% 
                  mutate(A = A + fixefect[['A']],
                         B = B + fixefect[['B']]) %>%
                  mutate(a=10/(1+exp(-A)),b=500/(1+exp(-B))))
  
}

# Save final Weibull parameters
write_delim(final_weibull_parameters, paste(spatial_data_directory, 'global_weib_coeff_oct.tsv', sep=''))
