#' @author Alberto Gil-Jimenez
#' Manuscript: Spatial relationships in the urothelial and head and neck tumor microenvironment predict response to combination immune checkpoint inhibitors


#' STEP 1 
#' 
#' Script to estimate initial Weibull parameters for the 1-NN distribution
#' 
#' Output is saved as a tab-delimited file (specified by output_file)
#'
#' @param spatial_data input spatial dataset, with at least columns: 
#'        "sample_id', "analysisregion", "Xcenter", "Ycenter", "phenotype"
#' @param spatial_areas input tissue areas for each sample with at least columns: 
#'        "sample_id', "Total.Area"
#' @param output_file output file name 
#' @examples
#' # Example code
#' Rscript --vanilla ./R/2__estimateInitialParameters.R results/step1_1nn_output.tsv ./results/step2_weibull_initial_params.tsv
#' @export

# Load packages
library(ComplexHeatmap)
library(tidyverse)
library(plyr)
library(cold) #function fixeff()
library(lme4)
library(glmmTMB)
library(fitdistrplus)
library(tidyr)
library(cold)
library(broom)
library(ggpubr) # stat_compare_means


# Get input system arguments
args = commandArgs(trailingOnly=TRUE)

# Double-check there is at least one argument (input file name)
if (length(args)==0) {
  stop("You must parse an input file", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  message("Output file not specified. Using default: ~/results/step2_weibull_initial_params.tsv")
  args[2] = "./results/step2_weibull_initial_params.tsv"
}
firstNNhistogram_coordinates      <- args[1]  # it's the output from 1__get1NNdistances.R
output_file                       <- args[2]

# Load 1-NN X/Y histogram coordinates dataframe (output from script 1__get1NNdistances.R)
all_distances_data <- read_delim(firstNNhistogram_coordinates, delim='\t',
                                 col_types=list(tumor_stroma=col_character()))

all_distances_data <- all_distances_data %>% as_tibble %>% 
  mutate(distance_window = WinMean) %>% 
  mutate(phenotype_combo = paste(phenotype_from, phenotype_to, sep='_to_')) %>% 
  dplyr::select(sample_id, phenotype_combo, `N.per.mm2.scaled`, distance_window) %>%
  mutate(new = `N.per.mm2.scaled` * 1000) %>%
  mutate(new =round(new))

# Recreate 1-NN histogram from coordinates (this is required for function "fitdistrplus")
all_combos_dists <- tibble()
for(tn in unique(all_distances_data %>% pull(sample_id))){
  for(combo in unique(all_distances_data %>% pull(phenotype_combo))){
    all_dists_tnum <- list()
    # dists <- all_distances_data %>% filter(sample_id == tn) %>% filter(phenotype_combo == combo) %>% pull(distance_window)
    # times <- all_distances_data %>% filter(sample_id == tn) %>% filter(phenotype_combo == combo) %>% pull(new)
    # 
    dists <- c(1,all_distances_data %>% filter(sample_id == tn) %>% filter(phenotype_combo == combo) %>% pull(distance_window) )
    times <- c(1, all_distances_data %>% filter(sample_id == tn) %>% filter(phenotype_combo == combo) %>% pull(new))
    
    if(length(times) == 2){
      dists <- c(1:298)
      times <- rep(0,length(c(1:298)))
    }
    for(i in 1:length(times)){
      all_dists_tnum <- c(all_dists_tnum, rep(dists[i], times[i]))
    } 
    aldists <- unlist(all_dists_tnum %>% as.numeric)
    aldistsdf <- as_tibble(data.frame(dists = aldists, 
                                      sample_id = rep(tn, length(aldists)),
                                      combo = rep(combo, length(aldists))))
    all_combos_dists <- bind_rows(all_combos_dists, aldistsdf)
  }
  
}


# Fit a Weibull distribution by Maximum likelihood MLE [this is an initial estimation that will be optimized in step 3)]
initial_params <- data.frame(matrix(ncol=5, nrow=0))
colnames(initial_params) <- c('term','estimate','std.error', 'sample_id','combo')
initial_params <- initial_params %>% 
  mutate(term=as.character(term), estimate=as.numeric(estimate), `std.error`=as.numeric(`std.error`),sample_id=as.character(sample_id),
         combo=as.character(combo))

for(combo in unique(all_distances_data %>%  pull(phenotype_combo))){
  # for different "sample_id"
  message(combo)
  
  for(tn in unique(all_combos_dists  %>% pull(sample_id))){
    # take into account only 1-NN dists < 100, because of the precision of fitdists (no of bins in histogram)
    # EXCEPT for From/To PanCK or negative (special case)
    # the parameters will be further optimized later
    if(str_detect(combo, 'PanCK\\+_') | str_detect(combo, 'negative_') | str_detect(combo, 'Cancer_') | str_detect(combo, 'Negative_')){
      params <- fitdist(all_combos_dists %>% 
                          filter(combo == combo) %>% filter(sample_id == tn)  %>%
                          pull(dists) ,"weibull") %>% summary 
    } else{
      params <- fitdist(all_combos_dists %>% 
                          filter(combo == combo) %>% filter(sample_id == tn) %>% filter(dists < 100) %>%
                          pull(dists) ,"weibull") %>% summary 
    }
    class(params) <- c("fitdist", "fitdistr") 
    params <- broom::tidy(params)
    params <- params %>% mutate(sample_id=rep(tn,nrow(params))) %>%
      mutate(combo=rep(combo,nrow(params)))
    
    initial_params <- bind_rows(initial_params, params)
  }
}

# Return a vector with Weibull distribution parameters estimated by maximum likelihood (will be called as 'initial parameters')
# write_delim(initial_params_weib, '~/nabucco_spatial/data/IMCISION/spatial_data/initial_params_weib_imcision.tsv',delim='\t')

write_delim(initial_params, 
              output_file,
              delim='\t')