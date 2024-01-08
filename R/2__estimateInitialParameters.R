# STEP 2 - Script to estimate initial Weibull parameters for the 1-NN distribution
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


spatial_data_directory <- '~/data/' # data directory
spatial_data_file      <- 'xy_data_exc_oct.tsv'
firstNNhistogram_coordinates <- 'AUCScaledSlides_300_tumorandstroma__thresholdStromaMargin_300.tsv' # output from 1__get1NNdistances.R


# Load 1-NN X/Y histogram coordinates dataframe (output from script 1__get1NNdistances.R)
all_distances_data <- read_delim(paste(spatial_data_directory, firstNNhistogram_coordinates, sep=''), delim='\t',
                                 col_types=list(tumor_stroma=col_character()))


all_distances_data <- all_distances_data %>% as_tibble %>% 
  mutate(distance_window = WinMean) %>% 
  mutate(phenotype_combo = paste(phenotype_from, phenotype_to, sep='_to_')) %>% 
  dplyr::select(tnumber, phenotype_combo, `N.per.mm2.scaled`, distance_window) %>%
  mutate(new = `N.per.mm2.scaled` * 1000) %>%
  mutate(new =round(new))

# Recreate 1-NN histogram from coordinates (this is required for function "fitdistrplus")
all_combos_dists <- tibble()
for(tn in unique(all_distances_data %>% pull(tnumber))){
  for(combo in unique(all_distances_data %>% pull(phenotype_combo))){
    all_dists_tnum <- list()
    # dists <- all_distances_data %>% filter(tnumber == tn) %>% filter(phenotype_combo == combo) %>% pull(distance_window)
    # times <- all_distances_data %>% filter(tnumber == tn) %>% filter(phenotype_combo == combo) %>% pull(new)
    # 
    dists <- c(1,all_distances_data %>% filter(tnumber == tn) %>% filter(phenotype_combo == combo) %>% pull(distance_window) )
    times <- c(1, all_distances_data %>% filter(tnumber == tn) %>% filter(phenotype_combo == combo) %>% pull(new))
    
    if(length(times) == 0){
      dists <- c(1:298)
      times <- rep(0,length(c(1:298)))
    }
    for(i in 1:length(times)){
      all_dists_tnum <- c(all_dists_tnum, rep(dists[i], times[i]))
    } 
    aldists <- unlist(all_dists_tnum %>% as.numeric)
    aldistsdf <- as_tibble(data.frame(dists = aldists, 
                                      tnumber = rep(tn, length(aldists)),
                                      combo = rep(combo, length(aldists))))
    all_combos_dists <- bind_rows(all_combos_dists, aldistsdf)
  }
  
}


# Fit a Weibull distribution by Maximum likelihood MLE [this is an initial estimation that will be optimized in step 3)]
initial_params <- data.frame(matrix(ncol=5, nrow=0))
colnames(initial_params) <- c('term','estimate','std.error', 'tnumber','combo')
initial_params <- initial_params %>% 
  mutate(term=as.character(term), estimate=as.numeric(estimate), `std.error`=as.numeric(`std.error`),tnumber=as.character(tnumber),
         combo=as.character(combo))

for(combo in unique(all_distances_data %>%  pull(phenotype_combo))){
  # for different "tnumber"
  message(combo)
  
  for(tn in unique(all_combos_dists  %>% pull(tnumber))){
    # take into account only 1-NN dists < 100, because of the precision of fitdists (no of bins in histogram)
    # EXCEPT for From/To PanCK or negative
    # the parameters will be further optimized later
    if(str_detect(combo, 'PanCK\\+_') | str_detect(combo, 'negative_')){
      params <- fitdist(all_combos_dists %>% 
                          filter(combo == combo) %>% filter(tnumber == tn)  %>%
                          pull(dists) ,"weibull") %>% summary 
    } else{
      params <- fitdist(all_combos_dists %>% 
                          filter(combo == combo) %>% filter(tnumber == tn) %>% filter(dists < 100) %>%
                          pull(dists) ,"weibull") %>% summary 
    }
    class(params) <- c("fitdist", "fitdistr") 
    params <- broom::tidy(params)
    params <- params %>% mutate(tnumber=rep(tn,nrow(params))) %>%
      mutate(combo=rep(combo,nrow(params)))
    
    initial_params <- bind_rows(initial_params, params)
  }
}

initial_params_weib <- initial_params

# Return a vector with Weibull distribution parameters estimated by maximum likelihood (will be called as 'initial parameters')
# write_delim(initial_params_weib, '~/nabucco_spatial/data/IMCISION/spatial_data/initial_params_weib_imcision.tsv',delim='\t')

write_delim(initial_params_weib, paste(spatial_data_directory,'initial_params_weib_oct21.tsv', sep=''),delim='\t')