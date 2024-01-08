# STEP 1 - Takes as an input the spatial coordinates of mIF data and computes the 
#  1-NN distances histogram for each sample / cell FROM / cell TO combination
library(tidyverse)
library(ggpubr)
library(readxl)
library(reshape2)
library(dbscan)
library(broom)
library(zoo)
library(knitr)
library(NMF,quietly = T)
library(MASS)
library(car)
library(doParallel)
library(foreach)
library(spatstat)
library(raster)
library(gridExtra)
library(data.table)

spatial_data_directory <- '~/data/' # data directory
spatial_data_file      <- 'xy_data_exc_oct.tsv' # File with x/y coordinates from the multiplex immunofluorescence experiment

# Load Multiplex coordinates data
# Dataframe with columns: "tnumber', "analysisregion", "layer", "Xmin", "Xmax", "Ymin", "Ymax", "Xcenter", "Ycenter", "phenotype", "cellarea"
dumdat <- read_delim(paste(spatial_data_directory, spatial_data_file,sep=''),delim='\t')
dumdat$analysisregion <- dumdat$`Analysis Region`


# Load dataframe with areas
# Dataframe with columns: "tnumber', "Tumor.Area.mm2", "Stroma.Area.mm2", "Total.Area"
# areasumm <- readRDS('~/nabucco/rmd/areasumm_oct21.rds')
# area is in 'area' column

areas <- read_delim(paste(spatial_data_file, tissue_areas_file, sep=''), delim='\t')

# First NNdist all-against-all phenotypes
#' @param x Multiplex immunofluorescence coordinate dataframe (with Xcenter, Ycenter, phenotype, etc)
#' @param phenotype1 Phenotype From
#' @param phenotype2 Phenotype To
distances <- function(x, phenotype1,phenotype2) {
  # Create a point pattern
  spdat <- ppp(
    x = x$Xcenter,
    y = x$Ycenter,
    window = owin(
      xrange = c(min(x$Xcenter), max(x$Xcenter)),
      yrange = c(min(x$Ycenter), max(x$Ycenter))
    ),
    marks = x$phenotype
  )
  phenos <- unique(x$phenotype)
  
  #if phenotypes to compare are equal, then the 2nd nearest neighbour must be used, else the current point is compared to itself
  if(phenotype1 != phenotype2){
    res <- nncross(spdat[which(spdat$marks == phenotype1)], spdat[which(spdat$marks == phenotype2)], k = 1, what = 'dist')
  } else {
    res <- nncross(spdat[which(spdat$marks == phenotype1)], spdat[which(spdat$marks == phenotype2)], k = 2, what = 'dist')
    
  }
  return(res)
}

# Define which phenotypes will be studied (UNIQUE CELL TYPES IN YOUR DATASET)
phenotypes <- c("CD8+CD3+", 
                "CD3+",
                "CD3+FoxP3+",
                "CD68+",
                "CD20+",
                "negative",
                "PanCK+")

#call cores (parallelization)
# Comment this because of limitation on core usability in darwin
# nc <- detectCores()
# cl <- makePSOCKcluster(floor(nc / 4) )
# registerDoParallel(cl)


# Compute 1-NN for all combinations of samples (n=24) x phenotype FROM (n=7) x phenotype TO (n=7)
spatresall <-
  foreach(i = unique(dumdat %>% pull(tnumber)),
          .packages = c('sp', 'spatstat'),
          .final = function(x) setNames(x, unique(dumdat %>% pull(tnumber)))) %dopar% {
            x <- dumdat[which(dumdat$tnumber == i),] 
            phenotypes <- setNames(phenotypes, phenotypes)
            spres <-
              lapply(phenotypes, function(p1) {
                lapply(phenotypes, function(p2) {
                  distances(x, p1, p2)
                })
              })
            return(spres)
          }

# #run when done
# stopCluster(cl)


#Reshape lists of lists of lists into long data frame
dfl <- lapply(spatresall, function(tnr){
  lapply(tnr, function(p1){
    rbindlist(lapply(p1,function(Distance) {
      as.data.frame(Distance)
    }), idcol = "phenotype_to")
  })
})

dfl2 <- lapply(dfl, function(tnr){
  rbindlist(tnr, idcol = "phenotype_from")
})

dfl3 <- rbindlist(dfl2, idcol = "tnumber")
dfl4 <- dfl3

#add pseudocounts on distances for log transform AND division by 2 for conversion pixels to micrometers
dfl4$Distance <- (dfl4$Distance + 1) / 2

# Retrieve only biopsy data
# NNTUR <- subset(dfl4, Biopsy == "TUR")
NNTUR <- dfl4

# free up some memory. these temps are no longer needed
rm(dfl)
rm(dfl2)
rm(dfl3)
rm(dfl4)
gc() # do garbage collection to free up memory
#count 1-NN distances with sliding window of 5 microns (this is to create a smoothed histogram)

# Explore data from 0 to 300 microns 
SlidingBins_300_tumorandstroma <- rollapply(c(0:300),5,function(y) { return(c(min(y),max(y)))})

# Implement sliding distance bins
CountRows <- function(x) {
  cnt <- apply(SlidingBins_300_tumorandstroma,1, function(y) {nrow(x[which(x$Distance >= y[1] & x$Distance <= y[2])])})
  return(cnt)
}

# Group distance vector into distance bins to calculate coordinates of histogram
BinCounts_300_tumorandstroma <- apply(SlidingBins_300_tumorandstroma, 1, function(x) {
  CurrentWindow <- subset(NNTUR, Distance >= x[1] & Distance <= x[2])
  GroupCounts <- CurrentWindow %>% group_by(tnumber, phenotype_from, phenotype_to) %>%
    dplyr::summarise(N = dplyr::n())
  return(GroupCounts)
})

names(BinCounts_300_tumorandstroma) <- c(1:297)
dfslides_300_tumorandstroma <- as.data.frame.matrix(SlidingBins_300_tumorandstroma)
dfslides_300_tumorandstroma$V3 <- apply(dfslides_300_tumorandstroma,1,mean)
dfslides_300_tumorandstroma$ID <- row.names(dfslides_300_tumorandstroma)

LongBinCounts_300_tumorandstroma <- bind_rows(BinCounts_300_tumorandstroma, .id = "nth_window")
# WinMean is the X coordinate of your "calculated histogram"
LongBinCounts_300_tumorandstroma$WinMean <- dfslides_300_tumorandstroma$V3[match(LongBinCounts_300_tumorandstroma$nth_window, dfslides_300_tumorandstroma$ID)]

# Match total areas of the tissue 
LongBinCounts_300_tumorandstroma$Area <- LongBinCounts_300_tumorandstroma %>%
  left_join(data.frame(areas), by='tnumber') %>% 
  mutate(Area=total_area) %>% pull(Area)

# N.per.mm2 is the Y coordinate of your "calculated histogram"
# Normalized for tissue area
LongBinCounts_300_tumorandstroma$N.per.mm2 <- LongBinCounts_300_tumorandstroma$N / LongBinCounts_300_tumorandstroma$Area

# Function to compute Area under the curve 
funAUC <- function(x,y){
  id <- order(x)
  AUC <- sum(diff(x[id]) * rollmean(y[id],2))
  return(AUC)
}

AUCScaledSlides_300_tumorandstroma <- LongBinCounts_300_tumorandstroma %>%
  dplyr::group_by(phenotype_from, phenotype_to, tnumber) %>% 
  dplyr::mutate(N.per.mm2.scaled = N.per.mm2 / funAUC(WinMean, N.per.mm2))

#double check if the AUC of AUCs is 1
AUCcheck_panck <- AUCScaledSlides_300_tumorandstroma %>%
  dplyr::group_by(phenotype_from, phenotype_to, tnumber) %>%
  dplyr::mutate(AUC = funAUC(WinMean, N.per.mm2.scaled))

# Store output (coordinates of histograms of 1-NN distances per patient)
write_delim(AUCScaledSlides_300_tumorandstroma %>% mutate(phenotype_combo=paste(phenotype_from,phenotype_to,sep='_to_')), 
            paste(spatial_data_directory, 'AUCScaledSlides_300_tumorandstroma__thresholdStromaMargin_', as.character(threshold_pixels/2), '.tsv', sep=''),delim='\t')

