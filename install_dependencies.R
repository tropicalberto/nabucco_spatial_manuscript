# For a full list of package versions used in the original manuscript, see 'Materials and Methods' at the published manuscript 
packages_dependencies = c('dplyr','tidyverse','ggpubr', 'RColorBrewer','tidyr',
                          'spatstat', 'fitdistrplus','patchwork','complexHeatmap','circlize','glmnet','pROC','nlme','ggrepel','plyr','ggrastr')

uninstalled_packages = setdiff(packages_dependencies, rownames(installed_packages()))

if(length(missing.packages) > 0){
  install.packages(uninstalled_packages)
}
