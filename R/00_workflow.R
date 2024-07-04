library(tidyverse)
library(purrr) #map, map2, pmap
library(emmeans) #emmans, contrast
library(ggrepel)
library(ggpubr)
library(patchwork) #?
library(miceRanger) #miceRanger
library(lme4) #lmer
library(scales) #rescale
library(vegan) #prc, scores
library(MuMIn) #r.squaredGLMM
library(here)
library(factoextra)
library(FactoMineR)
library(merTools) #predictInterval
library(paran)
library(metafor) #rma
library(broom)
library(ggh4x) #facet_nested
library(sf)
library(maps)
library(mapproj)

# Create subfolders
if(!dir.exists(here("plot"))){dir.create(here("plot"))}
if(!dir.exists(here("output"))){dir.create(here("output"))}
if(!dir.exists(here("data"))){dir.create(here("data"))}
# Download the processed_data/TransPlantNetwork.RData to "/data/" folder under working directory
## from OSF:https://osf.io/9874t/?view_only=90e7bafd8eed40eb960a22ce3098c613


nseed = 123456 # Seed number for setting seed for the miceRanger and predictInterval functions

ndataset = 100 # Number of datasets to be generated for the miceRanger
niteration  = 10 # Number of iterations for the miceRanger
nsim = 1000 # Number of iterations for predictInterval
set.nperm = 999 # set the number of permutations for the permutation tests for PRC

# Import the data
source(here("R","01_import.R"))

# Organize the data
## Get the nested data set, get the pools per experiment, get functional traits, compute community weighted means
source(here("R","02_organize.R"))

# Build principal response curves (PRC)
## Run permutation tests, produce PRC plots per experiment
source(here("R", "03_runPRC.R"))

# Rates of divergence and convergence at community level
## Get slopes per experiment, build linear models, build figures for the paper
source(here("R", "04_rates.R"))

emm_options(pbkrtest.limit = 6000)
emm_options(lmerTest.limit = 6000)
# Species pool weights 
source(here("R", "05_species_weights.R"))
## You might have memory usage issues for model predictions. 
## In that case, you can edit the memory usage of R.
## For macOS: library(usethis) usethis::edit_r_environ() 
## add R_MAX_VSIZE=32Gb in the open file and restart R. 

# Assess experimental variations
source(here("R", "06_experimental_variation.R"))

# Additionnal descriptive figures and tables
source(here("R", "07_descriptive_figures_tables.R"))


