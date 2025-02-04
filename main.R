#------------------------------------------------------------------------------#
#   Preliminaries
#------------------------------------------------------------------------------#
{library(rstan)
  library(rstantools)
  library(hmde)
  library(tidyverse)
  library(vioplot)
  library(ggridges)
  library(cowplot)
  library(parallel)
  library(grid)
  library(plotly)
  library(GGally)
  library(htmlwidgets)
}

sp_codes <- c("alsebl",
              "beilpe",
              "cordbi",
              "faraoc",
              "gar2in",
              "hirttr",
              "jac1co",
              "pri2co",
              "protpa",
              "protte",
              "quaras",
              "swars1",
              "swars2",
              "simaam",
              "tachve",
              "tet2pa",
              "tri2tu")

col_vec <- c("#f8766d",
             "#e7861b",
             "#cf9400",
             "#afa100",
             "#72b000",
             "#00b81f",
             "#00bf7d",
             "#00c1aa",
             "#00bdd0",
             "#00b4ef",
             "#00a5ff",
             "#77a5ff",
             "#9590ff",
             "#cf78ff",
             "#f066ea",
             "#ff62bc",
             "#fc717f")

#Check if output directory structure exists, construct if not
if(!dir.exists("input/input_ignore/")){
  dir.create("input/input_ignore/")
}
if(!dir.exists("output/")){
  dir.create("output/")
  dir.create("output/data")
  dir.create('output/figures')
  dir.create("output/figures/diagnostic")
}
if(!dir.exists("output/figures")){
  dir.create("output/figures")
}
if(!dir.exists("output/figures/diagnostic")){
  dir.create("output/figures/diagnostic")
}

#------------------------------------------------------------------------------#
#   Get sample from cleaned data.
#------------------------------------------------------------------------------#
sample <- FALSE
if(sample){
  source("R/Sampling.R", local = TRUE)
  sample_data(sample_size = 300,
              data_path = "input/input_ignore/tree_data.csv",
              out_path = "input/")
}

#------------------------------------------------------------------------------#
#   Run models
#------------------------------------------------------------------------------#
run_hierarchical_models <- FALSE
if(run_hierarchical_models){
  source("R/Model_Fitting.R", local = TRUE)
  fit_hierarchical_models(sp_codes = sp_codes,
                          rstan_file_path = "input/")
}

run_species_level_models <- TRUE
if(run_species_level_models){
  source("R/Model_Fitting.R", local = TRUE)
  fit_species_level_models(sp_codes = sp_codes,
                           rstan_file_path = "input/")
}

#------------------------------------------------------------------------------#
#   Extract estimates from fits
#------------------------------------------------------------------------------#
extract_estimates <- TRUE
if(extract_estimates){
  source("R/Extract_Estimates.R", local = TRUE)
  extract_fit_estimates(sp_codes = sp_codes,
                        fit_file_path = "input/input_ignore/",
                        rstan_file_path = "input/",
                        plot_diagnostics = TRUE,
                        warmup = 1500)
}


#------------------------------------------------------------------------------#
#   Run analysis
#------------------------------------------------------------------------------#
analysis <- FALSE
if(analysis){
  source("R/Fit_Analysis.R", local = TRUE)
  run_analysis(sp_codes = sp_codes,
               est_file_path = "output/data/")
}
