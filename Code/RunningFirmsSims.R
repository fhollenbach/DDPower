#### file to run state/county sims
library(groundhog)

pkgs <- c("devtools", 
          "tidyverse", 
          "here", 
          "fixest", 
          "data.table", 
          "foreach", 
          "broom", 
          "did2s", 
          "did", 
          "didimputation", 
          "doRNG", 
          "doParallel", 
          "modelsummary", 
          "etwfe", 
          "tinytable", 
          "marginaleffects", 
          "tictoc", 
          "DIDmultiplegtDYN",
          "future", 
          "future.apply",
          "purrr")
groundhog.day <- "2026-01-01"
groundhog.library(pkgs, groundhog.day, cores=10)


### load arrow outside groundhog
library(arrow)

#source(here("Code", "download_data_firms.R"))

source(here("Code", "Firm_Sim_ModelBased_functions.R"))


writeLines(c(""), here("SessionInfo_firm_sim.txt"))
sink(here("SessionInfo_firm_sim.txt"))
devtools::session_info()
sink()

source(here("Code", "Firm_Sims_ModelBased.R"))
