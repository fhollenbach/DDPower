#### file to run state/county sims
library(devtools)
library(tidyverse)
library(fixest)
library(broom)
library(tictoc)
library(did2s)
library(did)
library(didimputation)
library(modelsummary) #### use modelsummary
library(etwfe)
library(arrow)
library(foreach)
library(here)
library(doRNG)
library(doParallel)
source(here("Code", "StateCounty_Sim_functions.R"))


writeLines(c(""), here("SessionInfo_StateCounty_sim.txt"))
sink(here("SessionInfo_StateCounty_sim.txt"))
devtools::session_info()
sink()

source(here("Code", "StateCounty_Sim_opioids.R"))
