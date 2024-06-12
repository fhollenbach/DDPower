#### file to run firm sims
library(devtools)
library(tidyverse)
library(here)
library(fixest)
library(data.table)
library(foreach)
library(broom)
library(did2s)
library(did)
library(didimputation)
library(doRNG)
library(doParallel)
library(modelsummary) #### use modelsummary
library(etwfe)
library(arrow)

source(here("Code", "download_data_firms.R"))

source(here("Code", "Firm_Sim_functions.R"))


writeLines(c(""), here("SessionInfo_firm_sim.txt"))
sink(here("SessionInfo_firm_sim.txt"))
devtools::session_info()
sink()

source(here("scripts", "server", "Firm_Sims_spring24.R"))