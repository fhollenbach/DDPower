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
library(here)
library(tidyverse)

source(here("Code", "Firm_Sim_functions_scrambled.R"))

set.seed(12345)
comp <-  readRDS(here::here("Data", "simulation_data_new_panel.rds"))

mod_roa <- feols(roa ~ 1 | gvkey + fyear, cluster = "incorp", data = comp)
mod_logrev <- feols(log_rev ~ 1 | gvkey + fyear, cluster = "incorp", data = comp)

#### simulation parameter
pct_effect <- c(0.025, 0.05, 0.10, 0.15)
units <- 250
prop_treated <- 0.4
iterations <- 500
DV = c("roa", "log_rev")

simulation_para <- bind_rows(replicate(iterations, expand_grid(pct_effect, units, prop_treated, DV), simplify = FALSE)) %>%
    arrange(units, pct_effect, prop_treated, DV)
#names(simulation_para) <- c("notreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]


cl <- makeCluster(7, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(123)
foreach(i = 1:dim(simulation_para)[1],
          .verbose = FALSE,
          .errorhandling = "stop",
          .combine = rbind,
          .export = c("firm_simulation_scrambled", "mod_logrev", "mod_roa", "comp", "simulation_para"),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
          
          firm_simulation_scrambled(i)
    return(NULL)
}

stopImplicitCluster()
stopCluster(cl)

