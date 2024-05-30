


data <-  read_delim(here("Data", "ExternalDeathsCounty.txt"))  %>% 
  filter(County != "District of Columbia, DC") %>% 
  group_by(County) %>% 
  mutate(Suppress = if_else(any(Deaths == "Suppressed") | any(Deaths == "Missing"), 1, 0)) %>% 
  ungroup() %>% 
  mutate(Deaths = if_else(Deaths == "Suppressed" | Deaths == "Missing", NA, Deaths)) %>%
  filter(Suppress == 0) %>% 
  select(-c(Notes, `Year Code`, `Crude Rate`)) %>% 
  rename(CountyCode = `County Code`) %>% 
  mutate(Year = as.double(Year),
         Deaths = as.double(Deaths),
         Population = as.double(Population),
         death_rate = (Deaths/Population)*100000) %>%
  # create StateCode from first two letters of CountyCode
  mutate(StateCode = str_sub(CountyCode, 1, 2))




##### so lets have treatment effects of 2%, 5%, and 10%, 20%
pct_effect <- c(0.025, 0.05, 0.10, 0.15)
treated <- c(6, 12, 24)
iterations <- 500

simulation_para <- bind_rows(replicate(iterations, expand_grid(treated, pct_effect), simplify = FALSE))
names(simulation_para) <- c("notreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]


cl <- makeCluster(10)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(123)
tic("sim runs")
sim_results <- foreach(i = 1:dim(simulation_para)[1],
                       .verbose = FALSE,
                       .errorhandling = "stop",
                       .combine = rbind,
                       .export = c("sim_fun_statecounty", "data"),
                       .packages = c("tidyverse", "did",
                                     "did2s", "fixest",
                                     "didimputation", "here",
                                     "etwfe", "arrow")
) %dopar% {
  
  sim_fun_statecounty(i)
  return(NULL)
}
toc()
stopCluster(cl)



##run 2
cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(123)
tic("sim runs")
sim_results <- foreach(i = 650:dim(simulation_para)[1],
                       .verbose = FALSE,
                       .errorhandling = "stop",
                       .combine = rbind,
                       .export = c("sim_fun_external_county", "data"),
                       .packages = c("tidyverse", "did",
                                     "did2s", "fixest",
                                     "didimputation", "here",
                                     "etwfe", "arrow")
) %dopar% {
  
  sim_fun_external_county(i)
  return(NULL)
}
toc()
stopCluster(cl)



##run 3
cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(123)
tic("sim runs")
sim_results <- foreach(i = 1500:dim(simulation_para)[1],
                       .verbose = FALSE,
                       .errorhandling = "stop",
                       .combine = rbind,
                       .export = c("sim_fun_external_county", "data"),
                       .packages = c("tidyverse", "did",
                                     "did2s", "fixest",
                                     "didimputation", "here",
                                     "etwfe", "arrow")
) %dopar% {
  
  sim_fun_external_county(i)
  return(NULL)
}
toc()
stopCluster(cl)



##run 4
cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(123)
tic("sim runs")
sim_results <- foreach(i = 1900:dim(simulation_para)[1],
                       .verbose = FALSE,
                       .errorhandling = "stop",
                       .combine = rbind,
                       .export = c("sim_fun_external_county", "data"),
                       .packages = c("tidyverse", "did",
                                     "did2s", "fixest",
                                     "didimputation", "here",
                                     "etwfe", "arrow")
) %dopar% {
  
  sim_fun_external_county(i)
  return(NULL)
}
toc()
stopCluster(cl)



##run 5
cl <- makeCluster(10)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(123)
tic("sim runs")
sim_results <- foreach(i = 2700:dim(simulation_para)[1],
                       .verbose = FALSE,
                       .errorhandling = "stop",
                       .combine = rbind,
                       .export = c("sim_fun_external_county", "data"),
                       .packages = c("tidyverse", "did",
                                     "did2s", "fixest",
                                     "didimputation", "here",
                                     "etwfe", "arrow")
) %dopar% {
  
  sim_fun_external_county(i)
  return(NULL)
}
toc()
stopCluster(cl)

