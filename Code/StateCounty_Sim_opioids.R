

registerDoRNG(123)

data <-  read_delim(here("Data", "DrugsAlcohol.txt"))  %>% 
  filter(County != "District of Columbia, DC") %>% 
  mutate(Population = as.double(Population)) %>%
  group_by(County) %>% 
  mutate(Missing = if_else(any(Deaths == "Missing"), 1, 0),
         MissingPop = if_else(any(Population == 0), 1, 0), 
         allSupp = if_else(all(Deaths == "Suppressed"), 1, 0),
         meanSupp = mean(Deaths == "Suppressed")) %>% 
  ungroup() %>% 
  mutate(Suppressed = if_else(Deaths == "Suppressed", 1, 0),
         Deaths = if_else(Deaths == "Suppressed" | Deaths == "Missing", NA, Deaths)) %>%
  filter(Missing == 0 & MissingPop == 0) %>% 
  select(-c(Notes, `Year Code`, `Crude Rate`)) %>% 
  rename(CountyCode = `County Code`) %>% 
  mutate(Year = as.double(Year),
         Deaths = as.double(Deaths),
         Population = as.double(Population),
         death_rate = (Deaths/(Population/100000)),
         max_death = max(death_rate, na.rm = TRUE)*(Population/100000),
         Deaths = if_else(Suppressed == 1, sample(c(1:9), 1, replace = FALSE), Deaths),
         Deaths = if_else(Suppressed == 1 & Deaths > max_death, round(max_death, 0), Deaths),
         StateCode = str_sub(CountyCode, 1, 2),
         death_rate = Deaths/(Population/100000)) %>% 
  group_by(CountyCode) %>% 
  mutate(Zeros = if_else(any(Deaths == 0), 1, 0),
         allZeros = if_else(all(Deaths ==0), 1, 0),
         meanZero = mean(Deaths == 0)) %>% 
  ungroup() %>% 
  filter(meanZero < 0.10 & allSupp ==0 & meanSupp < 0.50)





length(unique(data$CountyCode))
summary(data)

##### so lets have treatment effects of 2%, 5%, and 10%, 20%
pct_effect <- c(0.025, 0.05, 0.10, 0.15)
treated <- c(6, 12, 24)
iterations <- 500

simulation_para <- bind_rows(replicate(iterations, expand_grid(treated, pct_effect), simplify = FALSE))
names(simulation_para) <- c("notreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]


cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(123)
tic("sim county state run")
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
## 3300
cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(1233300)
tic("sim county state run")
sim_results <- foreach(i = 3300:dim(simulation_para)[1],
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



##run 3
## 5050
cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()

registerDoRNG(1235050)
tic("sim county state run")
sim_results <- foreach(i = 5050:dim(simulation_para)[1],
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



