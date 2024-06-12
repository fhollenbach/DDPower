set.seed(12345)

comp <-  readRDS(here::here("Data", "simulation_data_new_panel.rds"))

#### simulation parameter
pct_effect <- c(0.025, 0.05, 0.10, 0.15)
units <- c(100, 250, 1000, 2000)
prop_treated <- c(0.1, 0.2, 0.4, 0.6)
iterations <- 500
DV = c("roa", "log_rev")

simulation_para <- bind_rows(replicate(iterations, expand_grid(pct_effect, units, prop_treated, DV), simplify = FALSE)) %>%
    arrange(units, pct_effect, prop_treated, DV)
#names(simulation_para) <- c("notreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]

iterations_save <- simulation_para$iteration[simulation_para$units %in% c(100, 1000)]


cl <- makeCluster(20, outfile = here("clustlogs.txt"))
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
          .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
    
        firm_simulation_panel(i)
    return(NULL)
}

stopImplicitCluster()
stopCluster(cl)


### run 2

cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(1237801)
foreach(i = 7801:dim(simulation_para)[1],
        .verbose = FALSE,
        .errorhandling = "stop",
        .combine = rbind,
        .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
          
          firm_simulation_panel(i)
          return(NULL)
        }

stopImplicitCluster()
stopCluster(cl)



### run 3

cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(12337501)
foreach(i = 37501:dim(simulation_para)[1],
        .verbose = FALSE,
        .errorhandling = "stop",
        .combine = rbind,
        .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
          
          firm_simulation_panel(i)
          return(NULL)
        }

stopImplicitCluster()
stopCluster(cl)



### run 3

cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(12338001)
foreach(i = 38001:dim(simulation_para)[1],
        .verbose = FALSE,
        .errorhandling = "stop",
        .combine = rbind,
        .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
          
          firm_simulation_panel(i)
          return(NULL)
        }

stopImplicitCluster()
stopCluster(cl)

### run 3

cl <- makeCluster(20, outfile = here("clustlogs.txt"))
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()


#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(12341001)
foreach(i = 41001:dim(simulation_para)[1],
        .verbose = FALSE,
        .errorhandling = "stop",
        .combine = rbind,
        .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
          
          firm_simulation_panel(i)
          return(NULL)
        }

stopImplicitCluster()
stopCluster(cl)



cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(1234200)
mdopar <- microbenchmark::microbenchmark(times = 1,
  foreach(i = 4200:dim(simulation_para)[1],
          .verbose = FALSE,
          .errorhandling = "stop",
          .combine = rbind,
          .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
    
        firm_simulation_panel(i)
    return(NULL)
}
)
stopImplicitCluster()
stopCluster(cl)





cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(12312000)
mdopar <- microbenchmark::microbenchmark(times = 1,
  foreach(i = 12000:dim(simulation_para)[1],
          .verbose = FALSE,
          .errorhandling = "stop",
          .combine = rbind,
          .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
        .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
    
        firm_simulation_panel(i)
    return(NULL)
}
)
stopImplicitCluster()
stopCluster(cl)



cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(12348400)
mdopar <- microbenchmark::microbenchmark(times = 1,
                                         foreach(i = 48400:dim(simulation_para)[1],
                                                 .verbose = FALSE,
                                                 .errorhandling = "stop",
                                                 .combine = rbind,
                                                 .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
                                                 .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
                                                   
                                                   firm_simulation_panel(i)
                                                   return(NULL)
                                                 }
)
stopImplicitCluster()
stopCluster(cl)



cl <- makeCluster(15)
registerDoParallel(cl)
#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()
registerDoRNG(12352800)
mdopar <- microbenchmark::microbenchmark(times = 1,
                                         foreach(i = 52800:dim(simulation_para)[1],
                                                 .verbose = FALSE,
                                                 .errorhandling = "stop",
                                                 .combine = rbind,
                                                 .export = c("firm_simulation_panel", "comp", "simulation_para", "iterations_save" ),
                                                 .packages = c("tidyverse", "did","did2s", "fixest", "didimputation", "here", "etwfe", "arrow"))%dopar%{
                                                   
                                                   firm_simulation_panel(i)
                                                   return(NULL)
                                                 }
)
stopImplicitCluster()
stopCluster(cl)