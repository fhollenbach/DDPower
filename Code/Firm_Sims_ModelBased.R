library(groundhog)

pkgs <- c("devtools", 
          "tidyverse", 
          "here", 
          "fixest", 
          "data.table", 
          "foreach", 
          "broom", 
          "did2s", 
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
          "future", 
          "future.apply",
          "purrr",
          "e1071")
groundhog.day <- "2026-01-01"
groundhog.library(pkgs, groundhog.day, cores=10)

### load arrow outside groundhog
library(arrow)

source(here("Code", "Firm_Sim_ModelBased_functions.R"))

set.seed(52)### based on random.org on Nov 17

comp <-  readRDS(here::here("Data", "simulation_data_new_panel_092025.rds")) %>% 
  mutate(log_rev = log(revt),
         log_rev_win = log(revt_win)) %>% 
  arrange(gvkey, fyear) %>% 
  group_by(gvkey) %>%
  mutate(log_rev_l1 = lag(log_rev, 1),
         log_rev_win_l1 = lag(log_rev_win, 1)) %>%
  ungroup() %>% 
  mutate(log_rev_fd = log_rev - log_rev_l1,
         log_rev_win_fd = log_rev_win - log_rev_win_l1)

#### stats for data description
mean(comp$revt_win)
sd(comp$revt_win)
skewness(comp$revt_win)
### unwin
mean(comp$revt)
sd(comp$revt)
skewness(comp$revt)

#### simulation parameter
pct_effect <- c(0, 0.05, 0.1, 0.15)
DV <- c("rev_win", "rev")
units <- c(100, 250, 500, 1000, 1500)
prop_treated <- 0.4
iterations <- 500
autocor <- c(0, 0.35)

NoFirms <- length(unique(comp$gvkey))
simulation_para <- bind_rows(replicate(iterations, expand_grid(pct_effect, units, prop_treated, DV, autocor), simplify = FALSE)) %>%
    arrange(autocor, DV, pct_effect, units, prop_treated)
#names(simulation_para) <- c("notreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]

iterations_save <- simulation_para$iteration[simulation_para$autocor == 0]


#### initial values
comp_t1 = comp %>% 
  filter(fyear == 2001) %>%
  select(gvkey, time, log_rev, log_rev_win) 

comp_draws <- comp %>% 
  filter(fyear > 2001) %>%
  select(gvkey, time, log_rev_fd, log_rev_win_fd) %>% 
  mutate(log_rev_fd = log_rev_fd %>%
           pmin(quantile(., 0.99, na.rm = TRUE)) %>%  # cap at 99th percentile
           pmax(quantile(., 0.01, na.rm = TRUE)),
         log_rev_win_fd = log_rev_win_fd %>%
                  pmin(quantile(., 0.99, na.rm = TRUE)) %>%  # cap at 99th percentile
                  pmax(quantile(., 0.01, na.rm = TRUE))      # cap at 1st percentile
  ) %>% 
  rename(random_gvkey = gvkey)
### how many treated firms
treated_units <- NoFirms*0.4

#### group sizes
group_sizes <- ceiling(treated_units/4) #### group sizes
ratios <- c(1, 3, 5) #### effect size ratios for each group
ratios <- ratios / mean(ratios)  # normalize so average = 1
true_prop_treated <- group_sizes*4/NoFirms #### actual exact proportion treated
#### draw firms for each group and sample size
#### 
####
####
set.seed(651) # from random.org Nov 17, 2025
####
group_sizes_samples <- c(100, 250, 500, 1000, 1500)*0.4/4 #### group sizes
untreated <- c(100, 250, 500, 1000, 1500) - (group_sizes_samples*4)
g2_g3_size <- group_sizes_samples
g1_size <- group_sizes_samples*2

firms <- tibble(gvkey = sample(unique(comp$gvkey), NoFirms, replace = FALSE),
                firm_id = 1:NoFirms, ### new ids
                first_treated = c(rep(0, NoFirms - (group_sizes*4)), rep(10, group_sizes*2), rep(13, group_sizes), rep(16, group_sizes))) %>% 
  left_join(comp_t1, by = "gvkey") %>% 
  arrange(first_treated, firm_id) %>% 
  group_by(first_treated) %>%
  mutate(counter = 1:n()) %>% 
  ungroup() %>%
  mutate(sample = case_when(first_treated == 0 & counter <= untreated[1] ~ 100,
                           first_treated == 0 & counter > untreated[1] & counter <= untreated[2] ~ 250,
                          first_treated == 0 & counter > untreated[2] & counter <= untreated[3] ~ 500,
                          first_treated == 0 & counter > untreated[3] & counter <= untreated[4] ~ 1000,
                          first_treated == 0 & counter > untreated[4] & counter <= untreated[5] ~ 1500,
                          (first_treated == 16 | first_treated == 13) & counter <= g2_g3_size[1] ~ 100,
                          (first_treated == 16 | first_treated == 13) & counter > g2_g3_size[1] & counter <= g2_g3_size[2] ~ 250,
                          (first_treated == 16 | first_treated == 13) & counter > g2_g3_size[2] & counter <= g2_g3_size[3]  ~ 500,
                          (first_treated == 16 | first_treated == 13) & counter > g2_g3_size[3] & counter <= g2_g3_size[4]~ 1000,
                          (first_treated == 16 | first_treated == 13) & counter > g2_g3_size[4] & counter <= g2_g3_size[5]~ 1500,
                          first_treated == 10 & counter <= g1_size[1] ~ 100,
                          first_treated == 10 & counter > g1_size[1] & counter <= g1_size[2] ~ 250,
                          first_treated == 10 & counter > g1_size[2] & counter <= g1_size[3] ~ 500,
                          first_treated == 10 & counter > g1_size[3] & counter <= g1_size[4] ~ 1000,
                          first_treated == 10 & counter > g1_size[4] & counter <= g1_size[5] ~ 1500))
                          




n_workers <- 5L

cl <- makeCluster(n_workers, type = "PSOCK")
registerDoParallel(cl)

## Reproducible RNG across workers
registerDoRNG(5306)  # any seed you like

clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(did)
    library(did2s)
    library(fixest)
    library(didimputation)
    library(etwfe)
    library(arrow)
    library(here)
    library(broom)
    
  })
  
  Sys.setenv(
    OMP_NUM_THREADS       = "1",
    MKL_NUM_THREADS       = "1",
    OPENBLAS_NUM_THREADS  = "1",
    VECLIB_MAXIMUM_THREADS = "1",   # relevant for macOS Accelerate
    RCPP_PARALLEL_NUM_THREADS = "1",
    ARROW_NUM_THREADS     = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  fixest::setFixest_nthreads(1L)
  
  NULL
})


clusterExport(
  cl,
  c("firms",
    "comp_draws",
    "simulation_para",
    "iterations_save",
    "ratios",
    "true_prop_treated",
    "firm_simulation_modelbased"),
  envir = environment()
)


n_sims <- nrow(simulation_para)  # or 30000 if hard-coded

chunk_size <- 200L   # tune this: 200–1000 is usually fine
chunk_ids <- split(
  seq_len(n_sims),
  ceiling(seq_len(n_sims) / chunk_size)
)

# directory for chunk result files
results_dir <- here("..", "SimResults_firms_modelbased", "chunks")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

for (c_idx in seq_along(chunk_ids)) {
  ids <- chunk_ids[[c_idx]]
  
  message("Running chunk ", c_idx, " of ", length(chunk_ids),
          " (iterations ", min(ids), "–", max(ids), ")")
  
  # parallel over this chunk's indices
  res <- foreach(i = ids,
                        .errorhandling = "stop",
                        .export = character(0),
                        .packages = character(0),
                        .combine = "rbind") %dopar% {
                          # firm_simulation_modelbased:
                          # - writes CSV itself when i %in% iterations_save
                          # - returns list(result=..., data_payload=..., file_name=...)
                          firm_simulation_modelbased(i)
                        }
  
  
  # Write this chunk's results to Parquet
  chunk_file <- file.path(
    results_dir,
    sprintf("results_chunk_%05d.parquet", c_idx)
  )
  arrow::write_parquet(res, chunk_file, compression = "zstd")
  
  # Free memory for this chunk
  rm(res); gc()
}

# Done with parallel work
stopCluster(cl)