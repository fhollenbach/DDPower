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
          "purrr")
groundhog.day <- "2025-11-01"
groundhog.library(pkgs, groundhog.day, cores=10)

source(here("Code", "StateCounty_gdp_modelbased_functions.R"))


registerDoRNG(9784)

data <-  read_csv(here("Data", "CAGDP9__ALL_AREAS_2001_2023.csv"))  %>% 
  filter(LineCode == 1) %>%  # All industry total
  filter(nchar(trimws(GeoFIPS)) == 5 | nchar(trimws(GeoFIPS)) == 6) %>%  # FIPS codes (may have leading space)
  filter(GeoName != "District of Columbia, DC") %>% 
  mutate(
    GeoFIPS = trimws(GeoFIPS),
    # Exclude state totals (end in 000) and US total
    is_county = !grepl("000$", GeoFIPS) & GeoFIPS != "00000"
  ) %>%
  filter(is_county) %>% 
  select(GeoFIPS, GeoName, as.character(2001:2023)) %>% 
  pivot_longer(
    -c(GeoFIPS, GeoName),
    names_to = "Year",
    values_to = "GDP"
  ) %>% 
  mutate(
    Year = as.integer(Year),
    # Handle suppressed values "(D)" and convert to numeric
    GDP = case_when(
      GDP == "(D)" ~ NA_character_,
      GDP == "(L)" ~ NA_character_,  # Less than $50,000
      GDP == "(NA)" ~ NA_character_,
      TRUE ~ GDP
    ),
    GDP = as.numeric(gsub(",", "", GDP)),  # Remove commas and convert
    # GDP is in thousands of dollars, convert to millions
    GDP = GDP/1000,
    # Extract state and county codes
    CountyCode = GeoFIPS,
    StateCode = substr(GeoFIPS, 1, 2),
    County = GeoName
  ) %>% 
  #left_join(pop_data, by = c("CountyCode", "Year")) %>%
  select(CountyCode, StateCode, County, Year, GDP) %>% #, total_population
  filter(!is.na(GDP) & GDP > 0 ) %>% #& !is.na(total_population
  mutate(log_GDP = log(GDP)) %>% #GDPpC = GDP/total_population,
  #filter(GDP < quantile(GDP, prob = 0.95)[1]) %>% 
  ### count obs by county
  group_by(CountyCode) %>%
  mutate(obs_count = n()) %>%
  ungroup() %>% 
  filter(obs_count == 23) %>%  # keep only counties with full obs
  arrange(CountyCode, Year) %>%
  group_by(CountyCode) %>%
  mutate(
    log_GDP_l1 = lag(log_GDP, 1),
    GDP_l1 = lag(GDP, 1),
    #GDPpC_l1 = lag(GDPpC, 1),
    GDP_fd = GDP - GDP_l1,
    #GDPpC_fd = GDPpC - GDPpC_l1,
    log_GDP_fd = log_GDP - log_GDP_l1  # This is approximately growth rate
  ) %>%
  ungroup()
  
#### for effect size calculation
ratios <- c(1, 3, 5) #### effect size ratios for each group
ratios <- ratios / mean(ratios)  # normalize so average = 1
#### draw firms for each group and sample size

# Calculate state-level average initial log_GDP and number of counties
state_chars <- data %>%
  filter(Year == min(Year)) %>%
  group_by(StateCode) %>%
  summarise(
    mean_initial_log_GDP = mean(log_GDP),
    n_counties = n()) %>%
  mutate(
    # Rank each variable separately
    rank_gdp = rank(mean_initial_log_GDP),
    rank_counties = rank(n_counties),
    # Combined rank (average of the two)
    combined_rank = (rank_gdp + rank_counties) / 2
  ) %>%
  arrange(combined_rank) %>%
  mutate(sorted_rank = row_number())
# Shuffle within blocks of 3 to add randomness while maintaining balance
set.seed(2784)
state_chars <- state_chars %>%
  mutate(pool = ((sorted_rank - 1) %% 3) + 1) %>%
  group_by(pool) %>%
  mutate(rand_within_pool = sample(n())) %>%
  ungroup()


treatment_assign <- state_chars %>%
  mutate(
    treated24 = case_when(
      pool == 1 & rand_within_pool <= 8 ~ "g1",
      pool == 2 & rand_within_pool <= 8 ~ "g2",
      pool == 3 & rand_within_pool <= 8 ~ "g3",
      TRUE ~ NA_character_
    ),
    treated12 = case_when(
      pool == 1 & rand_within_pool <= 4 ~ "g1",
      pool == 2 & rand_within_pool <= 4 ~ "g2",
      pool == 3 & rand_within_pool <= 4 ~ "g3",
      TRUE ~ NA_character_
    ),
    treated6 = case_when(
      pool == 1 & rand_within_pool <= 2 ~ "g1",
      pool == 2 & rand_within_pool <= 2 ~ "g2",
      pool == 3 & rand_within_pool <= 2 ~ "g3",
      TRUE ~ NA_character_
    )
  ) %>%
  select(StateCode, treated6, treated12, treated24, mean_initial_log_GDP, n_counties)

# Verify balance - now never-treated (NA) should also be balanced
treatment_assign %>%
  pivot_longer(cols = c(treated6, treated12, treated24),
               names_to = "scheme", values_to = "group") %>%
  mutate(group = replace_na(group, "never")) %>%
  group_by(scheme, group) %>%
  summarise(
    n_states = n(),
    mean_init_gdp = mean(mean_initial_log_GDP),
    mean_n_counties = mean(n_counties),
    total_counties = sum(n_counties),
    .groups = "drop"
  ) %>%
  print(n = 20)



data <- data %>%
  left_join(treatment_assign, by = "StateCode")



dat_t1 <- data %>% 
  filter(Year == min(Year)) %>% 
  select(CountyCode, Year, StateCode, log_GDP, GDP, treated6, treated12, treated24)

dat_draws <- data %>% 
  filter(Year > min(Year)) %>%
  select(CountyCode, Year, log_GDP_fd, GDP_fd) %>% 
  rename(random_countycode = CountyCode)
rm(data)

##### so lets have treatment effects of 2%, 5%, and 10%, 20%
pct_effect <- c(0, 0.01, 0.025, 0.05)
treated <- c(6, 12, 24)
iterations <- 500

simulation_para <- bind_rows(replicate(iterations, expand_grid(treated, pct_effect), simplify = FALSE))
names(simulation_para) <- c("numtreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]


n_workers <- 7L

cl <- makeCluster(n_workers, type = "PSOCK")
registerDoParallel(cl)

## Reproducible RNG across workers
registerDoRNG(8280)  # any seed you like

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
  c("dat_t1",
    "dat_draws",
    "simulation_para",
    "ratios",
    "sim_fun_statecounty"),
  envir = environment()
)


n_sims <- nrow(simulation_para)  # or 30000 if hard-coded

chunk_size <- 100L   # tune this: 200–1000 is usually fine
chunk_ids <- split(
  seq_len(n_sims),
  ceiling(seq_len(n_sims) / chunk_size)
)

# directory for chunk result files
results_dir <- here("..", "SimResults_statecounty_modelbased", "chunks")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
tic()
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
                   sim_fun_statecounty(i)
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
toc()
# Done with parallel work
stopCluster(cl)


