
sim_fun_statecounty <- function(i){
  #### draw three treatment groups
  no_treated <- simulation_para$notreated[i]
  effect <- simulation_para$pct_effect[i]
  states <- sample(unique(data$StateCode), 50, replace = FALSE)
  g1 <- states[1: (no_treated/3)] 
  g2 <- states[((no_treated/3) + 1) : (2*(no_treated/3))]
  g3 <- states[(2*(no_treated/3) + 1) : (3*(no_treated/3))]
  
  dat <- data %>% 
    mutate(first_treated = case_when(StateCode %in% g1 ~ 2006, 
                                     StateCode %in% g2 ~ 2012, 
                                     StateCode %in% g3 ~ 2016,
                                     TRUE ~ 0),
           first_treated_sa = if_else(first_treated == 0, 1000, first_treated), 
           treated = if_else(first_treated == 0, 0, 1),
           post_period = if_else(treated == 1 & Year >= first_treated, 1, 0),
           time_since_treatment = if_else(first_treated == 0, -1000, Year - first_treated),
           group_specific_effect = case_when(first_treated == 2006 ~ (effect / 3),
                                             first_treated == 2012 ~ (effect),
                                             first_treated == 2016 ~ (effect * 3),
                                             TRUE ~ 0)) 
  ### county level data
  dat_c <- dat %>%
    group_by(CountyCode) %>%
    mutate(ate = case_when(post_period == 0 ~ 0,
                           post_period == 1 ~ group_specific_effect / ((max(time_since_treatment)/ 2) + 1)), ### divide group specific ate by number of post treatment periods specific to each treatment group
           ate = if_else(post_period == 1, cumsum(ate), 0)) %>% ### event time ate increases over time
    rowwise() %>%
    mutate(deaths_remove = if_else(post_period == 0, 0, sum(rbinom(Deaths, 1, prob = ate))),
           deaths_sim = Deaths - deaths_remove,
           true_ate_prob = deaths_remove/Deaths,
           deathrate_sim = (deaths_sim/Population)*100000,
           y = deathrate_sim,
           true_death_rate = (Deaths/Population)*100000,
           true_att_calc = (deathrate_sim - true_death_rate)) %>% 
    ungroup() %>%
    mutate(sd_deathrate = sd(true_death_rate),### sd
           pre_treat_deathrate = mean(true_death_rate[treated == 1 & Year <= first_treated])) %>% ### pre treatment mean roa for treated units
    mutate(time = as.integer(Year),
           id = as.character(CountyCode),
           y = as.numeric(y),
           time_to_treatment = time_since_treatment) 
  
  ### state level data
  dat_s <- dat_c %>%
    group_by(StateCode, Year) %>% 
    summarize(Deaths = sum(Deaths),
              Population = sum(Population),
              first_treated = unique(first_treated),
              first_treated_sa = unique(first_treated_sa),
              treated = unique(treated),
              post_period = unique(post_period),
              time_since_treatment = unique(time_since_treatment),
              group_specific_effect = unique(group_specific_effect),
              ate = mean(ate),
              deaths_remove = sum(deaths_remove)) %>% 
    rowwise() %>%
    mutate(deaths_sim = Deaths - deaths_remove,
           true_ate_prob = deaths_remove/Deaths,
           deathrate_sim = (deaths_sim/Population)*100000,
           y = deathrate_sim,
           true_death_rate = (Deaths/Population)*100000,
           true_att_calc = (deathrate_sim - true_death_rate)) %>% 
    ungroup() %>%
    mutate(sd_deathrate = sd(true_death_rate),### sd
           pre_treat_deathrate = mean(true_death_rate[treated == 1 & Year <= first_treated])) %>% ### pre treatment mean roa for treated units
    mutate(time = as.integer(Year),
           id = as.character(StateCode),
           y = as.numeric(y),
           time_to_treatment = time_since_treatment) 
  
  #### calculate the true treatment effects to save
  ### county level
  avg_att_c <- dat_c %>%
    filter(post_period == 1) %>%
    pull(ate) %>%
    mean()
  
  avg_att_deathrate_c <- dat_c %>%
    filter(post_period == 1) %>%
    pull(true_att_calc) %>%
    mean()

  ### state level
  avg_att_s <- dat_s %>%
    filter(post_period == 1) %>%
    pull(ate) %>%
    mean()
  
  avg_att_deathrate_s <- dat_s %>%
    filter(post_period == 1) %>%
    pull(true_att_calc) %>%
    mean()
  ### write county data
  write_csv(dat_c[, c("y", "StateCode", "CountyCode", "Year", "post_period")], here("SimData_County_opioids", paste0("Data_county_",i,".csv")))  
  #### write state data
  write_csv(dat_s[, c("y", "StateCode", "Year", "post_period")], here("SimData_State_opioids", paste0("Data_state_",i,".csv")))  
  
  ############# first county analysis
  dat_c <- na.omit(dat_c)
  #setFixest_nthreads(2, save = TRUE)
  
  ### twfe
  twfe_c <- feols(y ~ i(post_period, 0) | id + time, cluster = "StateCode", data = dat_c) %>%
    tidy(conf.int = TRUE) %>%
    mutate(t = -100, 
           method = "TWFE") %>%
    select(t, estimate, std.error, method)
  
  wooldridge_est_c <- etwfe(
    fml  = y ~ 0,           # outcome ~ controls
    tvar = time,                  # time variable
    gvar = first_treated,
    gref = 0, # group variable (with bespoke ref. level)
    data = dat_c,                 # dataset
    vcov = ~ StateCode)            # vcov adjustment (here: clustered) 
  
  wooldridge_c  <- wooldridge_est_c %>%
    emfx("simple") %>%
    as_tibble() %>%
    mutate(t = -100, 
           method = "Mundlak") %>%
    select(t, estimate, std.error, method)
  dat_c$time_impute <- as.numeric(dat_c$time)
  dat_c$id_impute <- as.numeric(dat_c$id)
  dat_c$dv_impute <- as.numeric(dat_c$y)
  
  # callaway and sant'anna
  suppressWarnings(
    csa_est_c <- att_gt(yname = "y",
                      tname = "time",
                      idname = "id_impute",
                      gname = "first_treated",
                      xformla = ~1,
                      data = dat_c,
                      clustervars = "StateCode",
                      control_group = "notyettreated")
  )
  suppressWarnings(
    csa_c <- csa_est_c %>%
      aggte(type = "simple") %>%
      tidy(conf.int = FALSE) %>%
      mutate(t = -100,
             method = "CSA") %>% 
      select(t, estimate, std.error, method) 
  )

  bjs_c <- did_imputation(data = as.data.frame(dat_c), 
                        yname = "dv_impute", 
                        gname = "first_treated", 
                        tname = "time_impute", 
                        idname = "id_impute", 
                        first_stage = NULL, 
                        cluster_var = "StateCode") %>%
    mutate(t = -100,
           method = "BJS") %>% 
    select(t, estimate, std.error, method) 
  
  suppressMessages(
    gardner_c <- did2s(data = dat_c,
                     yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
                     second_stage = ~i(post_period, ref=0),
                     treatment = "post_period",
                     cluster_var = "StateCode") %>%
      tidy(conf.int = FALSE) %>% 
      mutate(t = -100,
             method = "Gardner") %>% 
      select(t, estimate, std.error, method)
  )
  
  # Sun and Abraham
  sab_est_c <- feols(y ~ sunab(first_treated_sa, time) | id + time,
                   cluster ~ StateCode,
                   data = dat_c) 
  sab_c <- sab_est_c %>%
    summary(agg = "att") %>%
    tidy(conf.int = FALSE) %>%
    mutate(t = -100,
           method = "SA") %>% 
    select(t, estimate, std.error, method) 
  
  pct_counties_treated <- length(unique(dat$CountyCode[dat$post_period == 1]))/length(unique(dat$CountyCode)) 
  no_counties_treated <- length(unique(dat$CountyCode[dat$post_period == 1]))
  
  results_c <- bind_rows(twfe_c, wooldridge_c, csa_c, bjs_c, gardner_c, sab_c) %>% 
    mutate(pct_counties_treated = pct_counties_treated,
           no_counties_treated = no_counties_treated,
          true_att = avg_att_c,
          true_att_deathrate = avg_att_deathrate_c,
          pretreatment_deathrate = unique(dat_c$pre_treat_deathrate),
          sd_deathrate = unique(dat_c$sd_deathrate)) %>%  
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))
  write_parquet(results_c, here("SimResults_County_opioids", paste0("results_county",i,".parquet")))
  rm(dat_c, results_c)
  
  ############# now state analysis
  dat_s <- na.omit(dat_s)
  #setFixest_nthreads(2, save = TRUE)
  
  ### twfe
  twfe_s <- feols(y ~ i(post_period, 0) | id + time, cluster = "StateCode", data = dat_s) %>%
    tidy(conf.int = TRUE) %>%
    mutate(t = -100, 
           method = "TWFE") %>%
    select(t, estimate, std.error, method)
  
  wooldridge_est_s <- etwfe(
    fml  = y ~ 0,           # outcome ~ controls
    tvar = time,                  # time variable
    gvar = first_treated,
    gref = 0, # group variable (with bespoke ref. level)
    data = dat_s,                 # dataset
    vcov = ~ StateCode)            # vcov adjustment (here: clustered) 
  
  wooldridge_s  <- wooldridge_est_s %>%
    emfx("simple") %>%
    as_tibble() %>%
    mutate(t = -100, 
           method = "Mundlak") %>%
    select(t, estimate, std.error, method)
  
  dat_s$time_impute <- as.numeric(dat_s$time)
  dat_s$id_impute <- as.numeric(dat_s$id)
  dat_s$dv_impute <- as.numeric(dat_s$y)
  
  # callaway and sant'anna
  suppressWarnings(
    csa_est_s <- att_gt(yname = "y",
                        tname = "time",
                        idname = "id_impute",
                        gname = "first_treated",
                        xformla = ~1,
                        data = dat_s,
                        clustervars = "StateCode",
                        control_group = "notyettreated")
  )
  suppressWarnings(
    csa_s <- csa_est_s %>%
      aggte(type = "simple") %>%
      tidy(conf.int = FALSE) %>%
      mutate(t = -100,
             method = "CSA") %>% 
      select(t, estimate, std.error, method) 
  )
  
  bjs_s <- did_imputation(data = as.data.frame(dat_s), 
                          yname = "dv_impute", 
                          gname = "first_treated", 
                          tname = "time_impute", 
                          idname = "id_impute", 
                          first_stage = NULL, 
                          cluster_var = "StateCode") %>%
    mutate(t = -100,
           method = "BJS") %>% 
    select(t, estimate, std.error, method) 
  
   
  suppressMessages(
    gardner_s <- did2s(data = dat_s,
                       yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
                       second_stage = ~i(post_period, ref=0),
                       treatment = "post_period",
                       cluster_var = "StateCode") %>%
      tidy(conf.int = FALSE) %>% 
      mutate(t = -100,
             method = "Gardner") %>% 
      select(t, estimate, std.error, method)
  )
  
  # Sun and Abraham
  sab_est_s <- feols(y ~ sunab(first_treated_sa, time) | id + time,
                     cluster ~ StateCode,
                     data = dat_s) 
  sab_s <- sab_est_s %>%
    summary(agg = "att") %>%
    tidy(conf.int = FALSE) %>%
    mutate(t = -100,
           method = "SA") %>% 
    select(t, estimate, std.error, method) 
  
  results_s <- bind_rows(twfe_s, wooldridge_s, csa_s, bjs_s, gardner_s, sab_s) %>% 
    mutate(true_att = avg_att_s,
           true_att_deathrate = avg_att_deathrate_s,
           pretreatment_deathrate = unique(dat_s$pre_treat_deathrate),
           sd_deathrate = unique(dat_s$sd_deathrate)) %>%  
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))
  write_parquet(results_s, here("SimResults_State_opioids", paste0("results_state",i,".parquet")))
  rm(results_s, dat_s)
  return(NULL)
}


