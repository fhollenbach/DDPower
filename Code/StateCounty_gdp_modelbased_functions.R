
sim_fun_statecounty <- function(i){

  year1 <- dat_t1
  if(simulation_para$numtreated[i] == 6){
    year1$treatment_group <- year1$treated6
  }
  if(simulation_para$numtreated[i] == 12){
    year1$treatment_group <- year1$treated12
  }
  if(simulation_para$numtreated[i] == 24){
    year1$treatment_group <- year1$treated24
  }
  
  year1 <- year1 %>% 
    mutate(first_treated = case_when(treatment_group == "g1" ~ 2009, 
                                     treatment_group == "g2" ~ 2015, 
                                     treatment_group == "g3"~ 2019,
                                     TRUE ~ 0),
           first_treated_sa = if_else(first_treated == 0, 1000, first_treated), 
           treated = if_else(first_treated == 0, 0, 1),
           random_countycode = sample(CountyCode)) %>% ###### randomization
           select(CountyCode, Year, StateCode, log_GDP, GDP, first_treated, first_treated_sa, random_countycode)
  
  
  skeleton <- expand_grid(CountyCode = year1$CountyCode, Year = 2002:2023) %>% 
    left_join(select(year1, c("CountyCode", "StateCode", "first_treated", "first_treated_sa", "random_countycode")), by = c("CountyCode")) %>% 
    left_join(dat_draws, by = c("Year", "random_countycode")) %>% 
    select(CountyCode, Year, StateCode, log_GDP = log_GDP_fd, GDP = GDP_fd, first_treated, first_treated_sa, random_countycode)### draw from fd
  
  dat_c <- year1 %>%
    bind_rows(skeleton) %>% 
    arrange(CountyCode, Year) %>% 
    mutate(treated = if_else(first_treated == 0, 0, 1)) %>%  ###evertreated indicator
    ### next useful variable to be used later
    mutate(post_period = if_else(treated == 1 & Year >= first_treated, 1, 0), ## indicator for post treatment period for treated units
           time_since_treatment = if_else(first_treated != 0, Year - first_treated, -1000)) %>% ### time since treatment, 0 at time of treatment
    arrange(CountyCode, Year) %>%
    group_by(CountyCode) %>% 
    mutate(log_GDP = cumsum(log_GDP),
           max_time_since_treatment = max(time_since_treatment)) %>% 
    ungroup()
    ## filter out any county that has a negative gdp per capita value
    #mutate(negGDP = if_else(any(log_GDP < 0), 1, 0)) %>%
    #ungroup() %>% 
    #filter(negGDP == 0)
    
  treated_obs_counts <- dat_c %>% #### calculating number of treated obs for each treatment group
    filter(post_period==1) %>% ### necessary to get correct group specific event time effect such that overall treatment effect will be correct
    group_by(first_treated) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    filter(first_treated != 0)
  rm(year1, skeleton)
  #### calculate needed effect sizes
  
  k <- simulation_para$pct_effect[i] * sum(treated_obs_counts$n_obs) / sum(treated_obs_counts$n_obs * ratios)
  d_g <- k * ratios ### group specific average att
  
  dat_c <- dat_c %>% 
    mutate(group_specific_effect = case_when(first_treated == 2009 ~ d_g[1]/mean(1:15),
                                             first_treated == 2015 ~ d_g[2]/mean(1:9),
                                             first_treated == 2019 ~ d_g[3]/mean(1:5), ### divide group specific ate by average of the number of post treatment periods specific to each treatment group
                                             TRUE ~ 0)) ### create dynamic effects for the three groups, depending on length of post-treatment period
    
  
  dat_c <- dat_c %>% 
    group_by(CountyCode) %>%
    mutate(att = ifelse(post_period == 1, group_specific_effect * (max_time_since_treatment + 1 - time_since_treatment), 0)) %>%  ### ate decreasing in post-treatment period
    mutate(log_GDP_sim = log_GDP + log(att+1),
           # True ATT in log points (what estimators should recover)
           true_att_log = log_GDP_sim - log_GDP,
           # Outcome variable for estimation
           y = log_GDP_sim) %>% 
    mutate(time = as.integer(Year),
           id = as.character(CountyCode),
           y = as.numeric(y))
  
  
  dat_s <- dat_c %>%
    group_by(StateCode, Year) %>% 
    summarize(log_GDP = mean(log_GDP),
              first_treated = first(first_treated),
              first_treated_sa = first(first_treated_sa),
              treated = first(treated),
              post_period = first(post_period),
              time_since_treatment = first(time_since_treatment)) %>%
    ungroup() 
  
  
  #### state treatment effects
  treated_obs_counts_state <- dat_s %>% #### calculating number of treated obs for each treatment group
    filter(post_period==1) %>% ### necessary to get correct group specific event time effect such that overall treatment effect will be correct
    group_by(first_treated) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    filter(first_treated != 0)
  #### calculate needed effect sizes
  
  k <- simulation_para$pct_effect[i] * sum(treated_obs_counts_state$n_obs) / sum(treated_obs_counts_state$n_obs * ratios)
  d_g <- k * ratios ### group specific average att
  
  dat_s <- dat_s %>% 
    mutate(group_specific_effect = case_when(first_treated == 2009 ~ d_g[1]/mean(1:15),
                                             first_treated == 2015 ~ d_g[2]/mean(1:9),
                                             first_treated == 2019 ~ d_g[3]/mean(1:5), ### divide group specific ate by average of the number of post treatment periods specific to each treatment group
                                             TRUE ~ 0)) %>% ### create dynamic effects for the three groups, depending on length of post-treatment period
    group_by(StateCode) %>% 
    mutate(max_time_since_treatment = max(time_since_treatment)) %>% 
    ungroup() %>% 
    mutate(att = ifelse(post_period == 1, group_specific_effect * (max_time_since_treatment + 1 - time_since_treatment), 0), ## ate decreasing in post-treatment period
           log_GDP_sim = log_GDP + log(att+1),
           true_att_log = log_GDP_sim - log_GDP,
           # Outcome variable for estimation
           y = log_GDP_sim) %>% 
    ungroup() %>% 
    mutate(time = as.integer(Year),
           id = as.character(StateCode),
           y = as.numeric(y))

  dat_c <- na.omit(dat_c)
  
  dat_c_summary <- dat_c %>%
    ungroup() %>% 
    summarize(
      pre_treat_log_gdp = mean(log_GDP[treated == 1L & Year < first_treated], na.rm = TRUE),
      sd_log_gdp = sd(log_GDP, na.rm = TRUE),
      true_att_calc = mean(att[post_period == 1L], na.rm = TRUE), ### true att in levels
      true_att_log = mean(true_att_log[post_period == 1L & treated == 1L], na.rm = TRUE), ### true att as proportion
      n_counties = length(unique(CountyCode))
      ) 
  
  
  #### calculate the true treatment effects to save

  ### write county data
  write_csv(dat_c[, c("y", "StateCode", "CountyCode", "Year", "post_period")], here("..", "SimData_County_modelbased", paste0("Data_county_",i,".csv")))  
  #### write state data
  write_csv(dat_s[, c("y", "StateCode", "Year", "post_period")], here("..", "SimData_State_modelbased", paste0("Data_state_",i,".csv")))  
  
  ############# first county analysis
  #setFixest_nthreads(2, save = TRUE)
  
  ### twfe
  twfe_c <- feols(y ~ i(post_period, 0) | id + time, cluster = "StateCode", data = dat_c) %>%
    tidy(conf.int = TRUE) %>%
    mutate(method = "TWFE") %>%
    select(estimate, std.error, method)
  
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
    mutate(method = "Mundlak") %>%
    select(estimate, std.error, method)
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
      mutate(method = "CSA") %>% 
      select(estimate, std.error, method) 
  )

  bjs_c <- did_imputation(data = as.data.frame(dat_c), 
                        yname = "dv_impute", 
                        gname = "first_treated", 
                        tname = "time_impute", 
                        idname = "id_impute", 
                        first_stage = NULL, 
                        cluster_var = "StateCode") %>%
    mutate(method = "BJS") %>% 
    select(estimate, std.error, method) 
  
  suppressMessages(
    gardner_c <- did2s(data = dat_c,
                     yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
                     second_stage = ~i(post_period, ref=0),
                     treatment = "post_period",
                     cluster_var = "StateCode") %>%
      tidy(conf.int = FALSE) %>% 
      mutate(method = "Gardner") %>% 
      select(estimate, std.error, method)
  )
  
  # Sun and Abraham
  sab_est_c <- feols(y ~ sunab(first_treated_sa, time) | id + time,
                   cluster ~ StateCode,
                   data = dat_c) 
  sab_c <- sab_est_c %>%
    summary(agg = "att") %>%
    tidy(conf.int = FALSE) %>%
    mutate(method = "SA") %>% 
    select(estimate, std.error, method) 
  
  results_c <- bind_rows(twfe_c, wooldridge_c, csa_c, bjs_c, gardner_c, sab_c) %>% 
    mutate(true_att = dat_c_summary$true_att_calc,
           true_att_log_gdp = dat_c_summary$true_att_log,
           pretreatment_log_gdp = dat_c_summary$pre_treat_log_gdp,
           sd_log_gdp = dat_c_summary$sd_log_gdp,
           level = "County", 
           num_counties = dat_c_summary$n_counties) %>%  
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))
  #write_parquet(results_c, here("SimResults_County_opioids", paste0("results_county",i,".parquet")))
  rm(list=setdiff(ls(), c("dat_draws", "dat_t1", "simulation_para", "ratios", "i", "results_c","dat_s")))
  ############# now state analysis
  dat_s <- na.omit(dat_s)
  #setFixest_nthreads(2, save = TRUE)
  
  ### twfe
  twfe_s <- feols(y ~ i(post_period, 0) | id + time, cluster = "StateCode", data = dat_s) %>%
    tidy(conf.int = TRUE) %>%
    mutate(method = "TWFE") %>%
    select(estimate, std.error, method)
  
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
    mutate(method = "Mundlak") %>%
    select(estimate, std.error, method)
  
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
      mutate(method = "CSA") %>% 
      select(estimate, std.error, method) 
  )
  
  bjs_s <- did_imputation(data = as.data.frame(dat_s), 
                          yname = "dv_impute", 
                          gname = "first_treated", 
                          tname = "time_impute", 
                          idname = "id_impute", 
                          first_stage = NULL, 
                          cluster_var = "StateCode") %>%
    mutate(method = "BJS") %>% 
    select(estimate, std.error, method) 
  
   
  suppressMessages(
    gardner_s <- did2s(data = dat_s,
                       yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
                       second_stage = ~i(post_period, ref=0),
                       treatment = "post_period",
                       cluster_var = "StateCode") %>%
      tidy(conf.int = FALSE) %>% 
      mutate(method = "Gardner") %>% 
      select(estimate, std.error, method)
  )
  
  # Sun and Abraham
  sab_est_s <- feols(y ~ sunab(first_treated_sa, time) | id + time,
                     cluster ~ StateCode,
                     data = dat_s) 
  sab_s <- sab_est_s %>%
    summary(agg = "att") %>%
    tidy(conf.int = FALSE) %>%
    mutate(method = "SA") %>% 
    select(estimate, std.error, method) 
  
  dat_s_summary <- dat_s %>%
    ungroup() %>% 
    summarize(
      pre_treat_log_gdp = mean(log_GDP[treated == 1L & Year < first_treated], na.rm = TRUE),
      sd_log_gdp = sd(log_GDP, na.rm = TRUE),
      true_att_calc = mean(att[post_period == 1L], na.rm = TRUE), ### true att in levels
      true_att_log = mean(true_att_log[post_period == 1L & treated == 1L], na.rm = TRUE) ### true att as proportion
    ) 
  
  results_s <- bind_rows(twfe_s, wooldridge_s, csa_s, bjs_s, gardner_s, sab_s) %>% 
    mutate(true_att = dat_s_summary$true_att_calc,
           true_att_log_gdp = dat_s_summary$true_att_log,
           pretreatment_log_gdp = dat_s_summary$pre_treat_log_gdp,
           sd_log_gdp = dat_s_summary$sd_log_gdp,
           level = "State",
           n_counties = NA) %>%  
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))

  
  results <- bind_rows(results_c, results_s)
  rm(list=setdiff(ls(), c("dat_draws", "dat_t1", "simulation_para", "ratios", "i", "results")))
  
  return(results)
}


