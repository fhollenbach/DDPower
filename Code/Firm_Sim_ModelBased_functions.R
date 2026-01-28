
safe_model <- function(expr, method_label) {
  tryCatch(
    expr,
    error = function(e) {
      tibble::tibble(
        estimate  = NA_real_,
        std.error = NA_real_,
        method    = method_label
      )
    }
  )
}



firm_simulation_modelbased <- function(i){
  
  firms_sample <- filter(firms, sample <= simulation_para$units[i] & is.na(sample) == FALSE) %>%
    select(gvkey, time, first_treated, log_rev, log_rev_win) %>% 
    mutate(random_gvkey = sample(gvkey))### sample of firms)) 
  skeleton <- expand_grid(gvkey = firms_sample$gvkey, time = 2:20) %>% 
    left_join(select(firms_sample, c("gvkey", "first_treated", "random_gvkey")), by = c("gvkey"))
  ##### which DV, winsorized or not
  if(simulation_para$DV[i] == "rev"){
    skeleton <- skeleton %>% 
      left_join(select(comp_draws, c("random_gvkey", "time", "log_rev_fd")), by = c("random_gvkey", "time")) %>% 
      select(gvkey, time, first_treated, outcome = log_rev_fd)### draw from fd
    firms_sample$outcome <- firms_sample$log_rev
  }
  if(simulation_para$DV[i] == "rev_win"){
    skeleton <- skeleton %>% 
      left_join(select(comp_draws, c("random_gvkey", "time", "log_rev_win_fd")), by = c("random_gvkey", "time")) %>% 
      select(gvkey, time, first_treated, outcome = log_rev_win_fd)### draw from fd
    firms_sample$outcome <- firms_sample$log_rev_win
  }
  
  dat_staggered <- firms_sample %>%
    select(gvkey, time, first_treated, outcome) %>%
    bind_rows(skeleton) %>% 
    arrange(gvkey, time) %>% 
    mutate(treated = if_else(first_treated == 0, 0, 1), ###evertreated indicator
           first_treated_sa = if_else(first_treated == 0, 1000, first_treated)) %>% ### special variable for S&A
    ### next useful variable to be used later
    mutate(post_period = if_else(treated == 1 & time >= first_treated, 1, 0), ## indicator for post treatment period for treated units
           time_since_treatment = if_else(first_treated != 0, time - first_treated, 0)) ### time since treatment, 0 at time of treatment
  
  treated_obs_counts <- dat_staggered %>% #### calculating number of treated obs for each treatment group
    filter(post_period==1) %>% ### necessary to get correct group specific event time effect such that overall treatment effect will be correct
    group_by(first_treated) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    filter(first_treated != 0)

  
  k <- simulation_para$pct_effect[i] * sum(treated_obs_counts$n_obs) / sum(treated_obs_counts$n_obs * ratios)
  d_g <- k * ratios ### group specific average att
  
  dat_staggered <- dat_staggered %>% 
    mutate(group_specific_effect = case_when(first_treated == 10 ~ d_g[1]/mean(1:11),
                                             first_treated == 13 ~ d_g[2]/mean(1:8),
                                             first_treated == 16 ~ d_g[3]/mean(1:5), ### divide group specific ate by average of the number of post treatment periods specific to each treatment group
                                             TRUE ~ 0)) %>% ### create dynamic effects for the three groups, depending on length of post-treatment period
    arrange(gvkey, time) %>%
    group_by(gvkey) %>% 
    mutate(dv = cumsum(outcome),
           max_time_since_treatment = max(time_since_treatment)) %>% ### event time ate increases over linearly over time
    group_by(gvkey) %>% 
    arrange(time, .by_group = TRUE) %>%  # make sure years are in order per firm
    mutate(
      lag_dv = lag(dv, 1),
      dv_trended = if_else(time > 1, dv + lag_dv*simulation_para$autocor[i], dv)) %>%
    ungroup() %>% 
    mutate(att = ifelse(post_period == 1, group_specific_effect * (max_time_since_treatment + 1 - time_since_treatment), 0), ### multiplicative ate on untransformed DV for treated in post-treatment period
           y_log = dv_trended + log(att+1),
           true_att_log = y_log - dv_trended) %>% #### new outcome 
    mutate(sd_dv = sd(dv_trended),### sd of outcome, without treatment
           pre_treat_dv = mean(dv_trended[treated == 1 & time < first_treated])) ### pre treatment mean of dv for treated units
  
  dat_staggered_summary <- dat_staggered %>%
    ungroup() %>% 
    summarize(
      att = mean(att[post_period == 1L]),
      sd_dv = sd(dv_trended, na.rm = TRUE),
      pre_treat_dv = mean(dv_trended[treated == 1L & time < first_treated], na.rm = TRUE), ### true att in levels
      true_att_log = mean(true_att_log[post_period == 1L & treated == 1L], na.rm = TRUE))
  
  
  dat_staggered <- dat_staggered %>%
    mutate(time = as.integer(time),
           firm_id = as.integer(gvkey),
           y_log = as.numeric(y_log),
           time_to_treatment = case_when(first_treated == 0 ~ -1000, ### for etwfe
                                         TRUE ~ time - first_treated))
  
  if(i %in% iterations_save){ ### saving relevant data sets for stata runs
    write_csv(dat_staggered[, c("y_log", "firm_id", "time", "post_period")], here("..", "SimData_firms_modelbased", paste0("Data_firms_staggered",i,".csv")))  
  }
  dat_staggered <- na.omit(dat_staggered)
  ### twfe
  suppressMessages(
    twfe <- safe_model({
      feols(y_log ~ i(post_period, 0) | firm_id + time, cluster = "firm_id", data = dat_staggered) %>%
      tidy(conf.int = TRUE) %>%
      mutate(t = -100, 
             method = "TWFE") %>%
      select(estimate, std.error, method)},
      method_label = "TFWE")
  )
  # 
  suppressMessages(
    wooldridge <- safe_model({
      ww <- etwfe(
      fml  = y_log ~ 1,           # outcome ~ controls
      tvar = time,                  # time variable
      gvar = first_treated,
      data = dat_staggered,                 # dataset
      vcov = ~firm_id) %>%           # vcov adjustment (here: clustered) 
      emfx("simple") %>%
      as_tibble() %>%
      mutate(t = -100, 
             method = "Mundlak") %>%
      select(estimate, std.error, method)},
      method_label = 'Mundlak')
  )
    # 
  # callaway and sant'anna
  suppressMessages(
    csa <- safe_model({
      att_gt(yname = "y_log",
                      tname = "time",
                      idname = "firm_id",
                      gname = "first_treated",
                      xformla = ~1,
                      data = dat_staggered,
                      clustervars = "firm_id", 
                      control_group = "notyettreated") %>%
      aggte(type = "simple") %>%
      tidy(conf.int = FALSE) %>%
      mutate(t = -100,
             method = "CSA") %>% 
      select(estimate, std.error, method)},
      method_label = 'CSA')
  )
  
  
  dat_staggered$time_impute <- as.numeric(dat_staggered$time)
  dat_staggered$id_impute <- as.numeric(dat_staggered$firm_id)
  dat_staggered$dv_impute <- as.numeric(dat_staggered$y_log)
  suppressMessages(
    bjs <- safe_model({
      did_imputation(data = as.data.frame(dat_staggered), yname = "dv_impute", gname = "first_treated", tname = "time_impute", idname = "id_impute", first_stage = NULL) %>%
      mutate(method = "BJS") %>% 
      select(estimate, std.error, method)},
      method_label = 'BJS')
  )
  suppressMessages(
    gardner <- safe_model({
      did2s(data = dat_staggered,
                     yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
                     second_stage = ~i(post_period, ref=0),
                     treatment = "post_period",
                     cluster_var = "id_impute") %>%
      tidy(conf.int = FALSE) %>% 
      mutate(method = "Gardner") %>% 
      select(estimate, std.error, method)},
      method_label = 'Gardner')
  )
  
  # Sun and Abraham
  suppressMessages(
    sab <- safe_model({feols(y_log ~ sunab(first_treated_sa, time) | firm_id + time,
                     cluster ~ firm_id,
                     data = dat_staggered) %>%
        summary(agg = "att") %>%
        tidy(conf.int = FALSE) %>%
        mutate(t = -100,
              method = "SA") %>% 
        select(estimate, std.error, method)},
        method_label = 'SA')
  )

  results_staggered <- bind_rows(twfe, wooldridge, csa, bjs, gardner, sab) %>% 
    mutate(true_att = dat_staggered_summary$true_att_log,
           att = dat_staggered_summary$att,
           pretreatment_dv = dat_staggered_summary$pre_treat_dv,
           sd_dv = dat_staggered_summary$sd_dv,
           true_prop_treated = true_prop_treated,
           type = "staggered") %>% 
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))
  rm(list=setdiff(ls(), c("firms",
                          "comp_draws",
                          "simulation_para",
                          "iterations_save",
                          "ratios",
                          "safe_model", 
                          "firms_sample", 
                          "skeleton", 
                          "results_staggered",
                          "simulation_para", 
                          "ratios", 
                          "i", 
                          "true_prop_treated")))
  
  
  ### non staggered
  dat_std <- firms_sample %>%
    select(gvkey, time, first_treated, outcome) %>%
    bind_rows(skeleton) %>% 
    arrange(gvkey, time) %>% 
    mutate(first_treated = if_else(first_treated != 0, 13, first_treated),
           treated = if_else(first_treated == 0, 0, 1), ###evertreated indicator
           first_treated_sa = if_else(first_treated == 0, 1000, first_treated)) %>% ### special variable for S&A
    ### next useful variable to be used later
    mutate(post_period = if_else(treated == 1 & time >= first_treated, 1, 0), ## indicator for post treatment period for treated units
           time_since_treatment = if_else(first_treated != 0, time - first_treated, 0)) ### time since treatment, 0 at time of treatment
  
  dat_std <- dat_std %>% 
    mutate(group_specific_effect = simulation_para$pct_effect[i]) %>% ### create dynamic effects for the three groups, depending on length of post-treatment period
    arrange(gvkey, time) %>%
    group_by(gvkey) %>% 
    mutate(dv = cumsum(outcome)) %>% ### event time ate increases over linearly over time
    group_by(gvkey) %>% 
    arrange(time, .by_group = TRUE) %>%  # make sure years are in order per firm
    mutate(
      lag_dv = lag(dv, 1),
      dv_trended = if_else(time > 1, dv + lag_dv*simulation_para$autocor[i], dv)) %>% 
    ungroup() %>% 
    mutate(att = ifelse(post_period == 1, group_specific_effect, 0), ### multiplicative ate on untransformed DV for treated in post-treatment period
           y_log = dv_trended + log(att+1),
           true_att_log = y_log - dv_trended) %>% #### new outcome 
    mutate(sd_dv = sd(dv_trended),### sd of outcome, without treatment
           pre_treat_dv = mean(dv_trended[treated == 1 & time < first_treated])) ### pre treatment mean of dv for treated units
  
  dat_std_summary <- dat_std %>%
    ungroup() %>% 
    summarize(
      att = mean(att[post_period == 1L]),
      sd_dv = sd(dv_trended, na.rm = TRUE),
      pre_treat_dv = mean(dv_trended[treated == 1L & time < first_treated], na.rm = TRUE), ### true att in levels
      true_att_log = mean(true_att_log[post_period == 1L & treated == 1L], na.rm = TRUE))
  
  dat_std <- dat_std %>%
    mutate(time = as.integer(time),
           firm_id = as.integer(gvkey),
           y_log = as.numeric(y_log),
           time_to_treatment = case_when(first_treated == 0 ~ -1000, ### for etwfe
                                         TRUE ~ time - first_treated))

  dat <- na.omit(dat_std)
  ### twfe
  suppressMessages(
    twfe <- safe_model({
      feols(y_log ~ i(post_period, 0) | firm_id + time, cluster = "firm_id", data = dat_std) %>%
        tidy(conf.int = TRUE) %>%
        mutate(t = -100, 
               method = "TWFE") %>%
        select(estimate, std.error, method)},
      method_label = 'TWFE')
  )
  
  suppressMessages(
    wooldridge <- safe_model({
      etwfe(
        fml  = y_log ~ 1,           # outcome ~ controls
        tvar = time,                  # time variable
        gvar = first_treated,
        data = dat_std,                 # dataset
        vcov = ~firm_id) %>%             # vcov adjustment (here: clustered) 
        emfx("simple") %>%
        as_tibble() %>%
        mutate(method = "Mundlak") %>%
        select(estimate, std.error, method)},
      method_label = 'Mundlak')
  )
  # callaway and sant'anna
  suppressMessages(
    csa <- safe_model({att_gt(yname = "y_log",
                      tname = "time",
                      idname = "firm_id",
                      gname = "first_treated",
                      xformla = ~1,
                      data = dat_std,
                      clustervars = "firm_id", 
                      control_group = "notyettreated") %>% 
        aggte(type = "simple") %>%
        tidy(conf.int = FALSE) %>%
        mutate(t = -100,
               method = "CSA") %>% 
        select(estimate, std.error, method)},
        method_label = 'CSA')
  )

  dat_std$time_impute <- as.numeric(dat_std$time)
  dat_std$id_impute <- as.numeric(dat_std$firm_id)
  dat_std$dv_impute <- as.numeric(dat_std$y_log)
  suppressMessages(
    bjs <- safe_model({
      did_imputation(data = as.data.frame(dat_std), yname = "dv_impute", gname = "first_treated", tname = "time_impute", idname = "id_impute", first_stage = NULL) %>%
        mutate(method = "BJS") %>% 
        select(estimate, std.error, method)},
      method_label = 'BJS')
  )
  suppressMessages(
    gardner <- safe_model({
      did2s(data = dat_std,
            yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
            second_stage = ~i(post_period, ref=0),
            treatment = "post_period",
            cluster_var = "id_impute") %>%
            tidy(conf.int = FALSE) %>% 
            mutate(method = "Gardner") %>% 
            select(estimate, std.error, method)},
      method_label = 'Gardner')
  )
  
  # Sun and Abraham
  suppressMessages(
    sab_est <- safe_model({
      feols(y_log ~ sunab(first_treated_sa, time) | firm_id + time,
                     cluster ~ firm_id,
                     data = dat_std) %>%
        summary(agg = "att") %>%
        tidy(conf.int = FALSE) %>%
        mutate(method = "SA") %>% 
        select(estimate, std.error, method)},
      method_label = 'SA')
  )
 
  results_std <- bind_rows(twfe, wooldridge, csa, bjs, gardner, sab_est) %>% 
    mutate(true_att = dat_std_summary$true_att_log,
           att = dat_std_summary$att,
           pretreatment_dv = dat_std_summary$pre_treat_dv,
           sd_dv = dat_std_summary$sd_dv,
           true_prop_treated = true_prop_treated,
           type = "standard") %>% 
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))
  
  results <- bind_rows(results_staggered, results_std)
  rm(list=setdiff(ls(), c("firms",
                          "comp_draws",
                          "simulation_para",
                          "iterations_save",
                          "ratios",
                          "true_prop_treated",
                          "firm_simulation_modelbased", 
                          "safe_model", 
                          "results", 
                          "simulation_para", 
                          "ratios", 
                          "i")))
  
  #write_parquet(results, paste("/work/DDpower/SimResults_firms_modelbased/results_firms",i,".parquet", sep =""))
  return(results)
}


