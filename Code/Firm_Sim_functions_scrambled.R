



firm_simulation_scrambled <- function(i){
  
  NUnits <- simulation_para$units[i]
  EffectSize <- simulation_para$pct_effect[i]
  PropTreated <- simulation_para$prop_treated[i]

  
  if(simulation_para$DV[i] == "roa"){
    mod <- mod_roa
  } else if(simulation_para$DV[i] == "log_rev"){
    mod <- mod_logrev
  } ### select dependent var)
  
  firm_fes <- fixef(mod)$gvkey
  n_firm_fes <- length(fixef(mod)$gvkey)
  year_fes <- fixef(mod)$fyear
  n_year_fes <- length(fixef(mod)$fyear)
  
  ### how many treated firms
  treated_units <- ceiling(NUnits*(PropTreated))
  #### group sizes
  group_sizes <- ceiling(treated_units/3)
  
  ### draw unique firms based on number of units
  sim_firms <- tibble(
    id = sample.int(10000, NUnits, replace = FALSE),
    firm_fe = sample(firm_fes, NUnits, replace = FALSE),
    incorp = sample(state.abb, NUnits, replace = TRUE),
    first_treated = c(rep(0, NUnits - (group_sizes*3)), rep(6, group_sizes),  rep(11, group_sizes),  rep(16, group_sizes))) %>%
    mutate(treated = if_else(first_treated == 0, 0, 1),
           first_treated_sa = if_else(first_treated == 0, 1000, first_treated)) ### special variable for S&A
  
  # pull year FE from the empirical distribution with replacement
  sim_year_fe <- tibble(
    time = 1:20,
    year_fe = sample(year_fes, 20, replace = TRUE)
  )
  
  dat <- crossing(sim_firms, sim_year_fe) %>% 
    mutate(resid = sample(mod$residuals, NUnits*20, replace = TRUE),
           dv = firm_fe + year_fe + resid) %>%
    arrange(id, time) %>% 
    mutate(sd_dv = sd(dv),### sd
           pre_treat_dv = mean(dv[treated == 1 & time <= first_treated])) %>%  ### pre treatment mean roa for treated units
    
    ### next useful variable to be used later
    mutate(post_period = if_else(treated == 1 & time >= first_treated, 1, 0), ## indicator for post treatment for treated units
           time_since_treatment = if_else(first_treated != 0, time - first_treated, 0), #### relative time indicator
           group_specific_effect = case_when(first_treated == 6 ~ (EffectSize*pre_treat_dv / 3),
                                             first_treated == 11 ~ (EffectSize*pre_treat_dv),
                                             first_treated == 16 ~ (EffectSize*pre_treat_dv * 3),
                                             TRUE ~ 0)) ### create heterogeneous effects, for the three groups
  
  
  dat <- dat %>%
    group_by(id) %>%
    mutate(ate = case_when(post_period == 0 ~ 0,
                           post_period == 1 ~ group_specific_effect / ((max(time_since_treatment)/ 2) + 1)), ### divide group specific ate by number of post treatment periods specific to each treatment group
           ate = if_else(post_period == 1, cumsum(ate), 0)) %>% ### event time ate increases over time
    rowwise() %>%
    mutate(y = dv + ate) %>%
    ungroup()
  
  #### calculate the true treatment effects to save
  avg_att <- dat %>%
    filter(post_period == 1) %>%
    pull(ate) %>%
    mean()

  dat <- dat %>%
    mutate(time = as.integer(time),
           id = as.integer(id),
           y = as.numeric(y),
           time_to_treatment = case_when(first_treated == 0 ~ -1000,
                                         TRUE ~ time - first_treated))
  
  dat <- na.omit(dat)
  ### twfe
  twfe <- feols(y ~ i(post_period, 0) | id + time, data = dat) %>%
    tidy(conf.int = TRUE) %>%
    mutate(t = -100, 
           method = "TWFE") %>%
    select(t, estimate, std.error, method)
  
  wooldridge_est <- etwfe(
    fml  = y ~ 0,           # outcome ~ controls
    tvar = time,                  # time variable
    gvar = first_treated,
    gref = 0, # group variable (with bespoke ref. level)
    data = dat,                 # dataset
    vcov = ~id)            # vcov adjustment (here: clustered) 
  
  
  wooldridge  <- wooldridge_est %>%
    emfx("simple") %>%
    as_tibble() %>%
    mutate(t = -100, 
           method = "Mundlak") %>%
    select(t, estimate, std.error, method)
  
  # callaway and sant'anna
  suppressWarnings(
    csa_est <- att_gt(yname = "y",
                      tname = "time",
                      idname = "id",
                      gname = "first_treated",
                      xformla = ~1,
                      data = dat,
                      clustervars = "id", 
                      control_group = "notyettreated")
  )
  suppressWarnings(
    csa <- csa_est %>%
      aggte(type = "simple") %>%
      tidy(conf.int = FALSE) %>%
      mutate(t = -100,
             method = "CSA") %>% 
      select(t, estimate, std.error, method) 
  )
  
  dat$time_impute <- as.numeric(dat$time)
  dat$id_impute <- as.numeric(dat$id)
  dat$dv_impute <- as.numeric(dat$y)
  
  bjs <- did_imputation(data = as.data.frame(dat), yname = "dv_impute", gname = "first_treated", tname = "time_impute", idname = "id_impute", first_stage = NULL) %>%
    mutate(t = -100,
           method = "BJS") %>% 
    select(t, estimate, std.error, method) 
  
  suppressMessages(
    gardner <- did2s(data = dat,
                     yname = "dv_impute", first_stage = ~ 0 | id_impute + time_impute,
                     second_stage = ~i(post_period, ref=0),
                     treatment = "post_period",
                     cluster_var = "id") %>%
      tidy(conf.int = FALSE) %>% 
      mutate(t = -100,
             method = "Gardner") %>% 
      select(t, estimate, std.error, method)
  )
  
  # Sun and Abraham
  sab_est <- feols(y ~ sunab(first_treated_sa, time) | id + time,
                   cluster ~ id,
                   data = dat) 
  sab <- sab_est %>%
    summary(agg = "att") %>%
    tidy(conf.int = FALSE) %>%
    mutate(t = -100,
           method = "SA") %>% 
    select(t, estimate, std.error, method) 
  
  results <- bind_rows(twfe, wooldridge, csa, bjs, gardner, sab) %>% 
    mutate(true_att = avg_att,
           pretreatment_dv = unique(dat$pre_treat_dv),
           sd_dv = unique(dat$sd_dv))
  results <- results %>% 
    bind_cols(slice(simulation_para[i,], rep(1, each = 6)))
  write_parquet(results, here("..", "SimResults_firms_scrambled", paste0("results_firms_scrambled",i,".parquet")))
  return(NULL)
}


