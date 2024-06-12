library(RPostgres)
library(here)
# source my passwords

source(here("..", "pw-r.R"))

# Connect to WRDS Server --------------------------------------------------
wrds_con <- dbConnect(Postgres(),
                      host = 'wrds-pgdata.wharton.upenn.edu',
                      port = 9737,
                      user = "fhoegb",
                      password = wrds,
                      dbname = 'wrds',
                      sslmode = 'require')

# download compustat data
comp <- tbl(wrds_con, sql("SELECT gvkey, fyear,  at, revt, ni, datafmt, popsrc, consol,indfmt, fic, sich FROM comp.funda")) %>%
  # filter the data - between 1979 and 2015, non-missing assets, and in the US
  filter(indfmt == 'INDL' & datafmt == 'STD' & popsrc == 'D' & consol == 'C' & !is.na(fyear) & 
           fyear %>% between(1990, 2020) &  fic == "USA" & at !=0) %>% 
  # make ROA variable
  group_by(gvkey) %>% 
  mutate(gvkey = as.numeric(gvkey),
         lag_at = lag(at, 1, order = fyear),
         #roa = ni/lag_at,
         roa = revt/lag_at) %>%  # maximum data gap by firm
  ungroup() %>% 
    # drop missing ROA
  filter(is.na(roa) == FALSE) %>% 
  group_by(gvkey) %>% 
  mutate(min_year = min(fyear),
         year_dif = fyear - lag(fyear, 1, order = fyear), ### check whether time series consec
         time = fyear - min_year, ### time counter
         max_time = max(time), ### time series length by firm
         first_gap = min(time[year_dif > 1]),
         last_gap = max(time[year_dif > 1]),
         length_to_gap = first_gap - 1,
         length_after_gap = max_time - last_gap) %>% 
  collect() %>% 
  arrange(gvkey, fyear) %>% 
  filter(max_time >= 20 & (length_to_gap >= 20 | length_after_gap >= 20 | is.na(first_gap))) %>% ### everything less than 20 years can already go
  ##### remove data afte/before gap that has shorter time series
  group_by(gvkey) %>%
  mutate(to_keep = if_else(length_to_gap > length_after_gap, "first", "last")) %>%
  ungroup() %>%
  filter((to_keep == "first" & time < last_gap) | is.na(to_keep) | (to_keep == "last" & time > first_gap)) %>%
  #### recreate time variables
  group_by(gvkey) %>% 
  mutate(min_year = min(fyear), ## first year now in data by firm, t = 1
         time = fyear - min_year  + 1, ### time counter starts at 1 for each firm
         max_time = max(time), ### last year for each firm
         time_dif = time - lag(time, 1, order = time), ## check for gaps
         time_dif = if_else(is.na(time_dif), 1, time_dif), ### first period is always 1, but not lag available
         max_gap = max(time_dif)) %>%
  ungroup() %>%
  filter(max_gap == 1 & max_time >= 20) %>%
  mutate(start_time = max_time - 19)


length(unique(comp$gvkey))
### 2800 firms with 20 years
table(comp$time)
# download comp header which has more info on location etc
comp_header <- tbl(wrds_con, sql("SELECT * FROM crsp.comphead")) %>% 
  mutate(gvkey = as.numeric(gvkey)) %>% 
  collect()

# merge in state of incorporation and industry information
comp <- comp %>% 
  left_join(comp_header %>% select(gvkey, incorp, sic)) %>% 
  # clean up SIC code - use historical sic code if available, if not use header sic
  mutate(sich = coalesce(sich, sic))

# make sure that each firm has at least 25 observations and no missingness
comp <- comp %>% 
  group_by(gvkey) %>% 
  add_tally() %>% 
  filter(n >= 20) %>% 
  ungroup()


# winsorize ROA at 99, and censor at -1
wins <- function(x) {
  # winsorize and return
  case_when(
    is.na(x) ~ NA_real_,
    x < quantile(x, 0.01, na.rm = TRUE) ~ quantile(x, 0.01, na.rm = TRUE),
    x > quantile(x, 0.99, na.rm = TRUE) ~ quantile(x, 0.99, na.rm = TRUE),
    TRUE ~ x
  )
}

# winsorize ROA by year
comp <- comp %>% 
  mutate(roa = wins(roa),
  revt = wins(revt), 
  log_rev = log(revt+1)) %>% 
  arrange(gvkey, fyear)

# save
saveRDS(comp, here::here("Data", "simulation_data_new_panel.rds"))
