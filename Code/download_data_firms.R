
library(groundhog)

pkgs <- c("here", "RPostgres", "DBI", "tidyverse")
groundhog.day <- "2025-09-01"
groundhog.library(pkgs, groundhog.day)


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
comp <- tbl(wrds_con, sql("SELECT gvkey, fyear, revt, datafmt, popsrc, consol,indfmt, fic, sich FROM comp.funda")) %>%
  # filter the data - between 1979 and 2015, non-missing assets, and in the US
  filter(indfmt == 'INDL' & datafmt == 'STD' & popsrc == 'D' & consol == 'C' & !is.na(fyear) & 
           fyear %>% between(2001, 2020) &  fic == "USA") %>% 
  # drop if missing at or revt
  filter(is.na(revt) == FALSE) %>% 
  mutate(gvkey = as.numeric(gvkey)) %>% 
  group_by(gvkey) %>% 
  mutate(min_year = min(fyear),
         year_dif = fyear - lag(fyear, 1, order = fyear), ### check whether time series consec
         time = (fyear - min_year) + 1, ### time counter
         max_time = max(time), ### time series length by firm
         max_gap = max(year_dif, na.rm = TRUE),
         max_year = max(fyear)) %>% 
  collect() %>% 
  arrange(gvkey, fyear) %>% 
  filter(min_year == 2001 & max_gap == 1 & fyear <= 2020 & max_year >=2020) %>% 
  select(-c(datafmt, popsrc, consol, indfmt, fic, min_year, year_dif, max_time, max_gap, max_year)) ### don't need these anymore
  


length(unique(comp$gvkey))
### 2176 firms with 20 years of data
table(comp$time)
# download comp header which has more info on location etc
comp_header <- tbl(wrds_con, sql("SELECT * FROM crsp.comphead")) %>% 
  mutate(gvkey = as.numeric(gvkey)) %>% 
  collect()


# merge in state of incorporation and industry information
comp_dat <- comp %>% 
  left_join(comp_header %>% select(gvkey, incorp, sic, state, gsector)) %>% 
  # clean up SIC code - use historical sic code if available, if not use header sic
  mutate(sich = coalesce(sich, sic))


# winsorize 
comp_dat <- comp_dat %>% 
  ungroup() %>% 
  mutate(revt_win = case_when(is.na(revt) ~ NA_real_,
                          revt < quantile(revt, 0.15, na.rm = TRUE) ~ quantile(revt, 0.15, na.rm = TRUE),
                          revt > quantile(revt, 0.85, na.rm = TRUE) ~ quantile(revt, 0.85, na.rm = TRUE),
                          TRUE ~ revt),
         revt = case_when(is.na(revt) ~ NA_real_,
                          revt < quantile(revt, 0.025, na.rm = TRUE) ~ quantile(revt, 0.025, na.rm = TRUE),
                          revt > quantile(revt, 0.975, na.rm = TRUE) ~ quantile(revt, 0.975, na.rm = TRUE),
                          TRUE ~ revt),
         log_rev = log(revt),
         log_rev_win = log(revt_win)) %>% 
  arrange(gvkey, fyear) %>% 
  group_by(gvkey) %>% 
  mutate(gvkey = as.numeric(gvkey)) %>% 
  ungroup() %>% 
  filter(fyear > 2000) %>% 
  mutate(time = (fyear - 2000)) ### time counter

length(unique(comp_dat$gvkey))
table(comp_dat$time)


# save
# downloaded on Monday, Sept. 15, 2025
saveRDS(comp_dat, here::here("Data", "simulation_data_new_panel_092025.rds"))
