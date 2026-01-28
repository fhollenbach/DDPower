#### file to run firm sims
library(groundhog)

pkgs <- c("e1071", "tidyverse", "here", "fixest", "modelsummary", "tinytable")
groundhog.day <- "2025-11-01"
groundhog.library(pkgs, groundhog.day)

theme_mfx <- function() {
  theme_minimal(base_family = "IBM Plex Sans Condensed") +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(size=20),
          plot.subtitle = element_text(size=20),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size = 24, face = "bold"),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(angle = 90, size = 24, face = "bold"),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position="bottom",
          legend.text=element_text(size=20))
}


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


mean(comp$log_rev)
sqrt(var(comp$log_rev))
skewness(comp$log_rev)

mean(comp$log_rev_win)
sqrt(var(comp$log_rev_win))
skewness(comp$log_rev_win)


mean(data$log_GDP)
sqrt(var(data$log_GDP))
skewness(data$log_GDP)
####### doing some descriptive stats on a single of the simulated datasets in each sim

########
#### load data for a single simulation 
#### state level gdp

data <- read_csv(here("..", "SimData_State_modelbased", "Data_state_6000.csv"))
pdat = panel(data, ~ StateCode + Year)
### lag dv model
m1 <- feols(y ~ l(y, 1) | StateCode + Year, data = pdat)
summary(m1)
### coefficien of 0.81709 on lag
mean(data$y)#7.398619
sqrt(var(data$y))# 1.165633
skewness(data$y) #0.5304769

##### firm winsorized rev
data <- read_csv(here("..", "SimData_firms_modelbased", "Data_firms_staggered19501.csv"))
pdat = panel(data, ~ firm_id + time)
### lag dv model
m1 <- feols(y_log ~ l(y_log, 1) | firm_id + time, data = pdat)
summary(m1)
### coefficien of 0.887096 on lag
mean(data$y_log)#6.026546
sqrt(var(data$y_log))#2.043287
skewness(data$y_log)#0.1878448

##### firm rev
data <- read_csv(here("..", "SimData_firms_modelbased", "Data_firms_staggered9501.csv"))
pdat = panel(data, ~ firm_id + time)
### lag dv model
m1 <- feols(y_log ~ l(y_log, 1) | firm_id + time, data = pdat)
summary(m1)
### coefficien of 0.860507 on lag
mean(data$y_log)#6.026546
sqrt(var(data$y_log))#2.871776
skewness(data$y_log)#-0.3473807




#######
pct_effect <- c(0.025, 0.05, 0.10, 0.15)
units <- c(100, 250, 1000, 2000)
prop_treated <- c(0.1, 0.2, 0.4, 0.6)
iterations <- 500
DV = c("roa", "log_rev")

simulation_para <- bind_rows(replicate(iterations, expand_grid(pct_effect, units, prop_treated, DV), simplify = FALSE)) %>%
  arrange(units, pct_effect, prop_treated, DV)
#names(simulation_para) <- c("notreated", "pct_effect")
simulation_para$iteration <- 1:dim(simulation_para)[1]

simulation_para[simulation_para$DV == 'log_rev' & simulation_para$pct_effect == 0.15 & simulation_para$units == 1000 & simulation_para$prop_treated == 0.4, ] 

simulation_para[simulation_para$DV == 'roa' & simulation_para$pct_effect == 0.15 & simulation_para$units == 1000 & simulation_para$prop_treated == 0.4, ] 

####### logrev for a single simulation

data_logrev <- read_csv(here("..", "SimData_firms", "Data_firms_46001.csv"))
pdat = panel(data_logrev, ~ gvkey + time)
### lag dv model
m1 <- feols(y ~ l(y, 1) | gvkey + time, data = pdat)
summary(m1)
### coefficien of 0.90 on lag
summary(data_logrev)
mean(data_logrev$y)
sd(data_logrev$y)
#### density
p <- ggplot(data_logrev, aes(x = y)) +
  geom_density() + theme_mfx() +
  labs(subtitle = "", y="Density", x = "Simulated Revenue (ln)")
ggsave(here("..", "output", "Density_logrev"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




####### roa for a single simulation

data_roa <- read_csv(here("..", "SimData_firms", "Data_firms_46501.csv"))
pdat = panel(data_roa, ~ gvkey + time)
### lag dv model
m1 <- feols(y ~ l(y, 1) | gvkey + time, data = pdat)
summary(m1)
### coefficien of 0.494 on lag
summary(data_roa)
mean(data_roa$y)
sd(data_roa$y)
#### density
p <- ggplot(data_roa, aes(x = y)) +
  geom_density() + theme_mfx() +
  labs(subtitle = "", y="Density", x = "Simulated ROA")
ggsave(here("..", "output", "Density_roa"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

