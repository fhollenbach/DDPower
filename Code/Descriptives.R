library(tidyverse)
library(fixest)


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

####### doing some descriptive stats on a single of the simulated datasets in each sim

########
#### load data for a single simulation 
#### state level drugs

data <- read_csv(here("..", "SimData_State_opioids", "Data_state_6000.csv"))
pdat = panel(data, ~ StateCode + Year)
### lag dv model
m1 <- feols(y ~ l(y, 1) | StateCode + Year, data = pdat)
summary(m1)
### coefficien of 0.94 on lag

#### density
p <- ggplot(data, aes(x = y)) +
  geom_density() + theme_mfx() +
  labs(subtitle = "", y="Density", x = "Simulated Death Rate (per 100,000)")
ggsave(here("..", "output", "Density_deathrate_state"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


data <- read_csv(here("..", "SimData_County_opioids", "Data_county_6000.csv"))
pdat = panel(data, ~ CountyCode + Year)
### lag dv model
m1 <- feols(y ~ l(y, 1) | CountyCode + Year, data = pdat)
summary(m1)
### coefficien of 0.397 on lag

#### density
p <- ggplot(data, aes(x = y)) +
  geom_density() + theme_mfx() +
  labs(subtitle = "", y="Density", x = "Simulated Death Rate (per 100,000)")
ggsave(here("..", "output", "Density_deathrate_county"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



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

