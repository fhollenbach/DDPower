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




