#### file to run state/county sims
library(groundhog)

pkgs <- c("here", 
          "tidyverse", 
          "modelsummary", 
          "viridis", 
          "patchwork", 
          "extrafont")
groundhog.day <- "2026-01-01"
groundhog.library(pkgs, groundhog.day, cores=10)
####load arrow outside of groundhog, required newer version than on groundhog
###### arrow 23.0.0
library(arrow)
font_import()
loadfonts()

datasets <- open_dataset("~/Documents/GitHub/SimResults_statecounty_modelbased")


theme_mfx <- function() {
  theme_minimal(base_family = "Fira Sans") +
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



### Mean Squared Error
mean_sq = function(est, true){
  mse = (sum((est - true)^2))/length(est)
  return(mse)
}

my_rmse <- function(est, true) {
  rmse = sqrt((sum((est - true)^2))/length(est))
  return(rmse)
}




##### setting color palette to ensure same colors with and without cdh
color_all <-  viridis(9, option = "turbo")
estimators <- c("BJS", "BJS Loo", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE")


color_main <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise")))]
color_main_notwfe <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise", "TWFE")))]
color_all_notwfe <- color_all[-(which(estimators %in% c("TWFE")))]

color_nostata <- color_all[-(which(estimators %in% c("BJS Loo", "dCDH", "Stepwise")))]



########################################################################
########################################################################
##################################################################
################# 

results <- datasets %>% 
  collect()


for_cdh <- results %>%
  select(-c(estimate, std.error, method)) %>%
  group_by(iteration, level) %>%
  summarize_all(unique)
#

### read in all csv files from stata runs, name method, and level

path_state <- here("..", "StataSim_state_modelbased")
path_county <- here("..", "StataSim_county_modelbased")

# Get all CSV files with 'CdH' in filename
files_cdh_state <- list.files(
  path = path_state,
  pattern = "CdH.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
CdH_state <- files_cdh_state %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  mutate(method = "dCDH", 
         level = "State") %>%
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         estimate = `1...3`,
         std.error = `1...4`)
### county
files_cdh_county <- list.files(
  path = path_county,
  pattern = "CdH.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
CdH_county <- files_cdh_county %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  mutate(method = "dCDH", 
         level = "County") %>%
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         estimate = `1...3`,
         std.error = `1...4`)

cdh <- bind_rows(CdH_state, CdH_county)

cdh <- for_cdh %>%
  left_join(cdh, by = c("iteration", "level")) %>%
  select(estimate, std.error, method, true_att, true_att_log_gdp, pretreatment_log_gdp, sd_log_gdp, level, numtreated, pct_effect, iteration)

##### BJS Loo and BJS Stata (just to check results are the same)
# Get all CSV files with 'imputation' in filename
files_bjs_state <- list.files(
  path = path_state,
  pattern = "imputation.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
bjs_state <- files_bjs_state %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         type = `1...3`,
         estimate = `1...4`,
         var = `1...5`) %>% 
  mutate(std.error = sqrt(var),
         method = ifelse(type == 0, "BJS Stata", "BJS LOO"), 
         level = "State") %>%
  select(-c(var, type)) 
### county
files_bjs_county <- list.files(
  path = path_county,
  pattern = "imputation.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
bjs_county <- files_bjs_county %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         type = `1...3`,
         estimate = `1...4`,
         var = `1...5`) %>% 
  mutate(std.error = sqrt(var),
         method = ifelse(type == 0, "BJS Stata", "BJS LOO"), 
         level = "County") %>%
  select(-c(var, type)) 

bjs_stata <- bind_rows(bjs_state, bjs_county)

bjs_stata <- for_cdh %>%
  left_join(bjs_stata, by = c("iteration", "level")) %>%
  select(estimate, std.error, method, true_att, true_att_log_gdp, pretreatment_log_gdp, sd_log_gdp, level, numtreated, pct_effect, iteration)




############ stepwise --- stata
# Get all CSV files with 'stepwise' in filename
files_step_state <- list.files(
  path = path_state,
  pattern = "stepwise.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
step_state <- files_step_state %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         estimate = `1...3`,
         var = `1...4`) %>% 
  mutate(std.error = sqrt(var),
         method = "Stepwise", 
         level = "State") %>%
  select(-c(var)) 
### county
files_step_county <- list.files(
  path = path_county,
  pattern = "stepwise.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
step_county <- files_step_county %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         estimate = `1...3`,
         var = `1...4`) %>% 
  mutate(std.error = sqrt(var),
         method = "Stepwise", 
         level = "County") %>%
  select(-c(var))

stepwise <- bind_rows(step_state, step_county)

stepwise <- for_cdh %>%
  left_join(stepwise, by = c("iteration", "level")) %>%
  select(estimate, std.error, method, true_att, true_att_log_gdp, pretreatment_log_gdp, sd_log_gdp, level, numtreated, pct_effect, iteration)


results <- bind_rows(results, cdh, bjs_stata, stepwise) %>% 
  mutate(method = case_when(method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "BJS LOO", "BJS Stata", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE"), order = TRUE)) %>% 
  arrange(iteration)
### clear some workspace
rm(list=setdiff(ls(), c("results", "mean_sq", "my_rmse", "theme_mfx", "color_main", "color_main_notwfe", "color_nostata")))
gc()

results <- results %>%
  mutate(conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att_log_gdp) >= conf.low & (true_att_log_gdp) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att_log_gdp),
         error_pct = error/true_att_log_gdp,
         abs.error = abs(error), 
         true_att_rounded = round(true_att_log_gdp, 3))

##### overall average of pre-treatment 
mean(results$pretreatment_log_gdp[results$method == "TWFE"])
mean(results$sd_log_gdp[results$method == "TWFE"])


#### for calculation of size adj. power
results_null <- results %>% 
  filter(true_att_rounded == 0) 

size_null <- results_null %>% 
  group_by(method, numtreated, level) %>%
  summarize(size_null = mean(p.value <= 0.05)) %>%
  ungroup()


summary_effect <- results %>%
  filter(true_att_rounded != 0) %>% 
  left_join(size_null, by = c("method", "numtreated", "level")) %>%
  group_by(method, true_att_rounded, numtreated, level) %>%
  summarize(mean_true = mean(true_att_log_gdp),
            mean_pretreat = mean(pretreatment_log_gdp),
            coverage = mean(CI_coverage, na.rm = TRUE),
            CI_size = mean(CI_size, na.rm = TRUE),
            mean_abs_error = mean(abs.error, na.rm = TRUE),
            mean_error = mean(error, na.rm = TRUE),
            mean_error_pct = mean(error_pct, na.rm = TRUE),
            share_sign = mean(Significant, na.rm = TRUE),
            rmse = my_rmse(estimate, true_att_log_gdp), 
            size = unique(size_null)) %>%
  mutate(numtreated = factor(numtreated, levels = c(6, 12, 24), labels = c(6, 12, 24), order = TRUE),
         level = factor(level, levels = c("State", "County")),
         att = factor(true_att_rounded, levels = c(0.01, 0.024, 0.048), labels = c("1%", "2.5%", "5%"), order = TRUE),
         delta_hat = qnorm(share_sign) - qnorm(size),
         R_hat = pnorm(delta_hat + qnorm(0.05)))


## lowest power with 6 treated units
min(summary_effect$share_sign[summary_effect$numtreated == 6 & summary_effect$level == "State"])

summary_by_sign <- results %>%
  filter(p.value <= 0.05 & true_att_rounded != 0) %>% 
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att_log_gdp),
         exageration_ratio = if_else(estimate < 0, NA, exageration_ratio)) %>%
  group_by(method, true_att_rounded, numtreated, level) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(numtreated = factor(numtreated, levels = c(6, 12, 24), labels = c(6, 12, 24), order = TRUE),
         level = factor(level, levels = c("State", "County")),
         att = factor(true_att_rounded, levels = c(0.01, 0.024, 0.048), labels = c("1%", "2.5%", "5%"), order = TRUE))


###### level
summary_effect_main_state <- summary_effect %>%
  filter(level == "State" & !(method %in% c("BJS Stata", "BJS LOO", "Stepwise")))

summary_effect_main_county <- summary_effect %>%
  filter(level == "County" & !(method %in% c("BJS Stata", "BJS LOO", "Stepwise")))

##### states
pl <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "1%" , ], aes(y = share_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "2.5%", ], aes(y = share_sign, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "5%", ], aes(y = share_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_gdp_state_main.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


### county
pl <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "1%" , ], aes(y = share_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "2.5%", ], aes(y = share_sign, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "5%", ], aes(y = share_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_gdp_county_main.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



############ size adjusted power (no twfe)
### state
summary_effect_mainnotwfe_state <- summary_effect_main_state %>%
  filter(!(method %in% c("TWFE")))

pl <-  ggplot(data = summary_effect_mainnotwfe_state[summary_effect_mainnotwfe_state$att == "1%" , ], aes(y = R_hat, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y=expression(atop(bold("Size Adjusted Power"), widehat(R)(alpha))), x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_mainnotwfe_state[summary_effect_mainnotwfe_state$att == "2.5%", ], aes(y = R_hat, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_mainnotwfe_state[summary_effect_mainnotwfe_state$att == "5%", ], aes(y = R_hat, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_gdp_state_main.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


### county
summary_effect_mainnotwfe_county <- summary_effect_main_county %>%
  filter(!(method %in% c("TWFE")))

pl <-  ggplot(data = summary_effect_mainnotwfe_county[summary_effect_mainnotwfe_county$att == "1%" , ], aes(y = R_hat, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y=expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_mainnotwfe_county[summary_effect_mainnotwfe_county$att == "2.5%", ], aes(y = R_hat, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_mainnotwfe_county[summary_effect_mainnotwfe_county$att == "5%", ], aes(y = R_hat, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_gdp_county_main.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### coverage 
##state
pl <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "1%" , ], aes(y = coverage, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y= "Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "2.5%", ], aes(y = coverage, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "5%", ], aes(y = coverage, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_gdp_state_main.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

### county
pl <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "1%" , ], aes(y = coverage, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y= "Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "2.5%", ], aes(y = coverage, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "5%", ], aes(y = coverage, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_gdp_county_main.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### CI length
pl <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "1%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.08) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "2.5%", ], aes(y = CI_size, x = numtreated , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.08)   #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "5%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.08)  #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_gdp_state_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### CI length
pl <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "1%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.08) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "2.5%", ], aes(y = CI_size, x = numtreated , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.08)   #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "5%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.08)  #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_gdp_county_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



###bias 

pl <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "1%", ], aes(y = mean_error, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "2.5%", ], aes(y = mean_error, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "5%", ], aes(y = mean_error, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_gdp_state_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


#####county
pl <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "1%", ], aes(y = mean_error, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "2.5%", ], aes(y = mean_error, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "5%", ], aes(y = mean_error, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_gdp_county_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### rmse
pl <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "1%", ], aes(y = rmse, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "2.5%", ], aes(y = rmse, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_state[summary_effect_main_state$att == "5%", ], aes(y = rmse, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures",  "rmse_gdp_state_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### rmse
pl <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "1%", ], aes(y = rmse, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "2.5%", ], aes(y = rmse, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main_county[summary_effect_main_county$att == "5%", ], aes(y = rmse, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures",  "rmse_gdp_county_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration & wrong sign
###### level
summary_by_sign_state <- summary_by_sign %>%
  filter(level == "State" & !(method %in% c("BJS Stata", "BJS LOO", "Stepwise", "TWFE")))

summary_by_sign_county <- summary_by_sign %>%
  filter(level == "County" & !(method %in% c("BJS Stata", "BJS LOO", "Stepwise", "TWFE")))



pl <-  ggplot(data = summary_by_sign_state[summary_by_sign_state$att == "1%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 3.5))
pl

pm <-  ggplot(data = summary_by_sign_state[summary_by_sign_state$att == "2.5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0.9, 3.5))
pm

pr <-  ggplot(data = summary_by_sign_state[summary_by_sign_state$att == "5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 3.5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_gdp_state_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




pl <-  ggplot(data = summary_by_sign_county[summary_by_sign_county$att == "1%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 3.5))
pl

pm <-  ggplot(data = summary_by_sign_county[summary_by_sign_county$att == "2.5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0.9, 3.5))
pm

pr <-  ggplot(data = summary_by_sign_county[summary_by_sign_county$att == "5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 3.5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_gdp_county_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_state[summary_by_sign_state$att == "1%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 8.5))
pl

pm <-  ggplot(data = summary_by_sign_state[summary_by_sign_state$att == "2.5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0, 8.5))
pm

pr <-  ggplot(data = summary_by_sign_state[summary_by_sign_state$att == "5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 8.5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_gdp_state_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




pl <-  ggplot(data = summary_by_sign_county[summary_by_sign_county$att == "1%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 12.5))
pl

pm <-  ggplot(data = summary_by_sign_county[summary_by_sign_county$att == "2.5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_main_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0, 12.5))
pm

pr <-  ggplot(data = summary_by_sign_county[summary_by_sign_county$att == "5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 12.5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_gdp_county_main.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

###################################################
###################################################
###################################################
###################################################
###################################################
######## all estimators

##### setting color palette to ensure same colors with and without cdh
color_all <-  viridis(9, option = "turbo")
estimators <- c("BJS", "BJS Loo", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE")


color_main <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise")))]
color_main_notwfe <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise", "TWFE")))]
color_all_notwfe <- color_all[-(which(estimators %in% c("TWFE")))]

color_nostata <- color_all[-(which(estimators %in% c("BJS Loo", "dCDH", "Stepwise")))]



###### level
summary_effect_all_state <- summary_effect %>%
  filter(level == "State" & !(method %in% c("BJS Stata")))

summary_effect_all_county <- summary_effect %>%
  filter(level == "County" & !(method %in% c("BJS Stata")))

##### states
pl <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "1%" , ], aes(y = share_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Power", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "2.5%", ], aes(y = share_sign, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "5%", ], aes(y = share_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_gdp_state_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


### county
pl <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "1%" , ], aes(y = share_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Power", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "2.5%", ], aes(y = share_sign, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "5%", ], aes(y = share_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_gdp_county_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



############ size adjusted power (no twfe)
summary_effect_allnotwfe_state <- summary_effect_all_state %>% 
  filter(method != "TWFE")

### state
pl <-  ggplot(data = summary_effect_allnotwfe_state[summary_effect_allnotwfe_state$att == "1%" , ], aes(y = R_hat, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y=expression(atop(bold("Size Adjusted Power"), widehat(R)(alpha))), x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_allnotwfe_state[summary_effect_allnotwfe_state$att == "2.5%", ], aes(y = R_hat, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_allnotwfe_state[summary_effect_allnotwfe_state$att == "5%", ], aes(y = R_hat, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_gdp_state_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


### county
summary_effect_allnotwfe_county <- summary_effect_all_county %>% 
  filter(method != "TWFE")
pl <-  ggplot(data = summary_effect_allnotwfe_county[summary_effect_allnotwfe_county$att == "1%" , ], aes(y = R_hat, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y=expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_allnotwfe_county[summary_effect_allnotwfe_county$att == "2.5%", ], aes(y = R_hat, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_allnotwfe_county[summary_effect_allnotwfe_county$att == "5%", ], aes(y = R_hat, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_gdp_county_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### coverage 
##state
pl <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "1%" , ], aes(y = coverage, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y= "Coverage", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "2.5%", ], aes(y = coverage, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "5%", ], aes(y = coverage, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_gdp_state_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

### county
pl <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "1%" , ], aes(y = coverage, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y= "Coverage", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "2.5%", ], aes(y = coverage, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "5%", ], aes(y = coverage, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_gdp_county_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### CI length
pl <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "1%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 0.09) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "2.5%", ], aes(y = CI_size, x = numtreated , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all) + ylim(0, 0.09)   #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "5%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 0.09)  #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_gdp_state_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### CI length
pl <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "1%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 0.09) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "2.5%", ], aes(y = CI_size, x = numtreated , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all) + ylim(0, 0.09)   #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "5%", ], aes(y = CI_size, x = numtreated, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 0.09)  #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_gdp_county_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



###bias 

pl <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "1%", ], aes(y = mean_error, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "2.5%", ], aes(y = mean_error, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "5%", ], aes(y = mean_error, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_gdp_state_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


#####county
pl <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "1%", ], aes(y = mean_error, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "2.5%", ], aes(y = mean_error, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "5%", ], aes(y = mean_error, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_gdp_county_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### rmse
pl <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "1%", ], aes(y = rmse, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "2.5%", ], aes(y = rmse, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_state[summary_effect_all_state$att == "5%", ], aes(y = rmse, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "rmse_gdp_state_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### rmse
pl <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "1%", ], aes(y = rmse, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "2.5%", ], aes(y = rmse, x = numtreated , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all_county[summary_effect_all_county$att == "5%", ], aes(y = rmse, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "rmse_gdp_county_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration & wrong sign
###### level
summary_by_sign_all_state <- summary_by_sign %>%
  filter(level == "State" & !(method %in% c("BJS Stata", "TWFE")))

summary_by_sign_all_county <- summary_by_sign %>%
  filter(level == "County" & !(method %in% c("BJS Stata", "TWFE")))



pl <-  ggplot(data = summary_by_sign_all_state[summary_by_sign_all_state$att == "1%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 3.5))
pl

pm <-  ggplot(data = summary_by_sign_all_state[summary_by_sign_all_state$att == "2.5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0.9, 3.5))
pm

pr <-  ggplot(data = summary_by_sign_all_state[summary_by_sign_all_state$att == "5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 3.5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_gdp_state_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




pl <-  ggplot(data = summary_by_sign_all_county[summary_by_sign_all_county$att == "1%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 3.5))
pl

pm <-  ggplot(data = summary_by_sign_all_county[summary_by_sign_all_county$att == "2.5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0.9, 3.5))
pm

pr <-  ggplot(data = summary_by_sign_all_county[summary_by_sign_all_county$att == "5%", ], aes(y = exageration_ratio, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 3.5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_gdp_county_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### pct wrong sign  

pl <-  ggplot(data = summary_by_sign_all_state[summary_by_sign_all_state$att == "1%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 9))
pl

pm <-  ggplot(data = summary_by_sign_all_state[summary_by_sign_all_state$att == "2.5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0, 9))
pm

pr <-  ggplot(data = summary_by_sign_all_state[summary_by_sign_all_state$att == "5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 9))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_gdp_state_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




pl <-  ggplot(data = summary_by_sign_all_county[summary_by_sign_all_county$att == "1%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 1%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 13))
pl

pm <-  ggplot(data = summary_by_sign_all_county[summary_by_sign_all_county$att == "2.5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 2.5%", y="", x = "Number of Treated States") + scale_color_manual(name="Method", values = color_all_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0, 13))
pm

pr <-  ggplot(data = summary_by_sign_all_county[summary_by_sign_all_county$att == "5%", ], aes(y = Pct_wrong_sign, x = numtreated, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 5%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 13))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_gdp_county_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


