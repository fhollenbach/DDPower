
library(groundhog)

pkgs <- c("here", 
          "tidyverse", 
          "viridis", 
          "patchwork", 
          "extrafont")
groundhog.day <- "2026-01-01"
groundhog.library(pkgs, groundhog.day)
####load arrow outside of groundhog, required newer version than on groundhog
###### arrow 23.0.0
library(arrow)
font_import()
loadfonts()


datasets <- open_dataset(here("..", "SimResults_firms_modelbased"))

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
  mse = (sum((est - true)^2, na.rm = TRUE))/length(est)
  return(mse)
}

my_rmse <- function(est, true) {
  rmse = sqrt((sum((est - true)^2, na.rm = TRUE))/length(est))
  return(rmse)
}



##### setting color palette to ensure same colors with and without cdh
color_all <-  viridis(9, option = "turbo")
estimators <- c("BJS", "BJS Loo", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE")

color_main <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise")))]
color_main_notwfe <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise", "TWFE")))]
color_all_notwfe <- color_all[-(which(estimators %in% c("TWFE")))]

color_nostata <- color_all[-(which(estimators %in% c("BJS Loo", "dCDH", "Stepwise")))]

###### main results, staggered dgp, 0 autocor

results <- datasets %>%
  filter(autocor == 0 & type == "staggered") %>% 
  collect()

mean(results$pretreatment_dv[results$DV == "rev_win"])
mean(results$sd_dv[results$DV == "rev_win"])

mean(results$pretreatment_dv[results$DV == "rev"])
mean(results$sd_dv[results$DV == "rev"])


for_cdh <- results %>%
  select(-c(estimate, std.error, method)) %>%
  group_by(iteration) %>%
  summarize_all(unique)

gc()
################ bring in Stata results
### read in all csv files from stata runs, name method, and level

path_stata <- here("..", "StataSim_firms_modelbased")


##### dCDH
# Get all CSV files with 'CdH' in filename
files_cdh_stata <- list.files(
  path = path_stata,
  pattern = "CdH.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
CdH_stata <- files_cdh_stata %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  mutate(method = "dCDH") %>%
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         estimate = `1...3`,
         std.error = `1...4`) %>% 
  left_join(for_cdh, by = "iteration") %>%
  select(estimate, std.error, method, true_att, att, pretreatment_dv, sd_dv, true_prop_treated, type, pct_effect, units, prop_treated, DV, autocor, iteration)


##### BJS Loo and BJS Stata (just to check results are the same)
files_bjs_stata <- list.files(
  path = path_stata,
  pattern = "imputation.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
bjs_stata <- files_bjs_stata %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         type = `1...3`,
         estimate = `1...4`,
         var = `1...5`) %>% 
  mutate(std.error = sqrt(var),
         method = ifelse(type == 0, "BJS Stata", "BJS LOO")) %>%
  select(-c(var, type)) %>% 
  left_join(for_cdh, by = "iteration") %>%
  select(estimate, std.error, method, true_att, att, pretreatment_dv, sd_dv, true_prop_treated, type, pct_effect, units, prop_treated, DV, autocor, iteration)
############ stepwise --- stata
files_step_stata <- list.files(
  path = path_stata,
  pattern = "stepwise.*\\.csv$",
  full.names = TRUE
)

# Read and row-bind
step_stata <- files_step_stata %>%
  lapply(read_csv, skip = 1) %>%
  bind_rows() %>% 
  select(-c(r1)) %>%
  rename(iteration = `1...2`,
         estimate = `1...3`,
         var = `1...4`) %>% 
  mutate(std.error = sqrt(var),
         method = "Stepwise") %>%
  select(-c(var)) %>% 
  left_join(for_cdh, by = "iteration") %>%
  select(estimate, std.error, method, true_att, att, pretreatment_dv, sd_dv, true_prop_treated, type, pct_effect, units, prop_treated, DV, autocor, iteration)




results <- bind_rows(results, CdH_stata, bjs_stata, step_stata) %>%
         mutate(method = if_else(method == "Mundlak", "ETWFE", method),
                method = if_else(method == "CSA", "CS", method),
                method = factor(method, levels = c("BJS", "BJS LOO", "BJS Stata", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE"), order = TRUE)) %>% 
  arrange(iteration)
#### compare R BJS to Stata BJS

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate - true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))

results_null <- results %>% 
  filter(att == 0) 

size_null <- results_null %>% 
  group_by(method, units, DV, type) %>%
  summarize(size_null = mean(p.value <= 0.05, na.rm = TRUE)) %>%
  ungroup()
  
summary_effect <- results %>%
  filter(pct_effect != 0) %>% 
  left_join(size_null, by = c("method", "units", "DV", "type")) %>%
  group_by(method, att, units, type, DV) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage, na.rm = TRUE),
            CI_size = mean(CI_size, na.rm = TRUE),
            mean_abs_error = mean(abs.error, na.rm = TRUE),
            mean_error = mean(error, na.rm = TRUE),
            mean_error_pct = mean(error_pct, na.rm = TRUE),
            share_sign = mean(significant, na.rm = TRUE),
            rmse = my_rmse(estimate, true_att), 
            size = unique(size_null)) %>%
  mutate(DV = factor(DV, levels = c("rev", "rev_win"), labels = c("Revenue", "Winsorized Revenue"), order = TRUE),
         type = factor(type, levels = c("standard", "staggered")),
         units = factor(units, levels = c(100, 250, 500, 1000, 1500), labels = c("100", "250", "500", "1000", "1500"), order = TRUE),
         att = factor(att, levels = c(0.05, 0.10, 0.15), labels = c("5%", "10%", "15%"), order = TRUE),
         delta_hat = qnorm(share_sign) - qnorm(size),
         R_hat = pnorm(delta_hat + qnorm(0.05)))
#### check that BJS and BJS Stata are the same
summary_effect$R_hat[summary_effect$method == "BJS"] 
summary_effect$R_hat[summary_effect$method == "BJS Stata"]

summary_by_sign_effect <- results %>%
  filter(p.value <= 0.05 & att != 0) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, att, units, type, DV) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(DV = factor(DV, levels = c("rev", "rev_win"), labels = c("Revenue", "Winsorized Revenue"), order = TRUE),
         type = factor(type, levels = c("standard", "staggered")),
         units = factor(units, levels = c(100, 250, 500, 1000, 1500), labels = c("100", "250", "500", "1000", "1500"), order = TRUE),
         att = factor(att, levels = c(0.05, 0.10, 0.15), labels = c("5%", "10%", "15%"), order = TRUE))

summary_effect_main <- summary_effect %>%
  filter(!(method %in% c("BJS Stata", "BJS LOO", "Stepwise")) & DV == "Winsorized Revenue" & type == "staggered")
##### power
pl <-  ggplot(data = summary_effect_main[summary_effect_main$att == "5%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main[summary_effect_main$att == "10%", ], aes(y = share_sign, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main[summary_effect_main$att == "15%", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_modelbased_winsorized.pdf"),
plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
#ggsave(here("..", "output", "power_modelbased_winsorized.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

################ now same plot but size adjusted power
summary_effect_mainnotwfe <- summary_effect_main %>%
  filter(!(method %in% c("TWFE")))
pl <-  ggplot(data = summary_effect_mainnotwfe[summary_effect_mainnotwfe$att == "5%" , ], aes(y = R_hat, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y=expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_mainnotwfe[summary_effect_mainnotwfe$att == "10%", ], aes(y = R_hat, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_mainnotwfe[summary_effect_mainnotwfe$att == "15%", ], aes(y = R_hat, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
########################################
######
##### exaggeration
summary_by_sign_main <- summary_by_sign_effect %>% 
  filter(!(method %in% c("BJS Stata", "BJS LOO", "Stepwise", "TWFE")) & DV == "Winsorized Revenue" & type == "staggered")

pl <-  ggplot(data = summary_by_sign_main[summary_by_sign_main$att == "5%", ], aes(y = exageration_ratio, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 5.5))


pm <-  ggplot(data = summary_by_sign_main[summary_by_sign_main$att == "10%", ], aes(y = exageration_ratio, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 5.5))


pr <-  ggplot(data = summary_by_sign_main[summary_by_sign_main$att == "15%", ], aes(y = exageration_ratio, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 5.5))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_main[summary_by_sign_main$att == "5%", ], aes(y = Pct_wrong_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 18)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_main[summary_by_sign_main$att == "10%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 18) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_main[summary_by_sign_main$att == "15%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 18)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### bias
pl <-  ggplot(data = summary_effect_main[summary_effect_main$att == "5%", ], aes(y = mean_error, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.03)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main[summary_effect_main$att == "10%", ], aes(y = mean_error, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.03) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main[summary_effect_main$att == "15%", ], aes(y = mean_error, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.03)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_modelbased_winsorized.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### coverage
pl <-  ggplot(data = summary_effect_main[summary_effect_main$att == "5%", ], aes(y = coverage, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_main[summary_effect_main$att == "10%", ], aes(y = coverage, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main[summary_effect_main$att == "15%", ], aes(y = coverage, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_winsorized_modelbased.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_main[summary_effect_main$att == "5%", ], aes(y = CI_size, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_main[summary_effect_main$att == "10%", ], aes(y = CI_size, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_main[summary_effect_main$att == "15%", ], aes(y = CI_size, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_winsorized_modelbased.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### rmse
pl <-  ggplot(data = summary_effect_main[summary_effect_main$att == "5%" , ], aes(y = rmse, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_main[summary_effect_main$att == "10%", ], aes(y = rmse, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_main[summary_effect_main$att == "15%", ], aes(y = rmse, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "rmse_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)





########################################
########################################
########################################
###### now compare windsorized to non-winsorized
summary_effect_win <- summary_effect %>%
  filter(!(method %in% c("BJS Stata", "BJS LOO", "Stepwise")) & att == "10%" & type == "staggered")

pl <-  ggplot(data = summary_effect_win[summary_effect_win$DV == "Revenue" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "Revenue", y = "Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_win[summary_effect_win$DV == "Winsorized Revenue", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "Winsorized Revenue", y= "", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "IBM Plex Sans Condensed")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "winsorized_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

################# now plot for manuscript with ATT = 10% comparing size adjusted and regular
pl <-  ggplot(data = summary_effect_mainnotwfe[summary_effect_mainnotwfe$att == "10%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 10%", y = "Power", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_mainnotwfe[summary_effect_mainnotwfe$att == "10%", ], aes(y = R_hat, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "", y= expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "Fira Sans")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_vs_sizeadj_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


#############################################
#############################################
#############################################
##### autocorrelation, staggered vs. non-staggered, winsorized vs. non-winsorized

results_auto <- datasets %>%
  collect() %>% 
  mutate(method = if_else(method == "Mundlak", "ETWFE", method),
         method = if_else(method == "CSA", "CS", method),
         method = factor(method, levels = c("BJS", "BJS LOO", "BJS Stata", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE"), order = TRUE)) %>% 
  arrange(iteration) %>% 
  filter(round(att,2) == 0.10 | round(att,2) == 0.00)


results_auto <- results_auto %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate - true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2),
         pct_effect_trans = exp(true_att)-1)

results_null <- results_auto %>% 
  filter(pct_effect == 0) 

size_null <- results_auto %>% 
  group_by(method, units, DV, type, autocor) %>%
  summarize(size_null = mean(p.value <= 0.05, na.rm = TRUE)) %>%
  ungroup()

summary_effect_auto <- results_auto %>%
  filter(att != 0) %>% 
  left_join(size_null, by = c("method", "units", "DV", "type", "autocor")) %>%
  group_by(method, att, units, type, DV, autocor) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage, na.rm = TRUE),
            CI_size = mean(CI_size, na.rm = TRUE),
            mean_abs_error = mean(abs.error, na.rm = TRUE),
            mean_error = mean(error, na.rm = TRUE),
            mean_error_pct = mean(error_pct, na.rm = TRUE),
            share_sign = mean(significant, na.rm = TRUE),
            rmse = my_rmse(estimate, true_att), 
            size = unique(size_null)) %>%
  mutate(DV = factor(DV, levels = c("rev", "rev_win"), labels = c("Revenue", "Winsorized Revenue"), order = TRUE),
         type = factor(type, levels = c("standard", "staggered")),
         units = factor(units, levels = c(100, 250, 500, 1000, 1500), labels = c("100", "250", "500", "1000", "1500"), order = TRUE),
         att = factor(att, levels = c(0.05, 0.10, 0.15), labels = c("5%", "10%", "15%"), order = TRUE),
         autocor = factor(autocor, levels = c(0.00, 0.35), labels = c("0", "0.35")),
         delta_hat = qnorm(share_sign) - qnorm(size),
         R_hat = pnorm(delta_hat + qnorm(0.05)))

summary_by_sign_effect_auto <- results_auto %>%
  filter(p.value <= 0.05 & att != 0) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, att, units, type, DV, autocor) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(DV = factor(DV, levels = c("rev", "rev_win"), labels = c("Revenue", "Winsorized Revenue"), order = TRUE),
         type = factor(type, levels = c("standard", "staggered")),
         units = factor(units, levels = c(100, 250, 500, 1000, 1500), labels = c("100", "250", "500", "1000", "1500"), order = TRUE),
         autocor = factor(autocor, levels = c(0.00, 0.35, 0.70), labels = c("0", "0.35", "0.7")),
         att = factor(att, levels = c(0.05, 0.10, 0.15), labels = c("5%", "10%", "15%"), order = TRUE))

summary_effect_auto_main <- summary_effect_auto %>%
  filter(DV == "Winsorized Revenue" & type == "staggered" & att == "10%")

##### power
pl <-  ggplot(data = summary_effect_auto_main[summary_effect_auto_main$autocor == "0" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "No add. Autocorrelation", y="Power", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_auto_main[summary_effect_auto_main$autocor == "0.35", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = expression(bold(rho ==  ~0.35)), y="", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())



p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "Fira Sans")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "autocor_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



#########################################
########################################
############## staggered, no staggered
summary_effect_auto_staggered <- summary_effect_auto %>%
  filter(DV == "Winsorized Revenue" & autocor == "0")

##### power
pl <-  ggplot(data = summary_effect_auto_staggered[summary_effect_auto_staggered$type == "staggered" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "Staggered & Dynamic", y="Power", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_auto_staggered[summary_effect_auto_staggered$type == "standard", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "Single Treatment", y="", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())



p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "Fira Sans")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "staggered_modelbased_winsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


####################################################################  
####################################################################  
#################################################################### 
#################################################################### 
################################## same plots but all non-winsorized
summary_effect_nonwin <- summary_effect %>%
  filter(!(method %in% c("BJS Stata", "BJS LOO", "Stepwise")) & DV == "Revenue" & type == "staggered")

##### power
pl <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "5%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "10%", ], aes(y = share_sign, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "15%", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
#ggsave(here("..", "output", "power_modelbased_winsorized.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

################ now same plot but size adjusted power
##### power
summary_effect_nonwin_notwfe <- summary_effect_nonwin %>%
  filter(!(method %in% c("TWFE")))
pl <-  ggplot(data = summary_effect_nonwin_notwfe[summary_effect_nonwin_notwfe$att == "5%" , ], aes(y = R_hat, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y=expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_nonwin_notwfe[summary_effect_nonwin_notwfe$att == "10%", ], aes(y = R_hat, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_nonwin_notwfe[summary_effect_nonwin_notwfe$att == "15%", ], aes(y = R_hat, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
########################################
######
##### exaggeration
summary_by_sign_nonwin <- summary_by_sign_effect %>% 
  filter(!(method %in% c("BJS Stata", "BJS LOO", "Stepwise", "TWFE")) & DV == "Revenue" & type == "staggered")

pl <-  ggplot(data = summary_by_sign_nonwin[summary_by_sign_nonwin$att == "5%", ], aes(y = exageration_ratio, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 6.5))


pm <-  ggplot(data = summary_by_sign_nonwin[summary_by_sign_nonwin$att == "10%", ], aes(y = exageration_ratio, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 6.5))


pr <-  ggplot(data = summary_by_sign_nonwin[summary_by_sign_nonwin$att == "15%", ], aes(y = exageration_ratio, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 6.5))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_nonwin[summary_by_sign_nonwin$att == "5%", ], aes(y = Pct_wrong_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 14)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_nonwin[summary_by_sign_nonwin$att == "10%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 14) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_nonwin[summary_by_sign_nonwin$att == "15%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 14)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


################# now plot for manuscript with ATT = 10% comparing size adjusted and regular
pl <-  ggplot(data = summary_effect_nonwin_notwfe[summary_effect_nonwin_notwfe$att == "10%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 10%", y = "Power", x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_nonwin_notwfe[summary_effect_nonwin_notwfe$att == "10%", ], aes(y = R_hat, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "", y= expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_main_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "IBM Plex Sans Condensed")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### bias
pl <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "5%", ], aes(y = mean_error, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.035)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "10%", ], aes(y = mean_error, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.035) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "15%", ], aes(y = mean_error, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.01, 0.035)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_modelbased_nonwinsorized.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### coverage
pl <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "5%", ], aes(y = coverage, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "10%", ], aes(y = coverage, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "15%", ], aes(y = coverage, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_nonwinsorized_modelbased.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "5%", ], aes(y = CI_size, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.6)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "10%", ], aes(y = CI_size, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.6) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "15%", ], aes(y = CI_size, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 0.6) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_nonwinsorized_modelbased.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### rmse
pl <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "5%" , ], aes(y = rmse, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "10%", ], aes(y = rmse, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_nonwin[summary_effect_nonwin$att == "15%", ], aes(y = rmse, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "rmse_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


########################################

#############################################
#############################################
#############################################
##### autocorrelation

summary_effect_auto_nonwin <- summary_effect_auto %>%
  filter(DV == "Revenue" & type == "staggered" & att == "10%")

##### power
pl <-  ggplot(data = summary_effect_auto_nonwin[summary_effect_auto_nonwin$autocor == "0" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "No add. Autocorrelation", y="Power", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_auto_nonwin[summary_effect_auto_nonwin$autocor == "0.35", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = expression(bold(rho ==  ~0.35)), y="", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())



p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "Fira Sans")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "autocor_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
#########################################
########################################
############## staggered, no staggered
summary_effect_auto_staggered_nonwin <- summary_effect_auto %>%
  filter(DV == "Revenue" & autocor == "0")

##### power
pl <-  ggplot(data = summary_effect_auto_staggered_nonwin[summary_effect_auto_staggered_nonwin$type == "staggered" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "Staggered \n& Dynamic", y="Power", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pr <-  ggplot(data = summary_effect_auto_staggered_nonwin[summary_effect_auto_staggered_nonwin$type == "standard", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "Single Treatment", y="", x = "") + scale_color_manual(name="Method", values = color_nostata) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())



p <- (pl + pr)/ wrap_elements(grid::textGrob("Number of Firms", gp = grid::gpar(fontsize = 24,  fontface = "bold", fontfamily = "Fira Sans")))/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.70, 0.05, 0.2))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans")) 
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "staggered_modelbased_nonwinsorized.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

######################################
######################################
######################################
####################### all estimators, winsorized

summary_effect_all <- summary_effect %>%
  filter(!(method %in% c("BJS Stata")) & DV == "Winsorized Revenue" & type == "staggered")
##### power
pl <-  ggplot(data = summary_effect_all[summary_effect_all$att == "5%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Power", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all[summary_effect_all$att == "10%", ], aes(y = share_sign, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all[summary_effect_all$att == "15%", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "power_modelbased_winsorized_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
#ggsave(here("..", "output", "power_modelbased_winsorized.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

################ now same plot but size adjusted power
##### power
summary_effect_allnotwfe <- summary_effect_all %>%
  filter(!(method %in% c("TWFE")))
pl <-  ggplot(data = summary_effect_allnotwfe[summary_effect_allnotwfe$att == "5%" , ], aes(y = R_hat, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
#pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y=expression(atop(bold("Size Adjusted Power"), bold(widehat(R)(alpha)))), x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_allnotwfe[summary_effect_allnotwfe$att == "10%", ], aes(y = R_hat, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
#pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_allnotwfe[summary_effect_allnotwfe$att == "15%", ], aes(y = R_hat, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
#pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "sizeadjusted_modelbased_winsorized_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)
########################################
######
##### exaggeration
summary_by_sign_all <- summary_by_sign_effect %>% 
  filter(!(method %in% c("BJS Stata", "TWFE")) & DV == "Winsorized Revenue" & type == "staggered")

pl <-  ggplot(data = summary_by_sign_all[summary_by_sign_all$att == "5%", ], aes(y = exageration_ratio, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 6))


pm <-  ggplot(data = summary_by_sign_all[summary_by_sign_all$att == "10%", ], aes(y = exageration_ratio, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all_notwfe)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 6))


pr <-  ggplot(data = summary_by_sign_all[summary_by_sign_all$att == "15%", ], aes(y = exageration_ratio, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", linewidth = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 6))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "exaggeration_modelbased_winsorized_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_all[summary_by_sign_all$att == "5%", ], aes(y = Pct_wrong_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 22)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_all[summary_by_sign_all$att == "10%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 22) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_all[summary_by_sign_all$att == "15%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all_notwfe) + ylim(0, 22)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "wrongsign_modelbased_winsorized_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### bias
pl <-  ggplot(data = summary_effect_all[summary_effect_all$att == "5%", ], aes(y = mean_error, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.035)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all[summary_effect_all$att == "10%", ], aes(y = mean_error, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.035) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all[summary_effect_all$att == "15%", ], aes(y = mean_error, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.01, 0.035)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "bias_modelbased_winsorized_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### coverage
pl <-  ggplot(data = summary_effect_all[summary_effect_all$att == "5%", ], aes(y = coverage, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_all[summary_effect_all$att == "10%", ], aes(y = coverage, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all[summary_effect_all$att == "15%", ], aes(y = coverage, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "coverage_winsorized_modelbased_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_all[summary_effect_all$att == "5%", ], aes(y = CI_size, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_all[summary_effect_all$att == "10%", ], aes(y = CI_size, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_all[summary_effect_all$att == "15%", ], aes(y = CI_size, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "CILength_winsorized_modelbased_allest.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### rmse
pl <-  ggplot(data = summary_effect_all[summary_effect_all$att == "5%" , ], aes(y = rmse, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5%", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.13)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_all[summary_effect_all$att == "10%", ], aes(y = rmse, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10%", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.13) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_all[summary_effect_all$att == "15%", ], aes(y = rmse, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15%", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.13)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("Fira Sans"))
ggsave(here("..", "..", "..", "Dropbox-CBS", "Florian Hollenbach", "Apps", "Overleaf", "DD power", "Figures", "rmse_modelbased_winsorized_allest.pdf"),
       plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


