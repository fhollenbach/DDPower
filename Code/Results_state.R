library(here)
library(tidyverse)
library(modelsummary) #### use modelsummary
library(viridis)
library(patchwork)     # For combining plots
library(arrow)

datasets <- open_dataset("~/Documents/GitHub/SimResults_State_opioids")


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

color_main <- color_all[-(which(estimators %in% c("BJS Loo", "Stepwise", "TWFE")))]
color_nocdh <- color_all[-(which(estimators %in% c("BJS Loo", "dCDH", "Stepwise", "TWFE")))]




########################################################################
########################################################################
##################################################################
################# varying prop treated and effect size

results <- datasets %>% 
  filter(t == -100 & method != "TWFE") %>% 
  collect()


for_cdh <- results %>%
  select(-c(t, estimate, std.error, method)) %>%
  group_by(iteration) %>%
  summarize_all(unique)
#
cdh1 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata1.csv"))[-1,]
cdh2 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata2.csv"))[-1,]
cdh3 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata3.csv"))[-1,]
cdh4 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata4.csv"))[-1,]

cdh <- bind_rows(cdh1, cdh2, cdh3, cdh4) %>%
  mutate(method = "CdH") %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         std.error = c4) %>%
  filter(t == -100)

cdh <- for_cdh %>%
  left_join(cdh, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, true_att_deathrate, pretreatment_deathrate, sd_deathrate, notreated, pct_effect, iteration)

results <- bind_rows(results, cdh) %>% ### no imp stata, no imp loo for now
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA"), order = TRUE)) %>% 
  arrange(iteration)


results <- results %>%
  mutate(conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att_deathrate) >= conf.low & (true_att_deathrate) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att_deathrate),
         error_pct = error/true_att_deathrate,
         abs.error = abs(error), 
         true_att_rounded = round(true_att_deathrate, 2))


summary_by_effect_prop <- results %>%
  group_by(method, pct_effect, notreated) %>%
  summarize(mean_true = mean(true_att_deathrate),
            mean_true_pct = mean(true_att),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att_deathrate)) %>%
  mutate(notreated = factor(notreated, levels = c(6, 12, 24), labels = c(6, 12, 24), order = TRUE),
         ate = factor(round(pct_effect, 3), levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE),
         mean_pct = round(mean_true_pct, 2))

## lowest power with 6 treated units
min(summary_by_effect_prop$share_sign[summary_by_effect_prop$notreated == 6])

summary_by_sign_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate > 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att_deathrate),
         exageration_ratio = if_else(estimate > 0, NA, exageration_ratio)) %>%
  group_by(method, pct_effect, notreated) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(notreated = factor(notreated, levels = c(6, 12, 24), labels = c(6, 12, 24), order = TRUE),
         ate = factor(round(pct_effect, 3), levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

### max wrong sign
max(summary_by_sign_prop$Pct_wrong_sign)
######


##### pct wrong by method for 6 treated
print(summary_by_sign_prop[summary_by_sign_prop$notreated == 6, c("method", "Pct_wrong_sign", "ate")], n= 24)
##### exaggeration by method for 6 treated
print(summary_by_sign_prop[summary_by_sign_prop$notreated == 6, c("method", "exageration_ratio", "ate")], n= 24)
##### exaggeration by method for 12 treated
print(summary_by_sign_prop[summary_by_sign_prop$notreated == 12, c("method", "exageration_ratio", "ate")], n= 24)


##### coverage
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = coverage, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = coverage, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = coverage, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = CI_size, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated States", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 13)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = CI_size, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated States", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_main) + ylim(0, 13) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = CI_size, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated States", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 13) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = share_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = share_sign, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = share_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### rmse
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = rmse, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 4)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = rmse, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 4) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = rmse, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 4)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = mean_error, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.3, 0.3)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = mean_error, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "Mean Error") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.3, 0.3) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = mean_error, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.3, 0.3)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### bias pct
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = mean_error_pct, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.1, 0.3)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = mean_error_pct, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "Mean Error") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.1, 0.3) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = mean_error_pct, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.1, 0.3)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_pct_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### exaggeration
pl <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 6, ], aes(y = exageration_ratio, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 5))
pl

pm <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 12, ], aes(y = exageration_ratio, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_main)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0.9, 5))
pm

pr <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 24, ], aes(y = exageration_ratio, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 5))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 6, ], aes(y = Pct_wrong_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 19)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 
pl

pm <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 12,  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_main) + ylim(0, 19) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  
pm

pr <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 24,  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 19)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_prop_state_opioids.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)







########################################################################
########################################################################
##################################################################
################# Alll estimators appendix

results <- datasets %>% 
  filter(t == -100) %>% 
  collect()


for_cdh <- results %>%
  select(-c(t, estimate, std.error, method)) %>%
  group_by(iteration) %>%
  summarize_all(unique)
#
cdh1 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata1.csv"))[-1,]
cdh2 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata2.csv"))[-1,]
cdh3 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata3.csv"))[-1,]
cdh4 <- read_csv(here("..", "StataSim_state", "sim_CdH_stata4.csv"))[-1,]

cdh <- bind_rows(cdh1, cdh2, cdh3, cdh4) %>%
  mutate(method = "CdH") %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         std.error = c4) %>%
  filter(t == -100)

cdh <- for_cdh %>%
  left_join(cdh, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, true_att_deathrate, pretreatment_deathrate, sd_deathrate, notreated, pct_effect, iteration)


imp1 <- read_csv(here("..", "StataSim_state", "sim_imputation_stata1.csv"))[-1,]
imp2 <- read_csv(here("..", "StataSim_state", "sim_imputation_stata2.csv"))[-1,]
imp3 <- read_csv(here("..", "StataSim_state", "sim_imputation_stata3.csv"))[-1,]
imp4 <- read_csv(here("..", "StataSim_state", "sim_imputation_stata4.csv"))[-1,]



imp <- bind_rows(imp1, imp2, imp3, imp4) %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         type = c2,
         t = c3,
         estimate = c4,
         var = c5) %>%
  mutate(std.error = sqrt(var),
         method = ifelse(type == 0, "BJS Stata", "BJS Stata LOO")) %>%
  select(iteration, t, estimate, std.error, method)

imp_stata <- for_cdh %>%
  left_join(imp[imp$method == "BJS Stata",], by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, true_att_deathrate, pretreatment_deathrate, sd_deathrate, notreated, pct_effect, iteration)


imp_loo <- for_cdh %>%
  left_join(imp[imp$method == "BJS Stata LOO",], by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, true_att_deathrate, pretreatment_deathrate, sd_deathrate, notreated, pct_effect, iteration) %>% 
  mutate(method = ifelse(method == "BJS Stata LOO", "BJS Loo", method))



step1 <- read_csv(here("..", "StataSim_state", "sim_stepwise_stata1.csv"))[-1,]
step2 <- read_csv(here("..", "StataSim_state", "sim_stepwise_stata2.csv"))[-1,]
step3 <- read_csv(here("..", "StataSim_state", "sim_stepwise_stata3.csv"))[-1,]
step4 <- read_csv(here("..", "StataSim_state", "sim_stepwise_stata4.csv"))[-1,]



stepwise <- bind_rows(step1, step2, step3, step4) %>%
  select(-c(`_rowname`)) %>%
  mutate(method = "Stepwise") %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         var = c4) %>%
  mutate(std.error = sqrt(var)) %>% 
  select(iteration, t, estimate, std.error, method)

stepwise <- for_cdh %>% 
  left_join(stepwise, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, true_att_deathrate, pretreatment_deathrate, sd_deathrate, notreated, pct_effect, iteration)


dim(imp_stata)
dim(imp_loo)
dim(stepwise)

rm(imp1, imp2, imp3, imp4)
rm(cdh1, cdh2, cdh3, cdh4)
rm(step1, step2, step3, step4)



results <- bind_rows(results, cdh, imp_loo, stepwise) %>% ### no imp stata, no imp loo for now
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "BJS Loo", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE"), order = TRUE)) %>% 
  arrange(iteration)


results <- results %>%
  mutate(conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att_deathrate) >= conf.low & (true_att_deathrate) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att_deathrate),
         error_pct = error/true_att_deathrate,
         abs.error = abs(error), 
         true_att_rounded = round(true_att_deathrate, 2))


summary_by_effect_prop <- results %>%
  group_by(method, pct_effect, notreated) %>%
  summarize(mean_true = mean(true_att_deathrate),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att_deathrate)) %>%
  mutate(notreated = factor(notreated, levels = c(6, 12, 24), labels = c(6, 12, 24), order = TRUE),
         ate = factor(round(pct_effect, 3), levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

## lowest power with 6 treated units
min(summary_by_effect_prop$share_sign[summary_by_effect_prop$notreated == 6])

summary_by_sign_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate > 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att_deathrate), 
         exageration_ratio = if_else(estimate > 0, NA, exageration_ratio)) %>%
  group_by(method, pct_effect, notreated) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(notreated = factor(notreated, levels = c(6, 12, 24), labels = c(6, 12, 24), order = TRUE),
         ate = factor(round(pct_effect, 3), levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

### max wrong sign
max(summary_by_sign_prop$Pct_wrong_sign)
######


##### pct wrong by method for 6 treated
print(summary_by_sign_prop[summary_by_sign_prop$notreated == 6, c("method", "Pct_wrong_sign", "ate")], n= 24)
##### exaggeration by method for 6 treated
print(summary_by_sign_prop[summary_by_sign_prop$notreated == 6, c("method", "exageration_ratio", "ate")], n= 24)
##### exaggeration by method for 12 treated
print(summary_by_sign_prop[summary_by_sign_prop$notreated == 12, c("method", "exageration_ratio", "ate")], n= 24)


##### coverage
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = coverage, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = coverage, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = coverage, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = CI_size, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated States", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 14)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = CI_size, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated States", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_all) + ylim(0, 14) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = CI_size, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated States", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 14) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = share_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Power", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = share_sign, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = share_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### rmse
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = rmse, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 4)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = rmse, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 4) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = rmse, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 4)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 6 , ], aes(y = mean_error, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .35), size = 3.5)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-2, 0.25)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))
pl

pm <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 12, ], aes(y = mean_error, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .35), size = 3.5)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "Mean Error") + scale_color_manual(name="Method", values = color_all)+ ylim(-2, 0.25) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pm

pr <-  ggplot(data = summary_by_effect_prop[summary_by_effect_prop$notreated == 24, ], aes(y = mean_error, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .35), size = 3.5)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-2, 0.25)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 6, ], aes(y = exageration_ratio, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_all)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(breaks = c(1, 3, 6, 9, 12), limits = c(0.9, 6))
pl

pm <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 12, ], aes(y = exageration_ratio, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_all)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +scale_y_continuous(limits = c(0.9, 6))
pm

pr <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 24, ], aes(y = exageration_ratio, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_all)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 6))
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 6, ], aes(y = Pct_wrong_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .5), size = 3.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "6 Treated Units", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 19)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 
pl

pm <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 12,  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .5), size = 3.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "12 Treated Units", y="", x = "True ATE (Percent Decrease)") + scale_color_manual(name="Method", values = color_all) + ylim(0, 19) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  
pm

pr <-  ggplot(data = summary_by_sign_prop[summary_by_sign_prop$notreated == 24,  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .5), size = 3.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "24 Treated Units", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 19)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 
pr

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_prop_state_opioids_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



