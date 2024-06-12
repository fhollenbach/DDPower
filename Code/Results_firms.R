library(here)
library(tidyverse)
library(viridis)
library(patchwork)     # For combining plots
library(arrow)

datasets <- open_dataset(here("..", "SimResults_firms"))



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
################# varying prop treated and effect size, set units to 100 
################### 
#### 100 observations, DV = log_rev
results <- datasets %>%
  filter(t == -100 & units == 100 & DV == "log_rev" & method != "TWFE") %>% 
  collect()


for_cdh <- results %>%
   select(-c(t, estimate, std.error, method)) %>%
   group_by(iteration) %>%
   summarize_all(unique)
#
cdh1 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata1.csv"))[-1,]
cdh2 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata2.csv"))[-1,]
cdh3 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata3.csv"))[-1,]
cdh4 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata4.csv"))[-1,]
cdh5 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata5.csv"))[-1,]
cdh6 <- read_csv(here("..", "StataSim_firms","sim_CdH_firm_stata6.csv"))[-1,]
cdh7 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata7.csv"))[-1,]
cdh8 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata8.csv"))[-1,]


cdh <- bind_rows(cdh1, cdh2, cdh3, cdh4, cdh5, cdh6, cdh7, cdh8) %>%
  mutate(method = "CdH") %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         std.error = c4) %>%
  filter(t == -100)

cdh <- for_cdh %>%
  left_join(cdh, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)

results <- bind_rows(results, cdh) %>% ### no imp stata, no imp loo for now
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA"), order = TRUE)) %>% 
  arrange(iteration)


gc()

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))


summary_effect_prop <- results %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

## lowest power , woah?
min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

### max wrong sign
max(summary_by_sign_effect_prop$Pct_wrong_sign)
######


##### pct wrong by method for 10% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", c("method", "Pct_wrong_sign", "ate")], n= 24)
##### exaggeration by method for 20% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", c("method", "exageration_ratio", "ate")], n= 24)
##### exaggeration by method for 40% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", c("method", "exageration_ratio", "ate")], n= 24)


##### coverage
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = coverage, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = coverage, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = coverage, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%", ], aes(y = CI_size, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "5 Treated Firms", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = CI_size, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "10 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0, 2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = CI_size, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "20 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 2) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = share_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = share_sign, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = share_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### rmse
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = rmse, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.5)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = rmse, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.5) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = rmse, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.5)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = mean_error, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.1, 0.1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = mean_error, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.1, 0.1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = mean_error, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.1, 0.1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = exageration_ratio, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 7))


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", ], aes(y = exageration_ratio, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 7))


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", ], aes(y = exageration_ratio, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 7))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = Pct_wrong_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 27)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0, 27) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 27)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_prop_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
##################################################################
################# varying prop treated and effect size, set units to 100 
################### 
#### 100 observations, DV = roa
results <- datasets %>%
  filter(t == -100 & units == 100 & DV == "roa" & method != "TWFE") %>% 
  collect()


for_cdh <- results %>%
  select(-c(t, estimate, std.error, method)) %>%
  group_by(iteration) %>%
  summarize_all(unique)


#
cdh1 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata1.csv"))[-1,]
cdh2 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata2.csv"))[-1,]
cdh3 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata3.csv"))[-1,]
cdh4 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata4.csv"))[-1,]
cdh5 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata5.csv"))[-1,]
cdh6 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata6.csv"))[-1,]
cdh7 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata7.csv"))[-1,]
cdh8 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata8.csv"))[-1,]


cdh <- bind_rows(cdh1, cdh2, cdh3, cdh4, cdh5, cdh6, cdh7, cdh8) %>%
  mutate(method = "CdH") %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         std.error = c4) %>%
  filter(t == -100)

cdh <- for_cdh %>%
  left_join(cdh, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


results <- bind_rows(results, cdh) %>% ### no imp stata, no imp loo for now
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA"), order = TRUE)) %>% 
  arrange(iteration)
gc()

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))


summary_effect_prop <- results %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

## lowest power , woah?
min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

### max wrong sign
max(summary_by_sign_effect_prop$Pct_wrong_sign)
######


##### pct wrong by method for 10% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", c("method", "Pct_wrong_sign", "ate")], n= 24)
##### exaggeration by method for 20% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", c("method", "exageration_ratio", "ate")], n= 24)
##### exaggeration by method for 40% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", c("method", "exageration_ratio", "ate")], n= 24)


##### coverage
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = coverage, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = coverage, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = coverage, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%", ], aes(y = CI_size, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "5 Treated Firms", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = CI_size, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "10 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = CI_size, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "20 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = share_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Power", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = share_sign, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = share_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### rmse
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = rmse, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.25)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = rmse, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.25) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = rmse, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(0, 0.25)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = mean_error, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.05, 0.05)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = mean_error, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.05, 0.05) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = mean_error, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main)+ ylim(-0.05, 0.05)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = exageration_ratio, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_main)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 12))


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", ], aes(y = exageration_ratio, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0,12))


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", ], aes(y = exageration_ratio, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 12))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = Pct_wrong_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 22)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_main) + ylim(0, 22) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_main) + ylim(0, 22)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_prop_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
####### alll estimators for appendix


results <- datasets %>%
  filter(t == -100 & units == 100 & DV == "log_rev") %>% 
  collect()


for_cdh <- results %>%
  select(-c(t, estimate, std.error, method)) %>%
  group_by(iteration) %>%
  summarize_all(unique)
#
cdh1 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata1.csv"))[-1,]
cdh2 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata2.csv"))[-1,]
cdh3 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata3.csv"))[-1,]
cdh4 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata4.csv"))[-1,]
cdh5 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata5.csv"))[-1,]
cdh6 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata6.csv"))[-1,]
cdh7 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata7.csv"))[-1,]
cdh8 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata8.csv"))[-1,]


cdh <- bind_rows(cdh1, cdh2, cdh3, cdh4, cdh5, cdh6, cdh7, cdh8) %>%
  mutate(method = "CdH") %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         std.error = c4) %>%
  filter(t == -100)

cdh <- for_cdh %>%
  left_join(cdh, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


imp1 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata1.csv"))[-1,]
imp2 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata2.csv"))[-1,]
imp3 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata3.csv"))[-1,]
imp4 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata4.csv"))[-1,]
imp5 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata5.csv"))[-1,]
imp6 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata6.csv"))[-1,]
imp7 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata7.csv"))[-1,]
imp8 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata8.csv"))[-1,]


imp <- bind_rows(imp1, imp2, imp3, imp4, imp5, imp6, imp7, imp8) %>%
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
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


imp_loo <- for_cdh %>%
  left_join(imp[imp$method == "BJS Stata LOO",], by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)



step1 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata1.csv"))[-1,]
step2 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata2.csv"))[-1,]
step3 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata3.csv"))[-1,]
step4 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata4.csv"))[-1,]
step5 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata5.csv"))[-1,]
step6 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata6.csv"))[-1,]
step7 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata7.csv"))[-1,]
step8 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata8.csv"))[-1,]


stepwise <- bind_rows(step1, step2, step3, step4, step5, step6, step7, step8) %>%
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
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


dim(imp_stata)
dim(imp_loo)
dim(stepwise)

rm(imp1, imp2, imp3, imp4, imp5, imp6, imp7, imp8)
rm(cdh1, cdh2, cdh3, cdh4, cdh5, cdh6, cdh7, cdh8)
rm(step1, step2, step3, step4, step5, step6, step7, step8)

results <- bind_rows(results, cdh, stepwise) %>% ### no imp stata, no imp loo for now
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE"), order = TRUE)) %>% 
  arrange(iteration)


gc()

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))


summary_effect_prop <- results %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

## lowest power , woah?
min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

### max wrong sign
max(summary_by_sign_effect_prop$Pct_wrong_sign)
######


##### pct wrong by method for 10% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", c("method", "Pct_wrong_sign", "ate")], n= 24)
##### exaggeration by method for 20% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", c("method", "exageration_ratio", "ate")], n= 24)
##### exaggeration by method for 40% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", c("method", "exageration_ratio", "ate")], n= 24)


##### coverage
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = coverage, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = coverage, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = coverage, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_prop_firms_logrev_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%", ], aes(y = CI_size, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "5 Treated Firms", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = CI_size, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "10 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0, 2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = CI_size, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "20 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 2) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_prop_firms_logrev_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = share_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Power", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = share_sign, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = share_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_prop_firms_logrev_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### rmse
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = rmse, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.5)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = rmse, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.5) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = rmse, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.5)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_prop_firms_logrev_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = mean_error, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.1, 0.3)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = mean_error, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.1, 0.3) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = mean_error, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.1, 0.3)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_prop_firms_logrev_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = exageration_ratio, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_all)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 6))


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", ], aes(y = exageration_ratio, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 6))


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", ], aes(y = exageration_ratio, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 6))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_prop_firms_logrev_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = Pct_wrong_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 27)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0, 27) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 27)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_prop_firms_logrev_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
##################################################################
################# varying prop treated and effect size, set units to 100 
###################  all estimators
#### 100 observations, DV = roa
results <- datasets %>%
  filter(t == -100 & units == 100 & DV == "roa") %>% 
  collect()


for_cdh <- results %>%
  select(-c(t, estimate, std.error, method)) %>%
  group_by(iteration) %>%
  summarize_all(unique)


#
cdh1 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata1.csv"))[-1,]
cdh2 <- read_csv(here("..", "StataSim_firms",  "sim_CdH_firm_stata2.csv"))[-1,]
cdh3 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata3.csv"))[-1,]
cdh4 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata4.csv"))[-1,]
cdh5 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata5.csv"))[-1,]
cdh6 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata6.csv"))[-1,]
cdh7 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata7.csv"))[-1,]
cdh8 <- read_csv(here("..", "StataSim_firms", "sim_CdH_firm_stata8.csv"))[-1,]


cdh <- bind_rows(cdh1, cdh2, cdh3, cdh4, cdh5, cdh6, cdh7, cdh8) %>%
  mutate(method = "CdH") %>%
  select(-c(`_rowname`)) %>%
  rename(iteration = c1,
         t = c2,
         estimate = c3,
         std.error = c4) %>%
  filter(t == -100)

cdh <- for_cdh %>%
  left_join(cdh, by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


imp1 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata1.csv"))[-1,]
imp2 <- read_csv(here("..", "StataSim_firms",  "sim_imputation_firm_stata2.csv"))[-1,]
imp3 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata3.csv"))[-1,]
imp4 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata4.csv"))[-1,]
imp5 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata5.csv"))[-1,]
imp6 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata6.csv"))[-1,]
imp7 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata7.csv"))[-1,]
imp8 <- read_csv(here("..", "StataSim_firms", "sim_imputation_firm_stata8.csv"))[-1,]


imp <- bind_rows(imp1, imp2, imp3, imp4, imp5, imp6, imp7, imp8) %>%
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
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


imp_loo <- for_cdh %>%
  left_join(imp[imp$method == "BJS Stata LOO",], by = "iteration") %>%
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)



step1 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata1.csv"))[-1,]
step2 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata2.csv"))[-1,]
step3 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata3.csv"))[-1,]
step4 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata4.csv"))[-1,]
step5 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata5.csv"))[-1,]
step6 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata6.csv"))[-1,]
step7 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata7.csv"))[-1,]
step8 <- read_csv(here("..", "StataSim_firms", "sim_stepwise_firm_stata8.csv"))[-1,]


stepwise <- bind_rows(step1, step2, step3, step4, step5, step6, step7, step8) %>%
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
  select(t, estimate, std.error, method, true_att, pretreatment_dv, sd_dv, pct_effect, units, prop_treated, DV, iteration)


dim(imp_stata)
dim(imp_loo)
dim(stepwise)

rm(imp1, imp2, imp3, imp4, imp5, imp6, imp7, imp8)
rm(cdh1, cdh2, cdh3, cdh4, cdh5, cdh6, cdh7, cdh8)
rm(step1, step2, step3, step4, step5, step6, step7, step8)

results <- bind_rows(results, cdh, stepwise) %>% ### no imp stata, no imp loo for now
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA", "Stepwise", "TWFE"), order = TRUE)) %>% 
  arrange(iteration)


gc()

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))


summary_effect_prop <- results %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

## lowest power , woah?
min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, prop_treated) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(prop_treated = factor(prop_treated, levels = c(0.1, 0.2, 0.4, 0.6), labels = c("10%", "20%", "40%", "60%"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

### max wrong sign
max(summary_by_sign_effect_prop$Pct_wrong_sign)
######


##### pct wrong by method for 10% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", c("method", "Pct_wrong_sign", "ate")], n= 24)
##### exaggeration by method for 20% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", c("method", "exageration_ratio", "ate")], n= 24)
##### exaggeration by method for 40% treated
print(summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", c("method", "exageration_ratio", "ate")], n= 24)


##### coverage
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = coverage, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = coverage, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = coverage, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0.2, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_prop_firms_roa_appAll.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%", ], aes(y = CI_size, x = ate, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "5 Treated Firms", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = CI_size, x = ate , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "10 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = CI_size, x = ate, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "20 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_prop_firms_roa_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = share_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Power", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = share_sign, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = share_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_prop_firms_roa_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### rmse
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = rmse, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = rmse, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = rmse, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_prop_firms_roa_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "10%" , ], aes(y = mean_error, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.2, 0.75)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "20%", ], aes(y = mean_error, x = ate , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.2, 0.75) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_prop[summary_effect_prop$prop_treated == "40%", ], aes(y = mean_error, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all)+ ylim(-0.2, 0.75)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_prop_firms_roa_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = exageration_ratio, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_all)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0, 13))


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%", ], aes(y = exageration_ratio, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 13))


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%", ], aes(y = exageration_ratio, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0, 13))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_prop_firms_roa_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "10%", ], aes(y = Pct_wrong_sign, x = ate, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "10 Treated Firms", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 22)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "20%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "20 Treated Firms", y="", x = "ATT in % of Mean") + scale_color_manual(name="Method", values = color_all) + ylim(0, 22) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_effect_prop[summary_by_sign_effect_prop$prop_treated == "40%",  ], aes(y = Pct_wrong_sign, x = ate, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "40 Treated Firms", y="", x = "") + scale_color_manual(name="Method", values = color_all) + ylim(0, 22)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_prop_firms_roa_appAll_spring24.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##############################################################################################
##############################################################################################
##############################################################################################
########################## same plots but 40% treated and over effect size and N --- log rev



results <- datasets %>%
  filter(t == -100 & prop_treated == 0.4 & DV == "log_rev" & method != "TWFE") %>% 
  collect()


results <- results %>% 
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA", "TWFE"), order = TRUE)) #%>% 

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))


summary_effect_units <- results %>%
  group_by(method, pct_att_mean, units) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(units = factor(units, levels = c(100, 250, 1000, 2500), labels = c("100", "250", "1000", "2500"), order = TRUE),
         ate = factor(round(pct_att_mean, 3), levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))



summary_by_sign_effect_units <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, units) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(units = factor(units, levels = c(100, 250, 1000, 2500), labels = c("100", "250", "1000", "2500"), order = TRUE),
         ate = factor(round(pct_att_mean, 3), levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))



##### coverage
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%", ], aes(y = coverage, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = coverage, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = coverage, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_units_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%", ], aes(y = CI_size, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.9)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = CI_size, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.9) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = CI_size, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.9) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_units_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="Power", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = share_sign, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_units_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### rmse
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%" , ], aes(y = rmse, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = rmse, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(0, 0.2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = rmse, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_units_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%" , ], aes(y = mean_error, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATE = 5% of Mean", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(-0.02, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = mean_error, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATE = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(-0.02, 0.02) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = mean_error, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATE = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(-0.02, 0.02)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_units_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "5%", ], aes(y = exageration_ratio, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATE = 5% of Mean", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_nocdh)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 3.5))


pm <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "10%", ], aes(y = exageration_ratio, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATE = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 3.5))


pr <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "15%", ], aes(y = exageration_ratio, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATE = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 3.5))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_units_logrev_firms.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "5%", ], aes(y = Pct_wrong_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATE = 5% of Mean", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 5)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "10%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATE = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 5) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "15%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATE = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 5)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_units_firms_logrev.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##############################################################################################
##############################################################################################
##############################################################################################
########################## same plots but 40% treated and over effect size and N --- log rev



results <- datasets %>%
  filter(t == -100 & prop_treated == 0.4 & DV == "roa" & method != "TWFE") %>% 
  collect()


results <- results %>% 
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA", "TWFE"), order = TRUE)) #%>% 

results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2))


summary_effect_units <- results %>%
  group_by(method, pct_att_mean, units) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(units = factor(units, levels = c(100, 250, 1000, 2500), labels = c("100", "250", "1000", "2500"), order = TRUE),
         ate = factor(round(pct_att_mean, 3), levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))



summary_by_sign_effect_units <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, units) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(units = factor(units, levels = c(100, 250, 1000, 2500), labels = c("100", "250", "1000", "2500"), order = TRUE),
         ate = factor(round(pct_att_mean, 3), levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))



##### coverage
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%", ], aes(y = coverage, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="Coverage", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = coverage, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = coverage, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "coverage_units_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%", ], aes(y = CI_size, x = units, color = method, group = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="CI Length", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.5)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = CI_size, x = units , color = method, group = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.5) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = CI_size, x = units, color = method, group = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.5) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "CIlength_units_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)




##### power
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%" , ], aes(y = share_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.85)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="Power", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = share_sign, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = share_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0.8, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "power_units_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### rmse
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%" , ], aes(y = rmse, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATT = 5% of Mean", y="RMSE", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))


pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = rmse, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATT = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(0, 0.2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = rmse, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATT = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(0, 0.2)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "rmse_units_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pl <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "5%" , ], aes(y = mean_error, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATE = 5% of Mean", y="Average Error", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(-0.01, 0.1)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2))

pm <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "10%", ], aes(y = mean_error, x = units , color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATE = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(-0.01, 0.1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

pr <-  ggplot(data = summary_effect_units[summary_effect_units$ate == "15%", ], aes(y = mean_error, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 0, color = "red", size = .75, alpha = 0.75)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATE = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh)+ ylim(-0.01, 0.1)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "bias_units_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pl <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "5%", ], aes(y = exageration_ratio, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATE = 5% of Mean", y="Exaggeration \nRatio", x = "") + scale_color_manual(name="Method", values = color_nocdh)
pl <- pl + guides(col = guide_legend(nrow = 2)) + scale_y_continuous(limits = c(0.9, 10))


pm <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "10%", ], aes(y = exageration_ratio, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATE = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh)
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous( limits = c(0.9, 10))


pr <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "15%", ], aes(y = exageration_ratio, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATE = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh)
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0.9, 10))


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') +
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) &
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "exaggeration_units_roa_firms.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### pct wrong sign  
pl <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "5%", ], aes(y = Pct_wrong_sign, x = units, color = method))
pl <- pl + geom_point(position = position_dodge2(width = .65), size = 6)
pl <- pl + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pl <- pl + labs(subtitle = "ATE = 5% of Mean", y="% Wrong \nSign", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 18)#+ ylim(-1, 1) #+ coord_flip()
pl <- pl + guides(col = guide_legend(nrow = 2)) 


pm <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "10%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(subtitle = "ATE = 10% of Mean", y="", x = "Number of Firms") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 18) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())  


pr <-  ggplot(data = summary_by_sign_effect_units[summary_by_sign_effect_units$ate == "15%",  ], aes(y = Pct_wrong_sign, x = units, color = method))
pr <- pr + geom_point(position = position_dodge2(width = .65), size = 6)
pr <- pr + theme_mfx() + theme(legend.position= 'none',legend.key=element_rect(fill='white'))
pr <- pr + labs(subtitle = "ATE = 15% of Mean", y="", x = "") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 18)#+ ylim(-1, 1) #+ coord_flip()
pr <- pr + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) 


p <- (pl + pm + pr)/guide_area() + plot_layout(guides = 'collect') + 
  plot_layout(heights = c(0.85, 0.15))
p <- p + plot_annotation(
  title = "") &
  theme(plot.title = element_text(size = 25, face = "bold")) & 
  theme(text = element_text("IBM Plex Sans Condensed"))
ggsave(here("..", "output", "wrongsign_units_firms_roa.pdf"), plot = p, device = cairo_pdf,  height = 10, width = 10* 1.618)






##############################################################################################
##############################################################################################
##############################################################################################
########################## now compare N = 250, prop treated = 0.4
#### compare over effect & DV

results <- datasets %>%
  filter(t == -100 & units == 250 & prop_treated == 0.4  & method != "TWFE") %>% 
  collect()


results <- results %>% 
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA", "TWFE"), order = TRUE)) #%>% 


results <- results %>%
  mutate(estimate = estimate,
         std.error = std.error,
         conf.low = estimate - qnorm(0.975) * std.error,
         conf.high = estimate + qnorm(0.975) * std.error,
         p.value = 2 * (1 - pnorm(abs(estimate / std.error))),
         CI_coverage = if_else((true_att) >= conf.low & (true_att) <= conf.high, 1, 0),
         CI_size = conf.high - conf.low,
         Significant = if_else(p.value <= 0.05, 1, 0),
         error = (estimate) - (true_att),
         error_pct = error/true_att,
         abs.error = abs(error), 
         pct_att_mean = round(true_att/pretreatment_dv,3),
         pct_att_sd = round(true_att/sd_dv, 3),
         true_att_rounded = round(true_att, 2),
         mean_sd_ratio = pretreatment_dv/sd_dv)


summary_effect_prop <- results %>%
  group_by(method, pct_att_mean, DV) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att),
            mean_sd_ratio = mean(mean_sd_ratio)) %>%
  mutate(DV = factor(DV, levels = c("log_rev", "roa"), labels = c("Revenue (ln)", "ROA"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))

min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all positive
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, DV) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(DV = factor(DV, levels = c("log_rev", "roa"), labels = c("Revenue (ln)", "ROA"), order = TRUE),
         ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE))


##### coverage
pm <-  ggplot(data = summary_effect_prop, aes(y = coverage, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Coverage", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "coverage_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### CI length
pm <-  ggplot(data = summary_effect_prop, aes(y = CI_size, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="CI Length", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0.2, 0.75) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "CIlength_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)

##### power
pm <-  ggplot(data = summary_effect_prop, aes(y = share_sign, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Power", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "power_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### rmse
pm <-  ggplot(data = summary_effect_prop, aes(y = rmse, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="RMSE", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 0.2) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "rmse_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### bias
pm <-  ggplot(data = summary_effect_prop, aes(y = mean_error, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Mean Error", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(-.05, 0.05) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "bias_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### exaggeration
pm <-  ggplot(data = summary_by_sign_effect_prop, aes(y = exageration_ratio, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Exaggeration \nRatio", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 6) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "exaggeration_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### pct wrong sign  
pm <-  ggplot(data = summary_by_sign_effect_prop, aes(y = exageration_ratio, x = ate , color = method, group = method)) + facet_wrap(~DV, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="% Wrong \nSign", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 9) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "wrongsign_att_dv.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)

