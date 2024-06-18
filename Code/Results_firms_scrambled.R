library(here)
library(tidyverse)
library(viridis)
library(patchwork)     # For combining plots
library(arrow)

datasets_panel <- open_dataset(here("..", "SimResults_firms"))
datasets_scrambled <- open_dataset(here("..", "SimResults_firms_scrambled"))



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
################# log rev
#### 
results_scrambled <- datasets_scrambled %>%
  filter(DV == "log_rev" & method != "TWFE") %>% 
  collect() %>% 
  mutate(Data = "Reconstructed Time Series")


results_panel <- datasets_panel %>%
  filter(DV == "log_rev" & method != "TWFE" & units == 250 & prop_treated == 0.4) %>% 
  collect() %>% 
  mutate(Data = "Full Time Series")



results <- bind_rows(results_panel, results_scrambled) %>% 
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA"), order = TRUE)) %>% 
  arrange(iteration)

rm(results_panel, results_scrambled)
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
  group_by(method, pct_att_mean, Data) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE)) %>% 
  filter(ate != "15%")

## lowest power , woah?
min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, Data) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE)) %>% 
  filter(ate != "15%")

### max wrong sign
max(summary_by_sign_effect_prop$Pct_wrong_sign)
######


##### coverage
pm <-  ggplot(data = summary_effect_prop, aes(y = coverage, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Coverage", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "coverage_att_panel_logrev.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### power
pm <-  ggplot(data = summary_effect_prop, aes(y = share_sign, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Power", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "power_att_panel_logrev.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### exaggeration
pm <-  ggplot(data = summary_by_sign_effect_prop, aes(y = exageration_ratio, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Exaggeration \nRatio", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 4) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "exaggeration_att_panel_logrev.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### pct wrong sign  
pm <-  ggplot(data = summary_by_sign_effect_prop, aes(y = Pct_wrong_sign, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="% Wrong \nSign", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 5) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "wrongsign_att_panel_logrev.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)







########################################################################
########################################################################
##################################################################
################# roa
#### 
results_scrambled <- datasets_scrambled %>%
  filter(DV == "roa" & method != "TWFE") %>% 
  collect() %>% 
  mutate(Data = "Reconstructed Time Series")


results_panel <- datasets_panel %>%
  filter(DV == "roa" & method != "TWFE" & units == 250 & prop_treated == 0.4) %>% 
  collect() %>% 
  mutate(Data = "Full Time Series")



results <- bind_rows(results_panel, results_scrambled) %>% 
  mutate(method = case_when(method == "CdH" ~ "dCDH",
                            method == "CSA" ~ "CS",
                            method == "Mundlak" ~ "ETWFE",
                            TRUE ~ method),
         method = factor(method, levels = c("BJS", "CS", "dCDH", "ETWFE", "Gardner", "SA"), order = TRUE)) %>% 
  arrange(iteration)

rm(results_panel, results_scrambled)
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
  group_by(method, pct_att_mean, Data) %>%
  summarize(mean_true = mean(true_att),
            mean_pct_sd = mean(pct_att_sd),
            coverage = mean(CI_coverage),
            CI_size = mean(CI_size),
            mean_abs_error = mean(abs.error),
            mean_error = mean(error),
            mean_error_pct = mean(error_pct),
            share_sign = mean(Significant),
            rmse = my_rmse(estimate, true_att)) %>%
  mutate(ate = factor(pct_att_mean, levels = c(0.025, 0.050, 0.100, 0.150), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE)) %>% 
  filter(ate != "15%")

## lowest power , woah?
min(summary_effect_prop$share_sign, na.rm = T)

summary_by_sign_effect_prop <- results %>%
  filter(p.value <= 0.05) %>%
  mutate(wrong_sign = if_else(estimate < 0, 1, 0), ### true effects are all negative
         exageration_ratio = abs(estimate)/abs(true_att),
         exageration_ratio = if_else(exageration_ratio < 0, NA, exageration_ratio)) %>%
  group_by(method, pct_att_mean, Data) %>%
  summarize(exageration_ratio = mean(exageration_ratio, na.rm = TRUE),
            Pct_wrong_sign = mean(wrong_sign)*100)  %>%
  mutate(ate = factor(pct_att_mean, levels = c(0.025, 0.05, 0.10, 0.15), labels = c("2.5%", "5%", "10%", "15%"), order = TRUE)) %>% 
  filter(ate != "15%")

### max wrong sign
max(summary_by_sign_effect_prop$Pct_wrong_sign)
######


##### coverage
pm <-  ggplot(data = summary_effect_prop, aes(y = coverage, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.95, color = "red", linewidth = .75, alpha = 0.75)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Coverage", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0.2, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "coverage_att_panel_roa.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### power
pm <-  ggplot(data = summary_effect_prop, aes(y = share_sign, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 0.8, color = "red", linewidth = .75, alpha = 0.85)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Power", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 1) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "power_att_panel_roa.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)



##### exaggeration
pm <-  ggplot(data = summary_by_sign_effect_prop, aes(y = exageration_ratio, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + geom_hline(yintercept = 1, color = "red", size = .75, alpha = 0.5)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="Exaggeration \nRatio", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 11) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "exaggeration_att_panel_roa.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)


##### pct wrong sign  
pm <-  ggplot(data = summary_by_sign_effect_prop, aes(y = Pct_wrong_sign, x = ate , color = method, group = method)) + facet_wrap(~Data, scales = "free_y")
pm <- pm + geom_point(position = position_dodge2(width = .65), size = 6)
pm <- pm + theme_mfx() + theme(legend.position= 'bottom',legend.key=element_rect(fill='white'))
pm <- pm + labs(y="% Wrong \nSign", x = "True ATT as % Mean") + scale_color_manual(name="Method", values = color_nocdh) + ylim(0, 18) #+ coord_flip()
pm <- pm + guides(col = guide_legend(nrow = 2)) +  theme(strip.background = element_blank(), strip.text.x = element_text(size = 24))
pm
ggsave(here("..", "output", "wrongsign_att_panel_roa.pdf"), plot = pm, device = cairo_pdf,  height = 10, width = 10* 1.618)





