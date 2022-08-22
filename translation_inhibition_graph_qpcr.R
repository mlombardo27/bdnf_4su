#remaking matteo's translation inhibition graph
#code from matteo 

library(tidyverse)
library(rstatix)
library(ggplot2)

value_sso <- read.table("C:/Users/mlomb/Desktop/tracked_files_github/bdnf_4su/data/qpcr_results.csv", sep = ",", header = T, stringsAsFactors = T)
value_sso$treatment <- factor(value_sso$treatment, levels = unique(value_sso$treatment))


stattest <- value_sso %>%
  group_by(gene) %>%
  pairwise_t_test(value ~ treatment, ref.group = "No_treatment") %>%
  add_significance() %>%
  add_xy_position()
stattest


nmd_plot_chx <- ggplot(value_sso, aes(x = treatment, y = value, color = treatment)) + 
  geom_point(position = 'jitter') + 
  stat_summary(fun.data = mean_se,  geom = "errorbar", aes(width = 0.2)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = treatment, alpha = .1), show.legend = F) + 
  facet_grid(cols = vars(gene), switch = "x") +
  #stat_pvalue_manual(stattest, hide.ns = T) + #, y.position = c(1.1, 1.18, 1.26, 1.1, 1.18, 1.26, 1.34, 1.42, 1.5), hide.ns = T) +
  theme_classic() + 
  #scale_color_brewer(palette = "Dark2") + 
  #scale_fill_brewer(palette = "Dark2") +
  ggeasy::easy_rotate_x_labels() +
  xlab("Treatment") + ylab("Percent of transcript after treatment")
plot(nmd_plot_chx)
