#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Fig 3
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 13 June 2026



# packages ####

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

here()



# START ####

rm(list = ls())




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A LOAD DATA ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data <- read_csv(
  here("data", "processed", "data_all.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(.default = "?")) %>%
  mutate(
    obs.year = as.factor(obs.year),
    eco.name = as.factor(eco.name),
    region = as.factor(region),
    mngm.type = as.factor(mngm.type),
    hydrology = as.factor(hydrology),
    site.type = as.factor(site.type),
    rest.meth = as.factor(rest.meth),
    rest.meth.type = as.factor(rest.meth.type),
    land.use.hist = as.factor(land.use.hist),
    region = if_else(region == "north", "North",
                     if_else(region == "centre", "Centre", "South")),
    region = factor(region, levels = c("North", "Centre", "South")),
    site.type = fct_relevel(site.type, "negative", "restored", "positive"),
    rest.meth = fct_relevel(rest.meth, "cus", "res", "dih", "mga"),
    rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "dih", "mga", "positive")
  )
str(data)


spe_rich_tot <- read.csv(
  here("outputs", "tables", 
       "EMM1_site_type_spe_rich_tot.csv")) %>%
  dplyr::mutate(site.type = factor(site.type, levels = c("negative", "restored", "positive")))


spe_rich_target <- read.csv(
  here("outputs", "tables", 
       "EMM1_site_type_spe_rich_target.csv")) %>%
  dplyr::mutate(site.type = factor(site.type, levels = c("negative", "restored", "positive")))


abu_tot <- read.csv(
  here("outputs", "tables", 
       "EMM1_site_type_abu_tot.csv")) %>%
  dplyr::mutate(site.type = factor(site.type, levels = c("negative", "restored", "positive")))


abu_target <- read.csv(
  here("outputs", "tables", 
       "EMM1_site_type_abu_target.csv")) %>%
  dplyr::mutate(site.type = factor(site.type, levels = c("negative", "restored", "positive")))


Hill2_tot <- read.csv(
  here("outputs", "tables", 
       "EMM1_site_type_Hill2_tot.csv")) %>%
  dplyr::mutate(site.type = factor(site.type, levels = c("negative", "restored", "positive")))


Hill2_target <- read.csv(
  here("outputs", "tables", 
       "EMM1_site_type_Hill2_target.csv")) %>%
  dplyr::mutate(site.type = factor(site.type, levels = c("negative", "restored", "positive")))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B PLOT AND SAVE ##############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my_theme <- function() {
  theme_test() +
    theme(
      text = element_text(size = 18, color = "black"),
      axis.title = element_blank(),
      axis.text = element_text(size = 18, color = "black"),
      plot.title = element_text(size = 14,  color = "black", face = "bold"),
      plot.subtitle = element_text(size = 18,  color = "black", face = "bold"),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(margin = ggplot2::margin(l = 2, unit = "mm"), size = 18),
      strip.background = element_rect(fill = "white", colour = "black")
    )
  }



## (A) Total species richness --------------------------------------------------

my_lim <- max(data$butterfly.rich.total)

plot_site_type_spe_rich_tot <- ggplot() +
  geom_violin(data = data, aes(site.type, butterfly.rich.total, fill = site.type)) +
  geom_point(aes(site.type, predicted), data = spe_rich_tot, size = 3) +
  geom_linerange(aes(site.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = spe_rich_tot) +
  scale_fill_manual(values = c("negative" = "#1e9b8a", "restored" = "#52c569", "positive" = "#c2df23")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("a", "b", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Negative\nn = 33 ",
      "Restored\nn = 121",
      "Positive\nn = 33 ") 
    ) +
  labs(title = "", subtitle = "(A) Total species richness",
       y = "", x = "")
print(plot_site_type_spe_rich_tot)



## (B) Target species richness -------------------------------------------------

my_lim <- max(data$butterfly.rich.total)

plot_site_type_spe_rich_target <- ggplot() +
  geom_violin(data = data, aes(site.type, butterfly.rich.target, fill = site.type)) +
  geom_point(aes(site.type, predicted), data = spe_rich_target, size = 3) +
  geom_linerange(aes(site.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = spe_rich_target) +
  scale_fill_manual(values = c("negative" = "#1e9b8a", "restored" = "#52c569", "positive" = "#c2df23")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("a", "a", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Negative\nn = 33 ",
      "Restored\nn = 121",
      "Positive\nn = 33 ")
    ) +
  labs(title = "", subtitle = "(B) Target species richness",
       y = "", x = "")
print(plot_site_type_spe_rich_target)



## (C) Total abundance ---------------------------------------------------------

my_lim <- max(data$butterfly.abu)

plot_site_type_abu_tot <- ggplot() +
  geom_violin(data = data, aes(site.type, butterfly.abu, fill = site.type)) +
  geom_point(aes(site.type, predicted), data = abu_tot, size = 3) +
  geom_linerange(aes(site.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = abu_tot) +
  scale_fill_manual(values = c("negative" = "#1e9b8a", "restored" = "#52c569", "positive" = "#c2df23")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("a", "b", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Negative\nn = 33 ",
      "Restored\nn = 121",
      "Positive\nn = 33 ")
    ) +
  labs(title = "", subtitle = "(C) Total abundance",
       y = "", x = "")
print(plot_site_type_abu_tot)



## (D) Target abundance --------------------------------------------------------

my_lim <- max(data$butterfly.abu)

plot_site_type_abu_target <- ggplot() +
  geom_violin(data = data, aes(site.type, butterfly.target.abu, fill = site.type)) +
  geom_point(aes(site.type, predicted), data = abu_target, size = 3) +
  geom_linerange(aes(site.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = abu_target) +
  scale_fill_manual(values = c("negative" = "#1e9b8a", "restored" = "#52c569", "positive" = "#c2df23")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("a", "b", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Negative\nn = 33 ",
      "Restored\nn = 121",
      "Positive\nn = 33 ")
    ) +
  labs(title = "", subtitle = "(D) Target abundance",
       y = "", x = "")
print(plot_site_type_abu_target)



## (E) Total Hill-Simpson diversity --------------------------------------------

my_lim <- max(data$butterfly.Hill2.total)

plot_site_type_Hill2_tot <- ggplot() +
  geom_violin(data = data, aes(site.type, butterfly.Hill2.total, fill = site.type)) +
  geom_point(aes(site.type, predicted), data = Hill2_tot, size = 3) +
  geom_linerange(aes(site.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = Hill2_tot) +
  scale_fill_manual(values = c("negative" = "#1e9b8a", "restored" = "#52c569", "positive" = "#c2df23")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("a", "a", "a"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Negative\nn = 33 ",
      "Restored\nn = 121",
      "Positive\nn = 33 ")
    ) +
  labs(title = "", subtitle = "(E) Total Hill-Simpson diversity",
       y = "", x = "")
print(plot_site_type_Hill2_tot)



## (F) Target Hill-Simpson ------------------------------------------------------

my_lim <- max(data$butterfly.Hill2.total)

plot_site_type_Hill2_target <- ggplot() +
  geom_violin(data = data, aes(site.type, butterfly.Hill2.target, fill = site.type)) +
  geom_point(aes(site.type, predicted), data = Hill2_target, size = 3) +
  geom_linerange(aes(site.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = Hill2_target) +
  scale_fill_manual(values = c("negative" = "#1e9b8a", "restored" = "#52c569", "positive" = "#c2df23")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("a", "a", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Negative\nn = 33 ",
      "Restored\nn = 121",
      "Positive\nn = 33 ")
    ) +
  labs(title = "", subtitle = "(F) Target Hill-Simpson diversity",
       y = "", x = "")
print(plot_site_type_Hill2_target)



## combine ---------------------------------------------------------------------

Fig3 <- ggarrange(plot_site_type_spe_rich_tot, 
                  plot_site_type_spe_rich_target, 
                  plot_site_type_abu_tot, 
                  plot_site_type_abu_target,
                  plot_site_type_Hill2_tot,
                  plot_site_type_Hill2_target,
                     ncol = 2, nrow = 3)
ggsave(here("outputs", "figures","Fig3_site_type.png"),
       plot = Fig3, width = 14, height = 16, dpi = 800)




# END ####