#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Fig 4
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 14 June 2026



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
  ) %>%
  filter(site.type == "restored")
str(data)

spe_rich_tot <- read.csv(
  here("outputs", "tables", 
       "EMM1_RM_spe_rich_tot.csv")) %>%
  dplyr::mutate(site.type = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))


spe_rich_target <- read.csv(
  here("outputs", "tables", 
       "EMM1_RM_spe_rich_target.csv")) %>%
  dplyr::mutate(site.type = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))


abu_tot <- read.csv(
  here("outputs", "tables", 
       "EMM1_RM_abu_tot.csv")) %>%
  dplyr::mutate(site.type = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))


abu_target <- read.csv(
  here("outputs", "tables", 
       "EMM1_RM_abu_target.csv")) %>%
  dplyr::mutate(site.type = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))


Hill2_tot <- read.csv(
  here("outputs", "tables", 
       "EMM1_RM_Hill2_tot.csv")) %>%
  dplyr::mutate(site.type = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))


Hill2_target <- read.csv(
  here("outputs", "tables", 
       "EMM1_RM_Hill2_target.csv")) %>%
  dplyr::mutate(site.type = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B PLOT AND SAVE ##############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

table(data$rest.meth)
# cus res dih mga 
# 21  38  40  22

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

plot_RM_spe_rich_tot <- ggplot() +
  geom_violin(data = data, aes(rest.meth, butterfly.rich.total, fill = rest.meth)) +
  geom_point(aes(rest.meth, predicted), data = spe_rich_tot, size = 3) +
  geom_linerange(aes(rest.meth, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = spe_rich_tot) +
  scale_fill_manual(values = c("cus" = "#fbb61a", 
                               "res" = "#ed6925", 
                               "dih" = "#bc3754", 
                               "mga" = "#781c6d")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3, 4), y = my_lim * 1,
           label = c("a", "ab", "(b)", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 21",
      "Regional seed\nmixture\nn = 38",
      "Direct \nharvesting\nn = 40",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "(A) Total species richness",
       y = "", x = "")
print(plot_RM_spe_rich_tot)



## (B) Target species richness -------------------------------------------------

my_lim <- max(data$butterfly.rich.total)

plot_RM_spe_rich_target <- ggplot() +
  geom_violin(data = data, aes(rest.meth, butterfly.rich.target, fill = rest.meth)) +
  geom_point(aes(rest.meth, predicted), data = spe_rich_target, size = 3) +
  geom_linerange(aes(rest.meth, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = spe_rich_target) +
  scale_fill_manual(values = c("cus" = "#fbb61a", 
                               "res" = "#ed6925", 
                               "dih" = "#bc3754", 
                               "mga" = "#781c6d")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3, 4), y = my_lim * 1,
           label = c("a", "b", "b", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 21",
      "Regional seed\nmixture\nn = 38",
      "Direct \nharvesting\nn = 40",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "(B) Target species richness",
       y = "", x = "")
print(plot_RM_spe_rich_target)



## (C) Total abundance ---------------------------------------------------------

my_lim <- max(data$butterfly.abu)

plot_RM_abu <- ggplot() +
  geom_violin(data = data, aes(rest.meth, butterfly.abu, fill = rest.meth)) +
  geom_point(aes(rest.meth, predicted), data = abu_tot, size = 3) +
  geom_linerange(aes(rest.meth, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = abu_tot) +
  scale_fill_manual(values = c("cus" = "#fbb61a", 
                               "res" = "#ed6925", 
                               "dih" = "#bc3754", 
                               "mga" = "#781c6d")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3, 4), y = my_lim * 1,
           label = c("a", "b", "b", "ab"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 21",
      "Regional seed\nmixture\nn = 38",
      "Direct \nharvesting\nn = 40",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "(C) Total abundance",
       y = "", x = "")
print(plot_RM_abu)



## (D) Target abundance --------------------------------------------------------

my_lim <- max(data$butterfly.abu)

plot_RM_abu_target <- ggplot() +
  geom_violin(data = data, aes(rest.meth, butterfly.target.abu, fill = rest.meth)) +
  geom_point(aes(rest.meth, predicted), data = abu_target, size = 3) +
  geom_linerange(aes(rest.meth, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = abu_target) +
  scale_fill_manual(values = c("cus" = "#fbb61a", 
                               "res" = "#ed6925", 
                               "dih" = "#bc3754", 
                               "mga" = "#781c6d")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3, 4), y = my_lim * 1,
           label = c("a", "b", "(b)", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 21",
      "Regional seed\nmixture\nn = 38",
      "Direct \nharvesting\nn = 40",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "(D) Target abundance",
       y = "", x = "")
print(plot_RM_abu_target)



## (E) Total Hill-Simpson diversity --------------------------------------------

my_lim <- max(data$butterfly.Hill2.total)

plot_RM_Hill2_tot <- ggplot() +
  geom_violin(data = data, aes(rest.meth, butterfly.Hill2.total, fill = rest.meth)) +
  geom_point(aes(rest.meth, predicted), data = Hill2_tot, size = 3) +
  geom_linerange(aes(rest.meth, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = Hill2_tot) +
  scale_fill_manual(values = c("cus" = "#fbb61a", 
                               "res" = "#ed6925", 
                               "dih" = "#bc3754", 
                               "mga" = "#781c6d")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3, 4), y = my_lim * 1,
           label = c("a", "a", "a", "a"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 21",
      "Regional seed\nmixture\nn = 38",
      "Direct \nharvesting\nn = 40",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "(E) Total Hill-Simpson diversity",
       y = "", x = "")
print(plot_RM_Hill2_tot)



## (F) Target Hill-Simpson diversity -------------------------------------------

my_lim <- max(data$butterfly.Hill2.total)

plot_RM_Hill2_target <- ggplot() +
  geom_violin(data = data, aes(rest.meth, butterfly.Hill2.target, fill = rest.meth)) +
  geom_point(aes(rest.meth, predicted), data = Hill2_target, size = 3) +
  geom_linerange(aes(rest.meth, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = Hill2_target) +
  scale_fill_manual(values = c("cus" = "#fbb61a", 
                               "res" = "#ed6925", 
                               "dih" = "#bc3754", 
                               "mga" = "#781c6d")) +
  my_theme() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3, 4), y = my_lim * 1,
           label = c("a", "ab", "ab", "b"),
           vjust = -0.5, size = 6) +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 21",
      "Regional seed\nmixture\nn = 38",
      "Direct \nharvesting\nn = 40",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "(F) Target Hill-Simpson diversity",
       y = "", x = "")
print(plot_RM_Hill2_target)



## combine ---------------------------------------------------------------------

Fig4 <- ggarrange(plot_RM_spe_rich_tot, 
                  plot_RM_spe_rich_target, 
                  plot_RM_abu, 
                  plot_RM_abu_target,
                  plot_RM_Hill2_tot,
                  plot_RM_Hill2_target,
                     ncol = 2, nrow = 3)
print(Fig4)
ggsave(here("outputs", "figures","Fig4_restoration_method.png"),
       plot = Fig4, width = 14, height = 16, dpi = 800)




# END ####