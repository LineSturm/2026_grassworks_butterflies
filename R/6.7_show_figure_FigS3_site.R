#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Fig S3 Effect of site variables
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 18 June 2026




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

## data all ####
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
  filter(site.type == "restored",
         mang.app.NEU.MM.1.minus != "NA")
str(data)



## prediction for total species richness ####
# plant.target.richness ***
Pred_site_var1_spe_rich_tot <- read.csv(here("outputs", "tables",
                                        "Pred_site_var1_spe_rich_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))



## prediction for target species richness ####
# plant.target.richness ***
Pred_site_var1_spe_rich_target <- read.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_spe_rich_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# site.cwm.pres.oek.f *
Pred_site_var2_spe_rich_target <- read.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var2_spe_rich_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# mngm.type (.)
EMM1_MType_spe_rich_target <- read.csv(file = here::here("outputs", "tables", 
                                          "EMM1_MType_spe_rich_target.csv"))



## prediction for total abundance ####
# plant.target.richness *
Pred_site_var1_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_abu_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# site.cwm.pres.oek.f (.)
Pred_site_var2_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                                     "Pred_site_var2_abu_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# veg.heterogeneity (.)
Pred_site_var3_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                                     "Pred_site_var3_abu_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))



## prediction for target abundance ####
# site.cwm.pres.oek.f **
Pred_site_var1_abu_target <- read.csv(file = here::here("outputs", "tables",
                                                     "Pred_site_var1_abu_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))



## prediction for total Hill-Simpson ####
# plant.target.richness **
Pred_site_var1_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                        "Pred_site_var1_Hill2_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# cover.vegetation *
Pred_site_var2_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                       "Pred_site_var2_Hill2_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))



## prediction for target Hill-Simpson ####
# plant.target.richness ***
Pred_site_var1_Hill2_target <- read.csv(file = here::here("outputs", "tables",
                                                       "Pred_site_var1_Hill2_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# cover.vegetation *
Pred_site_var2_Hill2_target <- read.csv(file = here::here("outputs", "tables",
                                                          "Pred_site_var2_Hill2_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B PLOT AND SAVE ##############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my_theme <- function() {
  theme_test() +
    theme(
      text = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 18, color = "black"), 
      plot.title = element_text(size = 14,  color = "black", face = "bold"),
      plot.subtitle = element_text(size = 18,  color = "black", face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(margin = ggplot2::margin(l = 2, unit = "mm"), size = 14),
      strip.background = element_rect(fill = "white", colour = "black")
    )
  }

my_theme2 <- function() {
  theme_test() +
    theme(
      text = element_text(size = 18, color = "black"),
      axis.title = element_blank(),
      axis.text = element_text(size = 18, color = "black"),
      plot.title = element_text(size = 14,  color = "black", face = "bold"),
      plot.subtitle = element_text(size = 18,  color = "black", face = "bold"),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(margin = ggplot2::margin(l = 2, unit = "mm"), size = 14),
      strip.background = element_rect(fill = "white", colour = "black")
    )
  }



## (A) Char. plant species richness -> total species richness ####

FigS3A <- 
  ggplot(data, aes(plant.target.richness, butterfly.rich.total, color = region )) +
  geom_point() +
  facet_wrap(~obs.year) +
  my_theme() +
  geom_line(data = Pred_site_var1_spe_rich_tot, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(A) Total species richness",
       y = "", x = "Char. plant species richness")
print(FigS3A)



## (B) Char. plant species richness -> target species richness ####

FigS3B <- 
  ggplot(data, aes(plant.target.richness, butterfly.rich.target, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var1_spe_rich_target, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(B) Target species richness",
       y = "", x = "Char. plant species richness")
print(FigS3B)



## (C) Ellenberg moisture -> target species richness ####

FigS3C <-
  ggplot(data, aes(site.cwm.pres.oek.f, butterfly.rich.target, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var2_spe_rich_target, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(C) Target species richness",
       y = "", x = "Ellenberg moisture value")
print(FigS3C)



## (D) Management type -> target species richness ####

table(data$mngm.type)

my_lim <- max(data$butterfly.rich.target)

FigS3D <- ggplot() +
  geom_violin(data = data, aes(mngm.type, butterfly.rich.target, fill = mngm.type)) +
  geom_point(aes(mngm.type, predicted), data = EMM1_MType_spe_rich_target, size = 3) +
  geom_linerange(aes(mngm.type, ymin = conf.low, ymax = conf.high),
                 linewidth = 1, linetype = 1.5, data = EMM1_MType_spe_rich_target) +
  scale_fill_manual(values = c("grazing" = "#b28500",
                               "mowing" = "#f6d543",
                               "both" = "#d8c27f")) +
  my_theme2() +
  ylim(0, my_lim*1.1) +
  # Add text on plot for significance
  annotate("text", x = c(1, 2, 3), y = my_lim * 1,
           label = c("ab", "a", "(b)"),
           vjust = -0.5, size = 5) +
  scale_x_discrete(
    labels = c(
      "Both\nn = 14",
      "Grazing\nn = 23",
      "Mowing\nn = 75")
    ) +
  labs(title = "", subtitle = "(D) Target species richness",
       y = "", x = "")
print(FigS3D)



## (E) Char. plant species richness -> total abundance ####

FigS3E <- 
  ggplot(data, aes(plant.target.richness, butterfly.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var1_abu_tot, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(E) Total abundance",
       y = "", x = "Char. plant species richness")
print(FigS3E)



## (F) Ellenberg moisture -> total abundance ####

FigS3F <- 
  ggplot(data, aes(site.cwm.pres.oek.f, butterfly.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var2_abu_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(F) Total abundance",
       y = "", x = "Ellenberg moisture value")
print(FigS3F)



## (G) Veg. heterogeneity -> total abundance ####

FigS3G <- 
  ggplot(data, aes(veg.heterogeneity, butterfly.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var3_abu_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(G) Total abundance",
       y = "", x = "Vegetation heterogeneity")
print(FigS3G)



## (H) Ellenberg moisture -> target abundance ####

FigS3H <- 
  ggplot(data, aes(site.cwm.pres.oek.f, butterfly.target.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var1_abu_target, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(H) Target abundance",
       y = "", x = "Ellenberg moisture value")
print(FigS3H)



## (I) Char. plant species richness -> total Hill-Simpson diversity ####

FigS3I <- 
  ggplot(data, aes(plant.target.richness, butterfly.Hill2.total)) +
  geom_point(aes(color = "")) +  # Dummy-Variable
  scale_color_manual(values = c("white")) + # color of Dummy-Points
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var1_Hill2_tot, linewidth = 1) +
  labs(title = "", subtitle = "(I) Total Hill-Simpson diversity",
       y = "", x = "Char. plant species richness", color = "")
print(FigS3I)



## (J) Cover vegetation -> total Hill-Simpson diversity ####

FigS3J <- 
  ggplot(data, aes(cover.vegetation, butterfly.Hill2.total)) +
  geom_point(aes(color = "")) +  # Dummy-Variable
  scale_color_manual(values = c("white")) + # color of Dummy-Points
  geom_point() +
  my_theme() +
  geom_line(data = Pred_site_var2_Hill2_tot, linewidth = 1) +
  labs(title = "", subtitle = "(J) Total Hill-Simpson diversity",
       y = "", x = "Vegetation cover", color = "")
print(FigS3J)



## (K) Char. plant species richness -> target Hill-Simpson diversity ####

FigS3K <- 
  ggplot(data, aes(plant.target.richness, butterfly.Hill2.target, color = region)) +
  geom_point() +
  facet_wrap(~obs.year) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  my_theme() +
  geom_line(data = Pred_site_var1_Hill2_target, linewidth = 1) +
  labs(title = "", subtitle = "(K) Target Hill-Simpson diversity",
       y = "", x = "Char. plant species richness")
print(FigS3K)



## (L) Cover vegetation -> target Hill-Simpson diversity ####

FigS3L <- 
  ggplot(data, aes(cover.vegetation, butterfly.Hill2.target, color = region)) +
  geom_point() +
  facet_wrap(~obs.year) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  my_theme() +
  geom_line(data = Pred_site_var2_Hill2_target, linewidth = 1) +
  labs(title = "", subtitle = "(L) Target Hill-Simpson diversity",
       y = "", x = "Vegetation cover")
print(FigS3L)



## combine ####

FigS3_1 <- ggarrange(FigS3A,
                     FigS3B,
                     FigS3C,
                     FigS3D,
                     FigS3E,
                     FigS3F,
                     FigS3G,
                     FigS3H,
                     FigS3I, # part 1
                     ncol = 3, nrow = 3)
print(FigS3_1)
ggsave(here("outputs", "figures","FigS3_page1_site.png"),
       plot = FigS3_1, width = 21, height = 16, dpi = 800)

FigS3_2 <- ggarrange(FigS3J, 
                     FigS3K, 
                     FigS3L, # part 2
                     ncol = 3, nrow = 1)
print(FigS3_2)
ggsave(here("outputs", "figures","FigS3_page2_site.png"),
       plot = FigS3_2, width = 21, height = 5.33, dpi = 800)




# END ####