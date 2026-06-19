#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Fig S6 Effect of restoration age (excluding cultivar sites)
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

## data_x_cus ####
data_x_cus <- read_csv(
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
  filter(rest.meth != "cus",
         mang.app.NEU.MM.1.minus != "NA")
str(data_x_cus)



## prediction for total species richness ####
# rest.age (.)
Pred_noCus_spe_rich_tot <- read.csv(file = here::here("outputs", "tables",
                                        "Pred_noCus_spe_rich_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))



## prediction for total Hill-Simpson diversity ####
# rest.age (.)
Pred_noCus_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                      "Pred_noCus_Hill2_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))



## prediction for target Hill-Simpson diversity ####
# rest.age *
Pred_noCus_Hill2_target <- read.csv(file = here::here("outputs", "tables",
                                                   "Pred_noCus_Hill2_target.csv")) %>%
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



## (A) Total species richness ####
FigS6A <- 
  ggplot(data_x_cus, aes(rest.age, butterfly.rich.total, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_noCus_spe_rich_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(A) Total species richness",
       y = "", x = "Age of restored sites")
print(FigS6A)



## (B) Total Hill-Simpson diversity ####
FigS6B <- 
  ggplot(data_x_cus, aes(rest.age, butterfly.Hill2.total)) +
  geom_point(aes(color = "")) +  # Dummy-Variable
  scale_color_manual(values = c("white")) + # color of Dummy-Points
  geom_point() +
  my_theme() +
  geom_line(data = Pred_noCus_Hill2_tot, linewidth = 1, linetype = "dashed") +
  labs(title = "", subtitle = "(B) Total Hill-Simpson diversity",
       y = "", x = "Age of restored sites", color = "")
print(FigS6B)



## (C) Target Hill-Simpson diversity ####
FigS6C <- 
  ggplot(data_x_cus, aes(rest.age, butterfly.Hill2.target, color = region)) +
  geom_point() +
  facet_wrap(~obs.year) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  my_theme() +
  geom_line(data = Pred_noCus_Hill2_target, linewidth = 1) +
  labs(title = "", subtitle = "(C) Target Hill-Simpson diversity",
       y = "", x = "Age of restored sites")
print(FigS6C)



## combine ####
FigS6 <- ggarrange(FigS6A,
                   FigS6B,
                   FigS6C,
                     ncol = 3, nrow = 1)
print(FigS6)
ggsave(here("outputs", "figures","FigS6_rest.age_no-cultivar.png"),
       plot = FigS6, width = 21, height = 5.33, dpi = 800)




# END ####