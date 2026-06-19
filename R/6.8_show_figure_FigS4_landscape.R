#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Fig S4 Effect of landscape variables
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
# 300 m radius
# perm.grassland.r300 (.)
Pred_Lr300_spe_rich_tot <- read.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr300_spe_rich_tot.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# 600 m radius -> no significant landscape variables
# 1200 m radius -> no significant landscape variables



## prediction for target species richness ####
# 300 m radius -> no significant landscape variables
# 600 m radius -> no significant landscape variables
# 1200 m radius -> no significant landscape variables



## prediction for total abundance ####
# 300 m radius
# perm.grassland.r300 (.)
Pred_Lr300_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                                      "Pred_Lr300_abu_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))

# 600 m radius
# perm.grassland.r600 (.)
Pred_Lr600_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                                 "Pred_Lr600_abu_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))
# rest.age x lc.shannon.div.r600 (.)
Pred_Lr600_Int_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                                 "Pred_Lr600_Int_abu_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))

# 1200 m radius
# perm.grassland.r1200 *
Pred_Lr1200_abu_tot <- read.csv(file = here::here("outputs", "tables",
                                                 "Pred_Lr1200_abu_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))



## prediction for target abundance ####
# 300 m radius
# perm.grassland.r300 (.)
Pred_Lr300_abu_target <- read.csv(file = here::here("outputs", "tables",
                                                 "Pred_Lr300_abu_target.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))

# 600 m radius
# rest.age x lc.shannon.div.r600 (.)
Pred_Lr600_abu_target <- read.csv(file = here::here("outputs", "tables",
                                                    "Pred_Lr600_abu_target.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))

# 1200 m radius -> no significant landscape variables



## prediction for total Hill-Simpson diversity ####
# 300 m radius -> no significant landscape variables

# 600 m radius
# Grassland plant species richness *
Pred_Lr600_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                    "Pred_Lr600_Hill2_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))
# rest.age x lc.shannon.div.r600 (.)
Pred_Lr600_Int_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                   "Pred_Lr600_Int_Hill2_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))

# 1200 m radius
# lk.grassland.per.min *
Pred_Lr1200_1_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                   "Pred_Lr1200_1_Hill2_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))
# field.size.avg.r1200 (.)
Pred_Lr1200_2_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                      "Pred_Lr1200_2_Hill2_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))
# agrar.r1200 (.)
Pred_Lr1200_3_Hill2_tot <- read.csv(file = here::here("outputs", "tables",
                                                      "Pred_Lr1200_3_Hill2_tot.csv")) %>%
  dplyr::mutate(region = factor(region, levels = c("North", "Centre", "South")))



## prediction for target Hill-Simpson diversity ####
# 300 m radius
# rest.age x lc.shannon.r300 (.)
Pred_Lr300_Hill2_target <- read.csv(file = here::here("outputs", "tables",
                                                   "Pred_Lr300_Hill2_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# 600 m radius
# rest.age x lc.shannon.r600 (.)
Pred_Lr600_Hill2_target <- read.csv(file = here::here("outputs", "tables",
                                                      "Pred_Lr600_Hill2_target.csv")) %>%
  dplyr::mutate(obs.year = as.factor(obs.year), 
                region = factor(region, levels = c("North", "Centre", "South")))

# 1200 m radius -> no significant landscape variables




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



## (A) Perm. grassland (300 m) -> total species richness ####

FigS4A <-
  ggplot(data, aes(perm.grassland.r300, butterfly.rich.total, color = region )) +
  geom_point() +
  facet_wrap(~obs.year) +
  my_theme() +
  geom_line(data = Pred_Lr300_spe_rich_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(A) Total species richness",
       y = "", x = "Perm. grassland (300 m)")
print(FigS4A)



## (B) Perm. grassland (300 m) -> Total abundance ####

FigS4B <-
  ggplot(data, aes(perm.grassland.r300, butterfly.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr300_abu_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(B) Total abundance",
       y = "", x = "Perm. grassland (300 m)")
print(FigS4B)



## (C) Perm. grassland (600 m) -> Total abundance ####

FigS4C <-
  ggplot(data, aes(perm.grassland.r600, butterfly.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr600_abu_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(C) Total abundance",
       y = "", x = "Perm. grassland (600 m)")
print(FigS4C)



## (D) Age x Shannon (600 m) -> Total abundance ####

FigS4D <- 
  ggplot(data, aes(rest.age, butterfly.abu, colour = region)) +
  geom_point() +
  geom_line(data = Pred_Lr600_Int_abu_tot,
            aes(y = butterfly.abu,
                linetype = lc_level),
            linewidth = 1) +
  my_theme() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  scale_linetype_manual(breaks = c("Low", "Mean", "High Landscape Shannon (600 m)"),
                        values = c("Low"  = "dotted",
                                   "Mean" = "twodash",
                                   "High Landscape Shannon (600 m)" = "solid")) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  labs(title = "", subtitle = "(D) Total abundance",
       y = "", x = "Age of restored sites")
print(FigS4D)



## (E) Perm. grassland (1200 m) -> Total abundance ####

FigS4E <- 
  ggplot(data, aes(perm.grassland.r1200, butterfly.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr1200_abu_tot, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(E) Total abundance",
       y = "", x = "Perm. grassland (1200 m)")
print(FigS4E)



## (F) Perm. grassland (300 m) -> Target abundance ####

FigS4F <-
  ggplot(data, aes(perm.grassland.r300, butterfly.target.abu, color = region )) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr300_abu_target, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(F) Target abundance",
       y = "", x = "Perm. grassland (300 m)")
print(FigS4F)



## (G) Age x Shannon (600 m) -> Target abundance ####

FigS4G <- 
  ggplot(data, aes(rest.age, butterfly.target.abu, colour = region)) +
  geom_point() +
  geom_line(data = Pred_Lr600_abu_target,
            aes(y = butterfly.target.abu,
                linetype = lc_level),
            linewidth = 1) +
  my_theme() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  scale_linetype_manual(breaks = c("Low", "Mean", "High Landscape Shannon (600 m)"),
                        values = c("Low"  = "dotted",
                                   "Mean" = "twodash",
                                   "High Landscape Shannon (600 m)" = "solid")) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(title = "", subtitle = "(G) Target abundance",
       y = "", x = "Age of restored sites")
print(FigS4G)



## (H) Grassland plant species richness (600 m) -> Total Hill-Simpson diversity ####

FigS4H <- 
  ggplot(data, aes(lk.grassland.per.min, butterfly.Hill2.total, color = region)) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr600_Hill2_tot, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(H) Total Hill-Simpson diversity",
       y = "", x = "Grassland plant species richness (600 m)")
print(FigS4H)



## (I) Age x Shannon (600 m) -> Total Hill-Simpson diversity ####

FigS4I <- 
  ggplot(data, aes(rest.age, butterfly.Hill2.total, colour = region)) +
  geom_point() +
  geom_line(data = Pred_Lr600_Int_Hill2_tot,
            aes(y = butterfly.Hill2.total,
                linetype = lc_level),
            linewidth = 1) +
  my_theme() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  scale_linetype_manual(breaks = c("Low", "Mean", "High Landscape Shannon (600 m)"),
                        values = c("Low"  = "dotted",
                                   "Mean" = "twodash",
                                   "High Landscape Shannon (600 m)" = "solid")) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(title = "", subtitle = "(I) Total Hill-Simpson diversity",
       y = "", x = "Age of restored sites")
print(FigS4I)



## (J) Grassland plant species richness (1200 m) -> Total Hill-Simpson diversity ####

FigS4J <- 
  ggplot(data, aes(lk.grassland.per.min, butterfly.Hill2.total, color = region)) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr1200_1_Hill2_tot, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(J) Total Hill-Simpson diversity",
       y = "", x = "Grassland plant species richness (1200 m)")
print(FigS4J)



## (K) Arable land (1200 m) -> Total Hill-Simpson diversity ####

FigS4K <- 
  ggplot(data, aes(agrar.r1200, butterfly.Hill2.total, color = region)) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr1200_2_Hill2_tot, linewidth = 1) +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(K) Total Hill-Simpson diversity",
       y = "", x = "Arable land (1200 m)")
print(FigS4K)



## (L) Avg. arable field size (1200 m) -> Total Hill-Simpson diversity ####

FigS4L <- 
  ggplot(data, aes(field.size.avg.r1200, butterfly.Hill2.total, color = region)) +
  geom_point() +
  my_theme() +
  geom_line(data = Pred_Lr1200_3_Hill2_tot, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  labs(title = "", subtitle = "(L) Total Hill-Simpson diversity",
       y = "", x = "Avg. arable field size (1200 m)")
print(FigS4L)



## (M) Age x Shannon (300 m) -> Target Hill-Simpson diversity ####

FigS4M <- 
  ggplot(data, aes(rest.age, butterfly.Hill2.target, colour = region)) +
  geom_point() +
  geom_line(data = Pred_Lr300_Hill2_target, aes(y = butterfly.Hill2.target,
                                                linetype = lc_level), linewidth = 1) +
  facet_wrap(~obs.year) +
  my_theme() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  scale_linetype_manual(breaks = c("Low", "Mean", "High Landscape Shannon (300 m)"),
                        values = c("Low"  = "dotted",
                                   "Mean" = "twodash",
                                   "High Landscape Shannon (300 m)" = "solid")) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(title = "", subtitle = "(M) Target Hill-Simpson diversity",
       y = "", x = "Age of restored sites")
print(FigS4M)



## (N) Age x Shannon (600 m) -> Target Hill-Simpson diversity ####

FigS4N <- 
  ggplot(data, aes(rest.age, butterfly.Hill2.target, colour = region)) +
  geom_point() +
  geom_line(data = Pred_Lr600_Hill2_target,
            aes(y = butterfly.Hill2.target,
                linetype = lc_level),
            linewidth = 1) +
  facet_wrap(~obs.year) +
  my_theme() +
  theme(legend.position = "bottom", legend.box = "vertical") +
  scale_color_manual(values = c("North" = "#bbc1c6", "Centre" = "#8e98a1", "South" = "#555b60")) +
  scale_linetype_manual(breaks = c("Low", "Mean", "High Landscape Shannon (600 m)"),
                        values = c("Low"  = "dotted",
                                   "Mean" = "twodash",
                                   "High Landscape Shannon (600 m)" = "solid")) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(title = "", subtitle = "(N) Target Hill-Simpson diversity",
       y = "", x = "Age of restored sites")
print(FigS4N)



## combine ####

FigS4_1 <- ggarrange(FigS4A,
                     FigS4B,
                     FigS4C,
                     FigS4D,
                     FigS4E,
                     FigS4F,
                     FigS4G,
                     FigS4H,
                     FigS4I, # part 1
                     ncol = 3, nrow = 3)
print(FigS4_1)
ggsave(here("outputs", "figures","FigS4_page1_landscape.png"),
       plot = FigS4_1, width = 21, height = 16, dpi = 800)


FigS4_2 <- ggarrange(FigS4J, 
                     FigS4K, 
                     FigS4L,
                     FigS4M,
                     FigS4N, # part 2
                     ncol = 3, nrow = 2)
print(FigS4_2)
ggsave(here("outputs", "figures","FigS4_page2_landscape.png"),
       plot = FigS4_2, width = 21, height = 10.66, dpi = 800)



# END ####