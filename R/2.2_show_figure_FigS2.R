#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Fig SUPP restoration method and age of restored sites
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 22 April 2026




# packages ####

library(here)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)

here()



# START ####

rm(list = ls())




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A LOAD DATA ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## data of 112 restored sites ####
data_restored_112 <- read_csv(
  here("data", "processed", "data_all.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(.default = "?")) %>% 
  dplyr::filter(site.type == "restored", 
                mang.app.NEU.MM.1.minus != "NA") %>%
  dplyr::select(id.site, rest.age, rest.meth, mang.app.NEU.MM.1.minus) %>%
  dplyr::mutate(rest.meth = factor(rest.meth, levels = c("cus", "res", "dih", "mga")))
str(data_restored_112)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B CREATE PLOT ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

table(data_restored_112$rest.meth)

my_theme <- function() {
  theme_test() +
    theme(
      text = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 14,  color = "black", face = "bold"),
      plot.subtitle = element_text(size = 18,  color = "black", face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(margin = ggplot2::margin(l = 2, unit = "mm"), size = 14),
      strip.background = element_rect(fill = "white", colour = "black")
      )
  }

plot_restmeth_age<- ggplot(data_restored_112, aes(x = rest.meth, y = rest.age)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  my_theme() +
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture\nn = 19",
      "Regional seed\nmixture\nn = 34",
      "Direct \nharvesting\nn = 37",
      "Management \nadaptation\nn = 22")
  ) +
  labs(title = "", subtitle = "Age of restored sites",
       y = "", x = "")
print(plot_restmeth_age)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C EXPORT PLOT ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ggsave(here("outputs", "figures","FigS2_rest-meth-age.png"),
       plot = plot_restmeth_age, width = 6, height = 6, dpi = 800)




# END ####