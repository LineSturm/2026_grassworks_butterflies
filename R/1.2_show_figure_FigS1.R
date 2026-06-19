#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# figure mass abundance of P. coridon
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 11 June 2026


# packages ####

library(tidyverse)
library(here)
library(ggplot2)
library(ggsignif)
library(ggpubr)



# START ####

rm(list = ls())




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A LOAD DATA ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## site environment ####
site_data <- read_csv(
  here("data", "raw", "data_site_environment.csv"),
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
    region = fct_relevel(region, "north", "centre", "south"),
    site.type = fct_relevel(site.type, "negative", "restored", "positive"),
    rest.meth = fct_relevel(rest.meth, "cus", "res", "dih", "mga"),
    rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "dih", "mga", "positive")
  ) %>%
  dplyr::select(id.site, region, hydrology) %>%
  distinct()
str(site_data)



## butterfly data ####
butterfly_data <- read_csv(
  here("data", "raw", "data_butterfly_raw.csv"), col_names = TRUE) %>%
  filter(
    # remove april observation
    observation != "April",
    # remove additionally found species 
    !sub.transect %in% c("N1", "N2", "Zusatz")) %>%
  summarise(butterfly.no=sum(butterfly.no), .by = c(id.site, butterfly.name)) 



## butterfly rl & traits ####
butterfly_rl_traits <- read_csv(
  here("data", "raw", "data_butterfly_rl_traits_grassworks.csv"), col_names = TRUE) %>%
  mutate(
    rang = as.factor(rang),
    rl.ger = as.factor(rl.ger),
    larval.feeding.type = as.factor(larval.feeding.type),
    habitat.preference = as.factor(habitat.preference),
    disp = as.factor(disp),
    target.species = as.factor(target.species)
  )
str(butterfly_rl_traits)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B PREPARE DATA ###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# add traits to butterfly data
butterfly_t_sum <- butterfly_data %>%
  left_join(butterfly_rl_traits %>%
              dplyr::select(butterfly.name, target.species, rang),
            by = "butterfly.name") 




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C PLOT #######################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my_theme <- function() {
  theme_test() + 
    theme(
      text = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 14,  color = "black", face = "bold"), 
      plot.subtitle = element_text(size = 20,  color = "black", face = "bold"), 
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(5, "mm"),
      legend.spacing.x = unit(0.5, "cm"),
      legend.spacing.y = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.text = element_text(margin = ggplot2::margin(l = 2, unit = "mm"), size = 16),
      panel.spacing.x = unit(0.25, "cm"),
      axis.text.x = element_text(angle = 90, vjust = 1, hjust=0)
      )
}



## (A) Total abundance ####

# illustrate abundance of P. coridon
butterfly_abundance_coridon <- butterfly_t_sum %>%
  filter(butterfly.name == "Polyommatus coridon") %>%
  group_by(id.site) %>%
  summarise(category = "Polyommatus coridon",
            butterfly.abu = sum(butterfly.no))

butterfly_abundance_all_others <- butterfly_t_sum %>%
  filter(butterfly.name != "Polyommatus coridon") %>%
  group_by(id.site) %>%
  summarise(category = "all others",
            butterfly.abu = sum(butterfly.no))

vis_coridon <- butterfly_abundance_coridon %>%
  rbind(butterfly_abundance_all_others) %>%
  left_join(site_data %>% select(id.site, region, hydrology), by = "id.site") %>%
  mutate(region = as.factor(region))

plot_vis_coridon <- vis_coridon %>%
  filter(region == "centre") %>%
  ggplot(aes(x = id.site, y = butterfly.abu, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("all others", bquote(italic("Polyommatus coridon")))) +
  my_theme() +
  labs(subtitle = "(A) Total abundance", y = "", x = "Sites of centre region")
plot(plot_vis_coridon)



## (B) Target abundance ####
butterfly_target_abundance_coridon <- butterfly_t_sum %>%
  filter(target.species == 1,
         butterfly.name == "Polyommatus coridon") %>%
  group_by(id.site) %>%
  summarise(category = "Polyommatus coridon",
            butterfly.abu = sum(butterfly.no))

butterfly_target_abundance_all_others <- butterfly_t_sum %>%
  filter(target.species == 1,
         butterfly.name != "Polyommatus coridon") %>%
  group_by(id.site) %>%
  summarise(category = "all others",
            butterfly.abu = sum(butterfly.no))

vis_target_coridon <- butterfly_target_abundance_coridon %>%
  rbind(butterfly_target_abundance_all_others) %>%
  left_join(site_data %>% select(id.site, region), by = "id.site") %>%
  mutate(region = as.factor(region))

plot_vis_target_coridon <- vis_target_coridon %>%
  filter(region == "centre") %>%
  ggplot(aes(x = id.site, y = butterfly.abu, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey", "black"),
                    labels = c("all others", bquote(italic("Polyommatus coridon")))) +
  my_theme() +
  labs(title = "", subtitle = "(B) Target abundance", 
       y = "", x = "Sites of centre region")
plot(plot_vis_target_coridon)



## combine ####
plot_vis_coridon_combined <- ggarrange(plot_vis_coridon, plot_vis_target_coridon,
                               ncol = 1, nrow = 2)
print(plot_vis_coridon_combined)
ggsave(here("outputs", "figures","FigS1_Polyommatus_coridon.png"),
       plot = plot_vis_coridon_combined, width = 15, height = 10, dpi = 800)




# END ####