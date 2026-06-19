#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# data exploration
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 13 June 2026




# packages ####

library(here)
library(tidyverse)
library(visdat)
library(naniar)
library(lattice)
library(rstatix)
library(Hmisc)
library(ggbeeswarm)
library(ggpubr)
library(openxlsx)

here()



# START ####

rm(list = ls())




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A LOAD DATA ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## data all ####
data_all <- read_csv(
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
    region = fct_relevel(region, "north", "centre", "south"),
    site.type = fct_relevel(site.type, "negative", "restored", "positive"),
    rest.meth = fct_relevel(rest.meth, "cus", "res", "dih", "mga"),
    rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "dih", "mga", "positive")
    )
str(data_all)



## data restored ####
data_restored <- data_all %>% dplyr::filter(site.type == "restored")
str(data_restored)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B DATA EXPLORATION ###########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Zuur et al. (2010) Methods Ecol Evol 
# https://doi.org/10.1111/2041-210X.12577



## 1 missing values ####

colSums(is.na(data_all)) 
# 66 reference sites have NA in rest.meth, land.use.hist, rest.age
# 1 mngm.type
# 75 management appropriateness (mang.app.NEU.MM.1.minus) (66 due to reference sites)



## 2 data distribution, outliers, zero-inflation ####


### a) butterfly response variables --------------------------------------------

# data distribution
# -> often right skrewed
ggplot(data_all, aes(x = butterfly.rich.total)) + geom_histogram(binwidth = 1)

ggplot(data_all, aes(x = butterfly.rich.target)) + geom_histogram(binwidth = 1)

ggplot(data_all, aes(x = butterfly.abu)) + geom_histogram(binwidth = 1)

ggplot(data_all, aes(x = butterfly.target.abu)) + geom_histogram(binwidth = 1)

ggplot(data_all, aes(x = butterfly.Hill2.total)) + geom_histogram(binwidth = 1)

ggplot(data_all, aes(x = butterfly.Hill2.target)) + geom_histogram(binwidth = 1)

# test for outliers
data_all %>% 
  dplyr::select(id.site, butterfly.rich.total) %>% 
  identify_outliers(butterfly.rich.total)

data_all %>% 
  dplyr::select(id.site, butterfly.rich.target) %>% 
  identify_outliers(butterfly.rich.target) 
# extreme outlier: M_TTH

data_all %>% 
  dplyr::select(id.site, butterfly.abu) %>% 
  identify_outliers(butterfly.abu)

data_all %>% 
  dplyr::select(id.site, butterfly.target.abu) %>% 
  identify_outliers(butterfly.target.abu)
# extreme outliers:
# S_JAU: many P. icarus
# M_HIR: many P. icarus
# M_TTO: many target species in general
# M_TRE: many P. icarus
# M_TTH: many P. icarus

data_all %>% 
  dplyr::select(id.site, butterfly.Hill2.total) %>% 
  identify_outliers(butterfly.Hill2.total)
# extreme outlier: M_HNB

data_all %>% 
  dplyr::select(id.site, butterfly.Hill2.target) %>% 
  identify_outliers(butterfly.Hill2.target)

# handling:
# -> values are reasonable, do not change input data 


### b) landscape variables of restored sites at the 1200-m radius --------------
mean(data_restored$agrar.r1200 == 0)
ggplot(data_restored, aes(x = agrar.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$cereal.r1200 == 0)
ggplot(data_restored, aes(x = cereal.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$rape.r1200 == 0)
ggplot(data_restored, aes(x = rape.r1200)) + geom_histogram(binwidth = 1) 
# tendency to zero-inflation

mean(data_restored$legume.r1200 == 0)
ggplot(data_restored, aes(x = legume.r1200)) + geom_histogram(binwidth = 1) 

mean(data_restored$special.crop.r1200 == 0)
ggplot(data_restored, aes(x = special.crop.r1200)) + geom_histogram(binwidth = 1) 

mean(data_restored$cult.grassland.r1200 == 0)
ggplot(data_restored, aes(x = cult.grassland.r1200)) + geom_histogram(binwidth = 1)
# zero-inflated

mean(data_restored$perm.grassland.r1200 == 0)
ggplot(data_restored, aes(x = perm.grassland.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$anthrop.herb.area.r1200 == 0)
ggplot(data_restored, aes(x = anthrop.herb.area.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$heathland.r1200 == 0)
ggplot(data_restored, aes(x = heathland.r1200)) + geom_histogram(binwidth = 1) 
# zero-inflated

mean(data_restored$margin.r1200 == 0)
ggplot(data_restored, aes(x = margin.r1200)) + geom_histogram(binwidth = 1) 
# zero-inflated

mean(data_restored$wetland.r1200 == 0)
ggplot(data_restored, aes(x = wetland.r1200)) + geom_histogram(binwidth = 1) 
# zero-inflated

mean(data_restored$mire.r1200 == 0)
ggplot(data_restored, aes(x = mire.r1200)) + geom_histogram(binwidth = 1) 
# zero-inflated

mean(data_restored$tree.avenue.r1200 == 0)
ggplot(data_restored, aes(x = tree.avenue.r1200)) + geom_histogram(binwidth = 1) 
# zero-inflated

mean(data_restored$forest.r1200 == 0)
ggplot(data_restored, aes(x = forest.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$grove.r1200 == 0)
ggplot(data_restored, aes(x = grove.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$hedge.r1200 == 0)
ggplot(data_restored, aes(x = hedge.r1200)) + geom_histogram(binwidth = 1) 
# zero-inflated

mean(data_restored$urban.r1200 == 0)
ggplot(data_restored, aes(x = urban.r1200)) + geom_histogram(binwidth = 1)

mean(data_restored$water.r1200 == 0)
ggplot(data_restored, aes(x = water.r1200)) + geom_histogram(binwidth = 1)


### c) summary and normality test ----------------------------------------------

# all sites
data_summary_distribution <- data_all %>%
  dplyr::select(-rest.age) %>%
  # Compute required summary statistics
  summarise(across(where(is.numeric), list(
    min = ~ round(min(.x, na.rm = TRUE),2),
    max = ~ round(max(.x, na.rm = TRUE),2),
    med = ~ round(median(.x, na.rm = TRUE),2),
    mean = ~ round(mean(.x, na.rm = TRUE),2),
    sd = ~ round(sd(.x, na.rm = TRUE),2),
    shapiro.statistic = ~ round(shapiro.test(.x)$statistic,4),
    shapiro.p.value = ~ round(shapiro.test(.x)$p.value, 4)), 
    .names = '{col}_{fn}')) %>%
  # Convert the data to long format
  pivot_longer(cols = everything(), names_to = c('variable.name', '.value'), 
               names_sep = '_') %>%
  mutate(data.set = "all sites", .before = variable.name) %>%
  mutate(distribution = ifelse(shapiro.p.value < 0.05, "non-normal", "normal"))

# restored sites
restored_summary_distribution <- data_restored %>%
  # Compute required summary statistics
  summarise(across(where(is.numeric), list(
    min = ~ round(min(.x, na.rm = TRUE),2),
    max = ~ round(max(.x, na.rm = TRUE),2),
    med = ~ round(median(.x, na.rm = TRUE),2),
    mean = ~ round(mean(.x, na.rm = TRUE),2),
    sd = ~ round(sd(.x, na.rm = TRUE),2),
    shapiro.statistic = ~ round(shapiro.test(.x)$statistic,4),
    shapiro.p.value = ~ round(shapiro.test(.x)$p.value, 4)), 
    .names = '{col}_{fn}')) %>%
  # Convert the data to long format
  pivot_longer(cols = everything(), names_to = c('variable.name', '.value'), 
               names_sep = '_') %>%
  mutate(data.set = "restored sites", .before = variable.name) %>%
  mutate(distribution = ifelse(shapiro.p.value < 0.05, "non-normal", "normal"))

# restored sites without NA values of management appropriateness 
restored_112_summary_distribution <- data_all %>%
  filter(site.type == "restored",
         mang.app.NEU.MM.1.minus != "NA") %>% 
  # Compute required summary statistics
  summarise(across(where(is.numeric), list(
    min = ~ round(min(.x, na.rm = TRUE),2),
    max = ~ round(max(.x, na.rm = TRUE),2),
    med = ~ round(median(.x, na.rm = TRUE),2),
    mean = ~ round(mean(.x, na.rm = TRUE),2),
    sd = ~ round(sd(.x, na.rm = TRUE),2),
    shapiro.statistic = ~ round(shapiro.test(.x)$statistic,4),
    shapiro.p.value = ~ round(shapiro.test(.x)$p.value, 4)), 
    .names = '{col}_{fn}')) %>%
  # Convert the data to long format
  pivot_longer(cols = everything(), names_to = c('variable.name', '.value'),
               names_sep = '_') %>%
  mutate(data.set = "112 restored sites", .before = variable.name) %>%
  mutate(distribution = ifelse(shapiro.p.value < 0.05, "non-normal", "normal"))



## 3 inspect categorical covariates ####

# all sites
table(data_all$region)
table(data_all$obs.year)

# restored sites
table(data_restored$region)
table(data_restored$obs.year)
table(data_restored$rest.meth) # Unbalanced ...but enough observations per level.
table(data_restored$rest.age) # Unbalanced ...not enough observations per level.
table(data_restored$mngm.type) # Unbalanced ...but enough observations per level.


# Was each restoration methods measured in every region?
table(data_restored$rest.meth, data_restored$region)
# Unbalanced, due to regional differences
# South: no mga, less cus
# Centre: less cus, even more less dih
# North: most balanced

# Was each land use history measured in every region?
table(data_restored$land.use.hist, data_restored$region)

# Was each land use history measured in every restoration method?
table(data_restored$land.use.hist, data_restored$rest.meth)
# Unbalanced: sowing (cus, res) mostly on arable land, mga mostly on grassland
# pattern expected

# Was each management type measured in every region?
table(data_all$mngm.type, data_all$region)
# Unbalanced, due to regional differences

# Was each management type measured in every restoration method?
table(data_restored$mngm.type, data_restored$rest.meth)
# Unbalanced, pattern expected 



## 4 check collinearity between continuous covariates ####

# Dormann et al. 2013 Ecography https://doi.org/10.1111/j.1600-0587.2012.07348.x
# We chose r > 0.6 (more conservative) as threshold 


### a) corr for 112 restored sites ---------------------------------------------

# zero-inflated landscape variables:
# cult.grassland, heathland, margin, wetland, mire, tree.avenue, hedge 

corr_restored <- data_restored %>% 
  dplyr::select(
    # butterfly response variables
    18:23, 25, 27, 
    # site and vegetation variables 
    13, 30:40,
    # landscape variables 300-m radius (excluding zero-inflated variables)
    41:49, 51, 57, 58, 60:62,
    # landscape variables 600-m radius (excluding zero-inflated variables)
    63:71, 73, 79, 80, 82:84,
    # landscape variables 1200-m radius (excluding zero-inflated variables)
    85:93, 95, 101, 102, 104:106)
corr_restored <- na.omit(corr_restored)
str(corr_restored)
corr_restored_result <- rcorr(as.matrix(corr_restored), type = c("spearman"))


### b) between continuous variables and factors --------------------------------

# explore patterns for plant.richness per region
ggplot(data_all, aes(x = region, y = plant.richness)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent")
# plant diversity increases from North to South

# explore patterns for plant.richness per restoration methods
ggplot(data_restored, aes(x = rest.meth, y = plant.richness)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent")
# restoration methods vary in their effectiveness

# explore patterns for restoration age per restoration methods
ggplot(data_restored, aes(x = rest.meth, y = rest.age)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent")
# restoration method interlinked with restoration age


### c) relationships between response and explanatory variables ----------------

# Plot response variables versus region, obs.year, site.type, rest.meth, age

# total species richness
ggplot(data_all, aes(x = region, y = butterfly.rich.total)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total species richness and region", y = "species richness", x = "")
ggplot(data_all, aes(x = obs.year, y = butterfly.rich.total)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total species richness and observation year", y = "species richness", x = "")
ggplot(data_restored, aes(x = rest.meth, y = butterfly.rich.total)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total species richness and restoration method", y = "species richness", x="")
ggplot(data_restored, aes(x = rest.age, y = butterfly.rich.total)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Total species richness and age of restoration",  y = "species richness", x="")

# total abundance 
ggplot(data_all, aes(x = region, y = butterfly.abu)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total abundance and region", y = "abundance", x = "")
ggplot(data_all, aes(x = obs.year, y = butterfly.abu)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total abundance and observation year", y = "species richness", x = "")
ggplot(data_restored, aes(x = rest.meth, y = butterfly.abu)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total abundance and restoration method", y = "abundance", x="")
ggplot(data_restored, aes(x = rest.age, y = butterfly.abu)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Total abundance and age of restoration",  y = "abundance", x="")

# total Hill-Simpson diversity
ggplot(data_all, aes(x = region, y = butterfly.Hill2.total)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total Hill-Simpson diversity and region", y = "abundance", x = "")
ggplot(data_all, aes(x = obs.year, y = butterfly.Hill2.total)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total Hill-Simpson diversity and observation year", y = "species richness", x = "")
ggplot(data_restored, aes(x = rest.meth, y = butterfly.Hill2.total)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Total Hill-Simpson diversity and restoration method", y = "abundance", x="")
ggplot(data_restored, aes(x = rest.age, y = butterfly.Hill2.total)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Total Hill-Simpson diversity and age of restoration",  y = "abundance", x="")

# target species richness
ggplot(data_all, aes(x = region, y = butterfly.rich.target)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target species richness and region", y = "species richness", x = "")
ggplot(data_all, aes(x = obs.year, y = butterfly.rich.target)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target species richness and observation year", y = "species richness", x = "")
ggplot(data_restored, aes(x = rest.meth, y = butterfly.rich.target)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target species richness and restoration method", y = "species richness", x="")
ggplot(data_restored, aes(x = rest.age, y = butterfly.rich.target)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Target species richness and age of restoration",  y = "species richness", x="")

# target abundance 
ggplot(data_all, aes(x = region, y = butterfly.target.abu)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target abundance and region", y = "abundance", x = "")
ggplot(data_all, aes(x = obs.year, y = butterfly.target.abu)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target abundance and observation year", y = "species richness", x = "")
ggplot(data_restored, aes(x = rest.meth, y = butterfly.target.abu)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target abundance and restoration method", y = "abundance", x="")
ggplot(data_restored, aes(x = rest.age, y = butterfly.target.abu)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Target abundance and age of restoration",  y = "abundance", x="")

# target Hill-Simpson diversity
ggplot(data_all, aes(x = region, y = butterfly.Hill2.target)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target Hill-Simpson diversity and region", y = "abundance", x = "")
ggplot(data_all, aes(x = obs.year, y = butterfly.Hill2.target)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target Hill-Simpson diversity and observation year", y = "species richness", x = "")
ggplot(data_restored, aes(x = rest.meth, y = butterfly.Hill2.target)) +
  geom_quasirandom(color = "grey") + 
  geom_boxplot(fill = "transparent") +
  labs(title = "Target Hill-Simpson diversity and restoration method", y = "abundance", x="")
ggplot(data_restored, aes(x = rest.age, y = butterfly.Hill2.target)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Target Hill-Simpson diversity and age of restoration",  y = "abundance", x="")


### d) interaction between region and observation year -------------------------

# accounting for regional differences
# -> often higher butterfly species richness in 2023 in north, opposite in south

data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.rich.total)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.rich.target)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.abu)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.target.abu)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.Hill2.total)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.Hill2.target)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")

data_restored %>% 
  ggplot(aes(x = obs.year, y = butterfly.rich.total)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_restored %>% 
  ggplot(aes(x = obs.year, y = butterfly.rich.target)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_restored %>% 
  ggplot(aes(x = obs.year, y = butterfly.abu)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_restored %>% 
  ggplot(aes(x = obs.year, y = butterfly.target.abu)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_restored %>% 
  ggplot(aes(x = obs.year, y = butterfly.Hill2.total)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")
data_all %>% 
  ggplot(aes(x = obs.year, y = butterfly.Hill2.target)) +
  geom_boxplot(show.legend = F) +
  facet_grid(cols=vars(region)) +
  stat_compare_means(method = "kruskal.test")




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C EXPORT #####################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

spearman_res_r=data.frame(corr_restored_result$r)
spearman_res_r %>% write.xlsx(here("outputs", "tables",
                                    "butterfly_correlation_matrix_spearman_r.xlsx"))

spearman_res_p=data.frame(corr_restored_result$P)
spearman_res_p %>% write.xlsx(here("outputs", "tables",
                                    "butterfly_correlation_matrix_spearman_p.xlsx"))



# END ####