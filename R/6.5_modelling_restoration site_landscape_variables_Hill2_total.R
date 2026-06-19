#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Influence of site and landscape variables - total Hill-Simpson
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 14 June 2026



# packages ####

library(here)
library(tidyverse)
library(DHARMa)
library(performance)
library(MASS)
library(emmeans)
library(glmmTMB)
library(TMB)
library(mgcv)
library(mgcViz)

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
  filter(site.type == "restored",
         mang.app.NEU.MM.1.minus != "NA")
str(data)

data_x_cus <- data %>%
  filter(rest.meth != "cus")

# here only blurred coordinates are included in the data, but exact coordinates 
# were used to test spatial auto-correlation.




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B MODELLING SITE VARIABLES ###################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm1 ####
lm1 <- lm(butterfly.Hill2.total ~ region + obs.year + region:obs.year
            + rest.age + cover.vegetation + site.cwm.pres.oek.f 
            + plant.target.richness + flower.cover
            + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
            data = data)
lm1_summary <- summary(lm1)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "summary_final_model_site.txt"))
lm1_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(lm1)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm1 <- simulateResiduals(lm1)
testDispersion(sr_lm1)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm1, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm1
# observed = -0.0057427, expected = -0.0090090, sd = 0.0301329, p-value = 0.9137
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
         "validation_final_model_site.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm1)) # ok
plotQQunif(sr_lm1)      # ok
plotResiduals(sr_lm1)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm1_anova <- car::Anova(lm1, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "anova_final_model_site.txt"))
lm1_anova
sink()


### f) export for plotting -----------------------------------------------------
# plant.target.richness **
newdat1 <- expand_grid(rest.age = mean(data$rest.age),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = seq(min(data$plant.target.richness),
                                                   max(data$plant.target.richness),
                                                   length = 100),
                       flower.cover = mean(data$flower.cover),
                       cover.vegetation = mean(data$cover.vegetation),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = "Centre",
                       obs.year = "2022") 
newdat1$butterfly.Hill2.total <- predict(lm1, newdata = newdat1, type = "response")            
newdat1 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_Hill2_tot.csv"), row.names = F)

# cover.vegetation *
newdat2 <- expand_grid(rest.age = mean(data$rest.age),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       cover.vegetation = seq(min(data$cover.vegetation),
                                                   max(data$cover.vegetation),
                                                   length = 100),
                       flower.cover = mean(data$flower.cover),
                       plant.target.richness = mean(data$plant.target.richness),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = "Centre",
                       obs.year = "2022") 
newdat2$butterfly.Hill2.total <- predict(lm1, newdata = newdat2, type = "response")            
newdat2 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var2_Hill2_tot.csv"), row.names = F)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C MODELLING LANDSCAPE VARIABLES (300 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# keep region for control 

## 1 lm2 ####
lm2 <- lm(butterfly.Hill2.total ~ region + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r300
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
          + field.size.avg.r300,
          data = data)
summary(lm2)


### a) check multicollinearity -------------------------------------------------
car::vif(lm2)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm2 <- simulateResiduals(lm2)
testDispersion(sr_lm2)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm2, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm2
# observed = -0.037997, expected = -0.009009, sd = 0.030147, p-value = 0.3363
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
hist(residuals(sr_lm2)) # ok
plotQQunif(sr_lm2)      # ok
plotResiduals(sr_lm2)   # not ok



## 2 lm2a ----
lm2a <- lm(log(butterfly.Hill2.total) ~ region + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r300
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
          + field.size.avg.r300,
          data = data)
lm2a_summary <- summary(lm2a)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "summary_final_model_landscape_r300.txt"))
lm2a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(lm2a)
# ok


### b) check overdispersion ----------------------------------------------------
sr_lm2a <- simulateResiduals(lm2a)
testDispersion(sr_lm2a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm2a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm2a
# observed = -0.015321, expected = -0.009009, sd = 0.030191, p-value = 0.8344
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
         "validation_final_model_landscape_r300.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm2a)) # ok
plotQQunif(sr_lm2a)      # ok
plotResiduals(sr_lm2a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm2a_anova <- car::Anova(lm2a, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "anova_final_model_landscape_r300.txt"))
lm2a_anova
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D MODELLING LANDSCAPE VARIABLES (600 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm3 ####
lm3 <- lm(butterfly.Hill2.total ~ region + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r600
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
          + field.size.avg.r600,
          data = data)
summary(lm3)


### a) check multicollinearity -------------------------------------------------
car::vif(lm3)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm3 <- simulateResiduals(lm3)
testDispersion(sr_lm3)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm3, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm3
# observed = -0.035746, expected = -0.009009, sd = 0.030145, p-value = 0.3751
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
hist(residuals(sr_lm3)) # ok
plotQQunif(sr_lm3)      # ok
plotResiduals(sr_lm3)   # not ok --> transformation



## 2 lm3a ####
lm3a <- lm(log(butterfly.Hill2.total) ~ region + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r600
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
          + field.size.avg.r600,
          data = data)
lm3a_summary <- summary(lm3a)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "summary_final_model_landscape_r600.txt"))
lm3a_summary
sink()
 

### a) check multicollinearity -------------------------------------------------
car::vif(lm3a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm3a <- simulateResiduals(lm3a)
testDispersion(sr_lm3a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm3a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm3a
# observed = -0.010738, expected = -0.009009, sd = 0.030184, p-value = 0.9543
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
         "validation_final_model_landscape_r600.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm3a)) # ok
plotQQunif(sr_lm3a)      # ok
plotResiduals(sr_lm3a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm3a_anova <- car::Anova(lm3a, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "anova_final_model_landscape_r600.txt"))
lm3a_anova
sink()


### f) export for plotting -----------------------------------------------------
# lk.grassland.per.min *
newdat3 <- expand_grid(rest.age = mean(data$rest.age),
                       plant.target.richness = mean(data$plant.target.richness),
                       cover.vegetation = mean(data$cover.vegetation),
                       agrar.r600 = mean(data$agrar.r600),
                       lk.grassland.per.min = seq(min(data$lk.grassland.per.min),
                                                  max(data$lk.grassland.per.min),
                                                  length = 100),
                       lc.shannon.div.r600 = mean(data$lc.shannon.div.r600),
                       field.size.avg.r600 = mean(data$field.size.avg.r600),
                       perm.grassland.r600 = mean(data$perm.grassland.r600),
                       region = unique(data$region)) 
sigma2 <- summary(lm3a)$sigma^2
newdat3$butterfly.Hill2.total <- 
  exp(predict(lm3a, newdata = newdat3) + 0.5 * sigma2)
newdat3 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr600_Hill2_tot.csv"), row.names = F)

# rest.age x lc.shannon.div.r600 (.)
lc_mean  <- mean(data$lc.shannon.div.r600)
lc_sd    <- sd(data$lc.shannon.div.r600)
newdat_int <- expand_grid(
  rest.age = seq(min(data$rest.age),
                 max(data$rest.age),
                 length = 100),
  lc.shannon.div.r600 = c(lc_mean - lc_sd,
                          lc_mean,
                          lc_mean + lc_sd),
  plant.target.richness = mean(data$plant.target.richness),
  cover.vegetation = mean(data$cover.vegetation),
  agrar.r600 = mean(data$agrar.r600),
  site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
  region = unique(data$region),
  perm.grassland.r600 = mean(data$perm.grassland.r600),
  lk.grassland.per.min = mean(data$lk.grassland.per.min),
  field.size.avg.r600 = mean(data$field.size.avg.r600))
sigma2 <- summary(lm3a)$sigma^2
newdat_int$butterfly.Hill2.total <- 
  exp(predict(lm3a, newdata = newdat_int) + 0.5 * sigma2)
newdat_int$lc_level <- factor(newdat_int$lc.shannon.div.r600,
                              labels = c("Low",
                                         "Mean",
                                         "High Landscape Shannon (600 m)"))
newdat_int %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr600_Int_Hill2_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E MODELLING LANDSCAPE VARIABLES (1200 m) #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm4 ####
lm4 <- lm(butterfly.Hill2.total ~ region + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r1200
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
          + field.size.avg.r1200,
          data = data)
summary(lm4)


### a) check multicollinearity -------------------------------------------------
car::vif(lm4)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm4 <- simulateResiduals(lm4)
testDispersion(sr_lm4)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm4, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm4
# observed = -0.029576, expected = -0.009009, sd = 0.030144, p-value = 0.4951
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
hist(residuals(sr_lm4)) # ok
plotQQunif(sr_lm4)      # ok
plotResiduals(sr_lm4)   # not ok



## 2 lm4a ####
lm4a <- lm(log(butterfly.Hill2.total) ~ region + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r1200
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
          + field.size.avg.r1200,
          data = data)
lm4a_summary <- summary(lm4a)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "summary_final_model_landscape_r1200.txt"))
lm4a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(lm4a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm4a <- simulateResiduals(lm4a)
testDispersion(sr_lm4a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm4a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm4a
# observed = -0.0096123, expected = -0.0090090, sd = 0.0301789, p-value = 0.9841
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
         "validation_final_model_landscape_r1200.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm4a)) # ok
plotQQunif(sr_lm4a)      # ok
plotResiduals(sr_lm4a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm4a_anova <- car::Anova(lm4a, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "anova_final_model_landscape_r1200.txt"))
lm4a_anova
sink()


### f) export for plotting -----------------------------------------------------
# lk.grassland.per.min *
newdat4 <- expand_grid(rest.age = mean(data$rest.age),
                       plant.target.richness = mean(data$plant.target.richness),
                       cover.vegetation = mean(data$cover.vegetation),
                       lk.grassland.per.min = seq(min(data$lk.grassland.per.min),
                                                  max(data$lk.grassland.per.min),
                                                  length = 100),
                       agrar.r1200 = mean(data$agrar.r1200),
                       lc.shannon.div.r1200 = mean(data$lc.shannon.div.r1200),
                       field.size.avg.r1200 = mean(data$field.size.avg.r1200),
                       perm.grassland.r1200 = mean(data$perm.grassland.r1200),
                       region = unique(data$region)) 
sigma2 <- summary(lm4a)$sigma^2
newdat4$butterfly.Hill2.total <- 
  exp(predict(lm4a, newdata = newdat4) + 0.5 * sigma2)
newdat4 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr1200_1_Hill2_tot.csv"), row.names = F)

# agrar.r1200 (.)
newdat5 <- expand_grid(rest.age = mean(data$rest.age),
                       plant.target.richness = mean(data$plant.target.richness),
                       cover.vegetation = mean(data$cover.vegetation),
                       agrar.r1200 = seq(min(data$agrar.r1200),
                                                  max(data$agrar.r1200),
                                                  length = 100),
                       lk.grassland.per.min = mean(data$lk.grassland.per.min),
                       lc.shannon.div.r1200 = mean(data$lc.shannon.div.r1200),
                       field.size.avg.r1200 = mean(data$field.size.avg.r1200),
                       perm.grassland.r1200 = mean(data$perm.grassland.r1200),
                       region = unique(data$region)) 
sigma2 <- summary(lm4a)$sigma^2
newdat5$butterfly.Hill2.total <- 
  exp(predict(lm4a, newdata = newdat5) + 0.5 * sigma2)
newdat5 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr1200_2_Hill2_tot.csv"), row.names = F)

# field.size.avg.r1200 (.)
newdat6 <- expand_grid(rest.age = mean(data$rest.age),
                       plant.target.richness = mean(data$plant.target.richness),
                       cover.vegetation = mean(data$cover.vegetation),
                       field.size.avg.r1200 = seq(min(data$field.size.avg.r1200),
                                         max(data$field.size.avg.r1200),
                                         length = 100),
                       lk.grassland.per.min = mean(data$lk.grassland.per.min),
                       lc.shannon.div.r1200 = mean(data$lc.shannon.div.r1200),
                       agrar.r1200 = mean(data$agrar.r1200),
                       perm.grassland.r1200 = mean(data$perm.grassland.r1200),
                       region = unique(data$region)) 
sigma2 <- summary(lm4a)$sigma^2
newdat6$butterfly.Hill2.total <- 
  exp(predict(lm4a, newdata = newdat6) + 0.5 * sigma2)
newdat6 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr1200_3_Hill2_tot.csv"), row.names = F)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F MODELLING SITE VARIABLES (excluding cultivar seed mixtures) ################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm_no_cus ####
lm_no_cus <- lm(butterfly.Hill2.total ~ region + obs.year + region:obs.year
                + rest.age + site.cwm.pres.oek.f  + cover.vegetation
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data_x_cus)
lm_no_cus_summary <- summary(lm_no_cus)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "summary_final_model_no_cus.txt"))
lm_no_cus_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(lm_no_cus)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm_no_cus <- simulateResiduals(lm_no_cus)
testDispersion(sr_lm_no_cus)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm_no_cus, x = data_x_cus$x, y = data_x_cus$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm_no_cus
# observed = -0.018743, expected = -0.010870, sd = 0.035409, p-value = 0.824
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
         "validation_final_model_no_cus.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm_no_cus)) # ok
plotQQunif(sr_lm_no_cus)      # ok
plotResiduals(sr_lm_no_cus)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm_no_cus_anova <- car::Anova(lm_no_cus, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_tot", 
          "anova_final_model_no_cus.txt"))
lm_no_cus_anova
sink()


### f) export for plotting -----------------------------------------------------
# rest.age (.)
newdat7 <- expand_grid(cover.vegetation = mean(data_x_cus$cover.vegetation),
                       site.cwm.pres.oek.f = mean(data_x_cus$site.cwm.pres.oek.f),
                       rest.age = seq(min(data_x_cus$rest.age),
                                      max(data_x_cus$rest.age),
                                      length = 100),
                       flower.cover = mean(data_x_cus$flower.cover),
                       plant.target.richness = mean(data_x_cus$plant.target.richness),
                       veg.heterogeneity = mean(data_x_cus$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data_x_cus$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = "Centre",
                       obs.year = "2022") 
newdat7$butterfly.Hill2.total <- predict(lm_no_cus, newdata = newdat7, type = "response")            
newdat7 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_noCus_Hill2_tot.csv"), row.names = F)




# END ####