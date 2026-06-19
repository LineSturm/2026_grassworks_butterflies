#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Influence of site and landscape variables - target Hill-Simpson
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
lm1 <- lm(butterfly.Hill2.target ~ region + obs.year + region:obs.year
            + rest.age + cover.vegetation + site.cwm.pres.oek.f 
            + plant.target.richness + flower.cover
            + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
            data = data)
lm1_summary <- summary(lm1)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
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
# observed = -0.032455, expected = -0.009009, sd = 0.030119, p-value = 0.4363
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
         "validation_final_model_site.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm1)) # ok
plotQQunif(sr_lm1)      # ok
plotResiduals(sr_lm1)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm1_anova <- car::Anova(lm1, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
          "anova_final_model_site.txt"))
lm1_anova
sink()


### f) export for plotting -----------------------------------------------------
# plant.target.richness ***
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
                       region = unique(data$region),
                       obs.year = unique(data$obs.year)) 
newdat1$butterfly.Hill2.target <- predict(lm1, newdata = newdat1, type = "response")            
newdat1 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_Hill2_target.csv"), row.names = F)

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
                       region = unique(data$region),
                       obs.year = unique(data$obs.year)) 
newdat2$butterfly.Hill2.target <- predict(lm1, newdata = newdat2, type = "response")            
newdat2 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var2_Hill2_target.csv"), row.names = F)

# rest.age **
newdat3 <- expand_grid(rest.age = seq(min(data$rest.age),
                                      max(data$rest.age),
                                      length = 100),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = mean(data$plant.target.richness),
                       flower.cover = mean(data$flower.cover),
                       cover.vegetation = mean(data$cover.vegetation),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data$region),
                       obs.year = unique(data$obs.year)) 
newdat3$butterfly.Hill2.target <- predict(lm1, newdata = newdat3, type = "response")            
newdat3 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var3_Hill2_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C MODELLING LANDSCAPE VARIABLES (300 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm2 ####
lm2 <- lm(butterfly.Hill2.target ~ region + obs.year + region:obs.year
          + plant.target.richness + cover.vegetation
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
# observed = -0.052238, expected = -0.009009, sd = 0.030111, p-value = 0.1511
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
hist(residuals(sr_lm2)) # ok
plotQQunif(sr_lm2)      # not ok
plotResiduals(sr_lm2)   # not ok



## 2 lm2a ####
lm2a <- lm(log(butterfly.Hill2.target+1) ~ region + obs.year + region:obs.year
          + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r300
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
          + field.size.avg.r300,
          data = data)
lm2a_summary <- summary(lm2a)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
          "summary_final_model_landscape_r300.txt"))
lm2a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(lm2a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm2a <- simulateResiduals(lm2a)
testDispersion(sr_lm2a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm2a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm2a
# observed = -0.064887, expected = -0.009009, sd = 0.030152, p-value = 0.06385
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
         "validation_final_model_landscape_r300.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm2a)) # ok
plotQQunif(sr_lm2a)      # ok
plotResiduals(sr_lm2a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm2a_anova <- car::Anova(lm2a, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
          "anova_final_model_landscape_r300.txt"))
lm2a_anova
sink()


### f) export for plotting -----------------------------------------------------
# rest.age x lc.shannon.div.r300 (.)
lc_mean  <- mean(data$lc.shannon.div.r300)
lc_sd    <- sd(data$lc.shannon.div.r300)
newdat4 <- expand_grid(
  rest.age = seq(min(data$rest.age),
                 max(data$rest.age),
                 length = 100),
  lc.shannon.div.r300 = c(lc_mean - lc_sd,
                          lc_mean,
                          lc_mean + lc_sd),
  plant.target.richness = mean(data$plant.target.richness),
  cover.vegetation = mean(data$cover.vegetation),
  agrar.r300 = mean(data$agrar.r300),
  site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
  region = unique(data$region),
  perm.grassland.r300 = mean(data$perm.grassland.r300),
  lk.grassland.per.min = mean(data$lk.grassland.per.min),
  field.size.avg.r300 = mean(data$field.size.avg.r300),
  obs.year = unique(data$obs.year))
sigma2 <- summary(lm2a)$sigma^2
newdat4$butterfly.Hill2.target <- 
  exp(predict(lm2a, newdata = newdat4) + 0.5 * sigma2) - 1
newdat4$lc_level <- factor(newdat4$lc.shannon.div.r300,
                           labels = c("Low",
                                      "Mean",
                                      "High Landscape Shannon (300 m)"))
newdat4 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr300_Hill2_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D MODELLING LANDSCAPE VARIABLES (600 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm3 ####
lm3 <- lm(butterfly.Hill2.target ~ region + obs.year + region:obs.year
          + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r600
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
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
# observed = -0.047476, expected = -0.009009, sd = 0.030116, p-value = 0.2015
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
hist(residuals(sr_lm3)) # ok
plotQQunif(sr_lm3)      # not ok
plotResiduals(sr_lm3)   # not ok



## 2 lm3a ####
lm3a <- lm(log(butterfly.Hill2.target+1) ~ region + obs.year + region:obs.year
          + plant.target.richness + cover.vegetation
          + lk.grassland.per.min + agrar.r600
          + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
          + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
          + field.size.avg.r600,
          data = data)
lm3a_summary <- summary(lm3a)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
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
# observed = -0.058685, expected = -0.009009, sd = 0.030162, p-value = 0.09956
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
         "validation_final_model_landscape_r600.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm3a)) # ok
plotQQunif(sr_lm3a)      # ok
plotResiduals(sr_lm3a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm3a_anova <- car::Anova(lm3a, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
          "anova_final_model_landscape_r600.txt"))
lm3a_anova
sink()


### f) export for plotting -----------------------------------------------------
# rest.age x lc.shannon.div.r600 (.)
lc_mean  <- mean(data$lc.shannon.div.r600)
lc_sd    <- sd(data$lc.shannon.div.r600)
newdat5 <- expand_grid(
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
  field.size.avg.r600 = mean(data$field.size.avg.r600),
  obs.year = unique(data$obs.year))
sigma2 <- summary(lm3a)$sigma^2
newdat5$butterfly.Hill2.target <- 
  exp(predict(lm3a, newdata = newdat5) + 0.5 * sigma2) - 1
newdat5$lc_level <- factor(newdat5$lc.shannon.div.r600,
                              labels = c("Low",
                                         "Mean",
                                         "High Landscape Shannon (600 m)"))
newdat5 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr600_Hill2_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E MODELLING LANDSCAPE VARIABLES (1200 m) #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm4 ####
lm4 <- lm(butterfly.Hill2.target ~ region + obs.year + region:obs.year
          + plant.target.richness + cover.vegetation
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
# observed = -0.032959, expected = -0.009009, sd = 0.030129, p-value = 0.4267
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
         "validation_final_model_landscape_r1200.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm4)) # ok
plotQQunif(sr_lm4)      # ok
plotResiduals(sr_lm4)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm4_anova <- car::Anova(lm4, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
          "anova_final_model_landscape_r1200.txt"))
lm4_anova
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F MODELLING SITE VARIABLES (excluding cultivar seed mixtures) ################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm_no_cus ####
lm_no_cus <- lm(butterfly.Hill2.target ~ region + obs.year + region:obs.year
                + rest.age + site.cwm.pres.oek.f  + cover.vegetation
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data_x_cus)
lm_no_cus_summary <- summary(lm_no_cus)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
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
# observed = -0.017306, expected = -0.010870, sd = 0.035390, p-value = 0.8557
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
         "validation_final_model_no_cus.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm_no_cus)) # ok
plotQQunif(sr_lm_no_cus)      # ok
plotResiduals(sr_lm_no_cus)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm_no_cus_anova <- car::Anova(lm_no_cus, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
          "anova_final_model_no_cus.txt"))
lm_no_cus_anova
sink()


### f) export for plotting -----------------------------------------------------
# rest.age *
newdat6 <- expand_grid(rest.age = seq(min(data_x_cus$rest.age),
                                      max(data_x_cus$rest.age),
                                      length = 100),
                       site.cwm.pres.oek.f = mean(data_x_cus$site.cwm.pres.oek.f),
                       plant.target.richness = mean(data_x_cus$plant.target.richness),
                       flower.cover = mean(data_x_cus$flower.cover),
                       cover.vegetation = mean(data_x_cus$cover.vegetation),
                       veg.heterogeneity = mean(data_x_cus$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data_x_cus$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data_x_cus$region),
                       obs.year = unique(data_x_cus$obs.year)) 
newdat6$butterfly.Hill2.target <- predict(lm_no_cus, newdata = newdat6, type = "response")            
newdat6 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_noCus_Hill2_target.csv"), row.names = F)




# END ####