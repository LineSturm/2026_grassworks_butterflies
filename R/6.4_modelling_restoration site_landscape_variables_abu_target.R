#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Q2 Influence of site and landscape variables - target abundance
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

## 1 glm1 ####
glm1 <- glm(butterfly.target.abu ~ region + obs.year + region:obs.year
            + rest.age + site.cwm.pres.oek.f + plant.target.richness + flower.cover
            + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
            data = data , family = "poisson")
summary(glm1)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1)
# exclude region:obs.year interaction



## 2 glm1a ####
glm1a <- update(glm1, .~. -obs.year:region) 
summary(glm1a)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm1a <- simulateResiduals(glm1a)
testDispersion(sr_glm1a) 
# switch to nb



## 3 glm1b ####
glm1b <- glm.nb(butterfly.target.abu ~ region + obs.year
                + rest.age + site.cwm.pres.oek.f + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data)
glm1b_summary <- summary(glm1b)
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "summary_final_model_site.txt"))
glm1b_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm1b)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm1b <- simulateResiduals(glm1b)
testDispersion(sr_glm1b)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm1b, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm1b
# observed = -0.020514, expected = -0.009009, sd = 0.030163, p-value = 0.7029
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_target", 
         "validation_final_model_site.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm1b)) # ok
plotQQunif(sr_glm1b)      # ok
plotResiduals(sr_glm1b)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm1b)
# R2m       R2c
# delta     0.4270164 0.4270164
# lognormal 0.5027619 0.5027619
# trigamma  0.3335290 0.3335290


### f) test fixed effects ------------------------------------------------------
glm1b_anova <- car::Anova(glm1b, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "anova_final_model_site.txt"))
glm1b_anova
sink()


### g) export for plotting -----------------------------------------------------
# site.cwm.pres.oek.f **
newdat1 <- expand_grid(rest.age = mean(data$rest.age),
                       plant.target.richness = mean(data$plant.target.richness),
                       site.cwm.pres.oek.f = seq(min(data$site.cwm.pres.oek.f),
                                                   max(data$site.cwm.pres.oek.f),
                                                   length = 100),
                       flower.cover = mean(data$flower.cover),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data$region),
                       obs.year = "2022") 
newdat1$butterfly.target.abu <- predict(glm1b, newdata = newdat1, type = "response")            
newdat1 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_abu_target.csv"), row.names = F)

# rest.age *
newdat2 <- expand_grid(rest.age = seq(min(data$rest.age),
                                         max(data$rest.age),
                                         length = 100),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = mean(data$plant.target.richness),
                       flower.cover = mean(data$flower.cover),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data$region),
                       obs.year = "2022") 
newdat2$butterfly.target.abu <- predict(glm1b, newdata = newdat2, type = "response")            
newdat2 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var2_abu_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C MODELLING LANDSCAPE VARIABLES (300 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm2 ####
glm2 <- glm(butterfly.target.abu ~ region + site.cwm.pres.oek.f
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
            + field.size.avg.r300,
            data = data , family = "poisson")
summary(glm2)


### a) check multicollinearity -------------------------------------------------
car::vif(glm2)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm2 <- simulateResiduals(glm2)
testDispersion(sr_glm2)
# switch to nb



## 2 glm2a ####
glm2a <- glm.nb(butterfly.target.abu ~ region + site.cwm.pres.oek.f
                + lk.grassland.per.min
                + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
                + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
                + field.size.avg.r300,
                data = data)
glm2a_summary <- summary(glm2a)
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "summary_final_model_landscape_r300.txt"))
glm2a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm2a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm2a <- simulateResiduals(glm2a)
testDispersion(sr_glm2a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm2a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm2a
# observed = -0.0020346, expected = -0.0090090, sd = 0.0301578, p-value = 0.8171
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_target", 
         "validation_final_model_landscape_r300.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm2a)) # ok
plotQQunif(sr_glm2a)      # ok
plotResiduals(sr_glm2a)   # ok
dev.off()


### e) R2  ---------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm2a)
# R2m       R2c
# delta     0.3890790 0.3890790
# lognormal 0.4637931 0.4637931
# trigamma  0.2992167 0.2992167


### f) test fixed effects  -----------------------------------------------------
glm2a_anova <- car::Anova(glm2a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "anova_final_model_landscape_r300.txt"))
glm2a_anova
sink()


### g) export for plotting  ----------------------------------------------------
# perm.grassland.r300 (.)
newdat3 <- expand_grid(perm.grassland.r300 = seq(min(data$perm.grassland.r300),
                                      max(data$perm.grassland.r300),
                                      length = 100),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       rest.age = mean(data$rest.age),
                       lc.shannon.div.r300 = mean(data$lc.shannon.div.r300),
                       field.size.avg.r300 = mean(data$field.size.avg.r300),
                       lk.grassland.per.min  = mean(data$lk.grassland.per.min),
                       region = unique(data$region))
newdat3$butterfly.target.abu <- predict(glm2a, newdata = newdat3, type = "response")
newdat3 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr300_abu_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D MODELLING LANDSCAPE VARIABLES (600 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm3 ####
glm3 <- glm(butterfly.target.abu ~ region + + site.cwm.pres.oek.f
            + lk.grassland.per.min 
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
            + field.size.avg.r600,
            data = data , family = "poisson")
summary(glm3)


### a) check multicollinearity -------------------------------------------------
car::vif(glm3)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm3 <- simulateResiduals(glm3)
testDispersion(sr_glm3)
# switch to nb



## 2 glm3a ####
glm3a <- glm.nb(butterfly.target.abu ~ region + site.cwm.pres.oek.f
                + lk.grassland.per.min
                + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
                + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
                + field.size.avg.r600,
            data = data)
glm3a_summary <- summary(glm3a)
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "summary_final_model_landscape_r600.txt"))
glm3a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm3a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm3a <- simulateResiduals(glm3a)
testDispersion(sr_glm3a)
# ok


### c) check spatial auto-correlation - ----------------------------------------
testSpatialAutocorrelation(glm3a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm3a
# observed = 0.028015, expected = -0.009009, sd = 0.030156, p-value = 0.2195
# alternative hypothesis: Distance-based autocorrelation


### d) validate model  ---------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_target", 
         "validation_final_model_landscape_r600.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm3a)) # ok
plotQQunif(sr_glm3a)      # ok
plotResiduals(sr_glm3a)   # ok
dev.off()


### e) R2  ---------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm3a)
# R2m       R2c
# delta     0.4006799 0.4006799
# lognormal 0.4769874 0.4769874
# trigamma  0.3080040 0.3080040  


### f) test fixed effects ------------------------------------------------------
glm3a_anova <- car::Anova(glm3a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "anova_final_model_landscape_r600.txt"))
glm3a_anova
sink()


### g) export for plotting -----------------------------------------------------
# effect of interaction
lc_mean  <- mean(data$lc.shannon.div.r600)
lc_sd    <- sd(data$lc.shannon.div.r600)
newdat_int <- expand_grid(
  rest.age = seq(min(data$rest.age),
                 max(data$rest.age),
                 length = 100),
  lc.shannon.div.r600 = c(lc_mean - lc_sd,
                          lc_mean,
                          lc_mean + lc_sd),
  site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
  region = unique(data$region),
  perm.grassland.r600 = mean(data$perm.grassland.r600),
  lk.grassland.per.min = mean(data$lk.grassland.per.min),
  field.size.avg.r600 = mean(data$field.size.avg.r600))
newdat_int$butterfly.target.abu <- predict(glm3a, newdata = newdat_int, type = "response")
newdat_int$lc_level <- factor(newdat_int$lc.shannon.div.r600,
                              labels = c("Low",
                                         "Mean",
                                         "High Landscape Shannon (600 m)"))
newdat_int %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr600_abu_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E MODELLING LANDSCAPE VARIABLES (1200 m) #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm4 ####
glm4 <- glm(butterfly.target.abu ~ region + site.cwm.pres.oek.f
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + field.size.avg.r1200,
            data = data , family = "poisson")
summary(glm4)


### a) check multicollinearity -------------------------------------------------
car::vif(glm4)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm4 <- simulateResiduals(glm4)
testDispersion(sr_glm4)
# switch to nb



## 2 glm4a ####
glm4a <- glm.nb(butterfly.target.abu ~ region + site.cwm.pres.oek.f
                + lk.grassland.per.min
                + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
                + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
                + field.size.avg.r1200,
                data = data)
glm4a_summary <- summary(glm4a)
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "summary_final_model_landscape_r1200.txt"))
glm4a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm4a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm4a <- simulateResiduals(glm4a)
testDispersion(sr_glm4a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm4a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm4a
# observed = 0.013779, expected = -0.009009, sd = 0.030166, p-value = 0.45
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_target", 
         "validation_final_model_landscape_r1200.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm4a)) # ok
plotQQunif(sr_glm4a)      # ok
plotResiduals(sr_glm4a)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm4a)
# R2m       R2c
# delta     0.4059507 0.4059507
# lognormal 0.4827090 0.4827090
# trigamma  0.3123434 0.3123434 


### f) test fixed effects ------------------------------------------------------
glm4a_anova <- car::Anova(glm4a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "anova_final_model_landscape_r1200.txt"))
glm4a_anova
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F MODELLING SITE VARIABLES (excluding cultivar seed mixtures) ################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm_no_cus ####
glm_no_cus <- glm(butterfly.target.abu ~ region + obs.year + region:obs.year
                + rest.age  + site.cwm.pres.oek.f
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data_x_cus , family = "poisson")
summary(glm_no_cus)


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus)
# exclude region:obs.year interaction



## 2 glm_no_cus_a ####
glm_no_cus_a <- update(glm_no_cus, .~. -region:obs.year)
summary(glm_no_cus_a)


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus_a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm_no_cus_a <- simulateResiduals(glm_no_cus_a)
testDispersion(sr_glm_no_cus_a)
# switch to nb



## 3 glm_no_cus_b ####
glm_no_cus_b <- glm.nb(butterfly.target.abu ~ region + obs.year
                       + rest.age  + site.cwm.pres.oek.f
                       + plant.target.richness + flower.cover
                       + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                    data = data_x_cus)
glm_no_cus_b_summary <- summary(glm_no_cus_b)
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "summary_final_model_no_cus.txt"))
glm_no_cus_b_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus_b)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm_no_cus_b <- simulateResiduals(glm_no_cus_b)
testDispersion(sr_glm_no_cus_b)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm_no_cus_b, x = data_x_cus$x, y = data_x_cus$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm_no_cus_b
# observed = 0.0049454, expected = -0.0108696, sd = 0.0354330, p-value = 0.6554
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_target", 
         "validation_final_model_no_cus.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm_no_cus_b)) # ok
plotQQunif(sr_glm_no_cus_b)      # ok
plotResiduals(sr_glm_no_cus_b)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm_no_cus_b)
# R2m       R2c
# delta     0.4369737 0.4369737
# lognormal 0.5017119 0.5017119
# trigamma  0.3583662 0.3583662 


### f) test fixed effects ------------------------------------------------------
glm_no_cus_b_anova <- car::Anova(glm_no_cus_b, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_target", 
          "anova_final_model_no_cus.txt"))
glm_no_cus_b_anova
sink()




# END ####