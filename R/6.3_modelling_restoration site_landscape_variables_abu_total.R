#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Influence of site and landscape variables - total abundance
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
glm1 <- glm(butterfly.abu ~ region + obs.year + region:obs.year
            + rest.age + site.cwm.pres.oek.f
            + plant.target.richness + flower.cover
            + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
            data = data , family = "poisson")
summary(glm1)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1)
# exclude region:obs.year interaction



## 2 glm1a ####
glm1a <- update(glm1, .~. - obs.year:region) 
summary(glm1a)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm1a <- simulateResiduals(glm1a)
testDispersion(sr_glm1a)
# switch to nb



## 3 glm1b ####
glm1b <- glm.nb(butterfly.abu ~ region + obs.year
                + rest.age + site.cwm.pres.oek.f
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data)
glm1b_summary <- summary(glm1b)
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "summary_final_model_site.txt"))
glm1b_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm1b)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm1b <- simulateResiduals(glm1b)
testDispersion(sr_glm1b) 


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm1b, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm1b
# observed = -0.048009, expected = -0.009009, sd = 0.030182, p-value = 0.1963
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_tot", 
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
# delta     0.3887787 0.3887787
# lognormal 0.4333056 0.4333056
# trigamma  0.3382060 0.3382060


### f) test fixed effects ------------------------------------------------------
anova_site <- car::Anova(glm1b, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "anova_final_model_site.txt"))
anova_site
sink()


### g) export for plotting
# Char. target richness
newdat1 <- expand_grid(rest.age = mean(data$rest.age),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = seq(min(data$plant.target.richness),
                                                   max(data$plant.target.richness),
                                                   length = 100),
                       flower.cover = mean(data$flower.cover),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data$region),
                       obs.year = "2022") 
newdat1$butterfly.abu <- predict(glm1b, newdata = newdat1, type = "response")            
newdat1 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_abu_tot.csv"), row.names = F)

# Ellenberg moisture
newdat2 <- expand_grid(rest.age = mean(data$rest.age),
                       site.cwm.pres.oek.f = seq(min(data$site.cwm.pres.oek.f),
                                                 max(data$site.cwm.pres.oek.f),
                                                 length = 100),
                       plant.target.richness = mean(data$plant.target.richness),
                       flower.cover = mean(data$flower.cover),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data$region),
                       obs.year = "2022") 
newdat2$butterfly.abu <- predict(glm1b, newdata = newdat2, type = "response")            
newdat2 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var2_abu_tot.csv"), row.names = F)

# Vegetation heterogeneity
newdat3 <- expand_grid(rest.age = mean(data$rest.age),
                       veg.heterogeneity = seq(min(data$veg.heterogeneity),
                                                 max(data$veg.heterogeneity),
                                                 length = 100),
                       plant.target.richness = mean(data$plant.target.richness),
                       flower.cover = mean(data$flower.cover),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       mang.app.NEU.MM.1.minus  = mean(data$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data$region),
                       obs.year = "2022") 
newdat3$butterfly.abu <- predict(glm1b, newdata = newdat3, type = "response")            
newdat3 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var3_abu_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C MODELLING LANDSCAPE VARIABLES (300 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm2 ####
glm2 <- glm(butterfly.abu ~ region + plant.target.richness + site.cwm.pres.oek.f 
            + veg.heterogeneity 
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
glm2a <- glm.nb(butterfly.abu ~ region + plant.target.richness + site.cwm.pres.oek.f 
             + veg.heterogeneity
             + lk.grassland.per.min
             + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
             + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
             + field.size.avg.r300,
             data = data)
glm2a_summary <- summary(glm2a)
sink(here("outputs", "models", "restored_sites", "abu_tot", 
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
# observed = -0.041316, expected = -0.009009, sd = 0.030165, p-value = 0.2842
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_tot", 
         "validation_final_model_landscape_r300.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm2a)) # ok
plotQQunif(sr_glm2a)      # ok
plotResiduals(sr_glm2a)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm2a)
# R2m       R2c
# delta     0.4057178 0.4057178
# lognormal 0.4498252 0.4498252
# trigamma  0.3553971 0.3553971     


### f) test fixed effects ------------------------------------------------------
glm2a_anova <- car::Anova(glm2a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "anova_final_model_landscape_r300.txt"))
glm2a_anova
sink()


### g) export for plotting -----------------------------------------------------
# perm.grassland.r300 (.)
newdat4 <- expand_grid(rest.age = mean(data$rest.age),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = mean(data$plant.target.richness),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       perm.grassland.r300 = seq(min(data$perm.grassland.r300),
                                                   max(data$perm.grassland.r300),
                                                   length = 100),
                       lk.grassland.per.min  = mean(data$lk.grassland.per.min),
                       lc.shannon.div.r300  = mean(data$lc.shannon.div.r300),
                       field.size.avg.r300  = mean(data$field.size.avg.r300),
                       region = unique(data$region))
newdat4$butterfly.abu <- predict(glm2a, newdata = newdat4, type = "response")
newdat4 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr300_abu_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D MODELLING LANDSCAPE VARIABLES (600 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm3 ####
glm3 <- glm(butterfly.abu ~ region + plant.target.richness + site.cwm.pres.oek.f 
            + veg.heterogeneity
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
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
glm3a <- glm.nb(butterfly.abu ~ region + plant.target.richness + site.cwm.pres.oek.f 
                + veg.heterogeneity
                + lk.grassland.per.min
                + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
                + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
                + field.size.avg.r600,
                data = data)
glm3a_summary <- summary(glm3a)
sink(here("outputs", "models", "restored_sites", "abu_tot", 
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


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm3a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm3a
# observed = -0.029995, expected = -0.009009, sd = 0.030178, p-value = 0.4868
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_tot", 
         "validation_final_model_landscape_r600.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm3a)) # ok
plotQQunif(sr_glm3a)      # ok
plotResiduals(sr_glm3a)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm3a)
# R2m       R2c
# delta     0.4175940 0.4175940
# lognormal 0.4616767 0.4616767
# trigamma  0.3671079 0.3671079


### f) test fixed effects ------------------------------------------------------
glm3a_anova <- car::Anova(glm3a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "anova_final_model_landscape_r600.txt"))
glm3a_anova
sink()


### g) export for plotting
# perm.grassland.r600 (.)
newdat5 <- expand_grid(rest.age = mean(data$rest.age),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = mean(data$plant.target.richness),
                       region = unique(data$region),
                       perm.grassland.r600 = seq(min(data$perm.grassland.r600),
                                                 max(data$perm.grassland.r600),
                                                 length = 100),
                       lk.grassland.per.min = mean(data$lk.grassland.per.min),
                       lc.shannon.div.r600 = mean(data$lc.shannon.div.r600),
                       field.size.avg.r600 = mean(data$field.size.avg.r600))
newdat5$butterfly.abu <- predict(glm3a, newdata = newdat5, type = "response")
newdat5 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr600_abu_tot.csv"), row.names = F)

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
  veg.heterogeneity = mean(data$veg.heterogeneity),
  site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
  region = unique(data$region),
  perm.grassland.r600 = mean(data$perm.grassland.r600),
  lk.grassland.per.min = mean(data$lk.grassland.per.min),
  field.size.avg.r600 = mean(data$field.size.avg.r600))
newdat_int$butterfly.abu <- predict(glm3a, newdata = newdat_int, type = "response")
newdat_int$lc_level <- factor(newdat_int$lc.shannon.div.r600,
                              labels = c("Low",
                                         "Mean",
                                         "High Landscape Shannon (600 m)"))
newdat_int %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr600_Int_abu_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E MODELLING LANDSCAPE VARIABLES (1200 m) #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm4 ####
glm4 <- glm(butterfly.abu ~ region + plant.target.richness + site.cwm.pres.oek.f 
            + veg.heterogeneity
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
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
glm4a <- glm.nb(butterfly.abu ~ region + plant.target.richness + site.cwm.pres.oek.f 
                + veg.heterogeneity
                + lk.grassland.per.min
                + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
                + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
                + field.size.avg.r1200,
                data = data)
glm4a_summary <- summary(glm4a)
sink(here("outputs", "models", "restored_sites", "abu_tot", 
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
# observed = -0.039231, expected = -0.009009, sd = 0.030177, p-value = 0.3166
# alternative hypothesis: Distance-based autocorrelatio


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_tot", 
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
# delta     0.4293159 0.4293159
# lognormal 0.4732986 0.4732986
# trigamma  0.3787596 0.3787596   


### f) test fixed effects ------------------------------------------------------
glm4a_anvova <- car::Anova(glm4a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "anova_final_model_landscape_r1200.txt"))
glm4a_anvova
sink()


### g) export for plotting -----------------------------------------------------
# perm.grassland.r1200 *
newdat6 <- expand_grid(rest.age = mean(data$rest.age),
                       veg.heterogeneity = mean(data$veg.heterogeneity),
                       site.cwm.pres.oek.f = mean(data$site.cwm.pres.oek.f),
                       plant.target.richness = mean(data$plant.target.richness),
                       region = unique(data$region),
                       perm.grassland.r1200 = seq(min(data$perm.grassland.r1200),
                                                   max(data$perm.grassland.r1200),
                                                   length = 100),
                       agrar.r1200 = mean(data$agrar.r1200),
                       lk.grassland.per.min = mean(data$lk.grassland.per.min),
                       lc.shannon.div.r1200 = mean(data$lc.shannon.div.r1200),
                       field.size.avg.r1200 = mean(data$field.size.avg.r1200)) 
newdat6$butterfly.abu <- predict(glm4a, newdata = newdat6, type = "response")            
newdat6 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr1200_abu_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F MODELLING SITE VARIABLES (excluding cultivar seed mixtures) ################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm_no_cus ####
glm_no_cus <- glm(butterfly.abu ~ region + obs.year + region:obs.year
                + rest.age + site.cwm.pres.oek.f
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data_x_cus , family = "poisson")
summary(glm_no_cus)


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus)
# exclude region:obs.year interaction



## 2 glm_no_cus_a ####
glm_no_cus_a <- update(glm_no_cus, .~.-region:obs.year)
summary(glm_no_cus_a)


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus_a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm_no_cus_a <- simulateResiduals(glm_no_cus_a)
testDispersion(sr_glm_no_cus_a)
# switch to nb



## 3 glm_no_cus_b ####
glm_no_cus_b <- glm.nb(butterfly.abu ~ region + obs.year
                       + rest.age + site.cwm.pres.oek.f
                       + plant.target.richness + flower.cover
                       + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                    data = data_x_cus)
glm_no_cus_summary <- summary(glm_no_cus_b)
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "summary_final_model_no_cus.txt"))
glm_no_cus_summary
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
# observed = -0.015075, expected = -0.010870, sd = 0.035466, p-value = 0.9056
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "abu_tot", 
         "validation_final_model_no_cus.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm_no_cus_b)) # ok
plotQQunif(sr_glm_no_cus_b)      # ok
plotResiduals(sr_glm_no_cus_b)   # ok
dev.off()


### e) R2  ---------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm_no_cus_b)
# R2m       R2c
# delta     0.3783558 0.3783558
# lognormal 0.4171556 0.4171556
# trigamma  0.3350141 0.3350141


### f) test fixed effects  -----------------------------------------------------
glm_no_cus_anova <- car::Anova(glm_no_cus_b, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "abu_tot", 
          "anova_final_model_no_cus.txt"))
glm_no_cus_anova
sink()




# END ####