#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Influence of site and landscape variables - target species richness
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
glm1 <- glm(butterfly.rich.target ~ region + obs.year + region:obs.year
            + rest.age + site.cwm.pres.oek.f
            + plant.target.richness + flower.cover
            + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
            data = data , family = "poisson")
summary(glm1)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1)



## 2 glm1a ####
glm1a <- update(glm1, .~. - region:obs.year) 
glm1a_summary <- summary(glm1a)
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "summary_final_model_site.txt"))
glm1a_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm1a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm1a <- simulateResiduals(glm1a)
testDispersion(sr_glm1a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm1a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm1a
# observed = 0.020417, expected = -0.009009, sd = 0.030156, p-value = 0.3292
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_target", 
         "validation_final_model_site.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm1a)) # ok
plotQQunif(sr_glm1a)      # ok
plotResiduals(sr_glm1a)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm1a)
# R2m       R2c
# delta     0.4921283 0.4921283
# lognormal 0.5139242 0.5139242
# trigamma  0.4684267 0.4684267    


### f) test fixed effects ------------------------------------------------------
anova_site <- car::Anova(glm1a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "anova_final_model_site.txt"))
anova_site
sink()


### g) export for plotting -----------------------------------------------------

# Char. plant species richness
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
newdat1$butterfly.rich.target <- predict(glm1a, newdata = newdat1, type = "response")            
newdat1 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_spe_rich_target.csv"), row.names = F)

# Ellenberg moisture value
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
newdat2$butterfly.rich.target <- predict(glm1a, newdata = newdat2, type = "response")
newdat2 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var2_spe_rich_target.csv"), row.names = F)

# restoration age
newdat3 <- expand_grid(rest.age = seq(min(data$rest.age),
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
newdat3$butterfly.rich.target <- predict(glm1a, newdata = newdat3, type = "response")            
newdat3 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var3_spe_rich_target.csv"), row.names = F)


### h) pairwise test ----
# Effect of mngm.type over all regions and year
EMM1 <- emmeans(glm1a, ~ mngm.type, type = "response")
EMM1
pairs(EMM1)
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  mngm.type = EMM1$mngm.type,
  predicted = EMM1$rate,
  std.error = EMM1$SE,
  conf.low = EMM1$asymp.LCL,
  conf.high = EMM1$asymp.UCL)
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_MType_spe_rich_target.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C MODELLING LANDSCAPE VARIABLES (300 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm2 ####
glm2 <- glm(butterfly.rich.target ~ region + mngm.type
            + plant.target.richness + site.cwm.pres.oek.f
            + lk.grassland.per.min 
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
            + field.size.avg.r300,
            data = data , family = "poisson")
glm2_summary <- summary(glm2)
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "summary_final_model_landscape_r300.txt"))
glm2_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm2)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm2 <- simulateResiduals(glm2)
testDispersion(sr_glm2)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm2, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm2
# observed = -0.0011872, expected = -0.0090090, sd = 0.0301567, p-value = 0.7953
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_target", 
         "validation_final_model_landscape_r300.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm2)) # ok
plotQQunif(sr_glm2)      # ok
plotResiduals(sr_glm2)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm2)
# R2m       R2c
# delta     0.4899877 0.4899877
# lognormal 0.5117843 0.5117843
# trigamma  0.4662946 0.4662946      


### f) test fixed effects ------------------------------------------------------
anova_r300 <- car::Anova(glm2, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "anova_final_model_landscape_r300.txt"))
anova_r300
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D MODELLING LANDSCAPE VARIABLES (600 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm3 ####
glm3 <- glm(butterfly.rich.target ~ region + mngm.type
            + plant.target.richness + site.cwm.pres.oek.f
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
            + field.size.avg.r600,
            data = data , family = "poisson")
glm3_summary <- summary(glm3)
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "summary_final_model_landscape_r600.txt"))
glm3_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm3)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm3 <- simulateResiduals(glm3)
testDispersion(sr_glm3)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm3, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm3
# observed = -0.0055896, expected = -0.0090090, sd = 0.0301708, p-value = 0.9098
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_target", 
         "validation_final_model_landscape_r600.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm3)) # ok
plotQQunif(sr_glm3)      # ok
plotResiduals(sr_glm3)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm3)
# R2m       R2c
# delta     0.4912835 0.4912835
# lognormal 0.5130798 0.5130798
# trigamma  0.4675851 0.4675851  


### f) test fixed effects ------------------------------------------------------
anova_r600 <- car::Anova(glm3, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "anova_final_model_landscape_r600.txt"))
anova_r600
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E MODELLING LANDSCAPE VARIABLES (1200 m) #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm4 ####
glm4 <- glm(butterfly.rich.target ~ region + mngm.type
            + plant.target.richness + site.cwm.pres.oek.f
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
            + field.size.avg.r1200,
            data = data , family = "poisson")
glm4_summary <- summary(glm4)
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "summary_final_model_landscape_r1200.txt"))
glm4_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm4)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm4 <- simulateResiduals(glm4)
testDispersion(sr_glm4)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm4, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm4
# observed = 0.039523, expected = -0.009009, sd = 0.030161, p-value = 0.1076
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_target", 
         "validation_final_model_landscape_r1200.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm4)) # ok
plotQQunif(sr_glm4)      # ok
plotResiduals(sr_glm4)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm4)
# R2m       R2c
# delta     0.4896034 0.4896034
# lognormal 0.5114001 0.5114001
# trigamma  0.4659119 0.4659119  


### f) test fixed effects ------------------------------------------------------
anova_r1200 <- car::Anova(glm4, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "anova_final_model_landscape_r1200.txt"))
anova_r1200
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F MODELLING SITE VARIABLES (excluding cultivar seed mixtures) ################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm_no_cus ####
glm_no_cus <- glm(butterfly.rich.target ~ region + obs.year + region:obs.year
                + rest.age + site.cwm.pres.oek.f
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data_x_cus , family = "poisson")
summary(glm_no_cus)


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus)
# exclude region:obs.year interaction



## 2 glm_no_cus_a ####
glm_no_cus_a <- update(glm_no_cus, .~. -region:obs.year)
summary_glm_no_cus_a <- summary(glm_no_cus_a)
sink(here("outputs", "models", "restored_sites", "spe_rich_target", 
          "summary_final_model_no_cus.txt"))
summary_glm_no_cus_a
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus_a)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm_no_cus_a <- simulateResiduals(glm_no_cus_a)
testDispersion(sr_glm_no_cus_a)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm_no_cus_a, x = data_x_cus$x, y = data_x_cus$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm_no_cus_a
# observed = 0.045880, expected = -0.010870, sd = 0.035426, p-value = 0.1092
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_target",
         "validation_final_model_no_cus.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm_no_cus_a)) # ok
plotQQunif(sr_glm_no_cus_a)      # ok
plotResiduals(sr_glm_no_cus_a)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm_no_cus_a)
# R2m       R2c
# delta     0.4535483 0.4535483
# lognormal 0.4734604 0.4734604
# trigamma  0.4321764 0.4321764


### f) test fixed effects ------------------------------------------------------
anova_no_cus <- car::Anova(glm_no_cus_a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_target",
          "anova_final_model_no_cus.txt"))
anova_no_cus
sink()




# END ####