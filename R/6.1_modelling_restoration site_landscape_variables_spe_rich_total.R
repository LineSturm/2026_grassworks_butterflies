#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Influence of site and landscape variables - total species richness
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
glm1 <- glm(butterfly.rich.total ~ region + obs.year + region:obs.year
            + rest.age + site.cwm.pres.oek.f 
            + plant.target.richness + flower.cover
            + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
            data = data, family = "poisson")
summary(glm1)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm1 <- simulateResiduals(glm1)
testDispersion(sr_glm1)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(sr_glm1, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm1
# observed = 0.059767, expected = -0.009009, sd = 0.030151, p-value = 0.02254
# alternative hypothesis: Distance-based autocorrelation



## 2 gam1 #####
gam1 <- gam(butterfly.rich.total ~ region + obs.year + region:obs.year
             + rest.age + site.cwm.pres.oek.f
             + plant.target.richness + flower.cover
             + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type
             + s(x, y),
             family = poisson(), data = data, method = "REML")
summary.gam(gam1)
gam1_summary <- summary.gam(gam1)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot",
          "summary_final_model_site.txt"))
gam1_summary
sink()
plot.gam(gam1)


### a) check multicollinearity --------------------------------------------------
check_collinearity(gam1)
# high correlation only for region, otherwise ok
# keep region, because otherwise predicted values do not fit observed data
concurvity(gam1, full = F)


### b) check over/underdispersion ----------------------------------------------
sr_gam1 <- simulateResiduals(gam1)
testDispersion(sr_gam1)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(gam1, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  gam1
# observed = -0.021166, expected = -0.009009, sd = 0.030151, p-value = 0.6868
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_tot",
         "validation_final_model_site.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(gam1)) # ok
plotQQunif(gam1)      # ok
plotResiduals(gam1)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
gam1_anova <- anova.gam(gam1)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot",
          "anova_final_model_site.txt"))
gam1_anova
sink()


### f) export for plotting -----------------------------------------------------
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
                       obs.year = unique(data$obs.year),
                       x = mean(data$x),
                       y = mean(data$y))
newdat1$butterfly.rich.total <- predict(gam1, newdata = newdat1, type = "response", exclude = "s(x,y)")          
newdat1 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_site_var1_spe_rich_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C MODELLING LANDSCAPE VARIABLES (300 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm2 ####
glm2 <- glm(butterfly.rich.total ~ region + obs.year + region:obs.year
            + plant.target.richness 
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r300, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r300, scale = FALSE)
            + field.size.avg.r300,
            data = data , family = "poisson")
glm2_summary <- summary(glm2)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot", 
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
# observed = 0.044566, expected = -0.009009, sd = 0.030156, p-value = 0.07564
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_tot", 
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
# delta     0.4963210 0.4963210
# lognormal 0.5046835 0.5046835
# trigamma  0.4876781 0.4876781    


### f) test fixed effects ------------------------------------------------------
anova_r300 <- car::Anova(glm2, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_tot", 
          "anova_final_model_landscape_r300.txt"))
anova_r300
sink()


### g) export for plotting ----
newdat2 <- expand_grid(rest.age = mean(data$rest.age),
                       plant.target.richness = mean(data$plant.target.richness),
                       region = unique(data$region),
                       obs.year = unique(data$obs.year),
                       lk.grassland.per.min = mean(data$lk.grassland.per.min),
                       lc.shannon.div.r300 = mean(data$lc.shannon.div.r300),
                       field.size.avg.r300  = mean(data$field.size.avg.r300),
                       perm.grassland.r300 = seq(min(data$perm.grassland.r300),
                                                 max(data$perm.grassland.r300),
                                                 length = 100))
newdat2$butterfly.rich.total <- predict(glm2, newdata = newdat2, type = "response")
newdat2 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_Lr300_spe_rich_tot.csv"), row.names = F)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D MODELLING LANDSCAPE VARIABLES (600 m) ######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm3 ####
glm3 <- glm(butterfly.rich.total ~ region + obs.year + region:obs.year
            + plant.target.richness
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
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm3, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm3
# observed = 0.053960, expected = -0.009009, sd = 0.030153, p-value = 0.03677
# alternative hypothesis: Distance-based autocorrelation



## 2 gam3 ####
gam3 <- gam(butterfly.rich.total ~ region + obs.year + region:obs.year
            + plant.target.richness + lk.grassland.per.min
            + field.size.avg.r600
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
            + s(x, y),
            family = poisson(), data = data, method = "REML")
summary.gam(gam3)
plot.gam(gam3)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam3)
concurvity(gam3, full = F)
# exclude region:obs.year interaction



## 2 gam3a ####
gam3a <- gam(butterfly.rich.total ~ region + obs.year 
            + plant.target.richness + lk.grassland.per.min
            + field.size.avg.r600
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r600, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r600, scale = FALSE)
            + s(x, y),
            family = poisson(), data = data, method = "REML")
gam3a_summary <- summary.gam(gam3a)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot",
          "summary_final_model_landscape_r600.txt"))
gam3a_summary
sink()
plot.gam(gam3a)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam3a)
# high correlation only for region
# keep region, because otherwise predicted values do not fit observed data
concurvity(gam3a, full = F)


### b) check over/underdispersion -----------------------------------------------
sr_gam3a <- simulateResiduals(gam3a)
testDispersion(sr_gam3a)  
# tendency to underdispersion


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(gam3a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  gam3a
# observed = -0.026103, expected = -0.009009, sd = 0.030164, p-value = 0.5709
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_tot",
         "validation_final_model_landscape_r600.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_gam3a)) # ok
plotQQunif(sr_gam3a)      # ok
plotResiduals(sr_gam3a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
anova_r600 <- anova.gam(gam3a)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot",
          "anova_final_model_landscape_r600.txt"))
anova_r600
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E MODELLING LANDSCAPE VARIABLES (1200 m) #####################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm4 ####
glm4 <- glm(butterfly.rich.total ~ region + obs.year + region:obs.year
            + plant.target.richness
            + lk.grassland.per.min
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
            + field.size.avg.r1200,
            data = data , family = "poisson")
summary(glm4)


### a) check multicollinearity -------------------------------------------------
car::vif(glm4)
# ok


### b) check under/overdispersion -----------------------------------------------
sr_glm4 <- simulateResiduals(glm4)
testDispersion(sr_glm4)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm4, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm4
# observed = 0.068884, expected = -0.009009, sd = 0.030159, p-value = 0.009801
# alternative hypothesis: Distance-based autocorrelation



## 2 gam4 ####
gam4 <- gam(butterfly.rich.total ~ region + obs.year + region:obs.year
            + plant.target.richness + lk.grassland.per.min
            + field.size.avg.r1200
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + s(x, y), 
            family = poisson(), data = data, method = "REML")
summary.gam(gam4)
plot.gam(gam4)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam4) 
concurvity(gam4, full = F)
# exclude region:obs.year interaction



## 3 gam4a ####
gam4a <- gam(butterfly.rich.total ~ region + obs.year
            + plant.target.richness + lk.grassland.per.min
            + field.size.avg.r1200
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + s(x, y), 
            family = poisson(), data = data, method = "REML")
summary.gam(gam4a)
plot.gam(gam4a)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam4a) 
# high correlation only for region
# keep region, because otherwise predicted values do not fit observed data
concurvity(gam4a, full = F)
# ok


### b) check over/underdispersion ----------------------------------------------
sr_gam4a <- simulateResiduals(gam4a) 
testDispersion(sr_gam4a) 
# underdispersion



## 4 gam4b ####
gam4b <- gam(butterfly.rich.total ~ region + obs.year
            + plant.target.richness + lk.grassland.per.min
            + field.size.avg.r1200
            + scale(rest.age, scale = FALSE) * scale(perm.grassland.r1200, scale = FALSE)
            + scale(rest.age, scale = FALSE) * scale(lc.shannon.div.r1200, scale = FALSE)
            + s(x, y), 
            family = tw(), data = data, method = "REML")
summary_r1200 <- summary.gam(gam4b)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot", 
          "summary_final_model_landscape_r1200.txt"))
summary_r1200
sink()
plot.gam(gam4b)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam4b) 
# high correlation only for region, otherwise ok
# keep region, because otherwise predicted values do not fit observed data
concurvity(gam4b, full = F)


### b) check over/underdispersion ----------------------------------------------
sr_gam4b <- simulateResiduals(gam4b) 
testDispersion(sr_gam4b)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(gam4b, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  gam4b
# observed = -0.029280, expected = -0.009009, sd = 0.030174, p-value = 0.5017
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_tot", 
         "validation_final_model_landscape_r1200.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_gam4b)) # ok
plotQQunif(sr_gam4b)      # ok
plotResiduals(sr_gam4b)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
anova_r1200 <- anova.gam(gam4b)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot", 
          "anova_final_model_landscape_r1200.txt"))
anova_r1200
sink()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F MODELLING SITE VARIABLES (excluding cultivar seed mixtures) ################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm_no_cus #####
glm_no_cus <- glm(butterfly.rich.total ~ region + obs.year + region:obs.year
                + rest.age + site.cwm.pres.oek.f
                + plant.target.richness + flower.cover
                + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                data = data_x_cus , family = "poisson")
summary(glm_no_cus)



### a) check multicollinearity -------------------------------------------------
car::vif(glm_no_cus)
# exclude region:obs.year interaction



## 2 glm_no_cus_a ####
glm_no_cus_a <- glm(butterfly.rich.total ~ region + obs.year 
                  + rest.age + site.cwm.pres.oek.f 
                  + plant.target.richness + flower.cover
                  + veg.heterogeneity + mang.app.NEU.MM.1.minus + mngm.type,
                  data = data_x_cus , family = "poisson")
summary_no_cus <- summary(glm_no_cus_a)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot", 
          "summary_final_model_no_cus.txt"))
summary_no_cus
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
# observed = 0.043864, expected = -0.010870, sd = 0.035470, p-value = 0.1228
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_tot", 
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
# delta     0.4728169 0.4728169
# lognormal 0.4807931 0.4807931
# trigamma  0.4645978 0.4645978


### f) test fixed effects ------------------------------------------------------
anova_no_cus <- car::Anova(glm_no_cus_a, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_tot", 
          "anova_final_model_no_cus.txt"))
anova_no_cus
sink()


### g) export for plotting -----------------------------------------------------
newdat3 <- expand_grid(plant.target.richness = mean(data_x_cus$plant.target.richness),
                       site.cwm.pres.oek.f = mean(data_x_cus$site.cwm.pres.oek.f),
                       rest.age = seq(min(data_x_cus$rest.age),
                                      max(data_x_cus$rest.age),
                                      length = 100),
                       flower.cover = mean(data_x_cus$flower.cover),
                       veg.heterogeneity = mean(data_x_cus$veg.heterogeneity),
                       mang.app.NEU.MM.1.minus  = mean(data_x_cus$mang.app.NEU.MM.1.minus),
                       mngm.type = "mowing",
                       region = unique(data_x_cus$region),
                       obs.year = "2022")
newdat3$butterfly.rich.total <- predict(glm_no_cus_a, newdata = newdat3, type = "response")
newdat3 %>% write.csv(file = here::here("outputs", "tables",
                                        "Pred_noCus_spe_rich_tot.csv"), row.names = F)




# END ####