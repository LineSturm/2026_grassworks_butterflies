#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Analyses of site type - total abundance
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 13 June 2026



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
  )
str(data)

# here only blurred coordinates are included in the data, but exact coordinates 
# were used to test spatial auto-correlation.




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B MODELLING ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm1 ####
glm1 <- glm(butterfly.abu ~ site.type + obs.year + region + site.cwm.pres.oek.f
                + obs.year:region,
                data = data, family = "poisson")
summary(glm1)



### a) check multicollinearity -------------------------------------------------
car::vif(glm1)
# exclude region:obs.year interaction



## 2 glm1a ####
glm1a <- glm(butterfly.abu ~ site.type + obs.year + region + site.cwm.pres.oek.f,
             data = data, family = "poisson")
summary(glm1a)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1a)
# ok


### b) check over/underdisperion -----------------------------------------------
sr_glm1a <- simulateResiduals(glm1a)
testDispersion(sr_glm1a)
# --> switch to nb



## 3 glm1a_nb ####
glm1a_nb <- glm.nb(butterfly.abu ~ site.type + obs.year + region + site.cwm.pres.oek.f,
             data = data)
glm1a_nb_summary <- summary(glm1a_nb)
sink(here("outputs", "models", "all_sites", "abu_tot",
          "summary_final_model_site_type.txt"))
glm1a_nb_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm1a_nb)
# ok


### b) check over/underdisperion ----
sr_glm1a_nb <- simulateResiduals(glm1a_nb)
testDispersion(sr_glm1a_nb)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm1a_nb, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm1a_nb
# observed = 0.0440426, expected = -0.0053763, sd = 0.0284113, p-value = 0.08196
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "all_sites", "abu_tot",
         "validation_final_model.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm1a_nb)) # ok
plotQQunif(sr_glm1a_nb)      # ok
plotResiduals(sr_glm1a_nb)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm1a_nb)
# R2m       R2c
# delta     0.3363968 0.3363968
# lognormal 0.3843728 0.3843728
# trigamma  0.2824757 0.2824757


### f) test fixed effects ------------------------------------------------------
glm1a_nb_anova <- car::Anova(glm1a_nb, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "all_sites", "abu_tot",
          "anova_final_model_site_type.txt"))
glm1a_nb_anova
sink()


### g) pairwise test -----------------------------------------------------------
EMM1 <- emmeans(glm1a_nb, ~ site.type, type = "response")
EMM1
pairs(EMM1)
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  site.type = EMM1$site.type,
  predicted = EMM1$response,
  std.error = EMM1$SE,
  conf.low = EMM1$asymp.LCL,
  conf.high = EMM1$asymp.UCL)
# export for plotting
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_site_type_abu_tot.csv"), row.names = F)




# END ####