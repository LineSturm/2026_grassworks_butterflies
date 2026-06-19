#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Analyses of restoration method - target species richness
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  ) %>%
  filter(site.type == "restored")
str(data)

# here only blurred coordinates are included in the data, but exact coordinates 
# were used to test spatial auto-correlation.




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B MODELLING ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm.RM ####
glm.RM <- glm(butterfly.rich.target ~ region + obs.year
              + rest.meth + site.cwm.pres.oek.f 
              + region:obs.year,
              data = data , family = "poisson")
summary(glm.RM)


### a) check multicollinearity -------------------------------------------------
car::vif(glm.RM)
# exclude region:obs.year interaction



## 2 glm.RMa ####
glm.RMa <- glm(butterfly.rich.target ~ region + obs.year
              + rest.meth + site.cwm.pres.oek.f ,
              data = data, family = "poisson")
glm.RM_summary <- summary(glm.RMa)
sink(here("outputs", "models", "restored_sites", "spe_rich_target",
          "summary_final_model_RM.txt"))
glm.RM_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm.RMa)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm.RMa <- simulateResiduals(glm.RMa)
testDispersion(sr_glm.RMa)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm.RMa, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm.RMa
# observed = 0.0437407, expected = -0.0083333, sd = 0.0277318, p-value = 0.06041
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_target", 
         "validation_final_model_RM.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm.RMa)) # ok
plotQQunif(sr_glm.RMa)      # ok
plotResiduals(sr_glm.RMa)   # ok
dev.off()



### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm.RMa)
# R2m       R2c
# delta     0.4583448 0.4583448
# lognormal 0.4802163 0.4802163
# trigamma  0.4346941 0.4346941 


### f) test fixed effects ------------------------------------------------------
glm.RMa_anova <- car::Anova(glm.RMa, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_target",
          "anova_final_model_RM.txt"))
glm.RMa_anova
sink()


### g) pairwise test ----
EMM1 <- emmeans(glm.RMa, ~ rest.meth, type = "response")
EMM1
pairs(EMM1)
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  rest.meth = EMM1$rest.meth,
  predicted = EMM1$rate,
  std.error = EMM1$SE,
  conf.low = EMM1$asymp.LCL,
  conf.high = EMM1$asymp.UCL)
# export for plotting
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_RM_spe_rich_target.csv"), row.names = F)




# END ####