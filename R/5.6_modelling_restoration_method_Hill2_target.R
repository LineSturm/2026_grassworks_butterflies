#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Analyses of restoration method - target Hill-Simpson
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
  ) %>%
  filter(site.type == "restored")
str(data)

# here only blurred coordinates are included in the data, but exact coordinates 
# were used to test spatial auto-correlation.




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B MODELLING ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 lm.RM ####
lm.RM <- lm(butterfly.Hill2.target ~ region + obs.year
              + rest.meth + site.cwm.pres.oek.f 
              + region:obs.year,
              data = data)
lm.RM_summary <- summary(lm.RM)
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target",
          "summary_final_model_RM.txt"))
lm.RM_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(lm.RM)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_lm.RM <- simulateResiduals(lm.RM)
testDispersion(sr_lm.RM)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm.RM, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm.RM
# observed = -0.0033267, expected = -0.0083333, sd = 0.0277223, p-value = 0.8567
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "Hill-Simpson_target", 
         "validation_final_model_RM.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_lm.RM)) # ok
plotQQunif(sr_lm.RM)      # ok
plotResiduals(sr_lm.RM)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
lm_anova <- car::Anova(lm.RM, type = "II", test.statistic = "F")
sink(here("outputs", "models", "restored_sites", "Hill-Simpson_target",
          "anova_final_model_RM.txt"))
lm_anova
sink()


### f) pairwise test -----------------------------------------------------------
EMM1 <- emmeans(lm.RM, ~ rest.meth, type = "response")
EMM1
pairs(EMM1)
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  rest.meth = EMM1$rest.meth,
  predicted = EMM1$emmean,
  std.error = EMM1$SE,
  conf.low = EMM1$lower.CL,
  conf.high = EMM1$upper.CL)
# export for plotting
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_RM_Hill2_target.csv"), row.names = F)




# END ####