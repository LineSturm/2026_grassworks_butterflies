#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Analyses of restoration method - total species richness
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
  filter(site.type == "restored")
str(data)

# here only blurred coordinates are included in the data, but exact coordinates 
# were used to test spatial auto-correlation.




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B MODELLING ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 glm.RM ####
glm.RM <- glm(butterfly.rich.total ~ region + obs.year
              + rest.meth + site.cwm.pres.oek.f 
              + region:obs.year,
              data = data , family = "poisson")
glm.RM_summary <- summary(glm.RM)
sink(here("outputs", "models", "restored_sites", "spe_rich_tot",
          "summary_final_model_RM.txt"))
glm.RM_summary
sink()


### a) check multicollinearity -------------------------------------------------
car::vif(glm.RM)
# ok


### b) check under/overdispersion ----------------------------------------------
sr_glm.RM <- simulateResiduals(glm.RM)
testDispersion(sr_glm.RM)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(glm.RM, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  glm.RM
# observed = 0.0445930, expected = -0.0083333, sd = 0.0277411, p-value = 0.05641
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
png(here("outputs", "models", "restored_sites", "spe_rich_tot", 
         "validation_final_model_RM.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(sr_glm.RM)) # ok
plotQQunif(sr_glm.RM)      # ok
plotResiduals(sr_glm.RM)   # ok
dev.off()


### e) R2 ----------------------------------------------------------------------
MuMIn::r.squaredGLMM(glm.RM)
# R2m       R2c
# delta     0.4163267 0.4163267
# lognormal 0.4245771 0.4245771
# trigamma  0.4078432 0.4078432


### f) test fixed effects ------------------------------------------------------
glm.RM_anova <- car::Anova(glm.RM, type = "II", test.statistic = "LR")
sink(here("outputs", "models", "restored_sites", "spe_rich_tot",
          "anova_final_model_RM.txt"))
glm.RM_anova
sink()


### g) pairwise test -----------------------------------------------------------
EMM1 <- emmeans(glm.RM, ~ rest.meth, type = "response")
EMM1
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  rest.meth = EMM1$rest.meth,
  predicted = EMM1$rate,
  std.error = EMM1$SE,
  conf.low = EMM1$asymp.LCL,
  conf.high = EMM1$asymp.UCL)
# export for plotting
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_RM_spe_rich_tot.csv"), row.names = F)




# END ####