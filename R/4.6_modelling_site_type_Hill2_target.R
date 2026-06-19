#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Analyses of site type - target Hill-Simpson
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

## 1 lm1 ####
lm1 <- lm(butterfly.Hill2.target ~ site.type + obs.year + region + site.cwm.pres.oek.f
          + region:obs.year,
          data = data)
summary(lm1)


### a) check multicollinearity -------------------------------------------------
car::vif(lm1)
# ok


### b) check over/underdisperion -----------------------------------------------
sr_lm1 <- simulateResiduals(lm1)
testDispersion(sr_lm1)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(lm1, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  lm1
# observed = 0.0523257, expected = -0.0053763, sd = 0.0283937, p-value = 0.04213
# alternative hypothesis: Distance-based autocorrelation



## 2 gam1 ####
gam1 <- gam(
  butterfly.Hill2.target ~ site.type + obs.year + region + region:obs.year
  + site.cwm.pres.oek.f 
  + s(x, y), 
  family = gaussian(), 
  data = data, method = "REML")
gam1_summary <- summary.gam(gam1)
sink(here("outputs", "models", "all_sites", "Hill-Simpson_target",
          "summary_final_model_site_type.txt"))
gam1_summary
sink()
plot.gam(gam1)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam1)
check_concurvity(gam1) 
concurvity(gam1)
concurvity(gam1, full = F)
# ok


### b) check over/underdisperion -----------------------------------------------
sr_gam1 <- simulateResiduals(gam1)
testDispersion(sr_gam1)
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(gam1, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  gam1
# observed = 0.0037004, expected = -0.0053763, sd = 0.0283948, p-value = 0.7492
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----
gam.check(gam1) 
png(here("outputs", "models", "all_sites", "Hill-Simpson_target",
         "validation_final_model.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(gam1)) # ok
plotQQunif(gam1)      # ok
plotResiduals(gam1)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
gam1_anova <- anova.gam(gam1)
sink(here("outputs", "models", "all_sites", "Hill-Simpson_target",
          "anova_final_model_site_type.txt"))
gam1_anova
sink()


### f) pairwise test -----------------------------------------------------------
EMM1 <- emmeans(gam1, ~ site.type, type = "response")
EMM1
pairs(EMM1)
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  site.type = EMM1$site.type,
  predicted = EMM1$emmean,
  std.error = EMM1$SE,
  conf.low = EMM1$lower.CL,
  conf.high = EMM1$upper.CL)
# export for plotting
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_site_type_Hill2_target.csv"), row.names = F)




# END ####