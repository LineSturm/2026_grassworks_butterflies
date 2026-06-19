#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Analyses of site type - total species richness
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
glm1 <- glm(butterfly.rich.total ~ site.type + obs.year + region + site.cwm.pres.oek.f
            + region:obs.year,
            data = data , family = "poisson")
summary(glm1)


### a) check multicollinearity -------------------------------------------------
# since categorial variables are included the generalized vif should be used
# values > √5 (2.2) indicate problematic multicollinearity
car::vif(glm1) 
# ok


### b) check over/underdisperion -----------------------------------------------
sr_glm1 <- simulateResiduals(glm1)
testDispersion(sr_glm1)
# switch to negative binomial



## 2 glm1_nb ####
glm1_nb <- glm.nb(butterfly.rich.total ~ site.type + obs.year + region 
                  + site.cwm.pres.oek.f
                  + obs.year:region,
                  data = data)
summary(glm1_nb)


### a) check multicollinearity -------------------------------------------------
car::vif(glm1_nb) 
# ok


### b) check over/underdisperion -----------------------------------------------
sr_glm1_nb <- simulateResiduals(glm1_nb)
testDispersion(sr_glm1_nb) 
# ok


### c) check spatial auto-correlation -------------------------------------------
testSpatialAutocorrelation(glm1_nb, x = data$x, y = data$y)
# data:  glm1_nb
# observed = 0.1122388, expected = -0.0053763, sd = 0.0283903, p-value = 3.431e-05
# alternative hypothesis: Distance-based autocorrelatio
# --> spatial auto-correlation



## 3 test glmm ####
# glmm.nat.reg <- glmmTMB(butterfly.rich.total ~ site.type + obs.year + region 
#                         + site.cwm.pres.oek.f
#                         #+ obs.year:region # exclude due to moderate collinearity
#                         + (1|nat.reg.NR), # we tested different levels natural regions of Germany
#                         data = data, family = "poisson")
# check_collinearity(glmm.nat.reg)
# sr_glmm.nat.reg <- simulateResiduals(glmm.nat.reg)
# testDispersion(sr_glmm.nat.reg)
# testSpatialAutocorrelation(sr_glmm.nat.reg, x = data$x, y = data$y)
# nat.reg.NR:  p-value = 3.616e-08
# nat.reg.GL:  p-value = 8.79e-11
# nat.reg.OD2: p-value = 1.347e-06
# nat.reg.OD3: p-value = 4.437e-09
# --> still spatial auto-correlation
# --> therefore, a generalised additive model was used



## 4 gam1 ####
m_gam <- gam(
  butterfly.rich.total ~ site.type + obs.year + region + region:obs.year
    + site.cwm.pres.oek.f
    + s(x, y), 
  family = poisson(), 
  data = data, method = "REML")
summary.gam(m_gam)
plot.gam(m_gam)


### a) check multicollinearity --------------------------------------------------
check_collinearity(m_gam) 
check_concurvity(m_gam) 
concurvity(m_gam)
concurvity(m_gam, full = F)
# exclude region:obs.year interaction



## 5 gam1a ####
gam1a <- gam(
  butterfly.rich.total ~ site.type + obs.year + region
  + site.cwm.pres.oek.f
  + s(x, y), 
  family = poisson(), 
  data = data, method = "REML")
gam1a_summary <- summary.gam(gam1a)
sink(here("outputs", "models", "all_sites", "spe_rich_tot",
          "summary_final_model_site_type.txt"))
gam1a_summary
sink()
plot.gam(gam1a)


### a) check multicollinearity -------------------------------------------------
check_collinearity(gam1a) 
check_concurvity(gam1a) 
concurvity(gam1a)
concurvity(gam1a, full = F)
# ok


### b) check over/underdispersion -----------------------------------------------
sr_gam1a <- simulateResiduals(gam1a) 
testDispersion(sr_gam1a) 
# ok


### c) check spatial auto-correlation ------------------------------------------
testSpatialAutocorrelation(gam1a, x = data$x, y = data$y)
# DHARMa Moran's I test for distance-based autocorrelation
# data:  gam1a
# observed = 0.0346147, expected = -0.0053763, sd = 0.0284020, p-value = 0.1591
# alternative hypothesis: Distance-based autocorrelation


### d) validate model ----------------------------------------------------------
gam.check(gam1a) 
png(here("outputs", "models", "all_sites", "spe_rich_tot",
         "validation_final_model.png"),
    width = 12, height = 12 * 1.618, units = "cm", res = 300)
par(mfrow = c(3, 1))
hist(residuals(gam1a)) # ok
plotQQunif(gam1a)      # ok
plotResiduals(gam1a)   # ok
dev.off()


### e) test fixed effects ------------------------------------------------------
gam1a_anova <- anova.gam(gam1a)
sink(here("outputs", "models", "all_sites", "spe_rich_tot",
          "anova_final_model_site_type.txt"))
gam1a_anova
sink()


### f) pairwise test -----------------------------------------------------------
# Effect of site.type over all regions and year
EMM1 <- emmeans(gam1a, ~ site.type, type = "response")
EMM1
pairs(EMM1)
EMM1 <- EMM1 %>% as.data.frame
dat_pred2 <- data.frame(
  site.type = EMM1$site.type,
  predicted = EMM1$rate,
  std.error = EMM1$SE,
  conf.low = EMM1$lower.CL,
  conf.high = EMM1$upper.CL)
# export for plotting
dat_pred2 %>% write.csv(file = here::here("outputs", "tables", 
                                          "EMM1_site_type_spe_rich_tot.csv"), row.names = F)




# END ####