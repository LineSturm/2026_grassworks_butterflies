#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# cRF analyses - total species richness
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 13 June 2026




# packages ####

library(here)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(tibble)
library(cowplot)
library(magrittr)
library(randomForest)
library(party)
library(caret)

here()



# START ####

rm(list = ls())




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A LOAD DATA ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_restored <- read_csv(
  here("data", "processed", "data_all.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(.default = "?")) %>%
  filter(site.type == "restored",
         mang.app.NEU.MM.1.minus != "NA") %>%
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
    region = fct_relevel(region, "north", "centre", "south"),
    site.type = fct_relevel(site.type, "negative", "restored", "positive"),
    rest.meth = fct_relevel(rest.meth, "cus", "res", "dih", "mga"),
    rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "dih", "mga", "positive")
  )
str(data_restored)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B CONDITIONAL RANDOM FOREST ##################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

colSums(is.na(data_restored)) 
# Management appropriateness: 9 sites
# N_EIC, N_PEV, N_MAU, 
# M_BEN, M_GRO, M_KOL, 
# S_NOZ, S_SHH, S_TIE
# these sites are dropped!



## 1 define data ####

data_r1200 <- data_restored %>%
  dplyr::select(butterfly.rich.total,
                # restoration variables
                land.use.hist,
                rest.age,
                # site variables     
                region,
                obs.year,
                pH.site,
                site.cwm.pres.oek.f,
                # vegetation variables
                plant.target.richness, 
                cover.vegetation,
                cover.bare.soil,
                cover.litter,
                flower.cover,
                veg.heterogeneity,
                # management variables 
                mngm.type,
                mang.app.NEU.MM.1.minus,     
                # landscape variables
                lk.grassland.per.min,
                field.size.avg.r1200,
                lc.shannon.div.r1200, 
                agrar.r1200, 
                legume.r1200, 
                rape.r1200, 
                special.crop.r1200,
                perm.grassland.r1200,
                grove.r1200,
                forest.r1200,
                urban.r1200,
                water.r1200) %>%
  na.omit() 
# 112 sites

set.seed(42)
indices_r1200<-sample(1:nrow(data_r1200),size=nrow(data_r1200)*0.8)
train_r1200 = data_r1200[indices_r1200,] # 89
test_r1200 = data_r1200[-indices_r1200,] # 23



## 2 find ntree ####

RFmodel_1200 <- randomForest(butterfly.rich.total~., data = train_r1200, 
                                  ntree = 2000)
RFmodel_1200
plot(RFmodel_1200)
ntree_best <- 1000



## 3 find mtry best ####

# hyper_grid <- expand.grid(
#   mtry = seq(2, 26, by = 1), # mtry-max has to be smaller than max number of variables!
#   ntree = seq(ntree_best, ntree_best, by = ntree_best / 10),
#   OOB_RMSE = 0
# ) 
# 
# # r1200
# for(i in 1:nrow(hyper_grid)) {
#   model <- party::cforest(formula = butterfly.rich.total~.,
#                           data    = train_r1200,
#                           controls = cforest_unbiased(ntree = hyper_grid$ntree[i],
#                                                       mtry = hyper_grid$mtry[i]))
#   error_mtry_1200 <- caret:::cforestStats(model)
#   # add OOB error to grid
#   hyper_grid$OOB_RMSE[i] <- error_mtry_1200[1]
# }
# 
# # top 10 performing models
# tune_r1200 <- hyper_grid %>%
#   dplyr::arrange(OOB_RMSE) %>%
#   mutate(response.var = "rich_total_1200")
# 
# tune_r1200 %>%
#   head(10)
# 
# # save
# tune_r1200 %>%
#    write.xlsx(here("outputs", "models", "cRF", "cRF_mtry", "tune_mtry_rich_total_1200_2026-03-12.xlsx"))

# load 
tune_mtry <- readxl::read_excel(
  here("outputs", "models", "cRF", "cRF_mtry", "tune_mtry_rich_total_1200_2026-03-12.xlsx"))

my_mtry_theme <- function() {
  theme_test() + 
    theme(
      text = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 20,  color = "black", face = "bold"), 
      plot.subtitle = element_text(size = 18,  color = "black"), 
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_blank(),
      panel.spacing.x = unit(0.25, "cm")
    )
  }

plot_tune_mtry_rich_total <- tune_mtry %>% 
  ggplot(aes(y = OOB_RMSE, x = mtry)) +
  scale_x_continuous(breaks = seq(0, 23, by = 2)) +
  scale_y_continuous(breaks = seq(4, 5, by = 0.1)) +
  geom_point(shape = 21, size = 2) +
  geom_smooth(se = FALSE, fullrange = FALSE, linewidth = 0.5) +
  my_mtry_theme() +
  labs(x = "mtry", y = "OOB RSME", 
       subtitle = "best mtry = 12, ntree = 1000",
       title = "(A) Total species richness")
print(plot_tune_mtry_rich_total)
mtry_best <- 12 



## 4 final cRF ####

# cf_1200 <- party::cforest(formula = butterfly.rich.total ~ .,
#                                data    = train_r1200,
#                                controls = cforest_unbiased(ntree = ntree_best,
#                                                            mtry = mtry_best))
# 
# # save model
# saveRDS(cf_1200, here("outputs", "models", "cRF", "cRF_RData", "cf_1200_rich_total_2026-03-12.RData"))

# load model
cf_1200 <- readRDS(here("outputs", "models", "cRF", "cRF_RData", "cf_1200_rich_total_2026-03-12.RData"))

# get stats
stats_1200 <- caret:::cforestStats(cf_1200)
stats_1200 <- as.data.frame(stats_1200)
stats_1200

# save stats training
stats_1200 %>% write.xlsx(here("outputs", "models", "cRF", "cRF_training",
"stats_training_rich_total_1200_2026-03-12.xlsx"))



## 5 prediction ####

pred_1200 = predict(cf_1200, newdata = test_r1200)
pred_RMSE_1200 <- round(sqrt(mean((pred_1200-test_r1200$butterfly.rich.total)**2)),2)
pred_RMSE_1200
# 4.01
pred_R2_1200 <- round(cor(pred_1200, test_r1200$butterfly.rich.total)**2,2)
pred_R2_1200
# 0.58

# save
plot(pred_1200, test_r1200$butterfly.rich.total,
     xlab = "Predicted",
     ylab = "Observed")
title(main = "(A) Total species richness", cex.main = 1.5)



## 6 importance ####

# imp_r1200 <- party::varimp(cf_1200, conditional = TRUE) %>%
#   sort(decreasing=TRUE) %>%
#   enframe() %>%
#   dplyr::top_n(24) %>%
#   mutate(name = as_factor(name),
#          sum.value = sum(value[value > 0]),
#          rel.imp = round(value/sum.value,2),
#          value = round(value, 2))
# imp_r1200 %>% write.xlsx(here("outputs", "tables", "imp_r1200_rich_total.xlsx"))




# END ####