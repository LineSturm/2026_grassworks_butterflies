#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project 
# Prepare butterfly data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Line Sturm
# 11 June 2026



# packages ####

library(renv)
library(here)
library(tidyverse)
library(writexl)
library(vegan)
library(hillR)
library(iNEXT)

# here::set_here()
here()



# START ####

rm(list = ls())
# installr::updateR(
#   browse_news = FALSE,
#   install_R = TRUE,
#   copy_packages = TRUE,
#   copy_site_files = TRUE,
#   keep_old_packages = FALSE,
#   update_packages = FALSE,
#   start_new_R = FALSE,
#   quit_R = TRUE
#   )




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A LOAD DATA ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## site coordinates ####
site_coords <- read_csv(
  here("data", "raw", "data_site_coords_blurred.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(.default = "?")) 
site_coords <- as.data.frame(site_coords)
str(site_coords)



## site natural region ####
nature_regions_Germany <- read_csv(
  here("data", "raw", "data_nature_regions_Germany.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(.default = "?")) %>%
  mutate(
    nat.reg.OD2 = as.factor(nat.reg.OD2),
    nat.reg.OD3 = as.factor(nat.reg.OD3),
    nat.reg.GL = as.factor(nat.reg.GL),
    nat.reg.NR = as.factor(nat.reg.NR))
str(nature_regions_Germany)



## site environment ####
site_data <- read_csv(
  here("data", "raw", "data_site_environment.csv"),
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
    region = fct_relevel(region, "north", "centre", "south"),
    site.type = fct_relevel(site.type, "negative", "restored", "positive"),
    rest.meth = fct_relevel(rest.meth, "cus", "res", "dih", "mga"),
    rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "dih", "mga", "positive")
    ) %>%
  distinct()
str(site_data)



## butterfly data  ####
butterfly_data_raw <- read_csv(
  here("data", "raw", "data_butterfly_raw.csv"), col_names = TRUE)
str(butterfly_data_raw)



## butterfly rl & traits ####
butterfly_rl_traits <- read_csv(
  here("data", "raw", "data_butterfly_rl_traits_grassworks.csv"), col_names = TRUE) %>%
  mutate(
    rang = as.factor(rang),
    rl.ger = as.factor(rl.ger),
    larval.feeding.type = as.factor(larval.feeding.type),
    habitat.preference = as.factor(habitat.preference),
    disp = as.factor(disp),
    target.species = as.factor(target.species)
    )
str(butterfly_rl_traits)



## vegetation variables ####
vegetation_variables <- read_csv(
  here("data", "raw", "data_vegetation_variables.csv"))
str(vegetation_variables)



## landscape variables ####
landscape_variables <- read_csv(
  here("data", "raw", "data_landscape_variables.csv")) 
str(landscape_variables)



## management appropriateness ####
management_appropriateness <- read_csv(
  here("data", "raw", "data_management_appropriateness.csv"))
str(management_appropriateness)




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B PREPARE DATA ###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## sum butterfly data ####
butterfly_data <- butterfly_data_raw %>%
  filter(
    # remove april observation
    observation != "April",
    # remove additionally found species 
    sub.transect != "zusatz")

# create subset of active and transect catches
butterfly_tac_sum <- butterfly_data %>%
  # group_by(id.site, butterfly.name) %>%
  summarise(butterfly.no=sum(butterfly.no), .by = c(id.site, butterfly.name)) 

# create subset of transect catches
butterfly_t_sum <- butterfly_data %>%
  filter(!sub.transect %in% c("N1", "N2")) %>%
  summarise(butterfly.no=sum(butterfly.no), .by = c(id.site, butterfly.name)) 



## add traits to butterfly data ####
butterfly_tac_sum <- butterfly_tac_sum %>%
  left_join(butterfly_rl_traits %>%
              dplyr::select(butterfly.name, target.species, family, 
                            subfamily.SegererHausmann2011, genus, rang),
            by = "butterfly.name") 

butterfly_t_sum <- butterfly_t_sum %>%
  left_join(butterfly_rl_traits %>%
              dplyr::select(butterfly.name, target.species, family, 
                            subfamily.SegererHausmann2011, genus, rang),
            by = "butterfly.name") 




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C CLEAN DATA #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# workflow:
# handle aggregates
# --> species is in an aggregate and other species of that aggregate occur on site
# --> all get aggregate name

# handle genus level species
# --> keep them if they don't occur with species of the same genus on site

# handle subfamily level species
# --> keep them if they don't occur with species of the same subfamily on site

# handle family level species
# --> keep them if they don't occur with species of the same family on site



## 1 handle aggregates ####

## observed aggregates
# Aricia agestis/artaxerxes-Komplex
# --> species were identified by distribution
# Colias hyale/alfacariensis-Komplex
# --> check
# Leptidea sinapis/juvernica-Komplex
# --> keep aggregate name
# Plebejus argus/idas/argyrognomon-Komplex
# --> check
# Zygaena purpuralis/minos-Komplex
# --> keep aggregate name


# check agg - transect + active catches
precheck_agg_tac <- butterfly_tac_sum %>%
  filter(rang == "SPE" |
           rang == "AGG") %>%
  dplyr::select(id.site, butterfly.name, rang) %>%
  distinct() %>%
  group_by(id.site) %>%
  mutate(hyale.agg = as.integer(butterfly.name == "Colias hyale/alfacariensis-Komplex"),
         hyale = as.integer(butterfly.name == "Colias hyale"),
         alfacariensis = as.integer(butterfly.name == "Colias alfacariensis"),
         argus.agg = as.integer(butterfly.name == "Plebejus argus/idas/argyrognomon-Komplex"),
         argus = as.integer(butterfly.name == "Plebejus argus"),
         idas = as.integer(butterfly.name == "Plebejus idas"),
         argyrognomon = as.integer(butterfly.name == "Plebejus argyrognomon")) %>%
  dplyr::select(-butterfly.name, - rang) %>%
  distinct()

check_agg_tac <- precheck_agg_tac %>%
  group_by(id.site) %>%
  summarise(
    check.hyale.agg = if_else(
      sum(hyale.agg, hyale, alfacariensis) > 1,
      "check", "ok"),
    check.argus.trw = if_else(
      sum(argus.agg, argus, idas, argyrognomon) > 1,
      "check", "ok")) %>%
  distinct()

if (any(check_agg_tac$check.hyale.agg == "check" | check_agg_tac$check.argus.trw == "check")) {
  print("to be checked")
} else {
  print("ok")
}

# if ok for transect and active catchess, also ok for transect only



## 2 handle genus taxa ####

# are there genus-level and minimum one species-level species at the same site?


### transect and active catches -------------------------------------------------
# Create a data frame with only the genus taxa
data_genus <- butterfly_tac_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "GEN")
sum(data_genus$butterfly.no)
# 19 individuals
# 10 cases of genus: Argynnis, Melitaea, Thymelicus, Zygaena

# Create a data frame with only the species and aggregate taxa
species_agg <- butterfly_tac_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "SPE" |
           rang == "AGG")

# Join the two data frames on matching id.site and name.plant
check_genus <- data_genus %>%
  inner_join(species_agg, by = "id.site",
             suffix = c(".genus", ".species"),
             relationship = "many-to-many") %>%
  filter(butterfly.name.genus == word(butterfly.name.species, 1)) %>%
  dplyr::select(id.site, butterfly.name.genus, butterfly.no.genus,
         butterfly.name.species, butterfly.no.species) %>%
  distinct(id.site, .keep_all = TRUE)
sum(check_genus$butterfly.no.genus)
# 6 cases where species and only-genus at the same site (12 individuals in total)

# delete cases
# because unclear if same species
butterfly_tac_2 <- butterfly_tac_sum %>%
  anti_join(check_genus, by = c("id.site", "butterfly.name" = "butterfly.name.genus"))


### transect only ---------------------------------------------------------------
# Create a data frame with only the genus taxa
data_genus_t <- butterfly_t_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "GEN")
sum(data_genus_t$butterfly.no)
# 7 individuals
# 4 cases of genus: Thymelicus, Zygaena

# Create a data frame with only the species and aggregate taxa
species_agg_t <- butterfly_t_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "SPE" |
           rang == "AGG")

# Join the two data frames on matching id.site and name.plant
check_genus_t <- data_genus_t %>%
  inner_join(species_agg_t, by = "id.site",
             suffix = c(".genus", ".species"),
             relationship = "many-to-many") %>%
  filter(butterfly.name.genus == word(butterfly.name.species, 1)) %>%
  dplyr::select(id.site, butterfly.name.genus, butterfly.no.genus,
         butterfly.name.species, butterfly.no.species) %>%
  distinct(id.site, .keep_all = TRUE)
sum(check_genus_t$butterfly.no.genus)
unique(check_genus_t$butterfly.name.genus)
# 4 cases where species and only-genus at the same site (7 individuals in total)

# delete cases
# because unclear if same species
butterfly_t_2 <- butterfly_t_sum %>%
  anti_join(check_genus_t, by = c("id.site", "butterfly.name" = "butterfly.name.genus"))



## 3 handle subfamily taxa ####

### transect and active catches -------------------------------------------------
# Create a data frame with only the subfamily taxa
data_sf <- butterfly_tac_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "SF")
sum(data_sf$butterfly.no)
# 471 individuals
# 183 cases of 8 subfamilies: Coliadinae, Heliconiinae, Hesperiinae, Lycaeninae,
# Nymphalinae, Pierinae, Polyommatinae, Satyrinae

# Create a data frame with only the species, aggregate and genus taxa
species_gen <- butterfly_tac_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang, subfamily.SegererHausmann2011) %>%
  filter(rang %in% c("SPE", "AGG", "GEN"))

# Join the two data frames on matching id.site and name.plant
check_sf <- data_sf %>%
  inner_join(species_gen, by = "id.site",
             suffix = c(".sf", ".species"),
             relationship = "many-to-many") %>%
  filter(butterfly.name.sf == subfamily.SegererHausmann2011) %>%
  dplyr::select(id.site, butterfly.name.sf, butterfly.no.sf,
         butterfly.name.species, butterfly.no.species) %>%
  distinct(id.site, .keep_all = TRUE)
sum(check_sf$butterfly.no.sf)
unique(check_sf$butterfly.name.sf)
# 133 cases of 6 subfamilies where species and only-genus at the same site (350 individuals in total)

# delete cases
# because unclear if same species
# but in total a high number of individuals!
butterfly_tac_3 <- butterfly_tac_2 %>%
  anti_join(check_sf, by = c("id.site", "butterfly.name" = "butterfly.name.sf"))


### transect only ---------------------------------------------------------------
# Create a data frame with only the subfamily taxa
data_sf_t <- butterfly_t_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "SF")
sum(data_sf_t$butterfly.no)
# 207 individuals
# 84 cases of 6 subfamilies: Coliadinae, Heliconiinae, Lycaeninae, Pierinae,
# Polyommatinae, Satyrinae

# Create a data frame with only the species, aggregate and genus taxa
species_gen_t <- butterfly_t_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang, subfamily.SegererHausmann2011) %>%
  filter(rang %in% c("SPE", "AGG", "GEN"))

# Join the two data frames on matching id.site and name.plant
check_sf_t <- data_sf_t %>%
  inner_join(species_gen_t, by = "id.site",
             suffix = c(".sf", ".species"),
             relationship = "many-to-many") %>%
  filter(butterfly.name.sf == subfamily.SegererHausmann2011) %>%
  dplyr::select(id.site, butterfly.name.sf, butterfly.no.sf,
         butterfly.name.species, butterfly.no.species) %>%
  distinct(id.site, .keep_all = TRUE)
sum(check_sf_t$butterfly.no.sf)
unique(check_sf_t$butterfly.name.sf)
# 84 cases of 4 subfamilies where species and only-genus at the same site (163 individuals in total)

# delete cases
# because unclear if same species
# but in total a high number of individuals
butterfly_t_3 <- butterfly_t_2 %>%
  anti_join(check_sf_t, by = c("id.site", "butterfly.name" = "butterfly.name.sf"))



## 4 handle family taxa ####

### transect and active catchess ------------------------------------------------
# Create a data frame with only the family taxa
data_fam <- butterfly_tac_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "FAM")
sum(data_fam$butterfly.no)
# 25 individuals
# 11 cases of 2 families: Hesperiidae, Nymphalidae

# Create a data frame with only the species, aggregate, genus and subfamily taxa
species_sf <- butterfly_tac_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang, family) %>%
  filter(rang %in% c("SPE", "AGG", "GEN", "SF"))

# Join the two data frames on matching id.site and name.plant
check_fam <- data_fam %>%
  inner_join(species_sf, by = "id.site",
             suffix = c(".fam", ".species"),
             relationship = "many-to-many") %>%
  filter(butterfly.name.fam == family) %>%
  dplyr::select(id.site, butterfly.name.fam, butterfly.no.fam,
         butterfly.name.species, butterfly.no.species) %>%
  distinct(id.site, .keep_all = TRUE)
sum(check_fam$butterfly.no.fam)
unique(check_fam$butterfly.name.fam)
# 11 cases of 2 subfamilies where species and only-genus at the same site (22 individuals in total)

# delete cases
# because unclear if same species
butterfly_tac_4 <- butterfly_tac_3 %>%
  anti_join(check_fam, by = c("id.site", "butterfly.name" = "butterfly.name.fam"))


### transect only ---------------------------------------------------------------
# Create a data frame with only the family taxa
data_fam_t <- butterfly_t_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang) %>%
  filter(rang == "FAM")
sum(data_fam_t$butterfly.no)
# 8 individuals
# 5 cases of 2 families: Hesperiidae, Nymphalidae

# Create a data frame with only the species, aggregate, genus and subfamily taxa
species_sf_t <- butterfly_t_sum %>%
  dplyr::select(id.site, butterfly.name, butterfly.no, rang, family) %>%
  filter(rang %in% c("SPE", "AGG", "GEN", "SF"))

# Join the two data frames on matching id.site and name.plant
check_fam_t <- data_fam_t %>%
  inner_join(species_sf_t, by = "id.site",
             suffix = c(".fam", ".species"),
             relationship = "many-to-many") %>%
  filter(butterfly.name.fam == family) %>%
  dplyr::select(id.site, butterfly.name.fam, butterfly.no.fam,
         butterfly.name.species, butterfly.no.species) %>%
  distinct(id.site, .keep_all = TRUE)
sum(check_fam_t$butterfly.no.fam)
unique(check_fam_t$butterfly.name.fam)
# 5 cases of 2 families where species and only-genus at the same site (7 individuals in total)

# delete cases
# because unclear if same species
butterfly_t_4 <- butterfly_t_3 %>%
  anti_join(check_fam_t, by = c("id.site", "butterfly.name" = "butterfly.name.fam"))



## 5 cleaned butterfly data ####

# remove butterfly name = NA or unkown


### transect + active catches  --------------------------------------------------
butterfly_tac_c <- butterfly_tac_4 %>%
  filter(!str_detect(butterfly.name, "^Unbekannt"),
         !is.na(butterfly.name)) %>%
  as.data.frame()

# deleted cases
butterfly_tac_deleted <- butterfly_tac_sum %>%
  anti_join(butterfly_tac_c, by = c("id.site", "butterfly.name"))
sum(butterfly_tac_deleted$butterfly.no)
unique(butterfly_tac_deleted$butterfly.name)
# 183 cases (428 individuals)


### transect only ---------------------------------------------------------------
butterfly_t_c <- butterfly_t_4 %>%
  filter(!str_detect(butterfly.name, "^Unbekannt"),
         !is.na(butterfly.name)) %>%
  as.data.frame()

# deleted cases
butterfly_t_deleted <- butterfly_t_sum %>%
  anti_join(butterfly_t_c, by = c("id.site", "butterfly.name"))
sum(butterfly_t_deleted$butterfly.no)
unique(butterfly_t_deleted$butterfly.name)
# 115 cases (209 individuals)
# most cases deleted are from the north region




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C CALCULATE DIVERSITY ########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 species richness ####


### a) total richness -----------------------------------------------------------
total_richness <- butterfly_tac_c %>%
  group_by(id.site) %>%
  summarise(butterfly.rich.total = n_distinct(butterfly.name)) %>%
  ungroup()


### b) target richness ----------------------------------------------------------
target_richness <- butterfly_tac_c %>%
  filter(target.species == 1) %>%
  group_by(id.site) %>%
  summarise(butterfly.rich.target = n_distinct(butterfly.name)) %>%
  ungroup() %>%
  complete(id.site = unique(butterfly_tac_c$id.site), fill = list(. = 0))
target_richness[is.na(target_richness)] <- 0



## 2 abundance ####
# use of unprocessed data of transect catches only

### a) total abundance ----------------------------------------------------------
total_abundance <- butterfly_t_sum %>%
  group_by(id.site) %>%
  summarise(butterfly.abu.coridon = sum(butterfly.no))

# new total abundance 
# with P. coridon excluded in central region, since there was a mass occurence
total_abundance_new <- butterfly_t_sum %>%
  left_join(site_data %>% dplyr::select(id.site, region), by = "id.site") %>%
  filter(!(region == "centre" & butterfly.name == "Polyommatus coridon")) %>%
  group_by(id.site) %>%
  summarise(butterfly.abu = sum(butterfly.no)) %>%
  complete(id.site = unique(butterfly_t_c$id.site)) %>%
  replace(is.na(.), 0)


### b) target abundance ---------------------------------------------------------
target_abundance <- butterfly_t_sum %>%
  filter(target.species == 1) %>%
  group_by(id.site) %>%
  summarise(butterfly.target.abu.coridon = sum(butterfly.no)) %>%
  complete(id.site = unique(butterfly_t_c$id.site)) %>%
  replace(is.na(.), 0)


# new target abundance 
# with P. coridon excluded in central region, since there was a mass occurence
target_abundance_new <- butterfly_t_sum %>%
  left_join(site_data %>% dplyr::select(id.site, region), by = "id.site") %>%
  filter(!(region == "centre" & butterfly.name == "Polyommatus coridon"), 
         target.species == 1) %>%
  group_by(id.site) %>%
  summarise(butterfly.target.abu = sum(butterfly.no)) %>%
  complete(id.site = unique(butterfly_t_c$id.site)) %>%
  replace(is.na(.), 0)



## 3 Hill-Shannon & Simpson Index (Hill q = 1, 2) ####
# only transect data, excluding P. coridon

# prepare: converting long into wide data (columns = species, rows = sites)
abundances_wide_total <- pivot_wider(
  butterfly_t_c %>%
    left_join(site_data %>% dplyr::select(id.site, region), by = "id.site") %>%
    filter(!(region == "centre" & butterfly.name == "Polyommatus coridon")) %>%
    select(id.site, butterfly.name, butterfly.no),
  id_cols = id.site,
  names_from = butterfly.name,
  values_from = butterfly.no, 
  values_fill = 0
  )

abundances_wide_target <- pivot_wider(
  butterfly_t_c %>%
    left_join(site_data %>% dplyr::select(id.site, region), by = "id.site") %>%
    filter(!(region == "centre" & butterfly.name == "Polyommatus coridon"), 
           target.species == 1) %>%
    select(id.site, butterfly.name, butterfly.no),
  id_cols = id.site,
  names_from = butterfly.name,
  values_from = butterfly.no,
  values_fill = 0
  )

# Hill_total <- abundances_wide_total["id.site"]
# Hill_total$butterfly.Hill1.total <- hill_taxa(comm = abundances_wide_total[,2:84],  q = 1)
# Hill_total$butterfly.Hill2.total <- hill_taxa(comm = abundances_wide_total[,2:84],  q = 2)
# 
# Hill_target <- abundances_wide_target["id.site"]
# Hill_target$butterfly.Hill1.target <- hill_taxa(comm = abundances_wide_target[,2:49],  q = 1)
# Hill_target$butterfly.Hill2.target <- hill_taxa(comm = abundances_wide_target[,2:49],  q = 2)
# # 8 sites missing (no target butterflies)
# Hill_target <- Hill_target %>%
#   complete(id.site = unique(butterfly_t_c$id.site)) %>%
#   replace(is.na(.), 0)

# combine
# data_butterfly_Hill_Shannon_Simpson_diversity <- Hill_total %>%
#   left_join(Hill_target, by = "id.site")
# data_butterfly_Hill_Shannon_Simpson_diversity %>% 
#   write.csv(
#     here("data", "processed", "data_butterfly_Hill_Shannon_Simpson_diversity.csv"), 
#     row.names = F)
Hill_Shannon_Simpson_diversity <- read_csv(
  here("data", "processed", "data_butterfly_Hill_Shannon_Simpson_diversity.csv"), 
  col_names = TRUE) 




## 4 Sampling completeness / coverage ####
# only transect data, excluding P. coridon

# prepare:
abundances_inext_total <- pivot_wider(
  butterfly_t_c %>%
    left_join(site_data %>% dplyr::select(id.site, region), by = "id.site") %>%
    filter(!(region == "centre" & butterfly.name == "Polyommatus coridon")) %>%
    select(id.site, butterfly.name, butterfly.no),
  id_cols = butterfly.name,
  names_from = id.site,
  values_from = butterfly.no,
  values_fill = 0)
abundances_inext_total <- column_to_rownames(abundances_inext_total, var = names(abundances_inext_total)[1])


### a) Observed for q0 -------------------------------------------------------------------------
# out <- iNEXT(abundances_inext_total, q = c(0), datatype = "abundance", nboot = 100)
# df_inext <- bind_rows(out$iNextEst)
# df_clean <- df_inext %>%
#   filter(!is.na(SC)) %>%
#   mutate(has_CI = !is.na(SC.LCL) & !is.na(SC.UCL))
# Observed_sampling_cov <- df_clean %>%
#   filter(Method == "Observed") %>%
#   group_by(Assemblage) %>%
#   slice(1) %>%
#   ungroup() %>%
#   rename(qD.obs = qD,
#          id.site = Assemblage,
#          hill0.observedsampling.cov = SC) %>%
#   select(-qD.obs)
# Observed_sampling_cov %>%
#   write.csv(file = here::here("data", "processed", "data_butterfly_hill0_observed_SC.csv"),
#             row.names = F)
Observed_sampling_cov <- read.csv(
  here("data", "processed", "data_butterfly_hill0_observed_SC.csv")) 
mean(Observed_sampling_cov$hill0.observed.cov) # 0.9086212
sd(Observed_sampling_cov$hill0.observed.cov)   # 0.1198141


### b) Hill-Simpson diversity for equal sampling coverage ----------------------
# Hill2_equal_sampling_cov <- estimateD(abundances_inext_total, q=c(2), datatype = "abundance", nboot = 100,
#                    base = "coverage")
# Hill2_equal_sampling_cov %>%
#   write.csv(file = here::here("data", "processed", "data_butterfly_Hill2_equal_SC.csv"),
#             row.names = F) 
Hill2_equal_sampling_cov <-
  read.csv(file = here::here("data", "processed", "data_butterfly_Hill2_equal_SC.csv")) %>%
  select(Assemblage, qD) %>%
  rename(hill2.equal.cov = qD,
         id.site = Assemblage)



## 5 combine butterfly diversity ####
butterfly_diversity <- total_richness %>%
  left_join(target_richness, by = "id.site") %>%
  left_join(total_abundance_new, by = "id.site") %>%
  left_join(target_abundance_new, by = "id.site") %>%
  left_join(Hill_Shannon_Simpson_diversity, by = "id.site") %>%
  left_join(Observed_sampling_cov, by = "id.site") %>%
  left_join(Hill2_equal_sampling_cov, by = "id.site") 
str(butterfly_diversity)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D DESCRIPTIVE ANALYSES ######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("site_coords",
                          "nature_regions_Germany",
                          "site_data", 
                          "butterfly_rl_traits",
                          "butterfly_t_c", 
                          "butterfly_diversity",
                          "vegetation_variables",
                          "landscape_variables",
                          "management_appropriateness")))



## correlation Hill-Simpson at observed vs. equal sample coverage ####
cor.test(butterfly_diversity$butterfly.Hill2.total, butterfly_diversity$hill2.equal.cov, method = "spearman")
# rho = 0.9627726 , p < 2.2e-16



## mean, median of species per site type ####
spe_counts_site_type <- butterfly_t_c %>% 
  filter(rang == "SPE" | rang == "AGG") %>%
  dplyr::select(id.site, butterfly.name, butterfly.no) %>% 
  # change name to agg.name 
  mutate(
    butterfly.name = 
      if_else(butterfly.name == "Aricia agestis", "Aricia agestis/artaxerxes-Komplex", 
              if_else(butterfly.name == "Colias hyale", "Colias hyale/alfacariensis-Komplex",
                      if_else(butterfly.name == "Melitaea aurelia", 
                              "Melitaea athalia/aurelia/britomartis-Komplex", 
                              if_else(butterfly.name == "Plebejus argus", 
                                      "Plebejus argus/idas/argyrognomon-Komplex", 
                                      butterfly.name)
                              )
                      )
              )
    ) %>%
  left_join(site_data %>% 
              dplyr::select(id.site, site.type),
            by = "id.site") %>%
  left_join(butterfly_rl_traits %>% 
              dplyr::select(butterfly.name, rl.ger, target.species),
            by = "butterfly.name") %>%
  arrange(butterfly.name) %>%
  dplyr::summarise(no.mean = round(mean(butterfly.no),2), 
                   no.median = round(median(butterfly.no),2),
                   sd.butterflies = round(sd(butterfly.no),2), 
                   .by = c(butterfly.name, target.species, rl.ger, site.type)) %>%  
  pivot_wider(id_cols = c(butterfly.name, target.species, rl.ger),
              names_from = site.type,
              values_from = c(no.mean, no.median, sd.butterflies),
              values_fill = 0)

# for P. coridon, but when excluded in centre region
spe_counts_site_type_coridon_new <- butterfly_t_c %>% 
  filter(butterfly.name == "Polyommatus coridon") %>%
  dplyr::select(id.site, butterfly.name, butterfly.no) %>%
  left_join(site_data %>% 
              dplyr::select(id.site, site.type, region),
            by = "id.site") %>%
  filter(!(region == "centre" & butterfly.name == "Polyommatus coridon")) %>%
  drop_na() %>%
  dplyr::summarise(no.mean = round(mean(butterfly.no) ,2), 
                   no.median = round(median(butterfly.no),2),
                   sd.butterflies = round(sd(butterfly.no),2),
                   .by = c(butterfly.name, site.type)) %>% 
  pivot_wider(id_cols = butterfly.name,
              names_from = site.type,
              values_from = c(no.mean, no.median, sd.butterflies),
              values_fill = 0)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E COMBINE DATA ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_all <- site_coords %>%
  left_join(nature_regions_Germany, by = "id.site") %>%
  left_join(site_data, by = "id.site") %>%
  left_join(butterfly_diversity, by = "id.site") %>%
  left_join(vegetation_variables, by = "id.site") %>%
  left_join(management_appropriateness, by = "id.site") %>%
  left_join(landscape_variables, by = "id.site")
str(data_all)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F EXPORT ####################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all",
                          "spe_counts_site_type",
                          "spe_counts_site_type_coridon_new")))


data_all %>% 
  write.csv(file = here::here("data", "processed", "data_all.csv"), row.names = F)

spe_counts_site_type %>%
  write_xlsx(
    here("outputs", "tables", "spe_counts_site_type.xlsx"))

spe_counts_site_type_coridon_new %>%
  write_xlsx(
    here("outputs", "tables", "spe_counts_site_type_coridon_new.xlsx"))




# END ####

sessionInfo()