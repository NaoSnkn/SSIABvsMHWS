# title: "Experiments_Comparison"
# description : Code for performance comparison between experiments/traits
# authors: "Naomi Sananikone & Holland Elder"

library(readxl)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(irr)
library(PMCMRplus)
library(rstatix)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(colorspace)
library(dichromat)
library(Polychrome)
library(RColorBrewer)
library(recolorize)
library(ggplot2)
library(patchwork)  
library(writexl)
library(ggnewscale)

# If don't want to remake full tables for analyses:
      # Load Results_cortable_filtered.xlsx if needed for following steps
      # Load Results_heatmap.xlsx if needed for following steps
setwd('D:/Folder/Comparison')
Results_cortable_filtered <- read_excel('Data/Results_cortable_filtered.xlsx')
Results_heatmap <- read_excel('Data/Results_heatmap.xlsx')

# Load previously output medians files for all species and Experiment_Traits
setwd('D:/Folder/Species')

ak_ahyper.medians_ED50 <- read_excel("Acropora kenti/Tables/ak_ahyper.medians_ED50.xlsx")
ak_aipam.medians_ED50 <- read_excel("Acropora kenti/Tables/ak_aipam.medians_ED50.xlsx")
ak_mhyper.medians_ED50 <- read_excel("Acropora kenti/Tables/ak_mhyper.medians_ED50_DHW.xlsx")
ak_mipam.medians_ED25 <- read_excel("Acropora kenti/Tables/ak_mipam.medians_ED25_DHW.xlsx")
ak_surv.medians_ED50 <- read_excel("Acropora kenti/Tables/ak_surv.medians_ED50_DHW.xlsx")
ak_ahyper.medians_ED25 <- read_excel("Acropora kenti/Tables/ak_ahyper.medians_ED25.xlsx")
ak_aipam.medians_ED25 <- read_excel("Acropora kenti/Tables/ak_aipam.medians_ED25.xlsx")
ak_mhyper.medians_ED25 <- read_excel("Acropora kenti/Tables/ak_mhyper.medians_ED25_DHW.xlsx")
ak_surv.medians_ED25 <- read_excel("Acropora kenti/Tables/ak_surv.medians_ED25_DHW.xlsx")

ah_ahyper.medians_ED50 <- read_excel("Acropora hyacinthus/Tables/ah_ahyper.medians_ED50.xlsx")
ah_aipam.medians_ED25 <- read_excel("Acropora hyacinthus/Tables/ah_aipam.medians_ED25.xlsx")
ah_mhyper.medians_ED50 <- read_excel("Acropora hyacinthus/Tables/ah_mhyper.medians_ED50_DHW.xlsx")
ah_surv.medians_ED50 <- read_excel("Acropora hyacinthus/Tables/ah_surv.medians_ED50_DHW.xlsx")
ah_surv.medians_ED50$Genet[ah_surv.medians_ED50$Genet == '271q'] <- '271qm'
ah_surv.medians_ED50$Genet[ah_surv.medians_ED50$Genet == '3q2tageat'] <- '3qm2tageat'
ah_ahyper.medians_ED25 <- read_excel("Acropora hyacinthus/Tables/ah_ahyper.medians_ED25.xlsx")
ah_mhyper.medians_ED25 <- read_excel("Acropora hyacinthus/Tables/ah_mhyper.medians_ED25_DHW.xlsx")
ah_mipam.medians_ED25 <- read_excel("Acropora hyacinthus/Tables/ah_mipam.medians_ED25_DHW.xlsx")
ah_surv.medians_ED25 <- read_excel("Acropora hyacinthus/Tables/ah_surv.medians_ED25_DHW.xlsx")
ah_surv.medians_ED25$Genet[ah_surv.medians_ED25$Genet == '271q'] <- '271qm'
ah_surv.medians_ED25$Genet[ah_surv.medians_ED25$Genet == '3q2tageat'] <- '3qm2tageat'

as_ahyper.medians_ED50 <- read_excel("Acropora spathulata/Tables/as_ahyper.medians_ED50.xlsx")
as_aipam.medians_ED50 <- read_excel("Acropora spathulata/Tables/as_aipam.medians_ED50.xlsx")
as_mhyper.medians_ED50 <- read_excel("Acropora spathulata/Tables/as_mhyper.medians_ED50_DHW.xlsx")
as_mipam.medians_ED50 <- read_excel("Acropora spathulata/Tables/as_mipam.medians_ED50_DHW.xlsx")
as_surv.medians_ED50 <- read_excel("Acropora spathulata/Tables/as_surv.medians_ED50_DHW.xlsx")
as_surv.medians_ED50$Genet[as_surv.medians_ED50$Genet == '1502q'] <- '1502?'
as_surv.medians_ED50$Genet[as_surv.medians_ED50$Genet == '1522q'] <- '1522?'
as_surv.medians_ED50$Genet[as_surv.medians_ED50$Genet == '1518q'] <- '1518?'
as_surv.medians_ED50$Genet[as_surv.medians_ED50$Genet == '1524q'] <- '1524?'

#### Get right table formats for analysis ####

## Format for spearman correlations

# Add columns for experiment, trait and value type.
ak_ahyper.medians_ED50 <- ak_ahyper.medians_ED50 |> 
  mutate(Experiment = 'Acute', Trait = 'NDVI', Values = ED50, Species = 'Acropora kenti', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_aipam.medians_ED50 <- ak_aipam.medians_ED50 |> 
  mutate(Experiment = 'Acute', Trait = 'FvFm', Values = ED50, Species = 'Acropora kenti', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_mhyper.medians_ED50 <- ak_mhyper.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'NDVI', Values = ED50, Species = 'Acropora kenti', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_mipam.medians_ED25 <- ak_mipam.medians_ED25 |> 
  mutate(Experiment = 'Moderate', Trait = 'FvFm', Values = ED25, Species = 'Acropora kenti', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_surv.medians_ED50 <- ak_surv.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'Survival', Values = ED50, Species = 'Acropora kenti', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_ahyper.medians_ED25 <- ak_ahyper.medians_ED25 |> 
  mutate(Experiment = 'Acute', Trait = 'NDVI', Values = ED25, Species = 'Acropora kenti', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_aipam.medians_ED25 <- ak_aipam.medians_ED25 |> 
  mutate(Experiment = 'Acute', Trait = 'FvFm', Values = ED25, Species = 'Acropora kenti', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_mhyper.medians_ED25 <- ak_mhyper.medians_ED25 |> 
  mutate(Experiment = 'Moderate', Trait = 'NDVI', Values = ED25, Species = 'Acropora kenti', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ak_surv.medians_ED25 <- ak_surv.medians_ED25 |> 
  mutate(Experiment = 'Moderate', Trait = 'Survival', Values = ED25, Species = 'Acropora kenti', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))

ah_ahyper.medians_ED50 <- ah_ahyper.medians_ED50 |> 
  mutate(Experiment = 'Acute', Trait = 'NDVI', Values = ED50, Species = 'Acropora hyacinthus', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_aipam.medians_ED25 <- ah_aipam.medians_ED25 |> 
  mutate(Experiment = 'Acute', Trait = 'FvFm', Values = ED25, Species = 'Acropora hyacinthus', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_mhyper.medians_ED50 <- ah_mhyper.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'NDVI', Values = ED50, Species = 'Acropora hyacinthus', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_surv.medians_ED50 <- ah_surv.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'Survival', Values = ED50, Species = 'Acropora hyacinthus', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_ahyper.medians_ED25 <- ah_ahyper.medians_ED25 |> 
  mutate(Experiment = 'Acute', Trait = 'NDVI', Values = ED25, Species = 'Acropora hyacinthus', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_mhyper.medians_ED25 <- ah_mhyper.medians_ED25 |> 
  mutate(Experiment = 'Moderate', Trait = 'NDVI', Values = ED25, Species = 'Acropora hyacinthus', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_mipam.medians_ED25 <- ah_mipam.medians_ED25 |> 
  mutate(Experiment = 'Moderate', Trait = 'FvFm', Values = ED25, Species = 'Acropora hyacinthus', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
ah_surv.medians_ED25 <- ah_surv.medians_ED25 |> 
  mutate(Experiment = 'Moderate', Trait = 'Survival', Values = ED25, Species = 'Acropora hyacinthus', ED = 'ED25') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))

as_ahyper.medians_ED50 <- as_ahyper.medians_ED50 |> 
  mutate(Experiment = 'Acute', Trait = 'NDVI', Values = ED50, Species = 'Acropora spathulata', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
as_aipam.medians_ED50 <- as_aipam.medians_ED50 |> 
  mutate(Experiment = 'Acute', Trait = 'FvFm', Values = ED50, Species = 'Acropora spathulata', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
as_mhyper.medians_ED50 <- as_mhyper.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'NDVI', Values = ED50, Species = 'Acropora spathulata', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
as_mipam.medians_ED50 <- as_mipam.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'FvFm', Values = ED50, Species = 'Acropora spathulata', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))
as_surv.medians_ED50 <- as_surv.medians_ED50 |> 
  mutate(Experiment = 'Moderate', Trait = 'Survival', Values = ED50, Species = 'Acropora spathulata', ED = 'ED50') |> 
  subset(select=c(Species, Genet,Experiment,Trait,Values, Rank, Weighted_rank, ED))

# Now you have 4 data sets with medians of weighted ranks, ranks, ED50s, Genets, Experiments and Traits.
# Let's connect these tables according to Experiments and Traits. 
# For it to work, we need:
# - 1 column as Genets, 
# - 1 column of data (ranks, weighted ranks, ED50s per combination of traits and experiment)
# - 1 column indicating wether the indicated values in the row is Rank, Weighted Rank or ED50 
New <- function(df) {
  experiment <- unique(df$Experiment)
  trait <- unique(df$Trait)
  df <- df[, !(names(df) %in% c("Experiment", "Trait"))]
  df |> 
    pivot_longer(cols=c(Values, Rank, Weighted_rank),
                 names_to = "Type",
                 values_to = paste0(experiment, "_", trait))
}
ak_ahyper.medians.sp_ED50 <- New(ak_ahyper.medians_ED50)
ak_aipam.medians.sp_ED50 <- New(ak_aipam.medians_ED50)
ak_mhyper.medians.sp_ED50 <- New(ak_mhyper.medians_ED50)
ak_mipam.medians.sp_ED25 <- New(ak_mipam.medians_ED25)
ak_surv.medians.sp_ED50 <- New(ak_surv.medians_ED50)
ak_ahyper.medians.sp_ED25 <- New(ak_ahyper.medians_ED25)
ak_aipam.medians.sp_ED25 <- New(ak_aipam.medians_ED25)
ak_mhyper.medians.sp_ED25 <- New(ak_mhyper.medians_ED25)
ak_surv.medians.sp_ED25 <- New(ak_surv.medians_ED25)

ah_ahyper.medians.sp_ED50 <- New(ah_ahyper.medians_ED50)
ah_aipam.medians.sp_ED25 <- New(ah_aipam.medians_ED25)
ah_mhyper.medians.sp_ED50 <- New(ah_mhyper.medians_ED50)
ah_surv.medians.sp_ED50 <- New(ah_surv.medians_ED50)
ah_ahyper.medians.sp_ED25 <- New(ah_ahyper.medians_ED25)
ah_mhyper.medians.sp_ED25 <- New(ah_mhyper.medians_ED25)
ah_mipam.medians.sp_ED25 <- New(ah_mipam.medians_ED25)
ah_surv.medians.sp_ED25 <- New(ah_surv.medians_ED25)

as_ahyper.medians.sp_ED50 <- New(as_ahyper.medians_ED50)
as_aipam.medians.sp_ED50 <- New(as_aipam.medians_ED50)
as_mhyper.medians.sp_ED50 <- New(as_mhyper.medians_ED50)
as_mipam.medians.sp_ED50 <- New(as_mipam.medians_ED50)
as_surv.medians.sp_ED50 <- New(as_surv.medians_ED50)

# Now let's merge this tables:
ak_all_tables <- list(ak_ahyper.medians.sp_ED50, ak_aipam.medians.sp_ED50, ak_mhyper.medians.sp_ED50, 
                      ak_mipam.medians.sp_ED25, ak_surv.medians.sp_ED50, ak_ahyper.medians.sp_ED25, 
                      ak_aipam.medians.sp_ED25, ak_mhyper.medians.sp_ED25, ak_surv.medians.sp_ED25)
ak_Results_cortable <- reduce(ak_all_tables, full_join, by = c("Genet", "Type","Species","ED")) |> 
  mutate(across(.cols = tidyselect::matches("\\.x$"),
              .fns = ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
              .names = "{sub('\\\\.x$', '', .col)}")) |> 
  dplyr::select(-matches("\\.x$|\\.y$")) 

ah_all_tables <- list(ah_ahyper.medians.sp_ED50, ah_aipam.medians.sp_ED25, ah_mhyper.medians.sp_ED50, 
                      ah_surv.medians.sp_ED50, ah_ahyper.medians.sp_ED25, ah_mhyper.medians.sp_ED25, 
                      ah_mipam.medians.sp_ED25, ah_surv.medians.sp_ED25)
ah_Results_cortable <- reduce(ah_all_tables, full_join, by = c("Genet", "Type","Species","ED")) |> 
  mutate(across(.cols = tidyselect::matches("\\.x$"),
                .fns = ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column()))),
                .names = "{sub('\\\\.x$', '', .col)}")) |> 
  dplyr::select(-matches("\\.x$|\\.y$")) 

as_all_tables <- list(as_ahyper.medians.sp_ED50, as_aipam.medians.sp_ED50, as_mhyper.medians.sp_ED50, 
                      as_mipam.medians.sp_ED50, as_surv.medians.sp_ED50)
as_Results_cortable <- reduce(as_all_tables, full_join, by = c("Genet", "Type","Species","ED"))

# Get final table for correlation test
Results_cortable_filtered <- bind_rows(ak_Results_cortable, ah_Results_cortable, as_Results_cortable)

## Format for performance consistency comparison

# Change format of individual tables
hm <- function(df) {
  df <- df |> 
    mutate(Experiment_Trait = paste(df$Experiment, df$Trait, sep = "_"))
  df <- df[ , ! colnames(df) %in% c("Experiment","Trait")]
  
}
ak_ahyper.medians.hm_ED50 <- hm(ak_ahyper.medians_ED50)
ak_aipam.medians.hm_ED50 <- hm(ak_aipam.medians_ED50)
ak_mhyper.medians.hm_ED50 <- hm(ak_mhyper.medians_ED50)
ak_mipam.medians.hm_ED25 <- hm(ak_mipam.medians_ED25)
ak_surv.medians.hm_ED50 <- hm(ak_surv.medians_ED50)
ak_ahyper.medians.hm_ED25 <- hm(ak_ahyper.medians_ED25)
ak_aipam.medians.hm_ED25 <- hm(ak_aipam.medians_ED25)
ak_mhyper.medians.hm_ED25 <- hm(ak_mhyper.medians_ED25)
ak_surv.medians.hm_ED25 <- hm(ak_surv.medians_ED25)

ah_ahyper.medians.hm_ED50 <- hm(ah_ahyper.medians_ED50)
ah_aipam.medians.hm_ED25 <- hm(ah_aipam.medians_ED25)
ah_mhyper.medians.hm_ED50 <- hm(ah_mhyper.medians_ED50)
ah_mipam.medians.hm_ED25 <- hm(ah_mipam.medians_ED25)
ah_surv.medians.hm_ED50 <- hm(ah_surv.medians_ED50)
ah_ahyper.medians.hm_ED25 <- hm(ah_ahyper.medians_ED25)
ah_mhyper.medians.hm_ED25 <- hm(ah_mhyper.medians_ED25)
ah_surv.medians.hm_ED25 <- hm(ah_surv.medians_ED25)

as_ahyper.medians.hm_ED50 <- hm(as_ahyper.medians_ED50)
as_aipam.medians.hm_ED50 <- hm(as_aipam.medians_ED50)
as_mhyper.medians.hm_ED50 <- hm(as_mhyper.medians_ED50)
as_mipam.medians.hm_ED50 <- hm(as_mipam.medians_ED50)
as_surv.medians.hm_ED50 <- hm(as_surv.medians_ED50)

# Merge tables to get final table for comparison
Results_heatmap <- rbind(ak_ahyper.medians.hm_ED50, ak_aipam.medians.hm_ED50,
                         ak_mhyper.medians.hm_ED50, ak_mipam.medians.hm_ED25, 
                         ak_surv.medians.hm_ED50, ak_ahyper.medians.hm_ED25, 
                         ak_aipam.medians.hm_ED25, ak_mhyper.medians.hm_ED25, 
                         ak_surv.medians.hm_ED25, ah_ahyper.medians.hm_ED50,
                         ah_aipam.medians.hm_ED25, ah_mhyper.medians.hm_ED50,
                         ah_mipam.medians.hm_ED25, ah_surv.medians.hm_ED50,
                         ah_ahyper.medians.hm_ED25, ah_mhyper.medians.hm_ED25, 
                         ah_surv.medians.hm_ED25, as_ahyper.medians.hm_ED50,
                         as_aipam.medians.hm_ED50,as_mhyper.medians.hm_ED50,
                         as_mipam.medians.hm_ED50, as_surv.medians.hm_ED50)

#### Pairwise Spearman Correlation ####

# Functions that skips NAs results 
    # (i.e. ignoring instead of stopping when no ED50 available)
safe_spearman <- function(x, y) {
  tryCatch(cor.test(x, y, method = "spearman", exact = FALSE)$estimate, 
           error = function(e) NA)
} # Ignore NAs
safe_spearman_pval <- function(x, y) {
  tryCatch(cor.test(x, y, method = "spearman", exact = FALSE)$p.value, 
           error = function(e) NA)
} #Ignore NAs

# Spearman Correlation Table
SpearmanComparison <- 
  Results_cortable_filtered |> 
  group_by(Species, Type, ED) |> 
  summarise(
    mipam_mhyper       = safe_spearman(Moderate_FvFm, Moderate_NDVI),
    mipam_ahyper       = safe_spearman(Moderate_FvFm, Acute_NDVI),
    aipam_mhyper       = safe_spearman(Acute_FvFm, Moderate_NDVI),
    mipam_aipam        = safe_spearman(Moderate_FvFm, Acute_FvFm),
    ahyper_mhyper      = safe_spearman(Acute_NDVI, Moderate_NDVI),
    aipam_ahyper       = safe_spearman(Acute_FvFm, Acute_NDVI),
    surv_aipam         = safe_spearman(Moderate_Survival, Acute_FvFm),
    surv_ahyper        = safe_spearman(Moderate_Survival, Acute_NDVI),
    surv_mipam         = safe_spearman(Moderate_Survival, Moderate_FvFm),
    surv_mhyper        = safe_spearman(Moderate_Survival, Moderate_NDVI),
    mipam_mhyper_PVAL  = safe_spearman_pval(Moderate_FvFm, Moderate_NDVI),
    mipam_ahyper_PVAL  = safe_spearman_pval(Moderate_FvFm, Acute_NDVI),
    aipam_mhyper_PVAL  = safe_spearman_pval(Acute_FvFm, Moderate_NDVI),
    aipam_ahyper_PVAL  = safe_spearman_pval(Acute_FvFm, Acute_NDVI),
    mipam_aipam_PVAL   = safe_spearman_pval(Moderate_FvFm, Acute_FvFm),
    ahyper_mhyper_PVAL = safe_spearman_pval(Acute_NDVI, Moderate_NDVI),
    surv_aipam_PVAL    = safe_spearman_pval(Moderate_Survival, Acute_FvFm),
    surv_ahyper_PVAL   = safe_spearman_pval(Moderate_Survival, Acute_NDVI),
    surv_mipam_PVAL    = safe_spearman_pval(Moderate_Survival, Moderate_FvFm),
    surv_mhyper_PVAL   = safe_spearman_pval(Moderate_Survival, Moderate_NDVI),
    .groups = "drop"
  ) |> 
  as.data.frame() |> 
  pivot_longer(
    cols = -c(Type, Species, ED),
    names_to = "Compaired_experiment", 
    values_to = "Correlation"
  ) |> 
  mutate(
    TEST = ifelse(grepl("_PVAL$", Compaired_experiment), "P-value", "Estimate"), 
    Compaired_experiment = gsub("_PVAL$", "", Compaired_experiment)
  ) |> 
  pivot_wider(
    names_from = TEST, 
    values_from = Correlation
  ) |> 
  filter(!is.na(Estimate)) |> # Remove NAs do to absence of ED50 for some combination
  group_by(Species, Type, Compaired_experiment) |> 
  slice_max(ED, with_ties = FALSE) |>  # Only keep ED50 comparison when available
  ungroup() |> 
  arrange(Species, -Estimate)

### CORRELATION PLOTS ####

# Define colors
color.is <- c("#3db0bf","#FFC300","#8e73d0")
color.fill <- c("#8ed8e201","#f9dd8905","#ccc2e701")

# TRAITS WITHIN EXPERIMENT

# Figure to compare mhyper and mipam
mhyper_mipam <- Results_cortable_filtered |>  
  filter(Type == "Weighted_rank",
         (Species == "Acropora hyacinthus" & ED == "ED25") |
           (Species == "Acropora kenti" & ED == "ED25") |
           (Species == "Acropora spathulata" & ED == "ED50")) |> 
ggplot(aes(x=Moderate_NDVI,y=Moderate_FvFm, fill=Species)) + 
   geom_point(aes(color=Species))+
   geom_smooth(method=lm, se=TRUE, color = '#2a2a2a') +
   stat_cor(method = "spearman") + 
   scale_color_manual(values = color.is) +
   scale_fill_manual(values = color.fill) +
   theme_classic() +
   theme(legend.position = "none",
         strip.background = element_blank(),
         strip.text = element_text(size = 11)
         ) +
  facet_wrap(~Species, scales = "free")

# Figure to compare ahyper and aipam
ahyper_aipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank',
           (Species == "Acropora hyacinthus" & ED == "ED25") |
             (Species == "Acropora kenti" & ED == "ED50") |
             (Species == "Acropora spathulata" & ED == "ED50")) |> 
    ggplot(aes(x=Acute_NDVI,y=Acute_FvFm, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
# EXPERIMENTS PER TRAIT
  
# Figure to compare mhyper and ahyper
mhyper_ahyper <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank' & ED == "ED50") |> 
    ggplot(aes(x=Moderate_NDVI,y=Acute_NDVI, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
 # Figure to compare aipam and mipam
aipam_mipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank',
           (Species == "Acropora hyacinthus" & ED == "ED25") |
             (Species == "Acropora kenti" & ED == "ED25") |
             (Species == "Acropora spathulata" & ED == "ED50")) |> 
    ggplot(aes(x=Moderate_FvFm,y=Acute_FvFm, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
# CROSSED

  # Figure to compare mhyper and aipam
mhyper_aipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank',
           (Species == "Acropora hyacinthus" & ED == "ED25") |
             (Species == "Acropora kenti" & ED == "ED50") |
             (Species == "Acropora spathulata" & ED == "ED50")) |> 
    ggplot(aes(x=Moderate_NDVI,y=Acute_FvFm, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
  # Figure to compare ahyper and mipam
ahyper_mipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank',
           (Species == "Acropora hyacinthus" & ED == "ED25") |
             (Species == "Acropora kenti" & ED == "ED25") |
             (Species == "Acropora spathulata" & ED == "ED50")) |> 
    ggplot(aes(x=Acute_NDVI,y=Moderate_FvFm, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
# TO SURVIVAL
  
  # Figure to compare mhyper vs surv
surv_mhyper <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank' & ED == "ED50") |> 
    ggplot(aes(x=Moderate_NDVI,y=Moderate_Survival, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
# Figure to compare ahyper vs surv
surv_ahyper <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank' & ED == "ED50") |> 
    ggplot(aes(x=Acute_NDVI,y=Moderate_Survival, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
# Figure to compare surv vs aipam
surv_aipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank',
           (Species == "Acropora hyacinthus" & ED == "ED25") |
             (Species == "Acropora kenti" & ED == "ED50") |
             (Species == "Acropora spathulata" & ED == "ED50")) |> 
    ggplot(aes(x=Acute_FvFm,y=Moderate_Survival, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")
  
 # Figure to compare suvr vs mipam
surv_mipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank',
           (Species == "Acropora hyacinthus" & ED == "ED25") |
             (Species == "Acropora kenti" & ED == "ED25") |
             (Species == "Acropora spathulata" & ED == "ED50")) |> 
    ggplot(aes(x=Moderate_FvFm,y=Moderate_Survival, fill= Species)) + 
  geom_point(aes(color=Species))+
  geom_smooth(method=lm, se=TRUE, color = '#242424') +
  stat_cor(method = "spearman") + 
  scale_color_manual(values = color.is) +
  scale_fill_manual(values = color.fill) +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11)
  ) +
  facet_wrap(~Species, scales = "free")

## FINAL CORPLOTS ##

Surv_vs_Acute <- surv_aipam/surv_ahyper
Surv_vs_Simulation <- surv_mipam/surv_mhyper
Acute_vs_Simulation <- aipam_mipam/mhyper_ahyper
FvFm_vs_NDVI <- mhyper_mipam/ahyper_aipam

#### Performance consistency test ####

## COMPARE DHWs ##

# Get a table of Top and Worst 20% Genets for each assay and their corresponding Survival's DHW
Surv_20percent <- Results_heatmap |> 
  group_by(Species, ED, Experiment_Trait) |> 
  arrange(Weighted_rank, .by_group = TRUE) |> 
  left_join(Results_heatmap |> 
              filter(Experiment_Trait == "Moderate_Survival")  |> 
  dplyr::select(Genet, ED, DHW_survival = Values)) |> # Adds a DHW column with Survival's DHW values
  mutate(row_number = row_number(),
         n = n(),
         top_cutoff = ceiling(n * 0.2), # Gets top 20%
         bottom_cutoff = ceiling(n * 0.2))  |> # Gets bottom 20% 
  mutate(TwentyPercent = case_when(
    row_number <= top_cutoff ~ "Top",
    row_number > n - bottom_cutoff ~ "Bottom",
    TRUE ~ NA_character_))  |> 
  filter(!is.na(TwentyPercent)) |> 
  as.data.frame() |>
  dplyr::select(Species, Genet, Values, Rank, Weighted_rank, Experiment_Trait, DHW_survival, TwentyPercent, ED) |> 
  group_by(Species, Experiment_Trait, TwentyPercent, ED) |> 
  summarize(DHW_survival_mean = mean(DHW_survival), # Gets Mean and CI of DHW's
            DHW_survival_lower = mean(DHW_survival)-sd(DHW_survival),
            DHW_survival_upper = mean(DHW_survival)+sd(DHW_survival)) 

# Get means of Survival's EDs
PopMeans <- rbind(read_excel('Acropora kenti/Tables/ak_PopMean_ED50.xlsx'), 
                  read_excel('Acropora hyacinthus/Tables/ah_PopMean_ED25.xlsx'), 
                  read_excel('Acropora spathulata/Tables/as_PopMean_ED50.xlsx'))
# or load PopMeans <- read_excel('D:/Folder/Comparison/Data/PopMeans.xlsx')

# Define colors
species_hline_colors <- c("Acropora hyacinthus" = "#0f4d5b80",
                          "Acropora kenti" = "#7b561080",
                          "Acropora spathulata" = "#45348680")
species_palette_bottom <- c("Acropora hyacinthus" = "#007281",
                            "Acropora kenti" = "#c69400",
                            "Acropora spathulata" = "#453486")
species_palette_top <- c("Acropora hyacinthus" = "#3db0bf",
                         "Acropora kenti" = "#FFC300",
                         "Acropora spathulata" = "#8e73d0")
TP_colors <- c(
  setNames(species_palette_bottom, paste(names(species_palette_bottom), "Bottom", sep = "_")),
  setNames(species_palette_top, paste(names(species_palette_top), "Top", sep = "_"))
)
PopMeans <- PopMeans |> 
  mutate(color = case_when(
    Species == "Acropora hyacinthus" ~ "#00728105",
    Species == "Acropora kenti" ~ "#c6940005",
    Species == "Acropora spathulata" ~ "#5b479c05"
  ))
Surv_20percent <- Surv_20percent |> 
  mutate(TP_Group = paste(Species, TwentyPercent, sep = "_"))

## FINAL PERFORMANCE CONSISTENCY PLOTS ##

# Acute vs Survival
TopBottom_20percent_survDHW <- Surv_20percent |> 
  filter(Experiment_Trait %in% c('Acute_FvFm', 'Acute_NDVI','Moderate_Survival'),
         (Species == "Acropora hyacinthus" & ED == "ED25") |
           (Species == "Acropora kenti" & ED == "ED50") |
           (Species == "Acropora spathulata" & ED == "ED50")) |> 
  ggplot(aes(x = Experiment_Trait, y = DHW_survival_mean, color = TP_Group, shape = TwentyPercent)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_linerange(aes(ymin = DHW_survival_lower, ymax = DHW_survival_upper), 
                 position = position_dodge(width = 0.4)) +
  scale_color_manual(values = TP_colors) +
  new_scale("color") + 
  geom_hline(data = PopMeans, aes(yintercept = pop.mean, color = Species), linewidth=0.7) +
  scale_color_manual(values = species_hline_colors) +
  geom_rect(data = PopMeans,
            aes(ymin = pop.lower, ymax = pop.upper, xmin = -Inf, xmax = Inf, fill = Species),
            color = NA, alpha = 0.08, inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(PopMeans$color, PopMeans$Species)) +
  xlab("Heat Stress Assay") +
  ylab("Survival's Degree Heating Weeks") +  
  facet_wrap(~Species, scales = 'free') +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

# All experiments
TopBottom_20percent_survDHW_all <- Surv_20percent |> 
  filter((Species == "Acropora hyacinthus" & ED == "ED25") |
           (Species == "Acropora kenti" & ED == "ED50") |
           (Species == "Acropora spathulata" & ED == "ED50")) |> 
  ggplot(aes(x = Experiment_Trait, y = DHW_survival_mean, color = TP_Group, shape = TwentyPercent)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_linerange(aes(ymin = DHW_survival_lower, ymax = DHW_survival_upper), 
                 position = position_dodge(width = 0.4)) +
  scale_color_manual(values = TP_colors) +
  new_scale("color") + 
  geom_hline(data = PopMeans, aes(yintercept = pop.mean, color = Species), linewidth=0.7) +
  scale_color_manual(values = species_hline_colors) +
  geom_rect(data = PopMeans,
            aes(ymin = pop.lower, ymax = pop.upper, xmin = -Inf, xmax = Inf, fill = Species),
            color = NA, alpha = 0.08, inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(PopMeans$color, PopMeans$Species)) +
  xlab("Heat Stress Assay") +
  ylab("Survival's Degree Heating Weeks") +  
  facet_wrap(~Species, scales = 'free') +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

# Plots Survival's DHW's mean accross heat stress assays for top and bottom 20% genet
Mortality_Population_DHW  <- (ah_mortality_pop + # from Ahya_Survival.R code
                                ak_mortality_pop + # from Akenti_Survival.R code
                                as_mortality_pop) # from Aspat_Survival.R code

### HEATMAPS ###

colors <- c("#0474BA", "#00A7E1","#EBEBEB", "#FFA630","#F17720")
color_gradient <- colorRampPalette(colors)

# A. kenti
ak_heatmap <- Results_heatmap |>
  filter(Species == "Acropora kenti" &
           ED == 'ED25') |> 
  group_by(Experiment_Trait) |> 
  mutate(Genet = factor(Genet, 
                        levels = Results_heatmap |> 
                          filter(Experiment_Trait == 'Moderate_Survival' &
                                   ED == 'ED25') |> 
                          arrange(Weighted_rank) |> 
                          pull(Genet)),
         Experiment_Trait = factor(Experiment_Trait, 
                                   levels = c('Moderate_Survival', 
                                              'Moderate_NDVI', 
                                              'Moderate_FvFm', 
                                              'Acute_FvFm', 
                                              'Acute_NDVI')))  |>
  ungroup() |> 
  ggplot(aes(x = Experiment_Trait, y = Genet, fill = Weighted_rank)) +
  geom_tile() +
  geom_text(aes(label=Genet)) +
  scale_fill_gradientn(colors=color_gradient(55))+
  labs(x ="Experiment and trait") +
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

# A. hyacinthus
ah_heatmap <- Results_heatmap |>
  filter(Species == "Acropora hyacinthus" &
           ED == 'ED25') |> 
  group_by(Experiment_Trait) |> 
  mutate(Genet = factor(Genet, 
                        levels = Results_heatmap |> 
                          filter(Experiment_Trait == 'Moderate_Survival' &
                                   ED == 'ED25') |> 
                          arrange(Weighted_rank) |> 
                          pull(Genet)),
         Experiment_Trait = factor(Experiment_Trait, 
                                   levels = c('Moderate_Survival', 
                                              'Moderate_NDVI', 
                                              'Moderate_FvFm', 
                                              'Acute_FvFm', 
                                              'Acute_NDVI')))  |>
  ungroup() |> 
  ggplot(aes(x = Experiment_Trait, y = Genet, fill = Weighted_rank)) +
  geom_tile() +
  geom_text(aes(label=Genet)) +
  scale_fill_gradientn(colors=color_gradient(26))+
  labs(x ="Experiment and trait") +
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

# A. spathulata
as_heatmap <- Results_heatmap |>
  filter(Species == "Acropora spathulata" &
           ED == 'ED50') |> 
  group_by(Experiment_Trait) |> 
  mutate(Genet = factor(Genet, 
                        levels = Results_heatmap |> 
                          filter(Experiment_Trait == 'Moderate_Survival' &
                                   ED == 'ED50') |> 
                          arrange(Weighted_rank) |> 
                          pull(Genet)),
         Experiment_Trait = factor(Experiment_Trait, 
                                   levels = c('Moderate_Survival', 
                                              'Moderate_NDVI', 
                                              'Moderate_FvFm', 
                                              'Acute_FvFm', 
                                              'Acute_NDVI')))  |>
  ungroup() |> 
  ggplot(aes(x = Experiment_Trait, y = Genet, fill = Weighted_rank)) +
  geom_tile() +
  geom_text(aes(label=Genet)) +
  scale_fill_gradientn(colors=color_gradient(38))+
  labs(x ="Experiment and trait") +
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
