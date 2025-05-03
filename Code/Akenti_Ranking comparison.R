# title: "Akenti_Comparison"
# description : Code for performance comparison between experiments/traits
# authors: "Naomi Sananikone"

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

setwd('D:/AIMS/Experiment/Analysis/Acropora kenti')

# Load previously output medians files
ahyper.medians <- read_excel("Output/Table/ak_ahyper.medians_ED50.xlsx")
aipam.medians <- read_excel("Output/Table/ak_aipam.medians_ED25.xlsx")
mhyper.medians <- read_excel("Output/Table/ak_mhyper.medians_ED50_DHW.xlsx")
mipam.medians <- read_excel("Output/Table/ak_mipam.medians_ED25_DHW.xlsx")
surv.medians <- read_excel("Output/Table/ak_surv.medians_ED50_DHW.xlsx")

# Add columns for experiment, trait and value type.
ahyper.medians <- ahyper.medians |> 
  mutate(Experiment = 'Acute', Trait = 'NDVI', Values = ED50) |> 
  subset(select=c(Genet,Experiment,Trait,Values, Rank, Weighted_rank))
aipam.medians <- aipam.medians |> 
  mutate(Experiment = 'Acute', Trait = 'FvFm', Values = ED25) |> 
  subset(select=c(Genet,Experiment,Trait,Values, Rank, Weighted_rank))
mhyper.medians <- mhyper.medians |> 
  mutate(Experiment = 'Moderate', Trait = 'NDVI', Values = ED50) |> 
  subset(select=c(Genet,Experiment,Trait,Values, Rank, Weighted_rank))
mipam.medians <- mipam.medians |> 
  mutate(Experiment = 'Moderate', Trait = 'FvFm', Values = ED25) |> 
  subset(select=c(Genet,Experiment,Trait,Values, Rank, Weighted_rank))
surv.medians <- surv.medians |> 
  mutate(Experiment = 'Moderate', Trait = 'Survival', Values = ED50) |> 
  subset(select=c(Genet,Experiment,Trait,Values, Rank, Weighted_rank))

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
ahyper.medians.sp <- New(ahyper.medians)
aipam.medians.sp <- New(aipam.medians)
mhyper.medians.sp <- New(mhyper.medians)
mipam.medians.sp <- New(mipam.medians)
surv.medians.sp <- New(surv.medians)

# Now let's merge this tables:
all_tables <- list(ahyper.medians.sp, aipam.medians.sp, mhyper.medians.sp, mipam.medians.sp,surv.medians.sp)
Results_cortable <- reduce(all_tables, full_join, by = c("Genet", "Type"))

# Verify that all genets are common to each table and if not, which ones are missing
length(unique(Results_cortable[["Genet"]])) # If ran before the following line, tells you how many Genets total.
Results_cortable_filtered <- na.omit(Results_cortable)
length(unique(Results_cortable_filtered[["Genet"]])) # Tells you how many Genets are left after cleaning.
setdiff(Results_cortable$Genet, Results_cortable_filtered$Genet )

# Let's compare our data now:

# Pairwise Spearman Correlation ####

# Compare raw ranks
SpearmanComparison <- 
  Results_cortable_filtered |> 
  group_by(Type) |> 
  summarise(mipam_mhyper = cor.test(Moderate_FvFm , Moderate_NDVI, method='spearman', exact = FALSE)$estimate,
            mipam_ahyper = cor.test(Moderate_FvFm , Acute_NDVI, method='spearman', exact = FALSE)$estimate,
            aipam_mhyper = cor.test(Acute_FvFm , Moderate_NDVI, method='spearman', exact = FALSE)$estimate,
            mipam_aipam = cor.test(Moderate_FvFm , Acute_FvFm, method='spearman', exact = FALSE)$estimate,
            ahyper_mhyper = cor.test(Acute_NDVI , Moderate_NDVI, method='spearman', exact = FALSE)$estimate,
            aipam_ahyper = cor.test(Acute_FvFm , Acute_NDVI, method='spearman', exact = FALSE)$estimate,
            surv_aipam = cor.test(Moderate_Survival , Acute_FvFm, method='spearman', exact = FALSE)$estimate,
            surv_ahyper = cor.test(Moderate_Survival , Acute_NDVI, method='spearman', exact = FALSE)$estimate,
            surv_mipam = cor.test(Moderate_Survival , Moderate_FvFm, method='spearman', exact = FALSE)$estimate,
            surv_mhyper = cor.test(Moderate_Survival , Moderate_NDVI, method='spearman', exact = FALSE)$estimate,
            mipam_mhyper_PVAL = cor.test(Moderate_FvFm , Moderate_NDVI, method='spearman', exact = FALSE)$p.value,
            mipam_ahyper_PVAL = cor.test(Moderate_FvFm , Acute_NDVI, method='spearman', exact = FALSE)$p.value,
            aipam_mhyper_PVAL = cor.test(Acute_FvFm , Moderate_NDVI, method='spearman', exact = FALSE)$p.value,
            aipam_ahyper_PVAL = cor.test(Acute_FvFm , Acute_NDVI, method='spearman', exact = FALSE)$p.value,
            mipam_aipam_PVAL = cor.test(Moderate_FvFm , Acute_FvFm, method='spearman', exact = FALSE)$p.value,
            ahyper_mhyper_PVAL = cor.test(Acute_NDVI , Moderate_NDVI, method='spearman', exact = FALSE)$p.value,
            surv_aipam_PVAL = cor.test(Moderate_Survival , Acute_FvFm, method='spearman', exact = FALSE)$p.value,
            surv_ahyper_PVAL = cor.test(Moderate_Survival , Acute_NDVI, method='spearman', exact = FALSE)$p.value,
            surv_mipam_PVAL = cor.test(Moderate_Survival , Moderate_FvFm, method='spearman', exact = FALSE)$p.value,
            surv_mhyper_PVAL = cor.test(Moderate_Survival , Moderate_NDVI, method='spearman', exact = FALSE)$p.value
            ) |> 
  as.data.frame() |> 
  pivot_longer(
    cols = -Type,  # Keep "Type" column, pivot everything else
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
  arrange(desc(Estimate))

# Display results
print(n=31, SpearmanComparison)

# CORRELATION PLOTS ####

cor.colors <- c('#511849',"#cc3399",'#900C3F','#C70039','#FF5733','#FF8D1A','#FFC300',
                '#EDDD53','#669900','#57C785','#00BAAD','#2A7B9B',"#666a86",'#3D3D6B')

# TRAITS WITHIN EXPERIMENT
# Figure to compare mhyper and mipam
  mhyper_mipam <- Results_cortable_filtered |>  
  filter(Type == 'Weighted_rank') |> 
  ggplot(aes(x=Moderate_NDVI,y=Moderate_FvFm)) + 
   geom_point(aes(color=Genet))+
   geom_smooth(method=lm, se=TRUE, color="black")+
   stat_cor(method = "spearman") + 
   scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
   theme_classic() +
   theme(legend.position = "none")
  
# Figure to compare ahyper and aipam
  ahyper_aipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Acute_NDVI,y=Acute_FvFm)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
# EXPERIMENTS PER TRAIT
  
# Figure to compare mhyper and ahyper
  mhyper_ahyper <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Moderate_NDVI,y=Acute_NDVI)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
 # Figure to compare aipam and mipam
  aipam_mipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Moderate_FvFm,y=Acute_FvFm)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
# CROSSED

  # Figure to compare mhyper and aipam
  mhyper_aipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Moderate_NDVI,y=Acute_FvFm)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
  # Figure to compare ahyper and mipam
  ahyper_mipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Acute_NDVI,y=Moderate_FvFm)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
# TO SURVIVAL
  
  # Figure to compare mhyper vs surv
  surv_mhyper <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Moderate_NDVI,y=Moderate_Survival)) + 
   geom_point(aes(color=Genet))+
   geom_smooth(method=lm, se=TRUE, color="black")+
   stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
# Figure to compare ahyper vs surv
  surv_ahyper <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Acute_NDVI,y=Moderate_Survival)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
# Figure to compare surv vs aipam
  surv_aipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Acute_FvFm,y=Moderate_Survival)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")
  
 # Figure to compare suvr vs mipam
  surv_mipam <- Results_cortable_filtered |>  
    filter(Type == 'Weighted_rank') |> 
    ggplot(aes(x=Moderate_FvFm,y=Moderate_Survival)) + 
    geom_point(aes(color=Genet))+
    geom_smooth(method=lm, se=TRUE, color="black")+
    stat_cor(method = "spearman")+ 
    scale_color_manual(values = colorRampPalette(cor.colors)(55)) +
    theme_classic() +
    theme(legend.position = "none")

# FINAL CORPLOT
  
cor.plot <- mget(c('surv_mipam','surv_aipam','surv_mhyper','surv_ahyper','aipam_mipam',
              'mhyper_ahyper','mhyper_mipam','ahyper_aipam','mhyper_aipam','ahyper_mipam'))
names(cor.plot) = c('surv_mipam','surv_aipam','surv_mhyper','surv_ahyper','aipam_mipam',
                    'mhyper_ahyper','mhyper_mipam','ahyper_aipam','mhyper_aipam','ahyper_mipam')
cor.legend <- cor.plot[["surv_mipam"]] + 
  guides(fill = "none", linetype = "none", shape = "none", size = "none", color = guide_legend(title = "Genet", nrow = 7)) +
  theme(legend.position = "bottom")
cor.plot.only <- lapply(cor.plot, function(p) p + theme(legend.position = "none"))

cor.legend.only <- get_legend(cor.legend)
cor.legend.only <- as_ggplot(cor.legend.only)

ak_corplot <- (cor.plot.only[[1]] + cor.plot.only[[5]] + cor.plot.only[[6]])/
  (cor.plot.only[[2]] + cor.plot.only[[7]] + cor.plot.only[[8]])/
  (cor.plot.only[[3]] + cor.plot.only[[9]] + cor.plot.only[[10]])/
  (cor.plot.only[[4]] + plot_spacer() + plot_spacer()) + inset_element(cor.legend.only, left = 0.3, bottom = 0, right = 1, top = 0.8)

# Let's make a heatmap now

# Tbale has to be in a different format
hm <- function(df) {
 
  df <- df |> 
    mutate(Experiment_Trait = paste(df$Experiment, df$Trait, sep = "_"))
  df <- df[ , ! colnames(df) %in% c("Experiment","Trait")]

}
ahyper.medians.hm <- hm(ahyper.medians)
aipam.medians.hm <- hm(aipam.medians)
mhyper.medians.hm <- hm(mhyper.medians)
mipam.medians.hm <- hm(mipam.medians)
surv.medians.hm <- hm(surv.medians)

Results_heatmap <- rbind(ahyper.medians.hm,aipam.medians.hm,mhyper.medians.hm,mipam.medians.hm, surv.medians.hm)

Results_heatmap |> 
  group_by(Experiment_Trait) |> 
  summarize(shapiro_test(Weighted_rank))
# None of them is normally distributed

colors <- c("#0474BA", "#00A7E1","#EBEBEB", "#FFA630","#F17720")
color_gradient <- colorRampPalette(colors)

# Check if it's colorblind friendly
barplot(1:5,col=colors)
barplot(1:5,col=dichromat(colors, type = "protan"))
barplot(1:5,col=dichromat(colors, type = "deutan"))
barplot(1:5,col=dichromat(colors, type = "tritan"))

# With ranks
ak_heatmap_Rank <- Results_heatmap |> 
  group_by(Experiment_Trait) |> 
  mutate(Genet = factor(Genet, 
                        levels = Results_heatmap |> 
                          filter(Experiment_Trait == 'Moderate_Survival') |> 
                          arrange(Weighted_rank) |> 
                          pull(Genet)),
         Experiment_Trait = factor(Experiment_Trait, 
                                   levels = c('Moderate_Survival', 
                                              'Moderate_NDVI', 
                                              'Moderate_FvFm', 
                                              'Acute_FvFm', 
                                              'Acute_NDVI')))  |> # Convert Genet to an ordered factor
  ungroup() |> 
  ggplot(aes(x = Experiment_Trait, y = Genet, fill = Rank)) +
  geom_tile() +
  geom_text(aes(label=Genet)) +
  scale_fill_gradientn(colors=color_gradient(55))+
  labs(x ="Experiment and trait") +
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

# With Weighted Ranks
ak_heatmap_Weighted <- Results_heatmap |> 
  group_by(Experiment_Trait) |> 
  mutate(Genet = factor(Genet, 
                        levels = Results_heatmap |> 
                          filter(Experiment_Trait == 'Moderate_Survival') |> 
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