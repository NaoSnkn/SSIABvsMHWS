# title: "Akenti_Moderate_Hyp"
# description : Code for analysis of Mean NDVI in the bleaching simulation
# authors: "Naomi Sananikone and Holland Elder"

library(tidyverse)      
library(rstanarm)       
library(cmdstanr)       
library(brms)           
library(coda)           
library(bayesplot)      
library(DHARMa)         
library(rstan)          
library(emmeans)        
library(broom)          
library(tidybayes)      
library(ggeffects)      
library(broom.mixed)    
library(bayestestR)     
library(see)            
library(easystats)      
library(patchwork)      
library(modelsummary)
library(car)            
library(easystats)
library(modelsummary)
library(readxl)
library(ggplot2)

setwd('D:/AIMS/Experiment/Analysis/Acropora kenti')

# Run the helperFunctions.R
source('Code/helperFunctions.R')

# Read in the data
ak_mhyper <- read_xlsx('Data/Akenti_Moderate_Hyp.xlsx', trim_ws=TRUE)

#### Data exploration and vizualisation ####

# Turn into factors
ak_mhyper <- ak_mhyper |> 
  mutate(Treatment = factor(Treatment),
         Date = factor(Date),
         Genet = factor(Genet),
         Reef = factor(Reef),
         Tank = factor(Tank),
         DHW = factor(DHW))

# Data distribution
ak_mhyper |> ggplot(aes(x=MeanNDVI)) + 
    geom_histogram(aes(y=..density..), alpha = .6)+
    geom_density(alpha=.2)

# Look at the data
ak_mhyper |> ggplot(aes(y = MeanNDVI, x = DHW, color=Genet)) +
  geom_point(position = position_jitter())+
  geom_smooth(aes(x=DHW, y=MeanNDVI)) +
  facet_wrap(~Treatment) +
  theme(legend.position = 'none')

ak_mhyper |> ggplot(aes(y = MeanNDVI, x = Treatment, color = DHW)) +
  geom_point(aes(color = DHW), 
             position = position_jitter())+
  geom_smooth(aes(group = DHW, color = DHW), se = FALSE) 

# Verifying that we don't have less than 2 fragments per genet

print(n=30, ak_mhyper |>
        filter(!is.na(MeanNDVI)) |> 
        filter(MeanNDVI > 0.00) |> 
        group_by(Date, Treatment) |> 
        count(Genet) |> 
        filter(n < 2))

# Remove NAs
ak_mhyper <- ak_mhyper |>
  filter(!is.na(MeanNDVI))

#### Model ####
## Priors, model fitting and dignostics 

# Formula
ak_mhyper.form <- bf(MeanNDVI ~ Treatment*Date*Genet,
                      phi = ~ Treatment:Date,
                      family = Beta(link = "logit"))

# Help defining priors if needed
ak_mhyper |> group_by(Treatment:Date) |> 
  summarise(median(logit(MeanNDVI), na.rm = TRUE),
            mad(logit(MeanNDVI), na.rm = TRUE))

# Priors
priors <- prior(normal(0,1), class = 'Intercept') +
          prior(normal(0,0.9), class = 'b') +
          prior(student_t(3,0,2.5), class = 'Intercept', dpar = 'phi') +
          prior(normal(0,1), class = 'b', dpar = 'phi') 

# Run model with sample priors only
ak_mhyper.brm <- brm(ak_mhyper.form, 
                    data = ak_mhyper,
                    prior = priors,
                    sample_prior = 'only',
                    iter = 9000,
                    warmup = 1000,
                    chains = 3, cores = 3,
                    thin = 20,
                    refresh = 0,
                    seed = 123,
                    control = list(adapt_delta = 0.99, max_treedepth = 20),
                    backend = 'cmdstan')

# Look at conditional effects to see if priors are too informative
ak_mhyper.brm |>  
  conditional_effects() |> 
  plot(points = TRUE)

# Update model with actual data
ak_mhyper.brm2 <- update(ak_mhyper.brm,
                       sample_prior = 'yes',
                       refresh = 0,
                       cores = 3, 
                       seed = 123) 

# Look at conditional effects to see if priors are too wide
ak_mhyper.brm2 |>  
  conditional_effects() |> 
  plot(points = TRUE)

# Model diagnostics
ak_mhyper.brm2$fit |> stan_trace()
ak_mhyper.brm2$fit |> stan_ac()
ak_mhyper.brm2$fit |> stan_rhat()
ak_mhyper.brm2$fit |> stan_ess()
ak_mhyper.brm2 |> pp_check(type = 'dens_overlay', ndraws = 100)

ak_mhyper.resids <- make_brms_dharma_res(ak_mhyper.brm2, integerResponse = FALSE)
testUniformity(ak_mhyper.resids)
plotResiduals(ak_mhyper.resids, quantreg = FALSE)
testDispersion(ak_mhyper.resids)

#### ED50 & ED25 Calculation, Ranking and Visualization ####

# Create a lookup tables that maps Date to Day and DHW
ak_mhyper.lookup_Day <- tribble(
  ~Date, ~Day,
  "20240322", 0,
  "20240404", 13,
  "20240410", 19,
  "20240415", 24,
  "20240419", 28)

ak_mhyper.lookup_DHW <- tribble(
  ~Date, ~DHW,
  "20240322", 0,
  "20240404", 4.19,
  "20240410", 6.86,
  "20240415", 9.09,
  "20240419", 10.87)

# Define a function for calculating the ED50
ffed50 <- function(.x, .y) {
  m <- max(.y)
  m2 <- m/2
  i_T2 <- which(.y < m2)[1]
  if (is.na(i_T2)) {
    return(.x[4])
  }
  i_T1 <- rev(which(.y[1:i_T2] > m2))[1]
  ED50 <- .x[i_T1] + ((m2 - .y[i_T1]) / (.y[i_T2] - .y[i_T1]) * (.x[i_T2] - .x[i_T1]))
  ED50
}
 # ff(a$Treatment, a$.value)
 # ff(.x, .y)

# Define a function for calculating the ED25
ffed25 <- function(.x, .y) {
  m <- max(.y)
  m2 <- m*0.75
  i_T2 <- which(.y < m2)[1]
  if (is.na(i_T2)) {
    return(.x[4])
  }
  i_T1 <- rev(which(.y[1:i_T2] > m2))[1]
  ED25 <- .x[i_T1] + ((m2 - .y[i_T1]) / (.y[i_T2] - .y[i_T1]) * (.x[i_T2] - .x[i_T1]))
  ED25
}
 # ff(a$Treatment, a$.value)
 # ff(.x, .y)

# Only keep Genets that are common to ALL experiments and traits
common_genets <- ak_mhyper.brm2[[2]] |> 
  filter(!(Genet %in% c("1868","1869","1899","1961","1968","1972?","2008","2035","1965","2010") 
           & !is.na(Genet))) |> 
  pull(Genet)
genets_amount <- n_distinct(common_genets)

## Get the predicted values for the response

### With DHW ###

# Get the predicted values for the response
ak_mhyper.pred_DHW <- emmeans(ak_mhyper.brm2, ~ Date | Genet | Treatment, type = "response") |>
  summary() |>
  as_tibble() |>
  left_join(ak_mhyper.lookup_DHW, by = c("Date" = "Date"))

## ED50s ##

# Calculate the full posteriors of each ED50, including weighted rank and ranks
ak_mhyper.draws_ED50_DHW <- emmeans(ak_mhyper.brm2, ~ Date | Genet | Treatment, type = "response") |>
  gather_emmeans_draws() |>
  filter(Treatment == '31.5') |> 
  left_join(ak_mhyper.lookup_DHW, by = c("Date" = "Date")) |>
  mutate(.value = plogis(.value)) |>
  arrange(Date) |>
  ungroup() |>
  group_by(Genet, .draw) |>
  summarise(ED50 = ffed50(DHW, .value)) |>
  ungroup() |>
  group_by(.draw) |>
  arrange(desc(ED50)) |> 
  filter(ED50 != "NA") |>
  filter(Genet %in% common_genets) |> 
  summarise(Genet,
            ED50,
            Weighted_rank = (ED50[1]-ED50)*genets_amount/(max(ED50)-min(ED50))+1, 
            Rank = rank(desc(ED50))) |> 
  group_by(Genet, .draw)

# Calculate the median for ED50s, Weighted _ranks, and Ranks draws
ak_mhyper.medians_ED50_DHW <- ak_mhyper.draws_ED50_DHW |> 
  group_by(Genet) |> 
  summarise(ED50_median = median(ED50, na.rm=TRUE),
            ED50_lower = HDInterval::hdi(ED50)[1],
            ED50_upper = HDInterval::hdi(ED50)[2]) |> 
  arrange(desc(ED50_median)) |> 
  mutate(ED50 = ED50_median,
         Weighted_rank = (ED50[1]-ED50)*genets_amount/(max(ED50)-min(ED50))+1, 
         Rank = rank(desc(ED50)))

## ED25s ##

## Calculate the full posteriors of each ED25, including weighted rank and ranks
ak_mhyper.draws_ED25_DHW <- emmeans(ak_mhyper.brm2, ~ Date | Genet | Treatment, type = "response") |>
  gather_emmeans_draws() |>
  filter(Treatment == '31.5') |> 
  left_join(ak_mhyper.lookup_DHW, by = c("Date" = "Date")) |>
  mutate(.value = plogis(.value)) |>
  arrange(Date) |>
  ungroup() |>
  group_by(Genet, .draw) |>
  summarise(ED25 = ffed25(DHW, .value)) |>
  ungroup() |>
  group_by(.draw) |>
  arrange(desc(ED25)) |>
  filter(ED25 != "NA") |> 
  filter(Genet %in% common_genets) |> 
  summarise(Genet,
            ED25) |> 
  group_by(Genet, .draw)

# Calculate the median for ED25s, Weighted _ranks, and Ranks draws
ak_mhyper.medians_ED25_DHW <- ak_mhyper.draws_ED25_DHW |> 
  group_by(Genet) |> 
  summarise(ED25_median = median(ED25, na.rm=TRUE),
            ED25_lower = HDInterval::hdi(ED25)[1],
            ED25_upper = HDInterval::hdi(ED25)[2]) |> 
  arrange(desc(ED25_median)) |> 
  mutate(ED25 = ED25_median,
         Weighted_rank = (ED25[1]-ED25)*genets_amount/(max(ED25)-min(ED25))+1, 
         Rank = rank(desc(ED25)))

### Plots ###

# Plot medians of ED50s
ak_mhyper_ED50s_DHW <- ak_mhyper.medians_ED50_DHW |>
  mutate(Genet = factor(Genet, levels=ak_mhyper.medians_ED50_DHW |>  arrange(ED50) |> pull(Genet))) |> 
  ggplot(aes(x = ED50, y = Genet)) +
  geom_pointrange(aes(xmin = ED50_lower, xmax = ED50_upper)) +
  labs(x = "Effective Dose at 50% (C°)") +
  theme_light()

# Plot medians of ED25s
ak_mhyper_ED25s_DHW <- ak_mhyper.medians_ED25_DHW |>
  mutate(Genet = factor(Genet, levels=ak_mhyper.medians_ED25_DHW |>  arrange(ED25) |> pull(Genet))) |> 
  ggplot(aes(x = ED25, y = Genet)) +
  geom_pointrange(aes(xmin = ED25_lower, xmax = ED25_upper)) +
  labs(x = "Effective Dose at 25% (C°)") +
  theme_light()

### With days ###

# Get the predicted values for the response
ak_mhyper.pred_Day <- emmeans(ak_mhyper.brm2, ~ Date | Genet | Treatment, type = "response") |>
  summary() |>
  as_tibble() |>
  left_join(ak_mhyper.lookup_Day, by = c("Date" = "Date"))

## ED50s ##

## Calculate the full posteriors of each ED50, including weighted rank and ranks
ak_mhyper.draws_ED50_Day <- emmeans(ak_mhyper.brm2, ~ Date | Genet | Treatment, type = "response") |>
  gather_emmeans_draws() |>
  filter(Treatment == '31.5') |> 
  left_join(ak_mhyper.lookup_Day, by = c("Date" = "Date")) |>
  mutate(.value = plogis(.value)) |>
  arrange(Date) |>
  ungroup() |>
  group_by(Genet, .draw) |>
  summarise(ED50 = ffed50(Day, .value)) |>
  ungroup() |>
  group_by(.draw) |>
  arrange(desc(ED50)) |>
  filter(ED50 != "NA") |> 
  filter(Genet %in% common_genets) |> 
  summarise(Genet,
            ED50) |> 
  group_by(Genet, .draw)

# Calculate the median for ED50s, Weighted _ranks, and Ranks draws
ak_mhyper.medians_ED50_Day <- ak_mhyper.draws_ED50_Day |> 
  group_by(Genet) |> 
  summarise(ED50_median = median(ED50, na.rm=TRUE),
            ED50_lower = HDInterval::hdi(ED50)[1],
            ED50_upper = HDInterval::hdi(ED50)[2]) |> 
  arrange(desc(ED50_median)) |> 
  mutate(ED50 = ED50_median,
         Weighted_rank = (ED50[1]-ED50)*genets_amount/(max(ED50)-min(ED50))+1, 
         Rank = rank(desc(ED50)))

## ED25s ##

## Calculate the full posteriors of each ED25, including weighted rank and ranks
ak_mhyper.draws_ED25_Day <- emmeans(ak_mhyper.brm2, ~ Date | Genet | Treatment, type = "response") |>
  gather_emmeans_draws() |>
  filter(Treatment == '31.5') |> 
  left_join(ak_mhyper.lookup_Day, by = c("Date" = "Date")) |>
  mutate(.value = plogis(.value)) |>
  arrange(Date) |>
  ungroup() |>
  group_by(Genet, .draw) |>
  summarise(ED25 = ffed25(Day, .value)) |>
  ungroup() |>
  group_by(.draw) |>
  arrange(desc(ED25)) |>
  filter(ED25 != "NA") |> 
  filter(Genet %in% common_genets) |> 
  summarise(Genet,
            ED25) |> 
  group_by(Genet, .draw)

# Calculate the median for ED25s, Weighted _ranks, and Ranks draws
ak_mhyper.medians_ED25_Day <- ak_mhyper.draws_ED25_Day |> 
  group_by(Genet) |> 
  summarise(ED25_median = median(ED25, na.rm=TRUE),
            ED25_lower = HDInterval::hdi(ED25)[1],
            ED25_upper = HDInterval::hdi(ED25)[2]) |> 
  arrange(desc(ED25_median)) |> 
  mutate(ED25 = ED25_median,
         Weighted_rank = (ED25[1]-ED25)*genets_amount/(max(ED25)-min(ED25))+1, 
         Rank = rank(desc(ED25)))

### Plots ###

# Plot medians of ED50s
ak_mhyper_ED50s_Day <- ak_mhyper.medians_ED50_Day |>
  mutate(Genet = factor(Genet, levels=ak_mhyper.medians_ED50_Day |>  arrange(ED50) |> pull(Genet))) |> 
  ggplot(aes(x = ED50, y = Genet)) +
  geom_pointrange(aes(xmin = ED50_lower, xmax = ED50_upper)) +
  labs(x = "Effective Dose at 50% (C°)") +
  theme_light()

# Plot medians of ED25s
ak_mhyper_ED25s_Day <- ak_mhyper.medians_ED25_Day |>
  mutate(Genet = factor(Genet, levels=ak_mhyper.medians_ED25_Day |>  arrange(ED25) |> pull(Genet))) |> 
  ggplot(aes(x = ED25, y = Genet)) +
  geom_pointrange(aes(xmin = ED25_lower, xmax = ED25_upper)) +
  labs(x = "Effective Dose at 25% (C°)") +
  theme_light()

### FINAL PLOT ###

ak_mhyper_EDs_DHW_perGenet <- ak_mhyper.pred_DHW |> 
  mutate(Genet = factor(Genet, levels=ak_mhyper.medians_ED50_DHW |>  arrange(desc(ED50)) |> pull(Genet))) |> 
  filter(Genet %in% common_genets & Genet != 'NA' & Treatment == "31.5") |> 
  ggplot(aes(y = response, x = DHW)) +
  facet_wrap(~Genet)+
  labs(y = "Mean NDVI", x ="Degree Heating Weeks") +
  geom_line(linewidth=0.75, lineend = "round", color = "#696969") +
  geom_vline(data = ak_mhyper.medians_ED25_DHW, aes(xintercept=ED25, color="ED25"),linewidth=0.75) +
  geom_vline(data = ak_mhyper.medians_ED25_DHW, aes(xintercept=ED25_lower, color="ED25"), alpha = 0.3) +
  geom_vline(data = ak_mhyper.medians_ED25_DHW, aes(xintercept=ED25_upper, color="ED25"), alpha = 0.3) +
  geom_rect(data = ak_mhyper.medians_ED25_DHW, 
            aes(xmin = ED25_lower, xmax = ED25_upper, ymin = -Inf, ymax = Inf, fill="ED25"), 
            alpha = 0.2, inherit.aes = FALSE)+
  geom_vline(data = ak_mhyper.medians_ED50_DHW, aes(xintercept=ED50, color = "ED50"),linewidth=0.75) +
  geom_vline(data = ak_mhyper.medians_ED50_DHW, aes(xintercept=ED50_lower, color = "ED50"), alpha = 0.3) +
  geom_vline(data = ak_mhyper.medians_ED50_DHW, aes(xintercept=ED50_upper, color = "ED50"), alpha = 0.3) +
  geom_rect(data = ak_mhyper.medians_ED50_DHW, 
            aes(xmin = ED50_lower, xmax = ED50_upper, ymin = -Inf, ymax = Inf, fill="ED50"),
            alpha = 0.2, inherit.aes = FALSE)+
  scale_color_manual(name = "Effective dose response", values = c("ED25" = "#f0a009", "ED50" = "#e74c3c")) +
  scale_fill_manual(name = "Effective dose response", values = c("ED25" = "#f0a009", "ED50" = "#e74c3c"))+
  theme_light()+
  theme(legend.position = "bottom",
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=15, colour='#696969', margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(size = 1, fill="white"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))
