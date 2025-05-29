# title: "Aspat_Acute_Hyp"
# description : Code for analysis of Mean NDVI in the acute heat stress assay
# author: "Naomi Sananikone and Holland Elder"

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
library(here)
library(MASS)
library(knitr)

# Set your working directory
setwd('D:/Folder/Species/Acropora spathulata')

# Run the helperFunctions.R
source('Code/helperFunctions.R')

# Read in the data 
as_ahyper <- read_xlsx('Data/Aspath_Acute_Hyp.xlsx', trim_ws=TRUE)

#### Data exploration and vizualisation ####

# Turn into factors
as_ahyper <- as_ahyper |> 
  mutate(Treatment = factor(Treatment, levels = c('MMM', 'MMM4', 'MMM7','MMM10')),
         Genet = factor(Genet),
         Reef = factor(Reef),
         Tank = factor(Tank),
         User = factor(User))

# Data distribution
as_ahyper |> ggplot(aes(x=MeanNDVI)) + 
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2)

# Look at the data
as_ahyper |> ggplot(aes(y = MeanNDVI, x = Treatment)) +
  geom_point(aes(color = Genet), 
             position = position_jitter())+
  theme(legend.position = 'none')

# Verifying that we don't have less than 2 fragments per genet
as_ahyper |> 
  filter(!is.na(MeanNDVI)) |> 
  filter(MeanNDVI > 0.00) |> 
  group_by(Treatment) |> 
  count(Genet) |> 
  filter(n < 2)

#### Model ####
## Priors, model fitting and dignostics 

# Formula
as_ahyper.form <- bf(MeanNDVI ~ Treatment*Genet,
                     phi = ~Treatment,
                     family = Beta(link = "logit"))

# Help defining priors if needed
as_ahyper |> group_by(Treatment) |> 
  summarise(median(logit(MeanNDVI), na.rm = TRUE),
            mad(logit(MeanNDVI), na.rm = TRUE))

# Priors
priors <- prior(normal(-0.67,0.6), class = 'Intercept') +
  prior(normal(0,2), class = 'b') +
  prior(student_t(3,0,2.5), class = 'Intercept', dpar = 'phi') +
  prior(normal(0,1), class = 'b', dpar = 'phi') 

# Run model with sample priors only
as_ahyper.brm2 <- brm(as_ahyper.form, 
                     data = as_ahyper,
                     prior = priors,
                     sample_prior = 'yes',
                     iter = 15000,
                     warmup = 1000,
                     chains = 3, cores = 3,
                     thin = 30,
                     refresh = 0,
                     seed = 123,
                     control = list(adapt_delta = 0.99, max_treedepth = 20),
                     backend = 'cmdstan')

# Look at conditional effects to see if priors are too informative
as_ahyper.brm |>  
  conditional_effects() |> 
  plot(points = TRUE)

# Update model with actual data
as_ahyper.brm2 <- update(as_ahyper.brm,
                         sample_prior = 'yes',
                         refresh = 0,
                         cores = 3, 
                         seed = 123) 

# Look at conditional effects to see if priors are too wide
as_ahyper.brm2 |>  
  conditional_effects() |> 
  plot(points = TRUE)

# Model diagnostics
as_ahyper.brm2$fit |> stan_trace()
as_ahyper.brm2$fit |> stan_ac()
as_ahyper.brm2$fit |> stan_rhat()
as_ahyper.brm2$fit |> stan_ess()
as_ahyper.brm2 |> pp_check(type = 'dens_overlay', ndraws = 100)

as_ahyper.resids <- make_brms_dharma_res(as_ahyper.brm2, integerResponse = FALSE)
testUniformity(as_ahyper.resids)
plotResiduals(as_ahyper.resids, quantreg = FALSE)
testDispersion(as_ahyper.resids)


#### ED50 & ED25 Calculation, Ranking and Visualization ####

# Create a lookup table that maps Treatments to Temperature
as_ahyper.lookup <- tribble(
  ~Treatment, ~Temperature,
  "MMM", 29,
  "MMM4", 33,
  "MMM7", 36,
  "MMM10", 39,
)

# Only keep Genets that are common to ALL experiments and traits
common_genets <- as_ahyper |> 
  filter(!(Genet %in% c("1501", "1519", "1536", "1538", "1539", "1540", "1572", 
                        "799", "956", "992", "993",  "995", "997", "998", "999",
                        '1515', '1566', '801',"1504","1582","802","990","1567") 
           & !is.na(Genet))) |> 
  pull(Genet)
genets_amount <- n_distinct(common_genets)

# Get the predicted values for the response
as_ahyper.pred <- emmeans(as_ahyper.brm2, ~ Treatment | Genet, type = "response") |>
  summary() |>
  as_tibble() |>
  left_join(as_ahyper.lookup, by = c("Treatment" = "Treatment")) 


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
# ff(a$Temperature, a$.value)
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
# ff(a$Temperature, a$.value)
# ff(.x, .y)

## ED50s ##

# Calculate the full posteriors of each ED50, including weighted rank and ranks
as_ahyper.draws_ED50 <- emmeans(as_ahyper.brm2, ~ Treatment | Genet, type = "response") |>
  gather_emmeans_draws() |>
  left_join(as_ahyper.lookup, by = c("Treatment" = "Treatment")) |>
  mutate(.value = plogis(.value)) |>
  arrange(Temperature) |>
  ungroup() |>
  group_by(Genet, .draw) |>
  summarise(ED50 = ffed50(Temperature, .value)) |>
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
as_ahyper.medians_ED50 <- as_ahyper.draws_ED50 |> 
  group_by(Genet) |> 
  summarise(ED50_median = median(ED50, na.rm=TRUE),
            ED50_lower = HDInterval::hdi(ED50)[1],
            ED50_upper = HDInterval::hdi(ED50)[2]) |> 
  arrange(desc(ED50_median)) |> 
  mutate(ED50 = ED50_median,
         Weighted_rank = (ED50[1]-ED50)*genets_amount/(max(ED50)-min(ED50))+1, 
         Rank = rank(desc(ED50)))

## ED25s ##

# Calculate the full posteriors of each ED25, including weighted rank and ranks
as_ahyper.draws_ED25 <- emmeans(as_ahyper.brm2, ~ Treatment | Genet, type = "response") |>
  gather_emmeans_draws() |>
  left_join(as_ahyper.lookup, by = c("Treatment" = "Treatment")) |>
  mutate(.value = plogis(.value)) |>
  arrange(Temperature) |>
  ungroup() |>
  group_by(Genet, .draw) |>
  summarise(ED25 = ffed25(Temperature, .value)) |>
  ungroup() |>
  group_by(.draw) |>
  filter(ED25 != "NA") |> 
  filter(Genet %in% common_genets) |> 
  summarise(Genet,
            ED25) |> 
  group_by(Genet, .draw)

# Calculate the median for ED25s, Weighted _ranks, and Ranks draws
as_ahyper.medians_ED25 <- as_ahyper.draws_ED25 |> 
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
as_ahyper_ED50s <- as_ahyper.medians_ED50 |>
  mutate(Genet = factor(Genet, levels=as_ahyper.medians_ED50 |>  arrange(ED50) |> pull(Genet))) |> 
  ggplot(aes(x = ED50, y = Genet)) +
  geom_pointrange(aes(xmin = ED50_lower, xmax = ED50_upper)) +
  labs(x = "Effective Dose at 50% (C°)") +
  theme_light()

# Plot medians of ED25s
as_ahyper_ED25s <- as_ahyper.medians_ED25 |>
  mutate(Genet = factor(Genet, levels=as_ahyper.medians_ED25 |>  arrange(ED25) |> pull(Genet))) |> 
  ggplot(aes(x = ED25, y = Genet)) +
  geom_pointrange(aes(xmin = ED25_lower, xmax = ED25_upper)) +
  labs(x = "Effective Dose at 25% (C°)") +
  theme_light()

### FINAL PLOT ###

as_ahyper_EDs_perGenet <- as_ahyper.pred |> 
  mutate(Genet = factor(Genet, levels=as_ahyper.medians_ED50 |>  arrange(desc(ED50)) |> pull(Genet))) |> 
  filter(Genet %in% common_genets) |>  
  ggplot(aes(y = response, x = Temperature)) +
  facet_wrap(~Genet)+
  labs(y = "Mean NDVI", x ="Temperature (C°)") +
  geom_line(linewidth=0.75, lineend = "round", color = "#696969") +
  geom_vline(data = as_ahyper.medians_ED25, aes(xintercept=ED25, color="ED25"),linewidth=0.75) +
  geom_vline(data = as_ahyper.medians_ED25, aes(xintercept=ED25_lower, color="ED25"), alpha = 0.3) +
  geom_vline(data = as_ahyper.medians_ED25, aes(xintercept=ED25_upper, color="ED25"), alpha = 0.3) +
  geom_rect(data = as_ahyper.medians_ED25, 
            aes(xmin = ED25_lower, xmax = ED25_upper, ymin = -Inf, ymax = Inf, fill="ED25"), 
            alpha = 0.2, inherit.aes = FALSE)+
  geom_vline(data = as_ahyper.medians_ED50, aes(xintercept=ED50, color = "ED50"),linewidth=0.75) +
  geom_vline(data = as_ahyper.medians_ED50, aes(xintercept=ED50_lower, color = "ED50"), alpha = 0.3) +
  geom_vline(data = as_ahyper.medians_ED50, aes(xintercept=ED50_upper, color = "ED50"), alpha = 0.3) +
  geom_rect(data = as_ahyper.medians_ED50, 
            aes(xmin = ED50_lower, xmax = ED50_upper, ymin = -Inf, ymax = Inf, fill="ED50"),
            alpha = 0.2, inherit.aes = FALSE)+
  scale_color_manual(name = "Effective dose response", values = c("ED25" = "#f0a009", "ED50" = "#e74c3c")) +
  scale_fill_manual(name = "Effective dose response", values = c("ED25" = "#f0a009", "ED50" = "#e74c3c"))+
  theme_light()+
  theme(legend.position = "bottom",
        axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        strip.text = element_text(size=15, colour='#696969', margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(size = 1, fill="white"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))