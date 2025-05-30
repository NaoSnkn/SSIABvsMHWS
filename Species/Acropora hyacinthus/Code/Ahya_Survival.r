# title: "Ahya_Survival"
# description : Survival Analysis for Acropora hyacinthus
# authors: "Holland Elder and Naomi Sananikone"

library(tidyverse) 
library(dplyr)
library(car)      
library(broom)    
library(ggfortify) 
library(knitr)    
library(emmeans)   
library(MASS)     
library(brms)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(standist) 
library(rstanarm)
library(cmdstanr)
library(ggeffects)
library(rstan)
library(DHARMa)
library(ggridges)
library(easystats)
library(patchwork)
library(tinytable)
library(modelsummary)
library(glmmTMB)
library(survival) 
library(vegan) 
library(coxme) 
library(corrplot) 
library(summarytools) 
library(cowplot)
library(RColorBrewer)
library(survminer)
library(readxl)
library(writexl)
library(openxlsx)

# Set your working directory
setwd('D:/Folder/Species/Acropora hyacinthus')

# Run the helperFunctions.R
source('Code/helperFunctions.R')

# Read in the data 
ah_surv <- read_xlsx('Data/Ahya_Survival.xlsx', trim_ws=TRUE)

# Turn into numeric and filter to treatment
ah_surv <- ah_surv |> 
  mutate(DHW =as.numeric(as.character(DHW)),
         Day =as.numeric(as.character(DayMeasured))) |> 
  filter(Temperature == '31.5')

# Data distribution
ah_surv |> ggplot(aes(x=Survival)) + 
  geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2)

# Look at the data
ah_surv |> ggplot(aes(y = Survival, x = DHW)) +
  geom_point(aes(color = Genet))+
  geom_smooth(aes(color = Genet))+
  theme(legend.position = 'none')+
  facet_wrap(~Genet)

# Genets 14 and 329 are only in one tank
ah_surv <- ah_surv |> filter(Genet != "14")
ah_surv <- ah_surv |> filter(Genet != "329")

# Let's make a table per Genet
GenetTables <- function(data) {
  unique_genets <- unique(data$Genet)
  GenetTables_list <- list()
  
  for (g in unique_genets) {
    subset_data <- subset(data, Genet == g)
    var_name <- paste0("ah_Genet", g, "_surv")
    assign(var_name, subset_data, envir = .GlobalEnv)
  }
  
  GenetTables_list <- setNames(
    lapply(unique_genets, function(g) {
      subset(data, Genet == g)
    }),
    paste0("ah_Genet", unique_genets, "_surv")
  )
  
  return(GenetTables_list)
} # Function makes 1 table per genets then a list of those tables
ah_GenetTables <- GenetTables(ah_surv)

#### DHW Model ####

# Formula
ah_survDHW.form <- bf(Survival|trials(1) ~ DHW + (1|Tank),
                      family = binomial(link = logit))

# Priors
priors <- prior(student_t(3, 0, 1), class = 'Intercept') +
  prior(normal(0,1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

# Run model for each genet with sample priors only
fit_brm <- function(table_list, formula, priors) {
  model_list <- list()
  
  for (name in names(table_list)) {
    message("Fitting model for: ", name)
    model <- brm(
      formula = formula,
      data = table_list[[name]],
      prior = priors,
      sample_prior = 'only',
      iter = 15000,
      warmup = 1000,
      chains = 3,
      cores = 3,
      thin = 10,
      control = list(adapt_delta = 0.99, max_treedepth = 20),
      refresh = 0,
      backend = "cmdstanr"
    )
    
    model_name <- paste0(name, "DHW.brm")
    model_list[[model_name]] <- model
  }
  
  return(model_list)
}
brm_models_DHW <- fit_brm(ah_GenetTables, ah_survDHW.form, priors)

# Look at conditional effects to see if priors are too informative
list_conditional_effects <- function(model_list) {
  for (model_name in names(model_list)) {
    cat("Plotting:", model_name, "\n")
    
    model <- model_list[[model_name]]
    
    # Try plotting and skip if it fails
    tryCatch({
      plot <- model |> conditional_effects() |> plot(points = TRUE) 
      print(plot)
      # Pause and wait for user input
      readline(prompt = "Press [Enter] to see the next plot...")
      
    }, error = function(e) {
      cat("⚠️ Could not plot", model_name, ":", e$message, "\n")
    })
  }
}
list_conditional_effects(brm_models_DHW)

# Update each genet's model with actual data
update_model_brm2 <- function(model_list) {
  updated_list <- list()
  
  for (model_name in names(model_list)) {
    cat("Updating model:", model_name, "\n")
    
    model <- model_list[[model_name]]
    
    tryCatch({
      updated_model <- update(
        model,
        sample_prior = "yes",
        refresh = 0,
        recompile = FALSE  # speeds up if model is already compiled
      )
      
      new_name <- sub("\\.brm$", ".brm2", model_name)
      updated_list[[new_name]] <- updated_model
      cat("✓ Updated:", model_name, "to", new_name, "\n")
      
    }, error = function(e) {
      cat("⚠️ Failed to update", model_name, ":", e$message, "\n")
    })
  }
  
  return(updated_list)
}
brm2_models_DHW <- update_model_brm2(brm_models_DHW)

# Look at conditional effects to see if priors are too wide
list_conditional_effects(brm2_models_DHW)

list_SUYR_prior_and_posterior <- function(model_list) {
  for (model_name in names(model_list)) {
    cat("Plotting:", model_name, "\n")
    
    model <- model_list[[model_name]]
    
    # Try plotting and skip if it fails
    tryCatch({
      plot <- model |> SUYR_prior_and_posterior()
      print(plot)
      # Pause and wait for user input
      readline(prompt = "Press [Enter] to see the next plot...")
      
    }, error = function(e) {
      cat("⚠️ Could not plot", model_name, ":", e$message, "\n")
    })
  }
}
list_SUYR_prior_and_posterior(brm2_models_DHW)

# Model diagnostics
generate_and_store_diagnostics <- function(model_list) {
  plot_list <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    fit <- model$fit
    
    cat("\n========== Diagnostics for:", model_name, "==========\n")
    
    # Named functions and suffixes
    plot_functions <- list(
      trace   = function() stan_trace(fit),
      ac      = function() stan_ac(fit),
      rhat    = function() stan_rhat(fit),
      ess     = function() stan_ess(fit),
      dens    = function() stan_dens(fit, separate_chains = TRUE),
      ppcheck = function() pp_check(model, type = 'dens_overlay', ndraws = 100)
    )
    
    for (suffix in names(plot_functions)) {
      plot_label <- paste0(model_name, "_", suffix)
      cat("\n--", plot_label, "--\n")
      readline(prompt = "Press [Enter] to show and store this plot...")
      
      tryCatch({
        p <- plot_functions[[suffix]]()
        print(p)
        plot_list[[plot_label]] <- p
      }, error = function(e) {
        cat("⚠️ Error in", plot_label, ":", e$message, "\n")
      })
    }
  }
  
  return(plot_list)
}
generate_and_store_diagnostics(brm2_models_DHW)
diagnostics_list_DHW <- generate_and_store_diagnostics(brm2_models_DHW)

list_dharma_diagnostics <- function(model_list) {
  dharma_plot_list <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    cat("\n========== DHARMa Diagnostics for:", model_name, "==========\n")
    
    # Create DHARMa residuals object once per model
    tryCatch({
      dharma_res <- make_brms_dharma_res(model, integerResponse = FALSE)
    }, error = function(e) {
      cat("⚠️ Failed to make DHARMa residuals for", model_name, ":", e$message, "\n")
      next
    })
    
    # Define diagnostics and their corresponding functions
    dharma_functions <- list(
      Uniformity = function() testUniformity(dharma_res, plot = TRUE),
      Residuals  = function() plotResiduals(dharma_res),
      Dispersion = function() testDispersion(dharma_res, plot = TRUE)
    )
    
    # Generate and optionally store plots
    for (suffix in names(dharma_functions)) {
      plot_label <- paste0(model_name, "_", suffix)
      cat("\n--", plot_label, "--\n")
      readline(prompt = "Press [Enter] to show and store this plot...")
      
      tryCatch({
        result <- dharma_functions[[suffix]]()
        
        # For Residuals, we get a ggplot object back
        # For Uniformity and Dispersion, we just show the base R plot (invisible)
        if (inherits(result, "ggplot")) {
          print(result)
          dharma_plot_list[[plot_label]] <- result
        } else {
          dharma_plot_list[[plot_label]] <- recordPlot()
        }
      }, error = function(e) {
        cat("⚠️ Error in", plot_label, ":", e$message, "\n")
      })
    }
  }
  
  return(dharma_plot_list)
}
DHARMa_diagnostics_DHW <- list_dharma_diagnostics(brm2_models_DHW)

# If wanting to load the models later: brm2_models_DHW <- readRDS('ah_surv.brm2_list_DHW')

#### ED50 & ED25 Calculation, Ranking and Visualization ####

# Get the predicted values for the response
generate_mean_DHW_table <- function(model_list) {
  all_means <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet from model name (assumes it's between "ak_Genet" and "_")
    genet <- stringr::str_extract(model_name, "(?<=ah_Genet)[^_]+")
    
    cat("Processing model:", model_name, "for Genet:", genet, "\n")
    
    tryCatch({
      # Common DHW vector
      DHWs <- c(0.00,0.16,0.39,0.69,1.07,1.51,2.85,3.30,3.74,4.19,4.63,5.08,5.53,
                5.97,6.42,6.86,7.31,7.75,8.20,8.65,9.09,9.54,9.98,10.43)
      
      # Only keep mean estimates
      mean_DHW <- model |>
        emmeans(~DHW, at = list(DHW = DHWs), type = "response") |>
        summary() |> 
        as_tibble() |> 
        mutate(Genet = genet)
      
      all_means[[model_name]] <- mean_DHW
      
    }, error = function(e) {
      cat("⚠️ Error in", model_name, ":", e$message, "\n")
    })
  }
  
  # Combine into one dataframe
  ah_mean_DHW <- bind_rows(all_means)
  return(ah_mean_DHW)
}
ah_mean_DHW <- generate_mean_DHW_table(brm2_models_DHW)

# Add Mortality probability to the table
ah_mean_DHW <- ah_mean_DHW |>
  mutate(PercentMort = (1 - prob) * 100,
         LowerPectMort = (1 - lower.HPD) * 100,
         UpperPerctMort = (1 - upper.HPD) * 100)

# Only keep Genets that are common to ALL experiments and traits
common_genets <- ah_mean_DHW |> 
  filter(!(Genet %in% c("14","275","329","uk2","335","1839","1840","301","320","305"))) |> 
  pull(Genet)
genets_amount <- n_distinct(common_genets)
ah_mean_DHW <- ah_mean_DHW |> 
  filter(Genet %in% common_genets) 

## ED50s ##

# Calculate the full posteriors of each ED50, including weighted rank and ranks
generate_ED50_draws_only <- function(model_list) {
  all_ED50_draws <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet name
    genet <- stringr::str_extract(model_name, "(?<=ah_Genet)[^_]+")
    cat("Processing model:", model_name, "for Genet:", genet, "\n")
    
    tryCatch({
      ed50_draws <- model |>
        as_draws_df(variable = "^b.*", regex = TRUE) |>
        mutate(
          ED50 = -1 * b_Intercept / b_DHW,
          Genet = genet,
          .draw = row_number()
        ) |> 
        dplyr::select(.draw, Genet, ED50)
      
      all_ED50_draws[[model_name]] <- ed50_draws
      
    }, error = function(e) {
      cat("⚠️ Error in", model_name, ":", e$message, "\n")
    })
  }
  
  # Combine all draws
  ah_surv.medians_ED50_draws <- bind_rows(all_ED50_draws)
  
  # ⚠️ Clean up Genet names
  ah_surv.medians_ED50_draws <- ah_surv.medians_ED50_draws |>
    mutate(Genet = trimws(as.character(Genet)))  # trim leading/trailing spaces
  
  # Print unique genet values to check what's actually there
  cat("Genets found:\n")
  print(unique(ah_surv.medians_ED50_draws$Genet))
  
  # Genets to exclude (make sure they're all characters)
  excluded_genets <- c("14","329")
  
  # Now apply filter
  ah_surv.medians_ED50_draws <- ah_surv.medians_ED50_draws |>
    filter(!(Genet %in% excluded_genets) & !is.na(Genet))
  
  cat("Remaining Genets:\n")
  print(unique(ah_surv.medians_ED50_draws$Genet))
  
  # Compute ranks
  ah_surv.medians_ED50_draws <- ah_surv.medians_ED50_draws |> 
    group_by(.draw) |>
    mutate(
      Rank = rank(desc(ED50)),
      Weighted_rank = (max(ED50) - ED50) * 56 / (max(ED50) - min(ED50)) + 1
    ) |>
    ungroup()
  
  return(ah_surv.medians_ED50_draws)
}
ah_surv.medians_ED50_draws <- generate_ED50_draws_only(brm2_models_DHW) |> 
  filter(Genet %in% common_genets) 

# Get the population ED50's mean from posteriors
ah_PopMean_ED50 <- ah_surv.medians_ED50_draws |> 
  summarize(pop.mean = mean(ED50),
            pop.lower = mean(ED50) - sd(ED50),
            pop.upper = mean(ED50) + sd(ED50)) |> 
  mutate(Species = 'Acropora hyacinthus')

# Calculate the median for ED50s, Weighted _ranks, and Ranks draws
ah_surv.medians_ED50_DHW <- ah_surv.medians_ED50_draws |> 
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
generate_ED25_draws_only <- function(model_list) {
  all_ED25_draws <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet name
    genet <- stringr::str_extract(model_name, "(?<=ah_Genet)[^_]+")
    cat("Processing model:", model_name, "for Genet:", genet, "\n")
    
    tryCatch({
      ED25_draws <- model |>
        as_draws_df(variable = "^b.*", regex = TRUE) |>
        mutate(
          ED25 = (qlogis(0.75) - b_Intercept) / b_DHW,
          Genet = genet,
          .draw = row_number()
        ) |> 
        dplyr::select(.draw, Genet, ED25)
      
      all_ED25_draws[[model_name]] <- ED25_draws
      
    }, error = function(e) {
      cat("⚠️ Error in", model_name, ":", e$message, "\n")
    })
  }
  
  # Combine all draws
  ah_surv.medians_ED25_draws <- bind_rows(all_ED25_draws)
  
  # ⚠️ Clean up Genet names
  ah_surv.medians_ED25_draws <- ah_surv.medians_ED25_draws |>
    mutate(Genet = trimws(as.character(Genet)))  # trim leading/trailing spaces
  
  # Print unique genet values to check what's actually there
  cat("Genets found:\n")
  print(unique(ah_surv.medians_ED25_draws$Genet))
  
  # Genets to exclude (make sure they're all characters)
  excluded_genets <- c("14","329")
  
  # Now apply filter
  ah_surv.medians_ED25_draws <- ah_surv.medians_ED25_draws |>
    filter(!(Genet %in% excluded_genets) & !is.na(Genet))
  
  cat("Remaining Genets:\n")
  print(unique(ah_surv.medians_ED25_draws$Genet))
  
  # Compute ranks
  ah_surv.medians_ED25_draws <- ah_surv.medians_ED25_draws |> 
    group_by(.draw) |>
    mutate(
      Rank = rank(desc(ED25)),
      Weighted_rank = (max(ED25) - ED25) * 56 / (max(ED25) - min(ED25)) + 1
    ) |>
    ungroup()
  
  return(ah_surv.medians_ED25_draws)
}
ah_surv.medians_ED25_draws <- generate_ED25_draws_only(brm2_models_DHW) |> 
  filter(Genet %in% common_genets)

# Get the population ED25's mean from posteriors
ah_PopMean_ED25 <- ah_surv.medians_ED25_draws |> 
  summarize(pop.mean = median(ED25),
            pop.lower = mean(ED25) - sd(ED25),
            pop.upper = mean(ED25) + sd(ED25)) |> 
  mutate(Species = 'Acropora hyacinthus')

# Calculate the median for ED25s, Weighted _ranks, and Ranks draws
ah_surv.medians_ED25_DHW <- ah_surv.medians_ED25_draws |> 
  group_by(Genet) |> 
  summarise(ED25_median = median(ED25, na.rm=TRUE),
            ED25_lower = HDInterval::hdi(ED25)[1],
            ED25_upper = HDInterval::hdi(ED25)[2]) |> 
  arrange(desc(ED25_median)) |> 
  mutate(ED25 = ED25_median,
         Weighted_rank = (ED25[1]-ED25)*genets_amount/(max(ED25)-min(ED25))+1, 
         Rank = rank(desc(ED25)))

# Helps arrange data for final plot
ah_surv.medians_ED50_DHW <- ah_surv.medians_ED50_DHW |> 
  mutate(Genet = factor(Genet, levels= ah_surv.medians_ED50_DHW |> 
                          arrange(Weighted_rank) |>
                          pull(Genet)))
ah_mean_DHW <- ah_mean_DHW |> 
  mutate(Genet = factor(Genet, levels= ah_surv.medians_ED50_DHW |> pull(Genet))) 

ah_surv.medians_ED25_DHW <- ah_surv.medians_ED25_DHW |> 
  mutate(Genet = factor(Genet, levels= ah_surv.medians_ED50_DHW |> pull(Genet))) 

### Plots ###

# Plot medians of ED50s
ah_surv_ED50s_DHW <- ah_surv.medians_ED50_DHW|>
  mutate(Genet = factor(Genet, levels=ah_surv.medians_ED50_DHW |>  arrange(ED50) |> pull(Genet))) |> 
  ggplot(aes(x = ED50, y = Genet)) +
  geom_pointrange(aes(xmin = ED50_lower, xmax = ED50_upper)) +
  labs(x = "Effective Dose Response at 50% (DHW)") +
  theme_light() 

# Plot medians of ED25s
ah_surv_ED25s_DHW <- ah_surv.medians_ED25_DHW|>
  mutate(Genet = factor(Genet, levels=ah_surv.medians_ED25_DHW |>  arrange(ED25) |> pull(Genet))) |> 
  ggplot(aes(x = ED25, y = Genet)) +
  geom_pointrange(aes(xmin = ED25_lower, xmax = ED25_upper)) +
  labs(x = "Effective Dose Response at 50% (DHW)") +
  theme_light()  

# Plot medians combined
ah_surv_EDs_medians_DHW <- ah_surv_ED50s_DHW + ah_surv_ED25s_DHW + 
  plot_annotation(tag_levels = 'A')

#### FINAL PLOTS ####

# ED's per Genet
ah_surv_EDs_DHW_perGenet <- ah_mean_DHW |> 
  mutate(Genet = factor(Genet, levels= ah_surv.medians_ED50_DHW |> arrange(Weighted_rank) |>  pull(Genet))) |> 
  filter(Genet != 'NA') |> 
  ggplot(aes(y = prob, x = DHW)) +
  geom_line(linewidth=0.75, lineend = "round", color = "#696969") +
  labs(y = "Survival probability (%)", x ="DHW") +
  geom_vline(data = ah_surv.medians_ED50_DHW, aes(xintercept=ED50, color = "ED50"),linewidth=0.75) +
  geom_vline(data = ah_surv.medians_ED50_DHW, aes(xintercept=ED50_lower, color = "ED50"), alpha = 0.3) +
  geom_vline(data = ah_surv.medians_ED50_DHW, aes(xintercept=ED50_upper, color = "ED50"), alpha = 0.3) +
  geom_rect(data = ah_surv.medians_ED50_DHW, 
            aes(xmin = ED50_lower, xmax = ED50_upper, ymin = -Inf, ymax = Inf, fill="ED50"),
            alpha = 0.2, inherit.aes = FALSE)+
  geom_vline(data = ah_surv.medians_ED25_DHW, aes(xintercept=ED25, color = "ED25"),linewidth=0.75) +
  geom_vline(data = ah_surv.medians_ED25_DHW, aes(xintercept=ED25_lower, color = "ED25"), alpha = 0.3) +
  geom_vline(data = ah_surv.medians_ED25_DHW, aes(xintercept=ED25_upper, color = "ED25"), alpha = 0.3) +
  geom_rect(data = ah_surv.medians_ED25_DHW, 
            aes(xmin = ED25_lower, xmax = ED25_upper, ymin = -Inf, ymax = Inf, fill="ED25"),
            alpha = 0.2, inherit.aes = FALSE)+
  scale_color_manual(name = "Effective dose response", values = c("ED25" = "#f0a009", "ED50" = "#e74c3c")) +
  scale_fill_manual(name = "Effective dose response", values = c("ED25" = "#f0a009", "ED50" = "#e74c3c"))+
  facet_wrap(~Genet)+
  theme_light()+
  theme(legend.position = "bottom",
        axis.title = element_text(size=12),
        strip.text = element_text(size=12, colour='#696969', margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(linewidth = 1, fill="white"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

# Percent Mortality & Population's mean ED50
ah_mortality_pop <- ah_mean_DHW |> 
  filter(Genet != 'NA') |> 
  ggplot(aes(y = PercentMort, x = DHW, group = Genet)) +
  geom_line(color = "#3db0bf50") +
  labs(y = "Mortality (%)", x ="Degree Heating Weeks") +
  geom_vline(data = ah_PopMean_ED25, aes(xintercept=pop.mean), color = "#0f4d5b", linewidth=0.71) +
  geom_vline(data = ah_PopMean_ED25, aes(xintercept=pop.upper), color = "#007281", linewidth=0.7, alpha = 0.3) +
  geom_vline(data = ah_PopMean_ED25, aes(xintercept=pop.lower), color = "#007281", linewidth=0.7, alpha = 0.3) + 
  geom_rect(data = ah_PopMean_ED25, 
            aes(xmin = pop.lower, xmax = pop.upper, ymin = -Inf, ymax = Inf),
            fill = "#3db0bf", alpha = 0.25, inherit.aes = FALSE)+
  theme_classic()+
  facet_wrap(~Species)+
  theme(strip.background = element_blank()) +
  xlim(0, 20) 