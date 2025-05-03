# title: "Survival Analysis Acropora kenti"
# description : Code for analysis of survival
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
easystats::easystats_update()
library(patchwork)
library(tinytable)
library(survival)  
library(corrplot)  
library(summarytools)  
library(cowplot)
library(RColorBrewer)
library(survminer)
library(readxl)
library(writexl)

setwd('D:/AIMS/Experiment/Analysis/Acropora kenti')

# Run the helperFunctions.R
source('Code/helperFunctions.R')

# Read in the data
ak_surv <- read_excel("Data/Akenti_survival.xlsx")

# Turn into numeric and filter to treatment
ak_surv <- ak_surv |> mutate(DHW =as.numeric(as.character(DHW)),
                             Day =as.numeric(as.character(DayMeasured))) |> 
  filter(Temperature == '31.5')

# Let's make a table per Genet
GenetTables <- function(data) {
  unique_genets <- unique(data$Genet)
  GenetTables_list <- list()
  
  for (g in unique_genets) {
    subset_data <- subset(data, Genet == g)
    var_name <- paste0("ak_Genet", g, "_surv")
    assign(var_name, subset_data, envir = .GlobalEnv)
  }
  
  GenetTables_list <- setNames(
    lapply(unique_genets, function(g) {
      subset(data, Genet == g)
    }),
    paste0("ak_Genet", unique_genets, "_surv")
  )
  
  return(GenetTables_list)
} # Function makes 1 table per genets then a list of those tables
ak_GenetTables <- GenetTables(ak_surv)

#### Model 2: Treats Day as continuous and Full Interaction ####

# Make models ####

ak_survday.form <- bf(Survival|trials(1) ~ Day + (1|Tank),
                      family = binomial(link = logit))

priors <- prior(student_t(3, 0, 1), class = 'Intercept') +
          prior(normal(0,1), class = 'b') +
          prior(student_t(3,0,1), class = 'sd')

fit_brm <- function(table_list, formula, priors) {
  model_list <- list()
  
  for (name in names(table_list)) {
    message("Fitting model for: ", name)
    model <- brm(
      formula = formula,
      data = table_list[[name]],
      prior = priors,
      sample_prior = 'only',
      iter = 5000,
      warmup = 1000,
      chains = 3,
      cores = 3,
      thin = 5,
      control = list(adapt_delta = 0.99, max_treedepth = 20),
      refresh = 0,
      backend = "cmdstanr"
    )
    
    model_name <- paste0(name, "day.brm")
    model_list[[model_name]] <- model
  }
  
  return(model_list)
}
brm_models_Day <- fit_brm(ak_GenetTables, ak_survday.form, priors)

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
list_conditional_effects(brm_models_Day)

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
brm2_models_Day <- update_model_brm2(brm_models_Day)

list_conditional_effects(brm2_models_Day)
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
list_SUYR_prior_and_posterior(brm2_models_Day)

  # Adjust Genets that have to ####
#   1861, 1864, 1883, 1884, 1889, 1990, 2008, 2017, 2027, 2030
ak_Genets_toredo <- mget(c("ak_Genet1861_surv","ak_Genet1864_surv", "ak_Genet1883_surv", "ak_Genet1884_surv", 
                         "ak_Genet1889_surv","ak_Genet2027_surv", "ak_Genet2030_surv"))
fit_brm_2.0 <- function(table_list, formula, priors) {
  model_list <- list()
  
  for (name in names(table_list)) {
    message("Fitting model for: ", name)
    model <- brm(
      formula = formula,
      data = table_list[[name]],
      prior = priors,
      sample_prior = 'only',
      iter = 9000,
      warmup = 1000,
      chains = 3,
      cores = 3,
      thin = 15,
      control = list(adapt_delta = 0.99, max_treedepth = 25),
      refresh = 0,
      backend = "cmdstanr"
    )
    
    model_name <- paste0(name, "day.brm")
    model_list[[model_name]] <- model
  }
  
  return(model_list)
}
brm_models_Day_redo <- fit_brm_2.0(ak_Genets_toredo, ak_survday.form, priors)
brm2_models_Day_redo <- update_model_brm2(brm_models_Day_redo)

# Genet 1990
ak_Genet1990_surv |> 
  summarise(mean(Survival, na.rm = TRUE),
            2.5*sd(Survival, na.rm = TRUE),
            2.5*sd(Survival)/sd(Day))

priors <- prior(student_t(3, 0.7, 1), class = 'Intercept') +
  prior(normal(0,0.1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

ak_Genet1990_surv.brm.base <- brm(ak_survday.form, 
                     data = ak_Genet1990_surv,
                     prior = priors,
                     sample_prior = 'only',
                     iter = 10000,
                     warmup = 1000,
                     chains = 3, cores = 3,
                     thin = 15,
                     refresh = 0,
                     seed = 123,
                     control = list(adapt_delta = 0.99, max_treedepth = 20),
                     backend = 'cmdstan')

ak_Genet1990_surv.brm2 <- update(ak_Genet1990_surv.brm.base,
                         sample_prior = 'yes',
                         refresh = 0,
                         cores = 3, 
                         seed = 123) 

# Genet 2008
ak_Genet2008_surv |> 
  summarise(mean(Survival, na.rm = TRUE),
            2.5*sd(Survival, na.rm = TRUE),
            2.5*sd(Survival)/sd(Day))

priors <- prior(student_t(3, 0.28, 1.1), class = 'Intercept') +
  prior(normal(0,0.1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

ak_Genet2008_surv.brm.base <- brm(ak_survday.form, 
                                  data = ak_Genet2008_surv,
                                  prior = priors,
                                  sample_prior = 'only',
                                  iter = 10000,
                                  warmup = 1000,
                                  chains = 3, cores = 3,
                                  thin = 15,
                                  refresh = 0,
                                  seed = 123,
                                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                                  backend = 'cmdstan')

ak_Genet2008_surv.brm2 <- update(ak_Genet2008_surv.brm.base,
                                 sample_prior = 'yes',
                                 refresh = 0,
                                 cores = 3, 
                                 seed = 123) 

# Genet 2017
ak_Genet2017_surv |> 
  summarise(mean(Survival, na.rm = TRUE),
            2.5*sd(Survival, na.rm = TRUE),
            2.5*sd(Survival)/sd(Day))

priors <- prior(student_t(3, 0.7, 1), class = 'Intercept') +
  prior(normal(0,0.1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

ak_Genet2017_surv.brm.base <- brm(ak_survday.form, 
                                  data = ak_Genet2017_surv,
                                  prior = priors,
                                  sample_prior = 'only',
                                  iter = 10000,
                                  warmup = 1000,
                                  chains = 3, cores = 3,
                                  thin = 15,
                                  refresh = 0,
                                  seed = 123,
                                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                                  backend = 'cmdstan')

ak_Genet2017_surv.brm2 <- update(ak_Genet2017_surv.brm.base,
                                 sample_prior = 'yes',
                                 refresh = 0,
                                 cores = 3, 
                                 seed = 123) 


# Add and replace better fitted genet to the original model list
brm2_models_Day[['ak_Genet1861_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet1861_survDay.brm2']]
brm2_models_Day[['ak_Genet1864_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet1864_survDay.brm2']]
brm2_models_Day[['ak_Genet1883_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet1883_survDay.brm2']]
brm2_models_Day[['ak_Genet1884_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet1884_survDay.brm2']]
brm2_models_Day[['ak_Genet1889_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet1889_survDay.brm2']]
brm2_models_Day[['ak_Genet2027_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet2027_survDay.brm2']]
brm2_models_Day[['ak_Genet2030_survday.brm2']] <- brm2_models_Day_redo[['ak_Genet2030_survDay.brm2']]
brm2_models_Day[['ak_Genet1990_survday.brm2']] <- ak_Genet1990_survDay.brm2
brm2_models_Day[['ak_Genet2008_survday.brm2']] <- ak_Genet2008_survDay.brm2
brm2_models_Day[['ak_Genet2017_survDay.brm2']] <- ak_Genet2017_survDay.brm2

# Make diagnostics ####

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
generate_and_store_diagnostics(brm2_models_Day)
diagnostics_list_Day <- generate_and_store_diagnostics(brm2_models_Day)

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
DHARMa_diagnostics_Day <- list_dharma_diagnostics(brm2_models_Day)

save_plots_png_tiff_pdf <- function(plot_list, output_dir = "plots", width = 8, height = 6, dpi = 300) {
  formats <- c("pdf", "png", "tiff")
  
  # Create folders if they don't exist
  for (fmt in formats) {
    dir.create(file.path(output_dir, fmt), recursive = TRUE, showWarnings = FALSE)
  }
  
  for (plot_name in names(plot_list)) {
    cat("Saving:", plot_name, "\n")
    plot_obj <- plot_list[[plot_name]]
    
    for (fmt in formats) {
      file_path <- file.path(output_dir, fmt, paste0(plot_name, ".", fmt))
      
      tryCatch({
        if (inherits(plot_obj, "ggplot")) {
          ggsave(filename = file_path, plot = plot_obj, width = width, height = height, dpi = dpi)
        } else if (inherits(plot_obj, "recordedplot")) {
          if (fmt == "pdf") {
            pdf(file_path, width = width, height = height)
          } else if (fmt == "png") {
            png(file_path, width = width * dpi, height = height * dpi, res = dpi)
          } else if (fmt == "tiff") {
            tiff(file_path, width = width * dpi, height = height * dpi, res = dpi)
          }
          replayPlot(plot_obj)
          dev.off()
        } else {
          warning("❗ Unknown plot type for", plot_name)
        }
      }, error = function(e) {
        cat("⚠️ Error saving", plot_name, "in", fmt, "format:", e$message, "\n")
      })
    }
  }
}
save_plots_png_tiff_pdf(diagnostics_list_Day, output_dir = "Output/Figure/Survival Diagnostics")
save_plots_png_tiff_pdf(DHARMa_diagnostics_Day, output_dir = "Output/Figure/Survival Diagnostics")

# If satisfied with all models, let's save them all in a list:
saveRDS(brm2_models_Day, file = "Output/Model/aken_brm2_list_Day.rds")
  # If wanting to load them later: brm2_models_Day <- readRDS("Output/Model/ahya_brm2_list.rds")

# Produce EDs and Figure ####

# Get the predicted values for the response
generate_mean_Day_table <- function(model_list) {
  all_means <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet from model name (assumes it's between "ak_Genet" and "_")
    genet <- stringr::str_extract(model_name, "(?<=ak_Genet)[^_]+")
    
    cat("Processing model:", model_name, "for Genet:", genet, "\n")
    
    tryCatch({
      # Common Day vector
      Days <- c(1,4,5,6,7,8,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
                25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
                43,44,45,46,47)
      
      # Only keep mean estimates
      mean_Day <- model |>
        emmeans(~Day, at = list(Day = Days), type = "response") |>
        summary() |> 
        as_tibble() |> 
        mutate(Genet = genet)
      
      all_means[[model_name]] <- mean_Day
      
    }, error = function(e) {
      cat("⚠️ Error in", model_name, ":", e$message, "\n")
    })
  }
  
  # Combine into one dataframe
  ak_mean_Day <- bind_rows(all_means)
  return(ak_mean_Day)
}
ak_mean_Day <- generate_mean_Day_table(brm2_models_Day)

# Calculate the full posteriors of each ED50, including weighted rank and ranks
generate_ED50_draws_only <- function(model_list) {
  all_ED50_draws <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet name
    genet <- stringr::str_extract(model_name, "(?<=ak_Genet)[^_]+")
    cat("Processing model:", model_name, "for Genet:", genet, "\n")
    
    tryCatch({
      ed50_draws <- model |>
        as_draws_df(variable = "^b.*", regex = TRUE) |>
        mutate(
          ED50 = -1 * b_Intercept / b_Day,
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
  ak_surv.medians_ED50_draws <- bind_rows(all_ED50_draws)
  
  # ⚠️ Clean up Genet names
  ak_surv.medians_ED50_draws <- ak_surv.medians_ED50_draws |>
    mutate(Genet = trimws(as.character(Genet)))  # trim leading/trailing spaces
  
  # Print unique genet values to check what's actually there
  cat("Genets found:\n")
  print(unique(ak_surv.medians_ED50_draws$Genet))
  
  # Genets to exclude (make sure they're all characters)
  excluded_genets <- c("1868","1869","1899","1961","1968","1972q","2008","2035")
  
  # Now apply filter
  ak_surv.medians_ED50_draws <- ak_surv.medians_ED50_draws |>
    filter(!(Genet %in% excluded_genets) & !is.na(Genet))
  
  cat("Remaining Genets:\n")
  print(unique(ak_surv.medians_ED50_draws$Genet))
  
  # Compute ranks
  ak_surv.medians_ED50_draws <- ak_surv.medians_ED50_draws |> 
    group_by(.draw) |>
    mutate(
      Rank = rank(desc(ED50)),
      Weighted_rank = (max(ED50) - ED50) * 56 / (max(ED50) - min(ED50)) + 1
    ) |>
    ungroup()
  
  return(ak_surv.medians_ED50_draws)
}
ak_surv.medians_ED50_draws <- generate_ED50_draws_only(brm2_models_Day)

# Calculate the median for ED50s, Weighted _ranks, and Ranks draws
ak_surv.medians_ED50_Day <- ak_surv.medians_ED50_draws |> 
  group_by(Genet) |> 
  summarise(ED50_median = median(ED50, na.rm=TRUE),
            ED50_lower = HDInterval::hdi(ED50)[1],
            ED50_upper = HDInterval::hdi(ED50)[2]) |> 
  arrange(desc(ED50_median)) |> 
  mutate(ED50 = ED50_median,
         Weighted_rank = (ED50[1]-ED50)*56/(max(ED50)-min(ED50))+1, 
         Rank = rank(desc(ED50)))

### Plots

# Plot medians of ED50s
ak_surv_ED50_Day <- ak_surv.medians_ED50_day|>
  mutate(Genet = factor(Genet, levels=ak_surv.medians_ED50_day |>  arrange(ED50) |> pull(Genet))) |> 
  ggplot(aes(x = ED50, y = Genet)) +
  geom_pointrange(aes(xmin = ED50_lower, xmax = ED50_upper)) +
  labs(x = "Effective Dose Response at 50% (Day)") +
  theme_light()

### FINAL PLOT ### 

ak_surv_EDs_Day_perGenet <- ak_mean_day |> 
  mutate(Genet = factor(Genet, levels=ak_surv.medians_ED50_day |> pull(Genet))) |>
  filter(Genet != 'NA') |> 
  ggplot(aes(y = prob, x = Day)) +
  geom_line(linewidth=0.75, lineend = "round", color = "#696969") +
  labs(y = "Survival probability (%)", x ="Day") +
  geom_vline(data = ak_surv.medians_ED50_day, aes(xintercept=ED50, color = "ED50"),linewidth=0.75) +
  geom_vline(data = ak_surv.medians_ED50_day, aes(xintercept=ED50_lower, color = "ED50"), alpha = 0.3) +
  geom_vline(data = ak_surv.medians_ED50_day, aes(xintercept=ED50_upper, color = "ED50"), alpha = 0.3) +
  geom_rect(data = ak_surv.medians_ED50_day, 
            aes(xmin = ED50_lower, xmax = ED50_upper, ymin = -Inf, ymax = Inf, fill="ED50"),
            alpha = 0.2, inherit.aes = FALSE)+
  scale_color_manual(name = "Effective dose response", values = c("ED50" = "#e74c3c")) +
  scale_fill_manual(name = "Effective dose response", values = c("ED50" = "#e74c3c"))+
  facet_wrap(~Genet)+
  theme_light()+
  theme(legend.position = "bottom",
        axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        strip.text = element_text(size=15, colour='#696969', margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(linewidth = 1, fill="white"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

#### Model: Treats DHW as continuous and Full Interaction ####
# Make models ####
ak_survDHW.form <- bf(Survival|trials(1) ~ DHW + (1|Tank),
                      family = binomial(link = logit))

priors <- prior(student_t(3, 0, 1), class = 'Intercept') +
  prior(normal(0,1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

fit_brm <- function(table_list, formula, priors) {
  model_list <- list()
  
  for (name in names(table_list)) {
    message("Fitting model for: ", name)
    model <- brm(
      formula = formula,
      data = table_list[[name]],
      prior = priors,
      sample_prior = 'only',
      iter = 5000,
      warmup = 1000,
      chains = 3,
      cores = 3,
      thin = 5,
      control = list(adapt_delta = 0.99, max_treedepth = 20),
      refresh = 0,
      backend = "cmdstanr"
    )
    
    model_name <- paste0(name, "DHW.brm")
    model_list[[model_name]] <- model
  }
  
  return(model_list)
}
brm_models_DHW <- fit_brm(ak_GenetTables, ak_survDHW.form, priors)

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

# Adjust Genets that have to ####
#   1861, 1864, 1883, 1884, 1889, 1990, 2008, 2017, 2027, 2030
ak_Genets_toredo <- mget(c("ak_Genet1861_surv","ak_Genet1864_surv", "ak_Genet1883_surv", "ak_Genet1884_surv", 
                           "ak_Genet1889_surv","ak_Genet2027_surv", "ak_Genet2030_surv"))
fit_brm_2.0 <- function(table_list, formula, priors) {
  model_list <- list()
  
  for (name in names(table_list)) {
    message("Fitting model for: ", name)
    model <- brm(
      formula = formula,
      data = table_list[[name]],
      prior = priors,
      sample_prior = 'only',
      iter = 9000,
      warmup = 1000,
      chains = 3,
      cores = 3,
      thin = 15,
      control = list(adapt_delta = 0.99, max_treedepth = 25),
      refresh = 0,
      backend = "cmdstanr"
    )
    
    model_name <- paste0(name, "DHW.brm")
    model_list[[model_name]] <- model
  }
  
  return(model_list)
}
brm_models_DHW_redo <- fit_brm_2.0(ak_Genets_toredo, ak_survDHW.form, priors)
brm2_models_DHW_redo <- update_model_brm2(brm_models_DHW_redo)

# Genet 1990
ak_Genet1990_surv |> 
  summarise(mean(Survival, na.rm = TRUE),
            2.5*sd(Survival, na.rm = TRUE),
            2.5*sd(Survival)/sd(DHW))

priors <- prior(student_t(3, 0.7, 1), class = 'Intercept') +
  prior(normal(0,0.1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

ak_Genet1990_surv.brm.base <- brm(ak_survDHW.form, 
                                  data = ak_Genet1990_surv,
                                  prior = priors,
                                  sample_prior = 'only',
                                  iter = 10000,
                                  warmup = 1000,
                                  chains = 3, cores = 3,
                                  thin = 15,
                                  refresh = 0,
                                  seed = 123,
                                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                                  backend = 'cmdstan')

ak_Genet1990_surv.brm2 <- update(ak_Genet1990_surv.brm.base,
                                 sample_prior = 'yes',
                                 refresh = 0,
                                 cores = 3, 
                                 seed = 123) 

# Genet 2008
ak_Genet2008_surv |> 
  summarise(mean(Survival, na.rm = TRUE),
            2.5*sd(Survival, na.rm = TRUE),
            2.5*sd(Survival)/sd(DHW))

priors <- prior(student_t(3, 0.28, 1.1), class = 'Intercept') +
  prior(normal(0,0.1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

ak_Genet2008_surv.brm.base <- brm(ak_survDHW.form, 
                                  data = ak_Genet2008_surv,
                                  prior = priors,
                                  sample_prior = 'only',
                                  iter = 10000,
                                  warmup = 1000,
                                  chains = 3, cores = 3,
                                  thin = 15,
                                  refresh = 0,
                                  seed = 123,
                                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                                  backend = 'cmdstan')

ak_Genet2008_surv.brm2 <- update(ak_Genet2008_surv.brm.base,
                                 sample_prior = 'yes',
                                 refresh = 0,
                                 cores = 3, 
                                 seed = 123) 

# Genet 2017
ak_Genet2017_surv |> 
  summarise(mean(Survival, na.rm = TRUE),
            2.5*sd(Survival, na.rm = TRUE),
            2.5*sd(Survival)/sd(DHW))

priors <- prior(student_t(3, 0.7, 1), class = 'Intercept') +
  prior(normal(0,0.1), class = 'b') +
  prior(student_t(3,0,1), class = 'sd')

ak_Genet2017_surv.brm.base <- brm(ak_survDHW.form, 
                                  data = ak_Genet2017_surv,
                                  prior = priors,
                                  sample_prior = 'only',
                                  iter = 10000,
                                  warmup = 1000,
                                  chains = 3, cores = 3,
                                  thin = 15,
                                  refresh = 0,
                                  seed = 123,
                                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                                  backend = 'cmdstan')

ak_Genet2017_surv.brm2 <- update(ak_Genet2017_surv.brm.base,
                                 sample_prior = 'yes',
                                 refresh = 0,
                                 cores = 3, 
                                 seed = 123) 

# Add and replace better fitted genet to the original model list
brm2_models_DHW[['ak_Genet1861_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet1861_survDHW.brm2']]
brm2_models_DHW[['ak_Genet1864_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet1864_survDHW.brm2']]
brm2_models_DHW[['ak_Genet1883_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet1883_survDHW.brm2']]
brm2_models_DHW[['ak_Genet1884_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet1884_survDHW.brm2']]
brm2_models_DHW[['ak_Genet1889_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet1889_survDHW.brm2']]
brm2_models_DHW[['ak_Genet2027_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet2027_survDHW.brm2']]
brm2_models_DHW[['ak_Genet2030_survDHW.brm2']] <- brm2_models_DHW_redo[['ak_Genet2030_survDHW.brm2']]
brm2_models_DHW[['ak_Genet1990_survDHW.brm2']] <- ak_Genet1990_survDHW.brm2
brm2_models_DHW[['ak_Genet2008_survDHW.brm2']] <- ak_Genet2008_survDHW.brm2
brm2_models_DHW[['ak_Genet2017_survDHW.brm2']] <- ak_Genet2017_survDHW.brm2

# Make diagnostics ####

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

save_plots_png_tiff_pdf <- function(plot_list, output_dir = "plots", width = 8, height = 6, dpi = 300) {
  formats <- c("pdf", "png", "tiff")
  
  # Create folders if they don't exist
  for (fmt in formats) {
    dir.create(file.path(output_dir, fmt), recursive = TRUE, showWarnings = FALSE)
  }
  
  for (plot_name in names(plot_list)) {
    cat("Saving:", plot_name, "\n")
    plot_obj <- plot_list[[plot_name]]
    
    for (fmt in formats) {
      file_path <- file.path(output_dir, fmt, paste0(plot_name, ".", fmt))
      
      tryCatch({
        if (inherits(plot_obj, "ggplot")) {
          ggsave(filename = file_path, plot = plot_obj, width = width, height = height, dpi = dpi)
        } else if (inherits(plot_obj, "recordedplot")) {
          if (fmt == "pdf") {
            pdf(file_path, width = width, height = height)
          } else if (fmt == "png") {
            png(file_path, width = width * dpi, height = height * dpi, res = dpi)
          } else if (fmt == "tiff") {
            tiff(file_path, width = width * dpi, height = height * dpi, res = dpi)
          }
          replayPlot(plot_obj)
          dev.off()
        } else {
          warning("❗ Unknown plot type for", plot_name)
        }
      }, error = function(e) {
        cat("⚠️ Error saving", plot_name, "in", fmt, "format:", e$message, "\n")
      })
    }
  }
}
save_plots_png_tiff_pdf(diagnostics_list_DHW, output_dir = "Output/Figure/Survival Diagnostics")
save_plots_png_tiff_pdf(DHARMa_diagnostics_DHW, output_dir = "Output/Figure/Survival Diagnostics")

# If satisfied with all models, let's save them all in a list:
saveRDS(brm2_models_DHW, file = "Output/Model/aken_brm2_list_DHW.rds")
# If wanting to load them later: brm2_models_DHW <- readRDS("Output/Model/ahya_brm2_list.rds")

# Produce EDs and Figure ####

# Get the predicted values for the response
generate_mean_DHW_table <- function(model_list) {
  all_means <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet from model name (assumes it's between "ak_Genet" and "_")
    genet <- stringr::str_extract(model_name, "(?<=ak_Genet)[^_]+")
    
    cat("Processing model:", model_name, "for Genet:", genet, "\n")
    
    tryCatch({
      # Common DHW vector
      DHWs <- c(0.00,0.16,0.39,0.69,1.07,1.51,2.85,3.30,3.74,4.19,4.63,5.08,5.53,
                5.97,6.42,6.86,7.31,7.75,8.20,8.65,9.09,9.54,9.98,10.43,10.87,
                11.32,11.77,12.21,12.66,13.01,13.55,13.99,14.44,14.89,15.33,15.78,
                16.22,16.67,17.11,17.56,18.01,18.45,18.90)
      
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
  ak_mean_DHW <- bind_rows(all_means)
  return(ak_mean_DHW)
}
ak_mean_DHW <- generate_mean_DHW_table(brm2_models_DHW)

# Calculate the full posteriors of each ED50, including weighted rank and ranks
generate_ED50_draws_only <- function(model_list) {
  all_ED50_draws <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    
    # Extract Genet name
    genet <- stringr::str_extract(model_name, "(?<=ak_Genet)[^_]+")
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
  ak_surv.medians_ED50_draws <- bind_rows(all_ED50_draws)
  
  # ⚠️ Clean up Genet names
  ak_surv.medians_ED50_draws <- ak_surv.medians_ED50_draws |>
    mutate(Genet = trimws(as.character(Genet)))  # trim leading/trailing spaces
  
  # Print unique genet values to check what's actually there
  cat("Genets found:\n")
  print(unique(ak_surv.medians_ED50_draws$Genet))
  
  # Genets to exclude (make sure they're all characters)
  excluded_genets <- c("1868","1869","1899","1961","1968","1972q","2008","2035")
  
  # Now apply filter
  ak_surv.medians_ED50_draws <- ak_surv.medians_ED50_draws |>
    filter(!(Genet %in% excluded_genets) & !is.na(Genet))
  
  cat("Remaining Genets:\n")
  print(unique(ak_surv.medians_ED50_draws$Genet))
  
  # Compute ranks
  ak_surv.medians_ED50_draws <- ak_surv.medians_ED50_draws |> 
    group_by(.draw) |>
    mutate(
      Rank = rank(desc(ED50)),
      Weighted_rank = (max(ED50) - ED50) * 56 / (max(ED50) - min(ED50)) + 1
    ) |>
    ungroup()
  
  return(ak_surv.medians_ED50_draws)
}
ak_surv.medians_ED50_draws <- generate_ED50_draws_only(brm2_models_DHW)

# Calculate the median for ED50s, Weighted _ranks, and Ranks draws
ak_surv.medians_ED50_DHW <- ak_surv.medians_ED50_draws |> 
  group_by(Genet) |> 
  summarise(ED50_median = median(ED50, na.rm=TRUE),
            ED50_lower = HDInterval::hdi(ED50)[1],
            ED50_upper = HDInterval::hdi(ED50)[2]) |> 
  arrange(desc(ED50_median)) |> 
  mutate(ED50 = ED50_median,
         Weighted_rank = (ED50[1]-ED50)*56/(max(ED50)-min(ED50))+1, 
         Rank = rank(desc(ED50)))

# Help arrange data for final plot
ak_mean_DHW <- ak_mean_DHW |> 
  mutate(Genet = factor(Genet, levels=ak_surv.medians_ED50_DHW |> pull(Genet)))
ak_surv.medians_ED50_DHW <- ak_surv.medians_ED50_DHW |> 
  mutate(Genet = factor(Genet, levels=ak_surv.medians_ED50_DHW |> pull(Genet)))

### Plots

# Plot medians of ED50s
ak_surv_ED50_DHW <- ak_surv.medians_ED50_DHW|>
  mutate(Genet = factor(Genet, levels=ak_surv.medians_ED50_DHW |>  arrange(ED50) |> pull(Genet))) |> 
  ggplot(aes(x = ED50, y = Genet)) +
  geom_pointrange(aes(xmin = ED50_lower, xmax = ED50_upper)) +
  labs(x = "Effective Dose Response at 50% (DHW)")  

### FINAL PLOT ###

ak_surv_EDs_DHW_perGenet <- ak_mean_DHW |> 
  mutate(Genet = factor(Genet, levels=ak_surv.medians_ED50_DHW |> pull(Genet))) |>
  filter(Genet != 'NA') |> 
  ggplot(aes(y = prob, x = DHW)) +
  geom_line(linewidth=0.75, lineend = "round", color = "#696969") +
  labs(y = "Survival probability (%)", x ="DHW") +
  geom_vline(data = ak_surv.medians_ED50_DHW, aes(xintercept=ED50, color = "ED50"),linewidth=0.75) +
  geom_vline(data = ak_surv.medians_ED50_DHW, aes(xintercept=ED50_lower, color = "ED50"), alpha = 0.3) +
  geom_vline(data = ak_surv.medians_ED50_DHW, aes(xintercept=ED50_upper, color = "ED50"), alpha = 0.3) +
  geom_rect(data = ak_surv.medians_ED50_DHW, 
            aes(xmin = ED50_lower, xmax = ED50_upper, ymin = -Inf, ymax = Inf, fill="ED50"),
            alpha = 0.2, inherit.aes = FALSE)+
  scale_color_manual(name = "Effective dose response", values = c("ED50" = "#e74c3c")) +
  scale_fill_manual(name = "Effective dose response", values = c("ED50" = "#e74c3c"))+
  facet_wrap(~Genet)+
  theme_light()+
  theme(legend.position = "bottom",
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=15, colour='#696969', margin = margin(0.1,0.1,0.1,0.1, "cm")),
        strip.background = element_rect(linewidth = 1, fill="white"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))