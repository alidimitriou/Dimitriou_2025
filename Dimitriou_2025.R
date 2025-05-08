# ---------------------------------------------------------------
## Habitat Use Models
# ---------------------------------------------------------------
# Load Required Libraries
packages <- c(
  "bayesplot", "brms", "cmdstanr", "dplyr", "ggdist", "ggplot2",
  "iNEXT", "kableExtra", "knitr", "lubridate", "patchwork",
  "posterior", "purrr", "readr", "sf", "spdep",
  "tidybayes", "tidyr", "tidyverse"
)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

# Prepare input/output directories
input_dir <- "~/Desktop/MSc_Chapter_1/output4"
output_dir <- "~/Desktop/MSc_Chapter_1/model_output"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Read cleaned dataset
model_data <- read.csv( "~/Desktop/MSc_Chapter_1/clean_model_data_scaled.csv")

# Define function for brms model
options(brms.backend = "cmdstanr", brms.parallel_chains = 4)

run_and_save_brms_model <- function(response_var, species_name, iter_val, warmup_val) {
  model <- brm(
    formula = as.formula(
      paste(response_var, "~",
            "mean_ndvi_scaled + canopy_cover_scaled +",
            "rec_month_100_scaled * dist_to_trail_decay_scaled +",
            "(1 | Deployment.Location.ID) + offset(log_effort)")
    ),
    data = model_data,
    family = negbinomial(),
    prior = c(
      prior(normal(0, 5), class = "b"),
      prior(normal(0, 5), class = "Intercept"),
      prior(cauchy(0, 2), class = "sd")
    ),
    chains = 4,
    seed = 123,
    control = list(adapt_delta = 0.999, max_treedepth = 20),
    save_pars = save_pars(all = TRUE),
    init = 0
  )
  model_path <- file.path(output_dir, paste0(species_name, "_model.rds"))
  saveRDS(model, model_path)
  cat(paste0("Model saved: ", model_path, "\n"))
  
  summary_path <- file.path(output_dir, paste0(species_name, "_summary.txt"))
  sink(summary_path)
  print(summary(model))
  sink()
  cat(paste0("Summary saved: ", summary_path, "\n"))
  
  plot_and_save_effect_sizes(model, species_name)
  
  gc()
  
  return(model)
}

# Run and save each species model
species_models <- list(
  "Snowshoe Hare" = "lepus_detections",
  "Black Bear" = "ursus_detections",
  "Deer" = "odocoileus_detections",
  "Coyote" = "canis_detections",
  "Bobcat" = "lynx_detections",
  "Marten" = "martes_detections",
  "Marmot" = "marmota_detections",
  "Rare Carnivores" = "rare_detections"
)

for (species_name in names(species_models)) {
  response_var <- species_models[[species_name]]
  cat(paste0("\nRunning model for: ", species_name, "...\n"))
  
  if (species_name %in% c("Bobcat", "Rare Carnivores")) { # to improve convergence for more sparsely detected species
    iter_val <- 5000
    warmup_val <- 1000
  } else {
    iter_val <- 6000
    warmup_val <- 2000
  }
  
  run_and_save_brms_model(response_var, species_name, iter_val, warmup_val)
}

# ---------------------------------------------------------------
## Effect size plots
# ---------------------------------------------------------------
# Extract all significant effects across models
prepare_all_significant_effects <- function(model_path, species_name) {
  model <- readRDS(model_path)
  posterior_samples <- as_draws_df(model)
  
  all_vars <- c(
    "b_mean_ndvi_scaled",
    "b_dist_to_trail_decay",
    "b_elevation_scaled",
    "b_canopy_cover_scaled",
    "b_rec_month_100_scaled",
    "b_rec_month_100_scaled:dist_to_trail_decay"
  )
  
  existing_vars <- all_vars[all_vars %in% colnames(posterior_samples)]
  if (length(existing_vars) == 0) return(NULL)
  
  posterior_samples %>%
    select(all_of(existing_vars)) %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Estimate") %>%
    mutate(Species = species_name) %>%
    group_by(Species, Variable) %>%
    summarize(
      median   = median(Estimate),
      lower95  = quantile(Estimate, 0.025),
      upper95  = quantile(Estimate, 0.975),
      lower80  = quantile(Estimate, 0.10),
      upper80  = quantile(Estimate, 0.90),
      .groups = "drop"
    ) %>%
    filter(lower95 > 0 | upper95 < 0)
}

#  Process All Models and Aggregate Significant Effects
model_files <- list.files(input_dir, pattern = "_model\\.rds$", full.names = TRUE)

all_significant_effects <- list()
for (model_path in model_files) {
  species_name <- gsub("_model\\.rds", "", basename(model_path))
  effects <- prepare_all_significant_effects(model_path, species_name)
  if (!is.null(effects)) all_significant_effects[[species_name]] <- effects
}

all_significant_effects_df <- bind_rows(all_significant_effects)
if (nrow(all_significant_effects_df) == 0) stop("No significant variables found across all models.")


# Plot Combined Significant Effects with Axis Lines

var_labels <- c(
  "b_mean_ndvi_scaled" = "NDVI",
  "b_dist_to_trail_decay" = "Distance to Trail",
  "b_elevation_scaled" = "Elevation",
  "b_canopy_cover_scaled" = "Canopy Cover",
  "b_rec_month_100_scaled" = "Recreation Park-Month",
  "b_rec_month_100_scaled:dist_to_trail_decay" = "Recreation × Distance"
)

plot_df <- all_significant_effects_df %>%
  mutate(Variable = factor(var_labels[Variable], levels = rev(unique(var_labels))))

combined_effects_plot <- ggplot(plot_df, aes(y = Variable)) +
  geom_segment(aes(x = lower80, xend = upper80, yend = Variable),
               linewidth = 3, color = "#1f78b4") +
  geom_segment(aes(x = lower95, xend = upper95, yend = Variable),
               linewidth = 1.2, color = "#084594") +
  geom_point(aes(x = median), size = 4.5, color = "#084594") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  facet_wrap(~Species, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8)
  ) +
  labs(x = "Posterior Estimate", y = NULL)

# Save Plot
ggsave(filename = file.path(output_dir, "All_Significant_Effects_80_95CI.png"),
       plot = combined_effects_plot, width = 12, height = 14, dpi = 300)


# ---------------------------------------------------------------
# Recreation effect plot 
# ---------------------------------------------------------------
# Recreation variable of interest
target_variable <- "b_rec_month_100_scaled"

# Extract and summarize
recreation_effects <- lapply(model_files, function(path) {
  model <- readRDS(path)
  species <- gsub("_model\\.rds", "", basename(path))
  samples <- as_draws_df(model)
  
  if (!target_variable %in% colnames(samples)) return(NULL)
  
  samples %>%
    select(all_of(target_variable)) %>%
    rename(Estimate = all_of(target_variable)) %>%
    mutate(Species = species) %>%
    group_by(Species) %>%
    summarize(
      median = median(Estimate),
      lower95 = quantile(Estimate, 0.025),
      upper95 = quantile(Estimate, 0.975),
      lower80 = quantile(Estimate, 0.10),
      upper80 = quantile(Estimate, 0.90),
      .groups = "drop"
    )
}) %>% bind_rows()

# Custom species order 
custom_order <- c("Ursus.americanus", "Lynx.rufus", "Odocoileus.hemionus", "Martes.americana", "Lepus.americanus")
others <- setdiff(recreation_effects$Species, custom_order)
ordered_species <- rev(c(custom_order, sort(others)))
recreation_effects$Species <- factor(recreation_effects$Species, levels = ordered_species)

# Plot
recreation_plot <- ggplot(recreation_effects, aes(y = Species)) +
  geom_segment(aes(x = lower80, xend = upper80, yend = Species), linewidth = 3, color = "#1f78b4") +
  geom_segment(aes(x = lower95, xend = upper95, yend = Species), linewidth = 1.2, color = "#084594") +
  geom_point(aes(x = median), size = 4.5, color = "#084594") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal(base_size = 18) +
  labs(x = "Posterior Estimate for Park-Month Recreation", y = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.spacing = unit(1.5, "lines")
  )

# Save
ggsave(filename = file.path(output_dir, "Recreation_Park_Month_Effects.png"),
       plot = recreation_plot, width = 12, height = 10, dpi = 300)


# ---------------------------------------------------------------
# Autocorrelation tests
# ---------------------------------------------------------------

# Define species-to-column mapping
species_name_map <- list(
  "Snowshoe Hare" = "lepus_detections",
  "Black Bear" = "ursus_detections",
  "Deer" = "odocoileus_detections",
  "Coyote" = "canis_detections",
  "Bobcat" = "lynx_detections",
  "Marten" = "martes_detections",
  "Marmot" = "marmota_detections",
  "Rare Carnivores" = "rare_detections"
)

# Define parks
parks <- c("GARI", "JOFF")

# Initialize empty dataframe to store Moran’s I results
moran_results <- data.frame(
  Species = character(),
  Park = character(),
  Month = character(),  
  Moran_I = numeric(),
  Expectation = numeric(),
  Variance = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each species
for (species_name in names(species_name_map)) {
  
  response_var <- species_name_map[[species_name]]
  
  cat("\n Using detection column:", response_var, "for species:", species_name, "\n")
  
  model_path <- paste0("~/Desktop/MSc_Chapter_1/output2/", species_name, "_model.rds")
  
  if (!file.exists(model_path)) {
    cat(" Model file not found for", species_name, "- Skipping.\n")
    next
  }
  
  model <- readRDS(model_path)
  
  # Get fitted values (predicted detection rates)
  predicted_values <- fitted(model, summary = TRUE)[, "Estimate"]
  
  # Compute residuals
  model_data <- model_data %>%
    mutate(
      Residuals = !!sym(response_var) - predicted_values,
      Month = floor_date(start_date, "month")
    )

  model_data_summary <- model_data %>%
    select(Month, Park, Deployment.Location.ID, Latitude, Longitude, Residuals)
  
  # Loop through each park
  for (park_name in parks) {
    
    # Get available months for the park
    available_months <- unique(model_data_summary$Month[model_data_summary$Park == park_name])
    
    if (length(available_months) == 0) {
      cat("️ No available months found for", species_name, "in", park_name, "- Skipping.\n")
      next
    }
    
    cat("Available months for", species_name, "in", park_name, ":", available_months, "\n")
    
    # Loop through each month
    for (selected_month in available_months) {
      
      cat("\n Running Moran's I for", species_name, "in", park_name, "for Month:", selected_month, "...\n")
      
      # Extract residuals for this month
      spatial_residuals <- model_data_summary %>%
        filter(Park == park_name, Month == selected_month) %>%
        select(Deployment.Location.ID, Latitude, Longitude, Residuals) %>%
        filter(!is.na(Residuals))  # Remove any NA values
      
      if (nrow(spatial_residuals) == 0) {
        cat(" No residuals found for", species_name, "in", park_name, "Month:", selected_month, "- Skipping.\n")
        next
      }
      
      # Convert to spatial object
      spatial_residuals_sf <- st_as_sf(spatial_residuals, coords = c("Longitude", "Latitude"), crs = 4326)
      
      # Ensure enough locations
      num_stations <- nrow(spatial_residuals_sf)
      if (num_stations < 4) {
        cat(" Not enough spatial locations for", species_name, "in", park_name, "Month:", selected_month, "- Skipping.\n")
        next
      }
      
      # **k = 4 for all months**
      k_neighbors <- min(4, num_stations - 1)  
      cat("Using", k_neighbors, "nearest neighbors\n")
      
      # Compute nearest neighbors
      coords <- st_coordinates(spatial_residuals_sf)
      nb <- knearneigh(coords, k = k_neighbors) %>% knn2nb()
      
      # Create spatial weights matrix
      listw <- nb2listw(nb, style = "W")
      
      # Run Moran’s I test
      moran_test <- moran.test(spatial_residuals$Residuals, listw)
      
      # Store results in dataframe
      moran_results <- rbind(moran_results, data.frame(
        Species = species_name,
        Park = park_name,
        Month = as.character(selected_month),  
        Moran_I = moran_test$estimate[1],
        Expectation = moran_test$estimate[2],
        Variance = moran_test$estimate[3],
        P_Value = moran_test$p.value
      ))
      
      # Print results
      print(moran_test)
    }
  }
}

# Save results to CSV
output_path <- "~/Desktop/MSc_Chapter_1/model_output/moran_results_monthly.csv"
write.csv(moran_results, output_path, row.names = FALSE)

# Display final results
cat("\n Moran's I monthly analysis complete. Results saved to:", output_path, "\n")
print(moran_results)

model_data$start_date <- as.Date(data$start_date)
model_data$MonthIndex <- as.numeric(as.factor(format(data$start_date, "%Y-%m")))

# Species to test
species_list <- c("Black Bear", "Marmot", "Marten", "Bobcat", "Snowshoe Hare", "Coyote", "Deer")

# Function to compute Durbin-Watson
compute_dw <- function(residuals) {
  sum_diff_sq <- sum(diff(residuals)^2)
  sum_sq <- sum(residuals^2)
  dw_stat <- sum_diff_sq / sum_sq
  return(dw_stat)
}

# Store results
dw_results <- data.frame(
  Species = character(),
  DW = numeric(),
  stringsAsFactors = FALSE
)

# Loop through species
for (species in species_list) {
  model_path <- file.path("~/Desktop/MSc_Chapter_1/output4", paste0(species, "_model.rds"))
  
  if (!file.exists(model_path)) {
    warning(paste("Model not found for", species))
    next
  }
  
  cat("Running manual Durbin-Watson for:", species, "\n")
  
  model <- readRDS(model_path)
  res <- tryCatch({
    residuals(model)[, 1]
  }, error = function(e) {
    warning(paste("Could not extract residuals for", species))
    return(NULL)
  })
  
  if (is.null(res)) next
  
  dw_stat <- compute_dw(res)
  
  dw_results <- rbind(dw_results, data.frame(
    Species = species,
    DW = round(dw_stat, 3)
  ))
}

# Save and print results
write.csv(dw_results, "~/Desktop/MSc_Chapter_1/model_output/manual_durbin_watson_results.csv", row.names = FALSE)
print(dw_results)

# ---------------------------------------------------------------
# Model Diagnostics
# ---------------------------------------------------------------
# Create output folders
trace_dir <- "~/Desktop/MSc_Chapter_1/output4/trace_plots"
pp_dir <- "~/Desktop/MSc_Chapter_1/output4/pp_checks"
dir.create(trace_dir, showWarnings = FALSE)
dir.create(pp_dir, showWarnings = FALSE)

# Loop through models
for (f in model_files) {
  model <- readRDS(f)
  species_name <- gsub("_model.rds", "", basename(f))
  cat("\nSaving diagnostics for:", species_name, "\n")
  
# Save posterior predictive check (mean) 
  pp_plot <- pp_check(model, type = "stat", stat = "mean", ndraws = 100) +
    ggtitle(paste("Posterior Predictive Check:", species_name)) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = file.path(pp_dir, paste0(species_name, "_ppcheck.png")),
    plot = pp_plot,
    width = 6,
    height = 4,
    dpi = 300
  )
  

# Save Trace Plot for Fixed Effects
  
  posterior_array <- as.array(model)
  b_params <- grep("^b_", dimnames(posterior_array)[[3]], value = TRUE)
  
  if (length(b_params) > 0) {
    trace_plot <- mcmc_trace(posterior_array, pars = b_params) +
      ggtitle(paste("Trace Plots:", species_name)) +
      theme_bw(base_size = 12)
    
    ggsave(
      filename = file.path(trace_dir, paste0(species_name, "_trace_plot.png")),
      plot = trace_plot,
      width = 10,
      height = 6,
      dpi = 300
    )
  } else {
    cat("No fixed effects found to plot.\n")
  }
}


# Prepare empty lists
pp_plots <- list()
trace_plots <- list()

# Collect plots
for (f in model_files) {
  model <- readRDS(f)
  species_name <- gsub("_model.rds", "", basename(f))
  cat("Preparing plots for:", species_name, "\n")
  
# Posterior predictive plot
  pp_plot <- pp_check(model, type = "stat", stat = "mean", ndraws = 100) +
    ggtitle(species_name) +
    theme_bw(base_size = 10)
  pp_plots[[species_name]] <- pp_plot
  
# Trace plot 
  posterior_array <- as.array(model)
  b_params <- grep("^b_", dimnames(posterior_array)[[3]], value = TRUE)
  if (length(b_params) > 0) {
    trace_plot <- mcmc_trace(posterior_array, pars = b_params) +
      ggtitle(species_name) +
      theme_bw(base_size = 10)
    trace_plots[[species_name]] <- trace_plot
  }
}

# Combine plots into one multi-panel layout

pp_combined <- wrap_plots(pp_plots, ncol = 3)
trace_combined <- wrap_plots(trace_plots, ncol = 1)

# Save combined plots
ggsave("~/Desktop/MSc_Chapter_1/output4/Combined_PP_Checks.pdf", pp_combined, width = 12, height = 9)
ggsave("~/Desktop/MSc_Chapter_1/output4/Combined_Trace_Plots.pdf", trace_combined, width = 10, height = 15)

# Extract Bayesian p-values 

# Initialize results list
bayes_pvals <- data.frame(
  Species = character(),
  Bayes_p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through models
for (f in model_files) {
  model <- readRDS(f)
  species_name <- gsub("_model.rds", "", basename(f))
  
  # Observed response variable name
  y_var <- as.character(model$formula$formula[[2]])
  y_obs <- model$data[[y_var]]
  
  # Generate posterior predictive draws
  y_rep <- posterior_predict(model, ndraws = 1000)
  
  # Test statistic: mean
  T_obs <- mean(y_obs)
  T_rep <- apply(y_rep, 1, mean)
  
  # Bayesian p-value
  p_bayes <- mean(T_rep >= T_obs)
  
  # Store in results
  bayes_pvals <- rbind(
    bayes_pvals,
    data.frame(Species = species_name, Bayes_p_value = round(p_bayes, 3))
  )
}

# View table
print(bayes_pvals)


# ---------------------------------------------------------------
## Community analysis 
# ---------------------------------------------------------------
   
# Define species to exclude
exclude_species <- c(
  "Unknown.species", "Homo.sapiens", "Bird.spp.", "Dendragapus.fuliginosus", "Tamias.striatus",
  "Mus.musculus", "Canis.familiaris", "Lagopus.leucura", "Felis.catus", "Neotamias.townsendii",
  "Phenacomys.intermedius", "Tamiasciurus.hudsonicus", "Glaucomys.sabrinus", "Myodes.gapperi",
  "Ochotona.princeps", "Peromyscus.maniculatus", "Tamiasciurus.douglasii", "Neotamias.minimus",
  "Neotamias.amoenus", "Sciurus.carolinensis", "Bonasa.umbellus", "Dendragapus.obscurus",
  "Falcipennis.canadensis", "Mustela.erminea"
)

# Generate and Save Accumulation Curves
plot_species_accumulation_q <- function(data_path, park_name, q = 0, y_label = "Species Richness", y_limit = c(0, 17)) {
  park_obs <- read.csv(data_path) %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d"))
  
  periods <- list(
    "2020" = interval(dmy("06-08-2020"), dmy("03-11-2020")),
    "2021" = interval(dmy("22-06-2021"), dmy("19-09-2021")),
    "2022" = interval(dmy("10-05-2022"), dmy("07-08-2022"))
  )
  
  year_colours <- c("2020" = "#E63946", "2021" = "#457B9D", "2022" = "#2A9D8F")
  plots <- list()
  
  for (year in names(periods)) {
    park_period <- park_obs %>%
      filter(Date %within% periods[[year]]) %>%
      group_by(Date) %>%
      summarise(
        Days = 1,
        across(where(is.numeric) & !any_of(c("Effort", "Deployment.Location.ID", exclude_species)), ~ sum(. > 0, na.rm = TRUE))
      )
    
    species_cols <- setdiff(names(park_period), c("Date", "Days"))
    if (length(species_cols) == 0) next
    
    inc_dat <- park_period %>% mutate(across(all_of(species_cols), ~ as.integer(. > 0)))
    
    total_days <- sum(park_period$Days)
    species_counts <- colSums(inc_dat[, species_cols, drop = FALSE])
    
    project_level <- list(project_level = c(total_days, species_counts))
    names(project_level$project_level)[1] <- "SamplingUnits"
    
    if (total_days > 0 && sum(species_counts) > 0) {
      out <- iNEXT(project_level, q = q, datatype = "incidence_freq", knots = 40, se = TRUE, conf = 0.95, nboot = 10)
      plot_data <- fortify(out, type = 1) %>%
        mutate(LineType = ifelse(Method == "Extrapolation", "dashed", "solid"))
      
      plots[[year]] <- ggplot(plot_data, aes(x = x, y = y)) +
        geom_line(aes(linetype = LineType), colour = year_colours[year], linewidth = 1.2) +
        geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), fill = year_colours[year], alpha = 0.2) +
        scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
        geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
        geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
        theme_minimal(base_size = 16) +
        theme(
          panel.grid = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 18, face = "bold")
        ) +
        labs(
          title = paste0(park_name, " - ", year),
          x = "Sampling Effort (Days)",
          y = y_label
        ) +
        coord_cartesian(ylim = y_limit)
    }
  }
  
  return(plots)
}

# Generate accumulation curves from processed data
joffre_richness <- plot_species_accumulation_q(
  "~/Desktop/MSc_Chapter_1/data/processed_data/JOFF_30min_Independent_daily_observations.csv",
  "Joffre Lakes", q = 0, y_label = "Species Richness", y_limit = c(0, 15)
)

garibaldi_richness <- plot_species_accumulation_q(
  "~/Desktop/MSc_Chapter_1/data/processed_data/GARI_30min_Independent_daily_observations.csv",
  "Garibaldi", q = 0, y_label = "Species Richness", y_limit = c(0, 15)
)

joffre_shannon <- plot_species_accumulation_q(
  "~/Desktop/MSc_Chapter_1/data/processed_data/JOFF_30min_Independent_daily_observations.csv",
  "Joffre Lakes", q = 1, y_label = "Shannon Diversity", y_limit = c(0, 8)
)

garibaldi_shannon <- plot_species_accumulation_q(
  "~/Desktop/MSc_Chapter_1/data/processed_data/GARI_30min_Independent_daily_observations.csv",
  "Garibaldi", q = 1, y_label = "Shannon Diversity", y_limit = c(0, 8)
)


accumulation_grid <- (
  joffre_richness[["2020"]] | joffre_richness[["2021"]] | joffre_richness[["2022"]]
) /
  (
    garibaldi_richness[["2020"]] | garibaldi_richness[["2021"]] | garibaldi_richness[["2022"]]
  ) /
  (
    joffre_shannon[["2020"]] | joffre_shannon[["2021"]] | joffre_shannon[["2022"]]
  ) /
  (
    garibaldi_shannon[["2020"]] | garibaldi_shannon[["2021"]] | garibaldi_shannon[["2022"]]
  )
# Save plots
ggsave(file.path(output_dir, "species_accumulation_4x3_grid.png"),
       plot = accumulation_grid, width = 18, height = 12, dpi = 300)

  
process_park_data_periods <- function(data_path, park_name) {
  park_obs <- read.csv(data_path) %>%
    mutate(Date = as.Date(Date))
  
  periods <- list(
    "2020" = interval(dmy("06-08-2020"), dmy("03-11-2020")),
    "2021" = interval(dmy("22-06-2021"), dmy("19-09-2021")),
    "2022" = interval(dmy("10-05-2022"), dmy("07-08-2022"))
  )
  
  richness_estimates <- shannon_estimates <- data.frame()
  
  for (year in names(periods)) {
    park_period <- park_obs %>% filter(Date %within% periods[[year]])
    if (nrow(park_period) == 0) next
    
    park_period <- park_period %>%
      select(-any_of(exclude_species)) %>%
      group_by(Date) %>%
      summarise(
        TotalEffort = sum(Effort, na.rm = TRUE),
        across(where(is.numeric) & !any_of("Effort"), ~ sum(. > 0, na.rm = TRUE))
      )
    
    species_cols <- setdiff(names(park_period), c("Date", "TotalEffort"))
    if (length(species_cols) == 0) next
    
    inc_dat <- park_period %>% mutate(across(all_of(species_cols), ~ as.integer(. > 0)))
    
    project_level <- list(project_level = c(sum(park_period$TotalEffort), colSums(inc_dat[, species_cols])))
    names(project_level$project_level)[1] <- "SamplingUnits"
    
    if (sum(project_level$project_level[-1]) > 0) {
      out <- iNEXT(project_level, q = c(0, 1), datatype = "incidence_freq", knots = 40, se = TRUE, conf = 0.95, nboot = 50)
      asy_estimates <- out$AsyEst %>%
        filter(Diversity %in% c("Species richness", "Shannon diversity")) %>%
        mutate(Year = year, Park = park_name) %>%
        select(Year, Park, Diversity, Estimator, LCL, UCL)
      
      richness_estimates <- bind_rows(richness_estimates, filter(asy_estimates, Diversity == "Species richness"))
      shannon_estimates <- bind_rows(shannon_estimates, filter(asy_estimates, Diversity == "Shannon diversity"))
    }
  }
  
  list(richness = richness_estimates, shannon = shannon_estimates)
}

# Run Diversity Processing for Each Park

joffre <- process_park_data_periods("~/Desktop/MSc_Chapter_1/data/processed_data/JOFF_30min_Independent_daily_observations.csv", "Joffre Lakes")
garibaldi <- process_park_data_periods("~/Desktop/MSc_Chapter_1/data/processed_data/GARI_30min_Independent_daily_observations.csv", "Garibaldi")

richness_estimates <- bind_rows(joffre$richness, garibaldi$richness)
shannon_estimates <- bind_rows(joffre$shannon, garibaldi$shannon)


# Calculate Evenness
custom_colors <- c("Joffre Lakes" = "#E63946", "Garibaldi" = "#457B9D")

evenness_estimates <- richness_estimates %>%
  rename(Richness = Estimator, LCL_R = LCL, UCL_R = UCL) %>%
  left_join(
    shannon_estimates %>%
      rename(Shannon_Hill = Estimator, LCL_S = LCL, UCL_S = UCL),
    by = c("Year", "Park")
  ) %>%
  mutate(
    Shannon = log(Shannon_Hill),
    Evenness = Shannon / log(Richness),
    SE_Shannon = (log(UCL_S) - log(LCL_S)) / (2 * 1.96),
    SE_Richness = (UCL_R - LCL_R) / (2 * 1.96),
    SE_Evenness = sqrt((1 / log(Richness))^2 * SE_Shannon^2 + (Shannon / (Richness * log(Richness)^2))^2 * SE_Richness^2),
    LCL_Evenness = Evenness - 1.96 * SE_Evenness,
    UCL_Evenness = Evenness + 1.96 * SE_Evenness
  )


# Plotting Functions and Outputs
plot_metrics <- function(data, y_label, show_legend = FALSE) {
  ggplot(data, aes(x = factor(Year), y = Estimator, color = Park)) +
    geom_point(size = 4, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.2, position = position_dodge(0.5)) +
    scale_color_manual(values = custom_colors) +
    theme_minimal(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_blank()
    ) +
    labs(x = "Year", y = y_label, color = "Park")
}

# Apply to each metric
richness_plot <- plot_metrics(richness_estimates, "Species Richness", show_legend = FALSE)
shannon_plot  <- plot_metrics(shannon_estimates, "Shannon Diversity", show_legend = TRUE)

# Evenness plot
evenness_plot <- ggplot(evenness_estimates, aes(x = factor(Year), y = Evenness, color = Park)) +
  geom_point(size = 4, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = LCL_Evenness, ymax = UCL_Evenness), width = 0.2, position = position_dodge(0.5)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(x = "Year", y = "Pielou's Evenness", color = "Park")

# Combine plots
diversity_grid <- richness_plot / shannon_plot / evenness_plot

# Save plot
ggsave(filename = file.path(output_dir, "diversity_metrics_grid.png"),
       plot = diversity_grid, width = 10, height = 12, dpi = 300)

# Save estimates
saveRDS(richness_estimates, "~/Desktop/MSc_Chapter_1/model_output/richness_estimates.rds")
saveRDS(shannon_estimates, "~/Desktop/MSc_Chapter_1/model_output/shannon_estimates.rds")
saveRDS(evenness_estimates, "~/Desktop/MSc_Chapter_1/model_output/evenness_estimates.rds")
 

# ---------------------------------------------------------------
# Bootstrap hypothesis test
# ---------------------------------------------------------------
  
# Bootstrap simulation function
simulate_bootstrap <- function(mean, se, n = 10000) {
  rnorm(n, mean = mean, sd = se)
}

# Function to run bootstrapped comparison
compare_groups_bootstrap <- function(data, group_var, metric, level1, level2) {
  sim1 <- data %>% filter(!!sym(group_var) == level1) %>% pull(Simulated) %>% .[[1]]
  sim2 <- data %>% filter(!!sym(group_var) == level2) %>% pull(Simulated) %>% .[[1]]
  
  est1 <- data %>% filter(!!sym(group_var) == level1) %>% pull(.data[[metric]])
  se1  <- data %>% filter(!!sym(group_var) == level1) %>% pull(.data[[paste0("SE_", metric)]])
  est2 <- data %>% filter(!!sym(group_var) == level2) %>% pull(.data[[metric]])
  se2  <- data %>% filter(!!sym(group_var) == level2) %>% pull(.data[[paste0("SE_", metric)]])
  
  diff_sim <- sim2 - sim1
  boot_p_value <- 2 * min(mean(diff_sim >= 0), mean(diff_sim <= 0))
  
  data.frame(
    Metric = metric,
    Comparison = paste(level1, "vs", level2),
    Mean_Diff = mean(diff_sim),
    CI_Lower = quantile(diff_sim, 0.025),
    CI_Upper = quantile(diff_sim, 0.975),
    Boot_P_Value = boot_p_value,
    Estimate_1 = est1,
    SE_1 = se1,
    Estimate_2 = est2,
    SE_2 = se2
  )
}

# Run Bootstrapped Comparisons
n_boot <- 1000
metrics <- c("Richness", "Shannon", "Evenness")
all_results <- list()

for (metric in metrics) {
  evenness_estimates <- evenness_estimates %>%
    mutate(Simulated = map2(.data[[metric]], SE_Evenness, ~ simulate_bootstrap(.x, .y, n_boot)))
  
  # Between-Year Comparisons (Within Each Park)
  for (park in unique(evenness_estimates$Park)) {
    park_data <- evenness_estimates %>% filter(Park == park)
    years <- unique(park_data$Year)
    
    for (i in 1:(length(years) - 1)) {
      for (j in (i + 1):length(years)) {
        res <- compare_groups_bootstrap(park_data, "Year", metric, years[i], years[j])
        res$Park <- park
        res$Comparison_Type <- "Between-Year"
        all_results <- append(all_results, list(res))
      }
    }
  }
  
  #Between-Park Comparisons (Within Each Year)
  for (year in unique(evenness_estimates$Year)) {
    year_data <- evenness_estimates %>% filter(Year == year)
    parks <- unique(year_data$Park)
    
    if (length(parks) == 2) {
      res <- compare_groups_bootstrap(year_data, "Park", metric, parks[1], parks[2])
      res$Year <- year
      res$Comparison_Type <- "Between-Park"
      all_results <- append(all_results, list(res))
    }
  }
}

# Combine and Correct Results
results_df <- bind_rows(all_results)
results_df$Boot_BY_P_Value <- p.adjust(results_df$Boot_P_Value, method = "BY")

# Round numeric values
results_df <- results_df %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save Results
write.csv(results_df, "~/Desktop/MSc_Chapter_1/model_output/diversity_bootstrap_BY_results.csv", row.names = FALSE)
 

