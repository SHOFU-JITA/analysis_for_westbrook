# =============================================================================
# Westbrook Project: Growth Strategy Cluster Analysis (Full Script)
# Updated: Automated Model Generation (Season-Matched Temp & Flow + Abundance)
# =============================================================================

# ===== 0. Load Packages =====
rm(list = ls(all.names = TRUE)) 

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(nnet)
library(broom)
library(purrr)
library(stringr)
library(ggplot2)
library(mclogit)

# =============================================================================
# 1. Basic Settings and Data Loading
# =============================================================================

file_path         <- "../csv/WB_combined_bkt_rev2.csv"
temp_file_path    <- "../rawdata/envDataWB.csv"
density_file_path <- "../JS/Abundance_Estimates_S3_Total_bkt.csv"
flow_file_path    <- "../Flow/analysis/Seasonal_Flow.csv" # Specify flow data path

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}
cat(paste0("Save directory: ", SCRIPT_DIR, "\n"))

num_clusters <- 3
CLUSTER_COLORS <- c("steelblue", "aquamarine3", "salmon3")

pdf_filename <- paste0("logit_env_flow_bkt", num_clusters,".pdf")
pdf_full_path <- file.path(SCRIPT_DIR, pdf_filename)

# --- Manual Plot Settings (High/Mean/Low) ---
MANUAL_PLOT_SETTINGS <- list(
  temp_summer = c(13, 15, 17),  
  temp_winter = NULL,
  temp_spring = NULL,
  temp_autumn = NULL,
  
  density_s3        = c(100, 300, 600), 
  raw_cohort_count  = c(20, 50, 100)  
)

cat(paste0("\nCreating PDF: ", pdf_full_path, " ...\n"))

if(!file.exists(file_path)) stop("Main data file not found: ", file_path)
df <- read.csv(file_path)
df$date <- as.Date(df$date)

cluster_col_name <- paste0("cluster_k", num_clusters)
df <- df %>%
  rename(ocluster_rename_dynamic = !!sym(cluster_col_name)) %>%
  mutate(
    ocluster_rename_dynamic = factor(ocluster_rename_dynamic, levels = as.character(1:num_clusters)),
    river = case_when(
      tolower(river) == "wb jimmy"    ~ "ol",
      tolower(river) == "wb mitchell" ~ "os",
      tolower(river) == "wb obear"    ~ "il",
      tolower(river) == "west brook"  ~ "wb",
      TRUE ~ tolower(river)
    )
  )

# =============================================================================
# 2. Create Model Dataset (Base)
# =============================================================================

df_analysis <- df %>%
  mutate(cohort_num = as.numeric(as.character(cohort))) %>% 
  filter(
    (river == "wb" & cohort_num %in% 2000:2015) |
      (river != "wb" & cohort_num %in% 2002:2015)
  )
cat(paste0("Data filtered: ", nrow(df_analysis), " rows (Original: ", nrow(df), " rows)\n"))

river_tags <- df_analysis %>%
  distinct(tag, river)

multi_tags <- river_tags %>%
  count(tag) %>%
  filter(n > 1) %>%
  pull(tag)

df_model <- df_analysis %>%
  mutate(
    river_grp = if_else(tag %in% multi_tags, "Multi", river),
    cohort_numeric = cohort_num
  ) %>%
  mutate(river_grp = factor(river_grp))

cat("Base data ready. Processing environmental data...\n")

# =============================================================================
# 3. Process Environmental Data (Temp, Density, Cohort Count, Flow)
# =============================================================================

# --- A. Process Temperature Data ---
raw_temp <- read.csv(temp_file_path)

temp_processed <- raw_temp %>%
  mutate(
    date = as.Date(date),
    month = month(date),
    year = year(date),
    river = case_when(
      tolower(river) == "wb jimmy"    ~ "ol",
      tolower(river) == "wb mitchell" ~ "os",
      tolower(river) == "wb obear"    ~ "il",
      tolower(river) == "west brook"  ~ "wb",
      TRUE ~ tolower(river)
    )
  ) %>%
  mutate(
    season = case_when(
      month %in% 2:4 ~ "temp_spring",
      month %in% 5:7 ~ "temp_summer",
      month %in% 8:10 ~ "temp_autumn",
      month %in% c(11, 12, 1) ~ "temp_winter",
      TRUE ~ NA_character_
    ),
    season_year = if_else(month == 1, year - 1, year)
  ) %>%
  filter(!is.na(season)) %>%
  group_by(river, season_year, season) %>%
  summarise(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = season, values_from = mean_temp)


# --- B. Process Density Data ---
raw_dens <- read.csv(density_file_path)

dens_processed <- raw_dens %>%
  filter(grepl("S3", season_label)) %>%
  mutate(
    season_year = as.numeric(str_sub(season_label, 1, 4)),
    river = tolower(river)
  ) %>%
  dplyr::select(river, season_year, density_s3 = N_est) 


# --- C. Calculate Observed Cohort Counts ---
cohort_counts <- df_analysis %>%
  group_by(river, cohort_num) %>%
  summarise(raw_cohort_count = n_distinct(tag), .groups = "drop") %>%
  mutate(river = tolower(river), season_year = cohort_num) %>%
  dplyr::select(river, season_year, raw_cohort_count)

# --- D. Process Flow Data ---
raw_flow <- read.csv(flow_file_path)

flow_processed <- raw_flow %>%
  mutate(
    river = case_when(
      tolower(Site_Name) == "OL"  ~ "ol",
      tolower(Site_Name) == "OS"  ~ "os",
      tolower(Site_Name) == "IL"  ~ "il",
      tolower(Site_Name) == "WB"  ~ "wb",
      TRUE ~ tolower(Site_Name)
    ),
    season_label = case_when(
      season == 1 ~ "spring",
      season == 2 ~ "summer",
      season == 3 ~ "autumn",
      season == 4 ~ "winter",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(season_label)) %>%
  pivot_wider(
    id_cols = c(river, year),
    names_from = season_label,
    values_from = c(Mean_Flow, High_Flow_Days, Low_Flow_Days),
    names_glue = "{.value}_{season_label}"
  ) %>%
  rename(season_year = year)

# --- E. Merge with Main Data ---
df_model_env <- df_model %>%
  filter(river_grp != "Multi") %>% 
  mutate(
    join_river = tolower(river_grp),
    join_year = cohort_numeric
  ) %>%
  left_join(temp_processed, by = c("join_river" = "river", "join_year" = "season_year")) %>%
  left_join(dens_processed, by = c("join_river" = "river", "join_year" = "season_year")) %>%
  left_join(cohort_counts, by = c("join_river" = "river", "join_year" = "season_year")) %>%
  left_join(flow_processed, by = c("join_river" = "river", "join_year" = "season_year"))

# Filter to complete cases for valid AIC comparison across all models
cols_to_check <- c(
  "temp_summer", "temp_winter", "temp_spring", "temp_autumn", "density_s3", "raw_cohort_count",
  names(df_model_env)[grepl("Flow", names(df_model_env))]
)
cols_to_check <- cols_to_check[cols_to_check %in% names(df_model_env)]
df_model_env <- df_model_env %>% drop_na(all_of(cols_to_check))

cat(paste0("Data merging complete: Using ", nrow(df_model_env), " complete rows.\n"))


# =============================================================================
# 3.5 Correlation Check
# =============================================================================
cat("\n--- Checking Correlations ---\n")
cor_vars <- cols_to_check
cor_matrix <- cor(df_model_env[, cor_vars], use = "pairwise.complete.obs")

cat("Correlation Matrix:\n")
print(round(cor_matrix, 2))
write.csv(cor_matrix, file.path(SCRIPT_DIR, paste0("Correlation_Matrix_", num_clusters, ".csv")))


# =============================================================================
# 4. Model Analysis (Automated Formula Generation)
# =============================================================================

# Scale numeric variables
scale_cols <- df_model_env %>% select(where(is.numeric), -cohort_numeric, -join_year) %>% names()
df_model_env_scaled <- df_model_env %>%
  mutate(across(all_of(scale_cols), ~ as.numeric(scale(.))))

cat("\n--- Generating Formulas (Season-Matched Temp & Flow + Abundance) ---\n")

formulas_mixed <- list()
formulas_mixed[["Null Model"]] <- ocluster_rename_dynamic ~ 1

seasons <- c("spring", "summer", "autumn", "winter")
flow_metrics <- c("Mean_Flow", "High_Flow_Days", "Low_Flow_Days")
abund_vars <- c("density_s3", "raw_cohort_count")

for (s in seasons) {
  t_var <- paste0("temp_", s)
  
  if (t_var %in% names(df_model_env_scaled)) {
    # [1] Temperature only
    formulas_mixed[[paste0(t_var, " Only")]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var))
    
    # [2] Temperature and Abundance (Density or Cohort Count)
    for (a_var in abund_vars) {
      if (a_var %in% names(df_model_env_scaled)) {
        formulas_mixed[[paste0(t_var, " + ", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "+", a_var))
        formulas_mixed[[paste0(t_var, " * ", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", a_var))
      }
    }
    
    for (fm in flow_metrics) {
      f_var <- paste0(fm, "_", s)
      
      if (f_var %in% names(df_model_env_scaled)) {
        # [3] Temperature and Flow
        formulas_mixed[[paste0(t_var, " + ", f_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "+", f_var))
        formulas_mixed[[paste0(t_var, " * ", f_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", f_var))
        
        # [4] Temperature, Flow, and Abundance
        for (a_var in abund_vars) {
          if (a_var %in% names(df_model_env_scaled)) {
            formulas_mixed[[paste0(t_var, " + ", f_var, " + ", a_var)]] <- 
              as.formula(paste("ocluster_rename_dynamic ~", t_var, "+", f_var, "+", a_var))
            formulas_mixed[[paste0(t_var, " * ", f_var, " * ", a_var)]] <- 
              as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", f_var, "*", a_var))
          }
        }
      }
    }
  }
}

cat(paste0("Generated ", length(formulas_mixed), " model formulas.\n"))

# --- Step 2: Model Comparison (AIC & R2) ---
get_ll_safe <- function(mod) {
  if (is.null(mod)) return(NA)
  ll <- try(as.numeric(logLik(mod)), silent = TRUE)
  if (!inherits(ll, "try-error") && length(ll) > 0 && !is.na(ll)) return(ll)
  if (!is.null(mod$deviance)) return(-0.5 * mod$deviance)
  return(NA)
}

cat("Calculating Null Model...\n")
null_fit <- tryCatch(
  mclogit::mblogit(ocluster_rename_dynamic ~ 1, 
                   random = ~ 1 | river_grp, 
                   data = df_model_env_scaled, 
                   maxit = 50),
  error = function(e) return(NULL)
)

logLik_null <- get_ll_safe(null_fit)

model_comparison <- map_dfr(names(formulas_mixed), function(name) {
  # cat(paste0("Running: ", name, " ...\n")) # Silenced to reduce console spam
  f <- formulas_mixed[[name]]
  
  fit <- tryCatch(
    mclogit::mblogit(f, 
                     random = ~ 1 | river_grp, 
                     data = df_model_env_scaled, 
                     maxit = 50),
    error = function(e) return(NULL)
  )
  
  ll_model <- get_ll_safe(fit)
  
  if (!is.na(ll_model) && !is.na(logLik_null)) {
    r2_mcfadden <- 1 - (ll_model / logLik_null)
    current_aic <- tryCatch(AIC(fit), error = function(e) {
      k <- length(coef(fit)) 
      return(2 * k - 2 * ll_model)
    })
    
    tibble(
      Model_Name  = name, 
      AIC         = current_aic, 
      McFadden_R2 = round(r2_mcfadden, 4),
      Status      = "Success"
    )
  } else {
    tibble(Model_Name = name, AIC = NA, McFadden_R2 = NA, Status = "Failed")
  }
})

model_ranking <- model_comparison %>%
  filter(Status == "Success") %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

print(head(model_ranking, 10))

csv_filename <- paste0("Modelselection_bkt_flow", num_clusters, ".csv")
write.csv(model_ranking, file.path(SCRIPT_DIR, csv_filename), row.names = FALSE)

# =============================================================================
# 5. Top 3 Models Detailed Analysis & Visualization
# =============================================================================

pdf(pdf_full_path, width = 12, height = 8)

if (nrow(model_ranking) > 0) {
  
  top_n_models <- min(3, nrow(model_ranking))
  
  for (m in 1:top_n_models) {
    
    best_model_name <- model_ranking$Model_Name[m]
    best_formula    <- formulas_mixed[[best_model_name]]
    
    cat("\n--- Plotting Model Rank", m, "---\n")
    cat(best_model_name, "\n")
    
    best_model <- mclogit::mblogit(
      best_formula,
      random = ~ 1 | river_grp,
      data   = df_model_env_scaled,
      maxit  = 50
    )
    
    predictors <- all.vars(best_formula)[-1]
    
    if (length(predictors) == 0) {
      cat("Null model: Skipping visualization\n")
    } else {
      plot_combinations <- list()
      
      if (length(predictors) == 1) {
        plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = "Dummy")
      } else if (length(predictors) == 2) {
        plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = predictors[2])
        plot_combinations[[2]] <- list(x_var = predictors[2], facet_var = predictors[1])
      } else if (length(predictors) >= 3) {
        # If 3 or more variables, prioritize temp and abundance (density/cohort) for the plot
        temp_var <- predictors[grepl("temp", predictors)][1]
        abund_var <- predictors[predictors %in% abund_vars][1]
        
        if (!is.na(temp_var) && !is.na(abund_var)) {
          plot_combinations[[1]] <- list(x_var = temp_var, facet_var = abund_var)
          plot_combinations[[2]] <- list(x_var = abund_var, facet_var = temp_var)
        } else {
          # Fallback if specific variables aren't matched
          plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = predictors[2])
        }
      }
      
      get_stats <- function(v) {
        list(mean = mean(df_model_env[[v]], na.rm = TRUE),
             sd   = sd(df_model_env[[v]], na.rm = TRUE))
      }
      
      for (combo in plot_combinations) {
        x_var     <- combo$x_var
        facet_var <- combo$facet_var
        
        x_stats <- get_stats(x_var)
        x_seq_real <- seq(
          min(df_model_env[[x_var]], na.rm = TRUE),
          max(df_model_env[[x_var]], na.rm = TRUE),
          length.out = 100
        )
        
        if (facet_var != "Dummy") {
          facet_stats <- get_stats(facet_var)
          
          if (!is.null(MANUAL_PLOT_SETTINGS[[facet_var]])) {
            facet_vals_real <- MANUAL_PLOT_SETTINGS[[facet_var]]
          } else {
            facet_vals_real <- c(
              quantile(df_model_env[[facet_var]], 0.1, na.rm = TRUE),
              mean(df_model_env[[facet_var]], na.rm = TRUE),
              quantile(df_model_env[[facet_var]], 0.9, na.rm = TRUE)
            )
          }
          prefix_labels <- c("Low", "Mean", "High")
          if (length(facet_vals_real) == 3) {
            facet_labels <- paste0(prefix_labels, " (", round(facet_vals_real, 1), ")")
          } else {
            facet_labels <- paste0("Val (", round(facet_vals_real, 1), ")")
          }
        } else {
          facet_vals_real <- 0
          facet_labels    <- "All data"
          facet_stats     <- list(mean = 0, sd = 1)
        }
        
        pred_grid_real <- expand.grid(
          Var_X = x_seq_real,
          Var_F = facet_vals_real
        )
        colnames(pred_grid_real) <- c(x_var, facet_var)
        
        pred_grid_real$Facet_Label <- factor(
          pred_grid_real[[facet_var]],
          levels = facet_vals_real,
          labels = facet_labels
        )
        
        # If a third variable is in the model, hold it at its mean
        for (p in predictors) {
          if (p != x_var && p != facet_var) {
            pred_grid_real[[p]] <- get_stats(p)$mean
          }
        }
        
        # Scale for prediction
        pred_grid_scaled <- pred_grid_real
        for (p in predictors) {
          if (p == facet_var && facet_var == "Dummy") next
          p_stats <- get_stats(p)
          pred_grid_scaled[[p]] <- (pred_grid_real[[p]] - p_stats$mean) / p_stats$sd
        }
        pred_grid_scaled$river_grp <- df_model_env_scaled$river_grp[1]
        
        set.seed(123)
        beta_hat <- coef(best_model)
        vcov_mat <- vcov(best_model)
        
        n_sim <- 1000
        beta_sim <- MASS::mvrnorm(n_sim, beta_hat, vcov_mat)
        sim_probs <- array(NA, dim = c(nrow(pred_grid_scaled), num_clusters, n_sim))
        
        for (i in 1:n_sim) {
          best_model$coefficients <- beta_sim[i, ] 
          sim_probs[, , i] <- predict(best_model, newdata = pred_grid_scaled, type = "response")
        }
        best_model$coefficients <- beta_hat 
        
        probs_q <- apply(sim_probs, c(1, 2), quantile, probs = c(0.025, 0.5, 0.975))
        
        plot_df <- map_dfr(1:num_clusters, function(k) {
          tibble(
            !!x_var := pred_grid_real[[x_var]],
            Facet_Label = pred_grid_real$Facet_Label,
            cluster = factor(k),
            prob  = probs_q["50%", , k],
            lower = probs_q["2.5%", , k],
            upper = probs_q["97.5%", , k]
          )
        })
        
        p <- ggplot(plot_df,
                    aes(x = .data[[x_var]],
                        y = prob,
                        color = cluster,
                        fill  = cluster)) +
          geom_ribbon(aes(ymin = lower, ymax = upper),
                      alpha = 0.2,
                      color = NA) +
          geom_line(linewidth = 1.8) +
          scale_y_continuous(labels = scales::percent,
                             limits = c(0, 1)) +
          scale_color_manual(values = CLUSTER_COLORS) +
          scale_fill_manual(values = CLUSTER_COLORS) +
          labs(
            title = paste("Model Rank", m, ":", best_model_name),
            subtitle = "Mixed-effects multinomial logit (95% CI)",
            x = x_var,
            y = "Probability"
          ) +
          theme_bw(base_size = 24) + 
          theme(
            text = element_text(face = "bold", color = "black"), 
            axis.text = element_text(color = "black"),            
            legend.position = "bottom",
            panel.border = element_blank(),             
            axis.line = element_line(color = "black"),
            strip.background = element_blank(),        
            strip.text = element_text(face = "bold", color = "black") 
          )
        
        if (facet_var != "Dummy") {
          p <- p + facet_wrap(~Facet_Label)
        }
        
        print(p) 
      }
    }
  }
}

dev.off()
cat("PDF creation complete\n")

