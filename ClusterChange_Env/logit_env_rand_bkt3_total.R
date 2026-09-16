# =============================================================================
# Westbrook Project: Growth Strategy Cluster Analysis (Full Script)
# Updated: Includes Interaction Terms (Temp * Density)
# Modified: Dynamic Labels with Values (e.g., "High (22°C)")
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
density_file_path_sp2 <- "../JS/Abundance_Estimates_S3_Total_bnt.csv" 
density_file_path_sp3 <- "../JS/Abundance_Estimates_S3_Total_ats.csv" 

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}
cat(paste0("Save directory: ", SCRIPT_DIR, "\n"))

num_clusters <- 3
CLUSTER_COLORS <- c("steelblue", "aquamarine3", "salmon3")

pdf_filename <- paste0("logit_env_rand_bkt", num_clusters,".pdf")
pdf_full_path <- file.path(SCRIPT_DIR, pdf_filename)

# --- Manual Plot Settings (High/Mean/Low) ---
MANUAL_PLOT_SETTINGS <- list(
  temp_summer = c(13, 15, 17),  
  temp_winter = NULL,
  temp_spring = NULL,
  temp_autumn = NULL,
  
  density_s3        = c(100, 300, 600), 
  Dens_total       = c(200, 600, 1000),
  raw_cohort_count = NULL  
)

cat(paste0("\nCreating PDF: ", pdf_full_path, " ...\n"))

if(!file.exists(file_path)) stop("Main data file not found: ", file_path)
df <- read.csv(file_path)
df$date <- as.Date(df$date)

cluster_col_name <- paste0("cluster_k", num_clusters)

df <- df %>%
  mutate(
    !!sym(cluster_col_name) := factor(!!sym(cluster_col_name), levels = as.character(1:num_clusters)),
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
    # riverがすでに"wb"等に変換されているので、それに合わせてフィルター条件も修正
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
# 3. Process Environmental Data (Temperature, Density, Cohort Count)
# =============================================================================

# --- A. Process Temperature Data ---
if (file.exists(temp_file_path)) {
  raw_temp <- read.csv(temp_file_path)
  
  temp_processed <- raw_temp %>%
    mutate(
      date = as.Date(date),
      month = month(date),
      year = year(date),
      # 川の名前を変換
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
  
} else {
  cat("Warning: Temperature file not found. Using mock data.\n")
  temp_processed <- df_model %>% distinct(river_grp, cohort_numeric) %>%
    rename(river = river_grp, season_year = cohort_numeric) %>%
    mutate(temp_summer = 20, temp_winter = 2, temp_spring = 10, temp_autumn = 15)
}

# --- B. Process Density Data ---
# 同じ構造のデータを処理するための共通関数
load_density <- function(path, col_name) {
  read.csv(path) %>%
    filter(grepl("S3", season_label)) %>%
    mutate(
      season_year = as.numeric(str_sub(season_label, 1, 4)),
      river = tolower(river)
    ) %>%
    dplyr::select(river, season_year, !!sym(col_name) := N_est)
}

dens_bkt <- load_density(density_file_path, "density_s3")
dens_sp2 <- load_density(density_file_path_sp2, "density_sp2")
dens_sp3 <- load_density(density_file_path_sp3, "density_sp3")

# 結合して Dens_total を計算
dens_processed <- dens_bkt %>%
  full_join(dens_sp2, by = c("river", "season_year")) %>%
  full_join(dens_sp3, by = c("river", "season_year")) %>%
  mutate(
    density_s3  = replace_na(density_s3, 0),
    density_sp2 = replace_na(density_sp2, 0),
    density_sp3 = replace_na(density_sp3, 0),
    Dens_total  = density_s3 + density_sp2 + density_sp3
  ) %>%
  dplyr::select(river, season_year, density_s3, Dens_total)

# --- C. Calculate Observed Cohort Counts ---
cohort_counts <- df_analysis %>%
  group_by(river, cohort_num) %>%
  summarise(raw_cohort_count = n_distinct(tag), .groups = "drop") %>%
  mutate(
    river = tolower(river),
    season_year = cohort_num
  ) %>%
  dplyr::select(river, season_year, raw_cohort_count)

# --- D. Merge with Main Data ---
df_model_env <- df_model %>%
  filter(river_grp != "Multi") %>% 
  mutate(
    join_river = tolower(river_grp),
    join_year = cohort_numeric
  ) %>%
  left_join(temp_processed, by = c("join_river" = "river", "join_year" = "season_year")) %>%
  left_join(dens_processed, by = c("join_river" = "river", "join_year" = "season_year")) %>%
  left_join(cohort_counts, by = c("join_river" = "river", "join_year" = "season_year")) %>%
  drop_na(temp_summer, temp_winter, temp_spring, temp_autumn, density_s3, Dens_total)

cat(paste0("Data merging complete: Using ", nrow(df_model_env), " rows.\n"))

# =============================================================================
# 4. Model Analysis (AIC Comparison & Pseudo R-squared)
# =============================================================================

# --- Step 1: Check Correlations ---
cat("\n--- Correlations among Predictors ---\n")
check_cols <- c("temp_spring", "temp_summer", "temp_autumn", "temp_winter", "density_s3", "Dens_total", "raw_cohort_count")
check_cols <- check_cols[check_cols %in% names(df_model_env)]
cor_matrix <- cor(df_model_env[check_cols], use = "complete.obs")
print(round(cor_matrix, 2))

# --- Step 2: Model Comparison (Mixed Effects mblogit: AIC & R2) ---
cat("\n--- Model Comparison (Mixed Effects mblogit: AIC & R2) ---\n")

df_model_env_scaled <- df_model_env %>%
  mutate(across(c(temp_summer, temp_winter, temp_spring, temp_autumn, 
                  density_s3, Dens_total, raw_cohort_count), ~ as.numeric(scale(.))))

predictors_list <- list(
  "Null Model"      = "1",
  "Summer Temp"     = "temp_summer",
  "Winter Temp"     = "temp_winter",
  "Spring Temp"     = "temp_spring",
  "Autumn Temp"     = "temp_autumn",
  "Density S3"      = "density_s3",
  "Count"           = "raw_cohort_count",
  "Summer + Dens"   = "temp_summer + density_s3",
  "Winter + Dens"   = "temp_winter + density_s3",
  "Spring + Dens"   = "temp_spring + density_s3",
  "Autumn + Dens"   = "temp_autumn + density_s3",
  "Summer * Dens"   = "temp_summer * density_s3",
  "Winter * Dens"   = "temp_winter * density_s3",
  "Spring * Dens"   = "temp_spring * density_s3",
  "Autumn * Dens"   = "temp_autumn * density_s3",
  "Summer * Count"  = "temp_summer * raw_cohort_count",
  "Winter * Count"  = "temp_winter * raw_cohort_count",
  "Spring * Count"  = "temp_spring * raw_cohort_count",
  "Autumn * Count"  = "temp_autumn * raw_cohort_count",
  "Summer + Count"  = "temp_summer + raw_cohort_count",
  "Winter + Count"  = "temp_winter + raw_cohort_count",
  "Spring + Count"  = "temp_spring + raw_cohort_count",
  "Autumn + Count"  = "temp_autumn + raw_cohort_count",
  "Dens Total"      = "Dens_total",
  "Summer + Total"  = "temp_summer + Dens_total",
  "Winter + Total"  = "temp_winter + Dens_total",
  "Summer * Total"  = "temp_summer * Dens_total",
  "Winter * Total"  = "temp_winter * Dens_total",
  "Autumn + Total"  = "temp_autumn + Dens_total",
  "Spring + Total"  = "temp_spring + Dens_total",
  "Autumn * Total"  = "temp_autumn * Dens_total",
  "Spring * Total"  = "temp_spring * Dens_total"
)

formulas_mixed <- lapply(predictors_list, function(preds) {
  as.formula(paste(cluster_col_name, "~", preds))
})

get_ll_safe <- function(mod) {
  if (is.null(mod)) return(NA)
  
  ll <- try(as.numeric(logLik(mod)), silent = TRUE)
  if (!inherits(ll, "try-error") && length(ll) > 0 && !is.na(ll)) return(ll)
  
  if (!is.null(mod$deviance)) return(-0.5 * mod$deviance)
  
  return(NA)
}

cat("Calculating Null Model...\n")
null_fit <- tryCatch(
  mclogit::mblogit(formulas_mixed[["Null Model"]], 
                   random = ~ 1 | river_grp, 
                   data = df_model_env_scaled, 
                   maxit = 50),
  error = function(e) return(NULL)
)

logLik_null <- get_ll_safe(null_fit)

if (is.na(logLik_null)) {
  cat("Warning: logLik failed, using manual deviance extraction.\n")
  logLik_null <- -0.5 * null_fit$deviance
}

model_comparison <- map_dfr(names(formulas_mixed), function(name) {
  cat(paste0("Running: ", name, " ...\n"))
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


print(model_ranking)

csv_filename <- paste0("Modelselection_rand_bkt", num_clusters, ".csv")
write.csv(model_ranking, file.path(SCRIPT_DIR, csv_filename), row.names = FALSE)

# =============================================================================
# 5. Best Model Detailed Analysis & Visualization
# =============================================================================

pdf(pdf_full_path, width = 12, height = 8)

if (nrow(model_ranking) > 0) {
  
  best_model_name <- model_ranking$Model_Name[1]
  best_formula    <- formulas_mixed[[best_model_name]]
  
  cat("\n--- Best Model ---\n")
  cat(best_model_name, "\n")
  
  best_model <- mclogit::mblogit(
    best_formula,
    random = ~ 1 | river_grp,
    data   = df_model_env_scaled,
    maxit  = 50
  )
  
  # ---------------------------------------------------------------------------
  # 1. Plot Pattern Setup
  # ---------------------------------------------------------------------------
  
  predictors <- all.vars(best_formula)[-1]
  
  if (length(predictors) == 0) {
    cat("Null model: Skipping visualization\n")
  } else {
    plot_combinations <- list()
    if (length(predictors) == 1) {
      plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = NULL)
    } else {
      plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = predictors[2])
      plot_combinations[[2]] <- list(x_var = predictors[2], facet_var = predictors[1])
    }
    
    get_stats <- function(v) {
      list(mean = mean(df_model_env[[v]], na.rm = TRUE),
           sd   = sd(df_model_env[[v]], na.rm = TRUE))
    }
    
    # =========================================================================
    # Generate Plots for Each Combination
    # =========================================================================
    for (combo in plot_combinations) {
      x_var     <- combo$x_var
      facet_var <- combo$facet_var
      
      x_stats <- get_stats(x_var)
      
      x_seq_real <- seq(
        min(df_model_env[[x_var]], na.rm = TRUE),
        max(df_model_env[[x_var]], na.rm = TRUE),
        length.out = 100
      )
      
      if (!is.null(facet_var)) {
        
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
        facet_var       <- "Dummy"
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
      
      pred_grid_scaled <- pred_grid_real
      pred_grid_scaled[[x_var]] <-
        (pred_grid_real[[x_var]] - x_stats$mean) / x_stats$sd
      
      if (facet_var != "Dummy") {
        pred_grid_scaled[[facet_var]] <-
          (pred_grid_real[[facet_var]] - facet_stats$mean) / facet_stats$sd
      }
      
      pred_grid_scaled$river_grp <- df_model_env_scaled$river_grp[1]
      
      set.seed(123)
      
      beta_hat <- coef(best_model)
      vcov_mat <- vcov(best_model)
      
      n_sim <- 1000
      beta_sim <- MASS::mvrnorm(n_sim, beta_hat, vcov_mat)
      
      sim_probs <- array(
        NA,
        dim = c(nrow(pred_grid_scaled), num_clusters, n_sim)
      )
      
      for (i in 1:n_sim) {
        best_model$coefficients <- beta_sim[i, ] 
        sim_probs[, , i] <- predict(
          best_model,
          newdata = pred_grid_scaled,
          type = "response"
        )
      }
      
      best_model$coefficients <- beta_hat 
      
      probs_q <- apply(
        sim_probs,
        c(1, 2),
        quantile,
        probs = c(0.025, 0.5, 0.975)
      )
      
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
          title = paste("Model prediction:", best_model_name),
          subtitle = "Mixed-effects multinomial logit (95% CI)",
          x = x_var,
          y = "Probability"
        ) +
        theme_bw(base_size = 24) + 
        theme(
          text = element_text(face = "bold", color = "black"), 
          axis.text = element_text(color = "black"),           
          legend.position = "bottom"
        )
      
      if (facet_var != "Dummy") {
        p <- p + facet_wrap(~Facet_Label)
      }
      
      print(p) 
    }
  }
}

dev.off()
cat("PDF creation complete\n")
