# ===================================================================
# Full R Script: CJS Analysis (Integrated Model)
# - No Interaction (Main Effects Only)
# - Includes Season, Age, Cohort, and Cluster Effects
# - Emigration Correction Included
# ===================================================================

# ===== 0. Environment Reset =====
rm(list = ls(all.names = TRUE))
gc()

# ===== 0.1. Configuration: Set Parameters =====
# Set the number of clusters
K_CLUSTERS <- 4

# Set the cohort range for analysis
START_COHORT <- 2000
END_COHORT   <- 2015

# ===== 0.5. Set Working Directory to Script Location =====
tryCatch({
  if (!require("rstudioapi")) install.packages("rstudioapi")
  if (rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getSourceEditorContext()$path
    if (script_path != "") {
      script_dir <- dirname(script_path)
      setwd(script_dir)
      message(paste0("Working directory changed to: ", script_dir))
    } else {
      warning("Script is not saved. Please save the file to set the working directory.")
    }
  }
}, error = function(e) {
  message("Could not automatically set working directory. Please set it manually if needed.")
})

# ===== 1. Load Packages =====
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, ggplot2, RMark, lubridate, purrr, tibble, rlang)

# -------------------------------------------------------------------
# Auto-generate variables based on K_CLUSTERS and cohort range
# -------------------------------------------------------------------
CLUSTER_COL_NAME   <- paste0("cluster_k", K_CLUSTERS)          
FILE_SUFFIX        <- paste0("bkt", K_CLUSTERS) 
EMIGRATION_FILE_PATH <- paste0("../Movingout/moveout_by_cluster_k", K_CLUSTERS, "_bkt.csv")
START_SEASON_INDEX <- 2

# ===================================================================
# 1.5. Input: Known Emigration Rates
# ===================================================================
DF_EMIGRATION <- read_csv(EMIGRATION_FILE_PATH) %>%
  select(
    cluster_id = all_of(CLUSTER_COL_NAME),
    season,
    ratio_move_sum
  ) %>%
  mutate(
    !!sym(CLUSTER_COL_NAME) := as.character(cluster_id),
    season_label = case_when(
      season == 1 ~ "1:Sp-Su",
      season == 2 ~ "2:Su-Au",
      season == 3 ~ "3:Au-Wi",
      season == 4 ~ "4:Wi-Sp",
      TRUE ~ as.character(season)
    ),
    emigration_rate = ratio_move_sum,
    fidelity = 1 - emigration_rate
  ) %>%
  select(
    all_of(CLUSTER_COL_NAME),
    season_label,
    emigration_rate,
    fidelity
  )

# ===================================================================
# 2. Load & Preprocess Data
# ===================================================================

df_raw <- read_csv("../Clustering/WB_bkt_clustered_rev2.csv")
df_raw$date <- as.Date(df_raw$date)

# ===== 3. Calculate SEASONS (with Winter Correction & Cohort Filter) =====
df <- df_raw %>%
  # Filter data by the specified cohort range
  filter(cohort >= START_COHORT & cohort <= END_COHORT) %>%
  mutate(
    month_val = month(date),
    cal_year = year(date),
    bio_year = if_else(month_val <= 2 & season == 4, cal_year - 1, cal_year),
    season_from_cohort = (bio_year - cohort) * 4 + (season - 1)
  ) %>%
  select(tag, !!sym(CLUSTER_COL_NAME), cohort, season_from_cohort) %>%
  distinct()

# ===================================================================
# 4. Create Capture History (Left Truncation: s3+)
# ===================================================================

capture_history_full <- df %>%
  mutate(capture = 1) %>%
  pivot_wider(
    id_cols = c(tag, !!sym(CLUSTER_COL_NAME), cohort),
    names_from = season_from_cohort, 
    values_from = capture, 
    values_fill = list(capture = 0)
  ) %>%
  arrange(tag)

season_cols <- names(capture_history_full)[!names(capture_history_full) %in% c("tag", CLUSTER_COL_NAME, "cohort")]
season_cols_sorted <- sort(as.numeric(season_cols))

keep_cols_num <- season_cols_sorted[season_cols_sorted >= START_SEASON_INDEX]
if (length(keep_cols_num) < 2) stop("Error: Not enough seasons left after truncation.")

keep_cols_str <- as.character(keep_cols_num)
new_col_names <- paste0("s", keep_cols_num)

capture_history <- capture_history_full %>%
  select(tag, !!sym(CLUSTER_COL_NAME), cohort, all_of(keep_cols_str)) %>%
  rename_at(vars(keep_cols_str), ~new_col_names) %>%
  filter(rowSums(select(., starts_with("s"))) > 0)

message(paste0("Data prepared. Analysis starts from s", START_SEASON_INDEX, "."))
message(paste0("Cohorts included: ", START_COHORT, " to ", END_COHORT))
message(paste0("N individuals: ", nrow(capture_history)))

# ===================================================================
# 5. CJS Integrated Model Function
# ===================================================================

run_cjs_integrated <- function(data, phi_formula, p_formula, start_offset) {
  ch_data <- data
  
  cols_to_paste <- names(ch_data)[!names(ch_data) %in% c("tag", CLUSTER_COL_NAME, "cohort")]
  ch_data$ch <- apply(ch_data[, cols_to_paste], 1, paste0, collapse = "")
  
  # Convert cohort and cluster to factors
  ch_data <- ch_data %>% 
    select(tag, ch, cohort, !!sym(CLUSTER_COL_NAME)) %>%
    mutate(across(c(cohort, !!sym(CLUSTER_COL_NAME)), as.factor))
  
  if (nrow(ch_data) < 2 || nchar(ch_data$ch[1]) <= 1) return(NULL)
  
  # --- process.data (All clusters combined) ---
  processed <- process.data(ch_data, model = "CJS", 
                            groups = c("cohort", CLUSTER_COL_NAME), 
                            begin.time = 0, 
                            time.intervals = rep(1, nchar(ch_data$ch[1]) - 1))
  
  # --- make.design.data ---
  design <- make.design.data(processed)
  
  design$Phi$time <- as.numeric(as.character(design$Phi$time))
  design$p$time   <- as.numeric(as.character(design$p$time))
  
  real_time_phi <- design$Phi$time + start_offset
  real_time_p   <- design$p$time + start_offset
  
  # Add age group effect
  design$Phi$age_group <- factor(
    ifelse(real_time_phi < 4, "Age0",
           ifelse(real_time_phi < 8, "Age1", "Age2+")),
    levels = c("Age1", "Age2+", "Age0")
  )
  
  # Add season effect
  season_labels <- c("1:Sp-Su", "2:Su-Au", "3:Au-Wi", "4:Wi-Sp")
  design$Phi$season_effect <- factor((real_time_phi %% 4) + 1, labels = season_labels)
  design$p$season_effect   <- factor((real_time_p %% 4) + 1, labels = season_labels)
  
  # --- Model Execution ---
  model <- mark(processed, design,
                model.parameters = list(Phi = list(formula = phi_formula),
                                        p = list(formula = p_formula)),
                silent = TRUE, output = FALSE, delete = TRUE)
  
  return(list(model = model, processed = processed, design = design))
}

# ===================================================================
# 6. Run Integrated Analysis (No Interactions)
# ===================================================================

cl_var <- CLUSTER_COL_NAME

# Define candidate models (Main effects only)
model_defs <- list(
  # 1. Null model (No cluster effect)
  list(name = "Phi(Age+Season+Cohort)_p(Season)", 
       phi_formula = ~age_group + season_effect + cohort, 
       p_formula = ~season_effect),
  
  # 2. Main effect model (Survival varies by cluster)
  list(name = "Phi(Cluster+Age+Season+Cohort)_p(Season)", 
       phi_formula = as.formula(paste0("~", cl_var, " + age_group + season_effect + cohort")), 
       p_formula = ~season_effect)
)

results_list <- list()

cat("\nRunning Integrated Models (Main Effects Only)...\n")

for (model_def in model_defs) {
  cat(paste0("Running: ", model_def$name, " ... "))
  
  res <- tryCatch({
    run_cjs_integrated(capture_history, 
                       model_def$phi_formula, 
                       model_def$p_formula, 
                       START_SEASON_INDEX)
  }, error = function(e) { cat("Error\n"); return(NULL) })
  
  if (!is.null(res)) {
    res$name <- model_def$name
    results_list[[model_def$name]] <- res
    cat(paste0("Done (AICc: ", round(res$model$results$AICc, 2), ")\n"))
  }
}

# --- Model selection via AICc ---
df_aic <- bind_rows(lapply(results_list, function(x) {
  data.frame(
    Model = x$name,
    npar = x$model$results$npar,
    AICc = x$model$results$AICc
  )
})) %>%
  arrange(AICc) %>%
  mutate(DeltaAICc = AICc - min(AICc))

cat("\n=== Model Selection Table ===\n")
print(df_aic)

# Extract best model
best_model_name <- df_aic$Model[1]
best_model_obj <- results_list[[best_model_name]]

cat(paste0("\nBest Model Selected: ", best_model_name, "\n"))

# ===================================================================
# 7. Extract Best Model Results
# ===================================================================

season_labels <- c("1:Sp-Su", "2:Su-Au", "3:Au-Wi", "4:Wi-Sp")

# ----- Extract survival (Phi) -----
phi_real_df <- get.real(best_model_obj$model, "Phi", se = TRUE, vcv = FALSE)

final_phi_results <- phi_real_df %>%
  rename(cluster = !!sym(cl_var)) %>%
  mutate(
    time = as.numeric(as.character(time)),
    season_index = time + START_SEASON_INDEX,
    season_num = (season_index %% 4) + 1,
    season_effect = factor(season_num, levels = 1:4, labels = season_labels)
  ) %>%
  group_by(cluster, season_index, season_effect) %>%
  summarise(
    surv = mean(estimate, na.rm = TRUE),
    lower_surv = mean(lcl, na.rm = TRUE),
    upper_surv = mean(ucl, na.rm = TRUE),
    .groups = "drop"
  )

# ----- Extract capture probability (p) -----
p_real_df <- get.real(best_model_obj$model, "p", se = TRUE, vcv = FALSE)

final_p_results <- p_real_df %>%
  rename(cluster = !!sym(cl_var)) %>%
  mutate(
    time = as.numeric(as.character(time)),
    season_index = time + START_SEASON_INDEX,
    season_num = (season_index %% 4) + 1,
    season_effect = factor(season_num, levels = 1:4, labels = season_labels)
  ) %>%
  group_by(cluster, season_index, season_effect) %>%
  summarise(
    cap_prob = mean(estimate, na.rm = TRUE),
    lower_p = mean(lcl, na.rm = TRUE),
    upper_p = mean(ucl, na.rm = TRUE),
    .groups = "drop"
  )

# ===================================================================
# 8. Apply Emigration Correction (For Phi Only)
# ===================================================================

if (nrow(final_phi_results) > 0) {
  
  all_results_corrected <- final_phi_results %>%
    mutate(cluster = as.character(cluster)) %>% 
    left_join(DF_EMIGRATION, by = setNames(c(CLUSTER_COL_NAME, "season_label"), c("cluster", "season_effect")))
  
  if (any(is.na(all_results_corrected$fidelity))) {
    missing_count <- sum(is.na(all_results_corrected$fidelity))
    warning(paste0("Warning: Emigration rates missing for ", missing_count, " rows."))
    all_results_corrected$fidelity[is.na(all_results_corrected$fidelity)] <- 1
  }
  
  all_results_best_models <- all_results_corrected %>%
    mutate(
      phi_apparent = surv,
      surv = pmin(surv / fidelity, 1.0),
      lower_surv = pmin(lower_surv / fidelity, 1.0),
      upper_surv = pmin(upper_surv / fidelity, 1.0)
    ) %>%
    arrange(cluster, season_index) %>%
    group_by(cluster) %>%
    mutate(
      cum_surv = cumprod(surv),
      cum_lower = cumprod(lower_surv),
      cum_upper = cumprod(upper_surv)
    ) %>%
    ungroup()
  
  message("Correction for emigration applied.")
  
} else {
  all_results_best_models <- data.frame()
}

# ===================================================================
# 9. Prepare Data for Graphs & Generate Plots
# ===================================================================

if (nrow(all_results_best_models) > 0) {
  all_results_best_models$cluster <- as.factor(all_results_best_models$cluster)
}

# index 2: Au-Wi, 3: Wi-Sp, 4: Sp-Su, 5: Su-Au
age0_annual_data <- all_results_best_models %>%
  filter(season_index %in% c(2, 3, 4, 5)) %>%
  group_by(cluster) %>%
  summarise(
    annual_surv_rate = prod(surv),
    annual_lower = prod(lower_surv),
    annual_upper = prod(upper_surv),
    .groups = "drop"
  )

# index 6: Au-Wi, 7: Wi-Sp, 8: Sp-Su, 9: Su-Au
age1_annual_data <- all_results_best_models %>%
  filter(season_index %in% c(6, 7, 8, 9)) %>%
  group_by(cluster) %>%
  summarise(
    annual_surv_rate = prod(surv),
    annual_lower = prod(lower_surv),
    annual_upper = prod(upper_surv),
    .groups = "drop"
  )

my_palette <- "Set1"

p_age0_annual <- ggplot(age0_annual_data, aes(x = cluster, y = annual_surv_rate, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8, color = "black") + 
  geom_errorbar(aes(ymin = annual_lower, ymax = annual_upper), width = 0.2, linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", annual_surv_rate)), 
            vjust = -2.0, size = 10, fontface = "bold") +
  scale_fill_brewer(palette = my_palette) +
  labs(
    title = "Annual Survival: Age 0-Aut to Age 1-Aut",
    subtitle = paste0("Model: ", best_model_name),
    x = "Cluster",
    y = "Survival Probability"
  ) +
  coord_cartesian(ylim = c(0, max(age0_annual_data$annual_upper) * 1.2)) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 20)),
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
    axis.text.x = element_text(size = 30, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

p_age1_annual <- ggplot(age1_annual_data, aes(x = cluster, y = annual_surv_rate, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8, color = "black") + 
  geom_errorbar(aes(ymin = annual_lower, ymax = annual_upper), width = 0.2, linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", annual_surv_rate)), 
            vjust = -2.0, size = 10, fontface = "bold") +
  scale_fill_brewer(palette = my_palette) +
  labs(
    title = "Annual Survival: Age 1-Aut to Age 2-Aut",
    subtitle = paste0("Model: ", best_model_name),
    x = "Cluster",
    y = "Survival Probability"
  ) +
  coord_cartesian(ylim = c(0, max(age1_annual_data$annual_upper) * 1.2)) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 20)),
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
    axis.text.x = element_text(size = 30, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
# ===================================================================
# 平均寿命（Expected Lifespan）の算出とグラフ作成（エラーバー付き）
# ===================================================================
# 既存の cum_surv, cum_lower, cum_upper の総和 × 0.25 で期待生存年数と近似的な信頼区間を計算
lifespan_data <- all_results_best_models %>%
  group_by(cluster) %>%
  summarise(
    mean_lifespan_years = sum(cum_surv) / 4,
    lifespan_lower = sum(cum_lower) / 4,
    lifespan_upper = sum(cum_upper) / 4,
    .groups = "drop"
  )

# 平均寿命のグラフ作成
p_lifespan <- ggplot(lifespan_data, aes(x = cluster, y = mean_lifespan_years, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8, color = "black") + 
  # エラーバーの追加
  geom_errorbar(aes(ymin = lifespan_lower, ymax = lifespan_upper), width = 0.2, linewidth = 1.2) +
  # エラーバーと重ならないように vjust を -2.0 に調整して数値を上へ配置
  geom_text(aes(label = sprintf("%.2f", mean_lifespan_years)), 
            vjust = -2.0, size = 10, fontface = "bold") +
  scale_fill_brewer(palette = my_palette) +
  labs(
    title = "Mean Lifespan from First Capture",
    subtitle = paste0("Model: ", best_model_name),
    x = "Cluster",
    y = "Expected Lifespan (Years)"
  ) +
  # Y軸の最大値をエラーバーの上限（lifespan_upper）に合わせて拡張
  coord_cartesian(ylim = c(0, max(lifespan_data$lifespan_upper) * 1.2)) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 20, hjust = 0.5, margin = margin(b = 20)),
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
    axis.text.x = element_text(size = 30, color = "black"),
    axis.text.y = element_text(size = 30, color = "black"),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
# ===================================================================
# 10. Export
# ===================================================================

pdf_filename <- paste0("Survival_Analysis_", FILE_SUFFIX, ".pdf")
pdf(pdf_filename, width = 11, height = 8.5)
print(p_age0_annual) 
print(p_age1_annual)
print(p_lifespan)
dev.off()
message(paste0("PDF exported: ", getwd(), "/", pdf_filename))

csv_filename <- paste0("survival_data_", FILE_SUFFIX, ".csv")
write_csv(all_results_best_models, csv_filename)
message(paste0("CSV exported: ", getwd(), "/", csv_filename))

message("All files exported successfully.")

# Print result table to console
cat("\n=== Final Capture Probabilities (p) ===\n")
print(head(final_p_results, 10))

