# ===================================================================
# Full R Script: CJS Analysis (Survival & Capture Probability)
# with Season, Age, Cluster, and COHORT Effects
# Includes Correction for Known Emigration Rates
# ===================================================================

# ===== 0. Environment Reset =====
rm(list = ls(all.names = TRUE))
gc()

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

CLUSTER_COL_NAME   <- "cluster_k9"        
FILE_SUFFIX        <- "bkt9"  
START_SEASON_INDEX <- 2

library(tidyverse)

# ===================================================================
# 1.5. Input: Known Emigration Rates
# ===================================================================
EMIGRATION_FILE_PATH <- "../Movingout/moveout_by_cluster_k9_bkt.csv"

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

print(head(DF_EMIGRATION))

# ===================================================================
# 2. Load & Preprocess Data
# ===================================================================

df_raw <- read_csv("../Clustering/WB_bkt_clustered_rev2.csv")
df_raw$date <- as.Date(df_raw$date)

# ===== 3. Calculate SEASONS (with Winter Correction) =====
df <- df_raw %>%
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
  filter(rowSums(select(., starts_with("s"))) > 0) %>%
  mutate(cohort = as.factor(cohort))

message(paste0("Data prepared. Analysis starts from s", START_SEASON_INDEX, "."))
message(paste0("N individuals: ", nrow(capture_history)))

# ===================================================================
# 5. CJS Model Function (Extracts BOTH Phi and p)
# ===================================================================

run_and_extract_cjs_cohort <- function(data, cluster_name, phi_formula, p_formula, start_offset) {
  ch_data <- data
  
  cols_to_paste <- names(ch_data)[!names(ch_data) %in% c("tag", CLUSTER_COL_NAME, "cohort")]
  ch_data$ch <- apply(ch_data[, cols_to_paste], 1, paste0, collapse = "")
  ch_data <- ch_data %>% select(tag, ch, cohort)
  
  if (nrow(ch_data) < 2 || nchar(ch_data$ch[1]) <= 1) return(NULL)
  
  # --- process.data ---
  processed <- process.data(ch_data, model = "CJS", 
                            groups = "cohort", 
                            begin.time = 0, 
                            time.intervals = rep(1, nchar(ch_data$ch[1]) - 1))
  
  # --- make.design.data ---
  design <- make.design.data(processed)
  
  design$Phi$time <- as.numeric(as.character(design$Phi$time))
  design$p$time   <- as.numeric(as.character(design$p$time))
  
  real_time_phi <- design$Phi$time + start_offset
  real_time_p   <- design$p$time + start_offset
  
  design$Phi$age_group <- factor(
    ifelse(real_time_phi < 4, "Age0",
           ifelse(real_time_phi < 8, "Age1", "Age2+")),
    levels = c("Age1", "Age2+", "Age0")
  )
  
  season_labels <- c("1:Sp-Su", "2:Su-Au", "3:Au-Wi", "4:Wi-Sp")
  design$Phi$season_effect <- factor((real_time_phi %% 4) + 1, labels = season_labels)
  design$p$season_effect   <- factor((real_time_p %% 4) + 1, labels = season_labels)
  
  design$Phi[[CLUSTER_COL_NAME]] <- cluster_name
  design$p[[CLUSTER_COL_NAME]]   <- cluster_name
  
  # --- Model Execution ---
  model <- mark(processed, design,
                model.parameters = list(Phi = list(formula = phi_formula),
                                        p = list(formula = p_formula)),
                silent = TRUE, output = FALSE, delete = TRUE)
  
  # --- Extract Phi (Survival) ---
  phi_raw <- model$results$real %>%
    tibble::rownames_to_column("param") %>%
    filter(grepl("Phi", param)) %>%
    mutate(cluster = cluster_name)
  
  # --- Extract p (Capture Probability) ---
  p_raw <- model$results$real %>%
    tibble::rownames_to_column("param") %>%
    filter(grepl("p", param)) %>%
    mutate(cluster = cluster_name)
  
  phi_with_covs <- tryCatch({
    idx <- as.numeric(rownames(phi_raw))
    covs <- model$design.data$Phi[idx, c("time", "group", "cohort", "age_group", "season_effect")]
    res <- cbind(phi_raw, covs)
    res$time <- as.numeric(as.character(res$time))
    res$season_index <- res$time + start_offset
    res
  }, error = function(e) {
    phi_raw$season_index <- (1:nrow(phi_raw) - 1) %% (nchar(ch_data$ch[1]) - 1) + start_offset
    return(phi_raw)
  })
  
  p_with_covs <- tryCatch({
    idx_p <- as.numeric(rownames(p_raw))
    covs_p <- model$design.data$p[idx_p, c("time", "group", "cohort", "season_effect")]
    res_p <- cbind(p_raw, covs_p)
    res_p$time <- as.numeric(as.character(res_p$time))
    res_p$season_index <- res_p$time + start_offset
    res_p
  }, error = function(e) {
    p_raw$season_index <- (1:nrow(p_raw) - 1) %% (nchar(ch_data$ch[1]) - 1) + start_offset
    return(p_raw)
  })
  
  return(list(phi = phi_with_covs, p = p_with_covs, model = model))
}

# ===================================================================
# 6. Run Analysis
# ===================================================================

cluster_list <- sort(unique(capture_history[[CLUSTER_COL_NAME]]))

model_defs <- list(
  list(name = "Phi_Age_Season_p_Season",  phi_formula = ~age_group + season_effect, p_formula = ~season_effect),
  list(name = "Phi_Cluster_p_Season",      phi_formula = as.formula(paste0("~", CLUSTER_COL_NAME)), p_formula = ~season_effect),
  list(name = "Phi_Age_Season_Cohort_p_Season", phi_formula = ~age_group + season_effect + cohort, p_formula = ~season_effect)
)

best_models_phi_data <- list()
best_models_p_data   <- list() 

cat("\nProcessing Clusters...\n")

for (cl_id in cluster_list) {
  current_cluster_data <- capture_history %>% filter(!!sym(CLUSTER_COL_NAME) == cl_id)
  cat(paste0("\nCluster: ", cl_id, " (n=", nrow(current_cluster_data), ")\n"))
  
  cluster_results <- list()
  
  for (model_def in model_defs) {
    tryCatch({
      result <- run_and_extract_cjs_cohort(
        current_cluster_data, cl_id, 
        model_def$phi_formula, model_def$p_formula, 
        start_offset = START_SEASON_INDEX
      )
      
      if (!is.null(result)) {
        cluster_results[[length(cluster_results) + 1]] <- list(
          name = model_def$name,
          AICc = result$model$results$AICc,
          phi = result$phi,
          p = result$p
        )
        cat(".") 
      }
    }, error = function(e) cat("x"))
  }
  
  if (length(cluster_results) > 0) {
    df_aic <- bind_rows(lapply(cluster_results, function(x) data.frame(name = x$name, AICc = x$AICc))) %>% arrange(AICc)
    print(head(df_aic, 1))
    
    best_res <- Filter(function(x) x$name == df_aic$name[1], cluster_results)[[1]]
    
    season_labels <- c("1:Sp-Su", "2:Su-Au", "3:Au-Wi", "4:Wi-Sp")
    
    # --- Aggregate Phi ---
    phi_agg <- best_res$phi %>%
      mutate(
        season_num = (season_index %% 4) + 1,
        season_effect = factor(season_num, levels = 1:4, labels = season_labels)
      ) %>%
      group_by(cluster, season_index, season_effect) %>%
      summarise(
        surv = mean(if("estimate" %in% names(.)) estimate else surv, na.rm=TRUE),
        lower_surv = mean(if("lcl" %in% names(.)) lcl else lower_surv, na.rm=TRUE),
        upper_surv = mean(if("ucl" %in% names(.)) ucl else upper_surv, na.rm=TRUE),
        .groups = "drop"
      )
    best_models_phi_data[[length(best_models_phi_data) + 1]] <- phi_agg
    
    # --- Aggregate p (Capture Prob) ---
    p_agg <- best_res$p %>%
      mutate(
        season_num = (season_index %% 4) + 1,
        season_effect = factor(season_num, levels = 1:4, labels = season_labels)
      ) %>%
      group_by(cluster, season_index, season_effect) %>%
      summarise(
        cap_prob = mean(if("estimate" %in% names(.)) estimate else p, na.rm=TRUE),
        lower_p = mean(if("lcl" %in% names(.)) lcl else lower_p, na.rm=TRUE),
        upper_p = mean(if("ucl" %in% names(.)) ucl else upper_p, na.rm=TRUE),
        .groups = "drop"
      )
    best_models_p_data[[length(best_models_p_data) + 1]] <- p_agg
  }
}

final_phi_results <- bind_rows(best_models_phi_data)
final_p_results   <- bind_rows(best_models_p_data)

# ===================================================================
# 7. Apply Emigration Correction (For Phi Only)
# ===================================================================

if (nrow(final_phi_results) > 0) {
  
  all_results_corrected <- final_phi_results %>%
    mutate(cluster = as.character(cluster)) %>% 
    left_join(DF_EMIGRATION, by = setNames(c(CLUSTER_COL_NAME, "season_label"), c("cluster", "season_effect")))
  
  if (any(is.na(all_results_corrected$fidelity))) {
    missing_count <- sum(is.na(all_results_corrected$fidelity))
    warning(paste0("Warning: Emigration rates are missing for ", missing_count, " rows. Check cluster IDs or season labels."))
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
# 8. Prepare Data for Graphs
# ===================================================================

if (nrow(all_results_best_models) > 0) {
  all_results_best_models$cluster <- as.factor(all_results_best_models$cluster)
}
if (nrow(final_p_results) > 0) {
  final_p_results$cluster <- as.factor(final_p_results$cluster) 
}

# Mean lifespan data
mean_life_data <- all_results_best_models %>%
  group_by(cluster) %>%
  summarise(
    mean_life_years = sum(cum_surv, na.rm = TRUE) / 4,
    lower_ci_years = sum(cum_lower, na.rm = TRUE) / 4,
    upper_ci_years = sum(cum_upper, na.rm = TRUE) / 4,
    .groups = "drop"
  )

# Seasonal survival data
seasonal_survival <- all_results_best_models %>%
  mutate(
    season_transition = factor((season_index %% 4) + 1, 
                               labels = c("1:Sp-Su", "2:Su-Au", "3:Au-Wi", "4:Wi-Sp"))
  )

# ===================================================================
# 9. Generate Plots
# ===================================================================

my_palette <- "Set1"
common_theme <- theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "#f0f0f0"), legend.position = "bottom")

# ===================================================================
# 9.6. Specific Plot: Age 0 Autumn to Age 1 Autumn Survival
# ===================================================================

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

# index 2: Au-Wi, 3: Wi-Sp, 4: Sp-Su, 5: Su-Au
age1_annual_data <- all_results_best_models %>%
  filter(season_index %in% c(6, 7, 8, 9)) %>%
  group_by(cluster) %>%
  summarise(
    annual_surv_rate = prod(surv),
    annual_lower = prod(lower_surv),
    annual_upper = prod(upper_surv),
    .groups = "drop"
  )

# ===================================================================
# 9.6. Specific Plot: Age 0 Autumn to Age 1 Autumn (Large Text)
# ===================================================================

p_age0_annual <- ggplot(age0_annual_data, aes(x = cluster, y = annual_surv_rate, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8, color = "black") + 
  geom_errorbar(aes(ymin = annual_lower, ymax = annual_upper), width = 0.2, linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", annual_surv_rate)), 
            vjust = -2.0, size = 10, fontface = "bold") +
  scale_fill_brewer(palette = my_palette) +
  labs(
    title = "Annual Survival: Age 0-Aut to Age 1-Aut",
    subtitle = "Calculated from True Survival Rates",
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

print(p_age0_annual)

p_age1_annual <- ggplot(age1_annual_data, aes(x = cluster, y = annual_surv_rate, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8, color = "black") + 
  geom_errorbar(aes(ymin = annual_lower, ymax = annual_upper), width = 0.2, linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.3f", annual_surv_rate)), 
            vjust = -2.0, size = 10, fontface = "bold") +
  scale_fill_brewer(palette = my_palette) +
  labs(
    title = "Annual Survival: Age 1-Aut to Age 2-Aut",
    subtitle = "Calculated from True Survival Rates",
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

print(p_age1_annual)

# ===================================================================
# 10. Export
# ===================================================================

# PDF Export
pdf_filename <- paste0("Survival_Analysis_", FILE_SUFFIX, ".pdf")
pdf(pdf_filename, width = 11, height = 8.5)
print(p_age0_annual) 
print(p_age1_annual)
dev.off()
message(paste0("PDF exported: ", getwd(), "/", pdf_filename))

# CSV Exports
write_csv(all_results_best_models, paste0("survival_data_", FILE_SUFFIX, ".csv"))

message("All files exported successfully.")

