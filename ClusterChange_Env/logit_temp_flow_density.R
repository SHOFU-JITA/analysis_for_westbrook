# =============================================================================
# Westbrook Project: Growth Strategy Cluster Analysis (Full Script)
# Updated: Automated Model Generation with True Density (per m^2) & Total Density
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

# --- River Area Settings (m^2) ---
# Please update these values to the actual surface areas of each river/site
RIVER_AREAS <- c(
  wb = 7096,  # West Brook
  ol = 935,   # WB Jimmy (Open Large)
  os = 598,   # WB Mitchell (Open Small)
  il = 893    # WB Obear (Isolated Large)
)

file_path              <- "csv/WB_combined_bkt_rev2.csv"
temp_file_path         <- "rawdata/envDataWB.csv"
flow_file_path         <- "Flow/hindcast/GR5H_Simulated_Discharge_WB_2000_2025.csv"

# 追加: 各魚種の密度データのパス
density_file_path_bkt  <- "JS/Abundance_Estimates_S3_Total_bkt.csv"
density_file_path_bnt  <- "JS/Abundance_Estimates_S3_Total_bnt.csv"
density_file_path_ats  <- "JS/Abundance_Estimates_S3_Total_ats.csv"

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}
cat(paste0("Save directory: ", SCRIPT_DIR, "\n"))

num_clusters <- 3
CLUSTER_COLORS <- c("steelblue", "aquamarine3", "salmon3")

pdf_filename <- paste0("logit_env_flow_bkt_dens", num_clusters,".pdf")
pdf_full_path <- file.path(SCRIPT_DIR, pdf_filename)

# --- Manual Plot Settings (High/Mean/Low) ---
# Note: abundance values are now true densities (individuals / m^2)
MANUAL_PLOT_SETTINGS <- list(
  temp_summer = c(13, 15, 17),  
  temp_winter = NULL,
  temp_spring = NULL,
  temp_autumn = NULL,
  
  density_s3        = NULL, 
  cohort_density    = NULL,
  dens_total        = NULL # 追加
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
      month %in% 3:5 ~ "temp_spring",
      month %in% 6:8 ~ "temp_summer",
      month %in% 9:11 ~ "temp_autumn",
      month %in% c(12, 1, 2) ~ "temp_winter",
      TRUE ~ NA_character_
    ),
    season_year = if_else(month == 12, year + 1, year)
  ) %>%
  filter(!is.na(season)) %>%
  group_by(river, season_year, season) %>%
  summarise(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = season, values_from = mean_temp)


# --- B. Process Density Data (Converted to True Density per m^2 & Total Density) ---
# 個別の魚種のファイルを読み込む補助関数
read_dens <- function(path, col_name) {
  read.csv(path) %>%
    filter(grepl("S3", season_label)) %>%
    mutate(
      season_year = as.numeric(str_sub(season_label, 1, 4)),
      river = tolower(river)
    ) %>%
    dplyr::select(river, season_year, !!col_name := N_est)
}

# 3魚種の読み込み
raw_dens_bkt <- read_dens(density_file_path_bkt, "N_bkt")
raw_dens_bnt <- read_dens(density_file_path_bnt, "N_bnt")
raw_dens_ats <- read_dens(density_file_path_ats, "N_ats")

# 全魚種を結合し、合計密度を計算
dens_processed <- raw_dens_bkt %>%
  full_join(raw_dens_bnt, by = c("river", "season_year")) %>%
  full_join(raw_dens_ats, by = c("river", "season_year")) %>%
  mutate(
    # NAの場合は0とする
    N_bkt = replace_na(N_bkt, 0),
    N_bnt = replace_na(N_bnt, 0),
    N_ats = replace_na(N_ats, 0),
    # 全個体数
    N_total = N_bkt + N_bnt + N_ats,
    
    # それぞれの密度を計算 (per m^2)
    density_s3 = N_bkt / unname(RIVER_AREAS[river]),
    dens_total = N_total / unname(RIVER_AREAS[river])
  ) %>%
  dplyr::select(river, season_year, density_s3, dens_total)


# --- C. Calculate Observed Cohort Densities (Converted to True Density per m^2) ---
cohort_counts <- df_analysis %>%
  group_by(river, cohort_num) %>%
  summarise(raw_cohort_count = n_distinct(tag), .groups = "drop") %>%
  mutate(
    river = tolower(river), 
    season_year = cohort_num,
    cohort_density = raw_cohort_count / unname(RIVER_AREAS[river])
  ) %>%
  dplyr::select(river, season_year, cohort_density)

# --- D. Process Flow Data ---
raw_flow <- read.csv(flow_file_path)

# 1. 1時間ごとのデータを日次（Daily）に集計
daily_flow <- raw_flow %>%
  mutate(
    Date = as.Date(substr(Date_Hour, 1, 10))
  ) %>%
  group_by(Site_Name, Date) %>%
  summarise(Daily_Q = mean(Q_sim_m3s, na.rm = TRUE), .groups = "drop")

# 2. 河川ごとに全期間を通した流量の5%タイルと95%タイルの閾値を計算
flow_thresholds <- daily_flow %>%
  group_by(Site_Name) %>%
  summarise(
    Q_05 = quantile(Daily_Q, 0.05, na.rm = TRUE),
    Q_95 = quantile(Daily_Q, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# 3. 季節区分の追加と集計
flow_processed <- daily_flow %>%
  left_join(flow_thresholds, by = "Site_Name") %>%
  mutate(
    year = year(Date),
    month = month(Date),
    season_label = case_when(
      month %in% 3:5 ~ "spring",
      month %in% 6:8 ~ "summer",
      month %in% 9:11 ~ "autumn",
      month %in% c(12, 1, 2) ~ "winter",
      TRUE ~ NA_character_
    ),
    season_year = if_else(month == 12, year + 1, year),
    is_high = if_else(Daily_Q >= Q_95, 1, 0),
    is_low  = if_else(Daily_Q <= Q_05, 1, 0)
  ) %>%
  group_by(Site_Name, season_year, season_label) %>%
  summarise(
    Mean_Flow      = mean(Daily_Q, na.rm = TRUE),
    High_Flow_Days = sum(is_high, na.rm = TRUE),
    Low_Flow_Days  = sum(is_low, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    river = case_when(
      tolower(Site_Name) %in% c("wb jimmy", "ol") ~ "ol",
      tolower(Site_Name) %in% c("wb mitchell", "os") ~ "os",
      tolower(Site_Name) %in% c("wb obear", "il") ~ "il",
      tolower(Site_Name) %in% c("west brook", "wb") ~ "wb",
      TRUE ~ tolower(Site_Name)
    )
  ) %>%
  pivot_wider(
    id_cols = c(river, season_year),
    names_from = season_label,
    values_from = c(Mean_Flow, High_Flow_Days, Low_Flow_Days),
    names_glue = "{.value}_{season_label}"
  )

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
  "temp_summer", "temp_winter", "temp_spring", "temp_autumn", 
  "density_s3", "dens_total", "cohort_density", # dens_totalを追加
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
# 4. Model Analysis (Automated Formula Generation - All Pairwise Interactions)
# =============================================================================

# Scale numeric variables
scale_cols <- df_model_env %>% select(where(is.numeric), -cohort_numeric, -join_year) %>% names()
df_model_env_scaled <- df_model_env %>%
  mutate(across(all_of(scale_cols), ~ as.numeric(scale(.))))

cat("\n--- Generating Formulas (Full model with all pairwise interactions) ---\n")

formulas_mixed <- list()
formulas_mixed[["Null Model"]] <- ocluster_rename_dynamic ~ 1

seasons <- c("spring", "summer", "autumn", "winter")
flow_metrics <- c("High_Flow_Days", "Low_Flow_Days")
abund_vars <- c("density_s3", "cohort_density") 

for (s in seasons) {
  t_var <- paste0("temp_", s)
  
  if (t_var %in% names(df_model_env_scaled)) {
    
    for (fm in flow_metrics) {
      f_var <- paste0(fm, "_", s)
      
      if (f_var %in% names(df_model_env_scaled)) {
        
        for (a_var in abund_vars) {
          if (a_var %in% names(df_model_env_scaled)) {
            
            # 変数名のセット
            vars <- c(t_var, f_var, a_var)
            
            # --- 階層性(Marginality)を満たす全サブセットモデルの生成 ---
            # 1. 1変数モデル
            formulas_mixed[[paste0(t_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var))
            formulas_mixed[[paste0(f_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", f_var))
            formulas_mixed[[paste0(a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", a_var))
            
            # 2. 2変数モデル (主効果のみ)
            formulas_mixed[[paste0(t_var, "+", f_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "+", f_var))
            formulas_mixed[[paste0(t_var, "+", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "+", a_var))
            formulas_mixed[[paste0(f_var, "+", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", f_var, "+", a_var))
            
            # 3. 3変数モデル (主効果のみ)
            formulas_mixed[[paste0(t_var, "+", f_var, "+", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "+", f_var, "+", a_var))
            
            # 4. 2変数 + その交互作用
            formulas_mixed[[paste0(t_var, "*", f_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", f_var))
            formulas_mixed[[paste0(t_var, "*", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", a_var))
            formulas_mixed[[paste0(f_var, "*", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", f_var, "*", a_var))
            
            # 5. 3変数 + 1つの2元交互作用
            formulas_mixed[[paste0("(", t_var, "*", f_var, ")+", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", f_var, "+", a_var))
            formulas_mixed[[paste0("(", t_var, "*", a_var, ")+", f_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", a_var, "+", f_var))
            formulas_mixed[[paste0("(", f_var, "*", a_var, ")+", t_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", f_var, "*", a_var, "+", t_var))
            
            # 6. 3変数 + 2つの2元交互作用
            formulas_mixed[[paste0(t_var, "*", f_var, "+", t_var, "*", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", f_var, "+", t_var, "*", a_var))
            formulas_mixed[[paste0(t_var, "*", f_var, "+", f_var, "*", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", f_var, "+", f_var, "*", a_var))
            formulas_mixed[[paste0(t_var, "*", a_var, "+", f_var, "*", a_var)]] <- as.formula(paste("ocluster_rename_dynamic ~", t_var, "*", a_var, "+", f_var, "*", a_var))
            
            # 7. フルモデル (3変数 + 全ての2元交互作用 = pairwise interactions)
            formulas_mixed[[paste0("Full Pairwise (", s, "_", fm, "_", a_var, ")")]] <- as.formula(paste("ocluster_rename_dynamic ~ (", t_var, "+", f_var, "+", a_var, ")^2"))
          }
        }
      }
    }  
  }
}

# 重複するモデル（名前違いで同じ構成のもの）を削除
formulas_mixed <- formulas_mixed[!duplicated(sapply(formulas_mixed, function(x) paste(deparse(x), collapse = "")))]

cat(paste0("Generated ", length(formulas_mixed), " unique model formulas.\n"))

# --- Step 2: Model Comparison (AIC & R2) ---
get_ll_safe <- function(mod) {
  if (is.null(mod)) return(NA)
  ll <- try(as.numeric(logLik(mod)), silent = TRUE)
  if (!inherits(ll, "try-error") && length(ll) > 0 && !is.na(ll)) return(ll)
  if (!is.null(mod$deviance)) return(-0.5 * mod$deviance)
  return(NA)
}

# Nullモデル（ベースライン）の計算
cat("Calculating Null Model for McFadden R2...\n")
null_fit <- tryCatch(
  mclogit::mblogit(ocluster_rename_dynamic ~ 1, 
                   random = ~ 1 | river_grp, 
                   data = df_model_env_scaled, 
                   maxit = 50),
  error = function(e) return(NULL)
)
logLik_null <- get_ll_safe(null_fit)

# 各モデルの評価
model_comparison <- map_dfr(names(formulas_mixed), function(name) {
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
    # McFaddenの疑似R2を計算
    r2_mcfadden <- 1 - (ll_model / logLik_null)
    
    current_aic <- tryCatch(AIC(fit), error = function(e) {
      k <- length(coef(fit)) 
      return(2 * k - 2 * ll_model)
    })
    
    tibble(
      Model_Name  = name, 
      AIC         = current_aic, 
      McFadden_R2 = round(r2_mcfadden, 4), # ここを復活
      Status      = "Success"
    )
  } else {
    tibble(Model_Name = name, AIC = NA, McFadden_R2 = NA, Status = "Failed")
  }
})

model_ranking <- model_comparison %>%
  filter(Status == "Success") %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    is_best_by_2 = (lead(AIC) - AIC) >= 2
  )

cat("\n--- Top 10 Models by AIC ---\n")
print(head(model_ranking, 10))

# ΔAIC >= 2 のチェック結果の出力
if (nrow(model_ranking) > 1) {
  best_model_diff <- model_ranking$AIC[2] - model_ranking$AIC[1]
  if(best_model_diff >= 2) {
    cat(sprintf("\n[CONCLUSION] The top model is clearly the best (ΔAIC to 2nd best = %.2f >= 2).\n", best_model_diff))
  } else {
    cat(sprintf("\n[CONCLUSION] The top model is NOT strictly better by 2 AIC points (ΔAIC to 2nd best = %.2f < 2).\n", best_model_diff))
    cat("Multiple competing models may need to be considered (Model Averaging or reporting top models).\n")
  }
}

csv_filename <- paste0("Modelselection_bkt_flow_dens", num_clusters, ".csv")
write.csv(model_ranking, file.path(SCRIPT_DIR, csv_filename), row.names = FALSE)
# =============================================================================
# 5. Top 1 Model Detailed Analysis & Visualization
# =============================================================================

pdf(pdf_full_path, width = 12, height = 8)

if (nrow(model_ranking) > 0) {
  
  top_n_models <- min(1, nrow(model_ranking))
  
  for (m in 1:top_n_models) {
    
    best_model_name <- model_ranking$Model_Name[m]
    best_formula    <- formulas_mixed[[best_model_name]]
    
    cat("\n--- Plotting Best Model ---\n")
    cat(best_model_name, "\n")
    best_model <- mclogit::mblogit(
      best_formula,
      random = ~ 1 | river_grp,
      data   = df_model_env_scaled,
      maxit  = 50
    )
    
    # =========================================================================
    # 自動チェック: 過学習・多重共線性・収束異常
    # =========================================================================
    coef_vals <- coef(best_model)
    se_vals   <- tryCatch(sqrt(diag(vcov(best_model))), error = function(e) rep(NA, length(coef_vals)))
    
    warning_msgs <- c()
    
    if (any(abs(coef_vals) > 10, na.rm = TRUE)) {
      warning_msgs <- c(warning_msgs, "・係数 (Estimate) の絶対値が10を超えているパラメータがあります。過学習や完全分離の可能性が高いです。")
    }
    if (any(se_vals > 10, na.rm = TRUE)) {
      warning_msgs <- c(warning_msgs, "・標準誤差 (Std. Error) が非常に大きい（>10）パラメータがあります。強い多重共線性やモデルの不安定性が疑われます。")
    }
    if (is.null(vcov(best_model)) || any(is.na(se_vals))) {
      warning_msgs <- c(warning_msgs, "・標準誤差が計算できませんでした。モデルが正しく収束していない可能性が高いです。")
    }
    
    if (length(warning_msgs) > 0) {
      cat("\n!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!\n")
      cat("モデル [", best_model_name, "] に異常なパラメータ値が検出されました：\n")
      for (msg in warning_msgs) cat(msg, "\n")
      cat("=> 見かけ上AICが良くても、グラフの予測確率が極端になるなど信頼性に欠ける可能性があります。\n")
      cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
    }
    
    predictors <- all.vars(best_formula)[-1]
    
    if (length(predictors) == 0) {
      cat("Null model: Skipping visualization\n")
    } else {
      plot_combinations <- list()
      
      temp_var  <- predictors[grepl("temp", predictors)][1]
      abund_var <- predictors[predictors %in% abund_vars][1]
      flow_var  <- predictors[grepl("Flow", predictors)][1]
      
      if (length(predictors) == 1) {
        plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = "Dummy")
      } else {
        idx <- 1
        
        # Temp vs Density 
        if (!is.na(temp_var) && !is.na(abund_var)) {
          plot_combinations[[idx]] <- list(x_var = temp_var, facet_var = abund_var) 
          idx <- idx + 1
          plot_combinations[[idx]] <- list(x_var = abund_var, facet_var = temp_var) 
          idx <- idx + 1
        }
        
        # Flow vs Density / Temp
        if (!is.na(flow_var) && !is.na(abund_var)) {
          plot_combinations[[idx]] <- list(x_var = flow_var, facet_var = abund_var) 
          idx <- idx + 1
        }
        
        if (!is.na(flow_var) && !is.na(temp_var)) {
          plot_combinations[[idx]] <- list(x_var = flow_var, facet_var = temp_var)  
          idx <- idx + 1
        }
        
        if (length(plot_combinations) == 0 && length(predictors) >= 2) {
          plot_combinations[[1]] <- list(x_var = predictors[1], facet_var = predictors[2])
          plot_combinations[[2]] <- list(x_var = predictors[2], facet_var = predictors[1])
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
        
        # プロット範囲を5%〜95%に制限
        x_quantiles <- quantile(df_model_env[[x_var]], probs = c(0.05, 0.95), na.rm = TRUE)
        x_seq_real <- seq(
          x_quantiles[1],
          x_quantiles[2],
          length.out = 100
        )
        
        if (facet_var != "Dummy") {
          facet_stats <- get_stats(facet_var)
          
          if (!is.null(MANUAL_PLOT_SETTINGS[[facet_var]])) {
            facet_vals_real <- MANUAL_PLOT_SETTINGS[[facet_var]]
          } else {
            facet_vals_real <- c(
              quantile(df_model_env[[facet_var]], 0.2, na.rm = TRUE), 
              mean(df_model_env[[facet_var]], na.rm = TRUE),          
              quantile(df_model_env[[facet_var]], 0.8, na.rm = TRUE)  
            )
          }
          prefix_labels <- c("Low", "Mean", "High")
          if (length(facet_vals_real) == 3) {
            facet_labels <- paste0(prefix_labels, " (", round(facet_vals_real, 3), ")")
          } else {
            facet_labels <- paste0("Val (", round(facet_vals_real, 3), ")")
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
        
        for (p in predictors) {
          if (p != x_var && p != facet_var) {
            pred_grid_real[[p]] <- get_stats(p)$mean
          }
        }
        
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
        
        plot_title <- if (facet_var != "Dummy") {
          paste0(facet_var, "-", x_var)
        } else {
          x_var
        }
        
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
            title = plot_title,
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

