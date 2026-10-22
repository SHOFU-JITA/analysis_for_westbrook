# =============================================================================
# Data Preparation and GR5H + CemaNeige Hydrological Modeling (Hourly)
# Focus: 2020-2025 Calibration (KGE) / 2000-2025 Simulation (Hindcast)
# =============================================================================

# 1. Initialize environment, load packages, and set working directory
rm(list = ls(all.names = TRUE))
library(tidyverse)
library(lubridate)
library(rstudioapi)
library(airGR)
library(zoo)
library(patchwork)
library(scales)

if (rstudioapi::isAvailable()) {
  try({
    current_path <- rstudioapi::getSourceEditorContext()$path
    if (current_path != "") {
      setwd(dirname(current_path))
      message("Working directory set to: ", getwd())
    }
  }, silent = TRUE)
}

output_filename <- "GR5H_Simulated_Discharge_WB_2000_2025_summer.csv"
pdf_output <- "GR5H_Calibration_Plots_BySite.pdf"

# 2. Define site dictionaries, catchment areas (sq km), and site-specific parameters
dict_eco <- tibble(
  raw_name = c("Jimmy Brook", "Obear Brook  Low", "Mitchell Brook", "West Brook 1 (Lower)"),
  unified_name = c("OL", "IL", "OS", "WB") 
)

site_params <- tibble(
  Site_Name = c("OL", "IL", "OS", "WB"),
  Area_sqkm = c(2.46, 1.33, 1.06, 21.78),
  snow_catch_factor = c(1.1, 1.15, 1.1, 1.1),
  et_multiplier = c(0.95, 0.9, 0.85, 0.9)
)

# =============================================================================
# 3. Load Data & 4. Format/Standardize (サマータイム補正のみ適用)
# =============================================================================
eco_draught_data <- read_csv("EcoDrought_Continuous_MA.csv", show_col_types = FALSE)

eco_hourly <- eco_draught_data %>%
  left_join(dict_eco, by = c("Site_Name" = "raw_name")) %>%
  filter(!is.na(unified_name)) %>%
  mutate(
    Site_Name = unified_name,
    Date_Hour = floor_date(as.POSIXct(DateTime_EST, tz = "EST"), "hour") 
  ) %>%
  group_by(Site_Name, Date_Hour) %>%
  summarize(
    mean_discharge = mean(Discharge_Hobo_cfs, na.rm = TRUE),
    .groups = "drop"
  )

site_list <- c("OL", "IL", "OS", "WB")
weather_list <- list()

for (site in site_list) {
  file_name <- paste0(site, "_hourly_weather.csv")
  
  if (file.exists(file_name)) {
    temp_df <- read_csv(file_name, show_col_types = FALSE) %>%
      mutate(
        Site_Name = site,
        # 1. Open-Meteoの時刻(America/New_York: 夏時間あり)としてパース
        Time_NY = as.POSIXct(as.character(Time), format = "%Y-%m-%d %H:%M", tz = "America/New_York"),
        Time_NY = if_else(
          is.na(Time_NY),
          as.POSIXct(as.character(Time), tz = "America/New_York"),
          Time_NY
        ),
        # 2. サマータイムをESTに変換した上で、
        #    Open-Meteoの「前1時間積算」を流量データ(floor_date)のタイミングに合わせるために「- hours(1)」を適用
        Date_Hour = with_tz(Time_NY, tzone = "EST") - hours(1)
      ) %>%
      filter(!is.na(Date_Hour)) %>%
      rename(Temp_C = Temperature_C) %>%
      # サマータイム終了時の重複時刻を平均・合算
      group_by(Site_Name, Date_Hour) %>%
      summarize(
        Temp_C = mean(Temp_C, na.rm = TRUE),
        Precipitation_mm = sum(Precipitation_mm, na.rm = TRUE),
        .groups = "drop"
      )
    
    weather_list[[site]] <- temp_df
  } else {
    warning("File not found: ", file_name)
  }
}
weather_hourly <- bind_rows(weather_list)

# =============================================================================
# 5. Merge datasets, calculate Diurnal PET, and Apply Snow Undercatch
# =============================================================================
min_time <- min(weather_hourly$Date_Hour, eco_hourly$Date_Hour, na.rm = TRUE)
max_time <- max(weather_hourly$Date_Hour, eco_hourly$Date_Hour, na.rm = TRUE)
full_time_grid <- expand_grid(
  Site_Name = site_list,
  Date_Hour = seq(min_time, max_time, by = "hour")
)

master_data <- full_time_grid %>%
  left_join(weather_hourly, by = c("Site_Name", "Date_Hour")) %>%
  left_join(eco_hourly, by = c("Site_Name", "Date_Hour")) %>%
  arrange(Site_Name, Date_Hour) %>%
  group_by(Site_Name) %>%
  mutate(
    Temp_C = na.approx(Temp_C, na.rm = FALSE, rule = 2),
    Precipitation_mm = replace_na(Precipitation_mm, 0)
  ) %>%
  ungroup()

site_latitude <- 42.431868 

daily_temp_summary <- master_data %>%
  mutate(Date = as.Date(Date_Hour, tz = "EST")) %>%
  group_by(Site_Name, Date) %>%
  summarize(mean_daily_temp = mean(Temp_C, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Julian_Day = yday(Date),
    PET_daily = PE_Oudin(JD = Julian_Day, Temp = mean_daily_temp, Lat = site_latitude, LatUnit = "deg")
  )

master_data <- master_data %>%
  mutate(
    Date = as.Date(Date_Hour, tz = "EST"),
    Hour = hour(Date_Hour)
  ) %>%
  left_join(daily_temp_summary %>% select(Site_Name, Date, PET_daily), by = c("Site_Name", "Date")) %>%
  left_join(site_params, by = "Site_Name") %>%
  mutate(
    Precipitation_mm = if_else(Temp_C <= 0, Precipitation_mm * snow_catch_factor, Precipitation_mm),
    Diurnal_Weight = if_else(Hour >= 6 & Hour <= 18, sin((Hour - 6) * pi / 12) + 0.1, 0.1),
    PET_hourly = (PET_daily * (Diurnal_Weight / 8.839437)) * et_multiplier,
    Volume_m3_hour = mean_discharge * 0.0283168 * 3600,
    Area_m2 = Area_sqkm * 1000000,
    Q_mm_hr = (Volume_m3_hour / Area_m2) * 1000
  ) %>%
  select(-Date, -Hour, -PET_daily, -Diurnal_Weight, -snow_catch_factor, -et_multiplier)


# =============================================================================
# 6. ステップワイズ & 変換空間ハイブリッド・キャリブレーション
# =============================================================================

# --- airGRのGR5H(5変数)とCemaNeige(2変数)の変換関数を結合するラッパー関数 ---
Transfo_CemaGR5H <- function(ParamIn, Direction) {
  if (is.vector(ParamIn)) ParamIn <- matrix(ParamIn, nrow = 1)
  p_gr5h <- TransfoParam_GR5H(ParamIn = ParamIn[, 1:5, drop = FALSE], Direction = Direction)
  p_cema <- TransfoParam_CemaNeige(ParamIn = ParamIn[, 6:7, drop = FALSE], Direction = Direction)
  return(as.numeric(cbind(p_gr5h, p_cema)))
}

results_list <- list()

for (target_site in site_list) {
  
  message("==================================================")
  message("Processing site: ", target_site)
  
  df_model <- master_data %>%
    filter(Site_Name == target_site) %>%
    arrange(Date_Hour) %>%
    mutate(PET_hourly = na.approx(PET_hourly, na.rm = FALSE, rule = 2))
  
  Ind_WarmUp_Calib <- which(df_model$Date_Hour >= as.POSIXct("2018-01-01 00:00:00", tz = "EST") & 
                              df_model$Date_Hour <= as.POSIXct("2019-12-31 23:00:00", tz = "EST"))
  Ind_Run_Calib <- which(df_model$Date_Hour >= as.POSIXct("2020-01-01 00:00:00", tz = "EST") & 
                           df_model$Date_Hour <= as.POSIXct("2025-12-31 23:00:00", tz = "EST"))
  
  Ind_WarmUp_Sim <- which(df_model$Date_Hour >= as.POSIXct("1999-01-01 00:00:00", tz = "EST") & 
                            df_model$Date_Hour <= as.POSIXct("1999-12-31 23:00:00", tz = "EST"))
  Ind_Run_Sim <- which(df_model$Date_Hour >= as.POSIXct("2000-01-01 00:00:00", tz = "EST") & 
                         df_model$Date_Hour <= as.POSIXct("2025-12-31 23:00:00", tz = "EST"))
  
  if (target_site == "WB") {
    catchment_elevation <- 200
    hypso_curve <- seq(from = 150, to = 250, length.out = 101)
  } else {
    catchment_elevation <- 150 
    hypso_curve <- seq(from = 130, to = 170, length.out = 101) 
  }
  
  calib_months <- month(df_model$Date_Hour[Ind_Run_Calib])
  bool_summer_calib <- (calib_months >= 6 & calib_months <= 10)
  
  # -------------------------------------------------------------------
  # 【STEP 1】純粋なGR5H（雪なし）による夏期特化キャリブレーション (X1-X5)
  # -------------------------------------------------------------------
  InputsModel_GR5H <- CreateInputsModel(
    FUN_MOD = RunModel_GR5H, DatesR = df_model$Date_Hour, 
    Precip = df_model$Precipitation_mm, PotEvap = df_model$PET_hourly
  )
  RunOptions_Calib_GR5H <- CreateRunOptions(
    FUN_MOD = RunModel_GR5H, InputsModel = InputsModel_GR5H, 
    IndPeriod_WarmUp = Ind_WarmUp_Calib, IndPeriod_Run = Ind_Run_Calib        
  )
  InputsCrit_Step1 <- CreateInputsCrit(
    FUN_CRIT = list(ErrorCrit_KGE, ErrorCrit_KGE), InputsModel = InputsModel_GR5H, 
    RunOptions = RunOptions_Calib_GR5H,
    Obs = list(df_model$Q_mm_hr[Ind_Run_Calib], df_model$Q_mm_hr[Ind_Run_Calib]),
    transfo = list("", "sqrt"), Weights = list(0.6, 0.4), 
    BoolCrit = list(bool_summer_calib, bool_summer_calib) 
  )
  
  CalibOptions_Step1 <- CreateCalibOptions(FUN_MOD = RunModel_GR5H, FUN_CALIB = Calibration_Michel)
  OutputsCalib_Step1 <- Calibration_Michel(
    InputsModel = InputsModel_GR5H, RunOptions = RunOptions_Calib_GR5H, 
    InputsCrit = InputsCrit_Step1, CalibOptions = CalibOptions_Step1, FUN_MOD = RunModel_GR5H
  )
  Param_GR5H <- OutputsCalib_Step1$ParamFinalR
  
  # X4は0.25hを許容し、X1の極端な暴走のみ1500でクリップ
  Param_GR5H[1] <- min(max(Param_GR5H[1], 50), 1500)
  
  message("  -> Step 1 (Summer GR5H) X1-X5: ", paste(round(Param_GR5H, 3), collapse = ", "))
  
  # -------------------------------------------------------------------
  # 【STEP 2】CemaNeige+GR5H の探索（ルートA: 固定探索 vs ルートB: 通年自由探索）
  # -------------------------------------------------------------------
  InputsModel_Cema <- CreateInputsModel(
    FUN_MOD = RunModel_CemaNeigeGR5H, DatesR = df_model$Date_Hour, 
    Precip = df_model$Precipitation_mm, PotEvap = df_model$PET_hourly,
    TempMean = df_model$Temp_C, ZInputs = catchment_elevation, 
    HypsoData = hypso_curve, NLayers = 5                     
  )
  RunOptions_Calib_Cema <- CreateRunOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H, InputsModel = InputsModel_Cema, 
    IndPeriod_WarmUp = Ind_WarmUp_Calib, IndPeriod_Run = Ind_Run_Calib        
  )
  InputsCrit_Step2 <- CreateInputsCrit(
    FUN_CRIT = list(ErrorCrit_KGE, ErrorCrit_KGE), InputsModel = InputsModel_Cema, 
    RunOptions = RunOptions_Calib_Cema,
    Obs = list(df_model$Q_mm_hr[Ind_Run_Calib], df_model$Q_mm_hr[Ind_Run_Calib]),
    transfo = list("", "sqrt"), Weights = list(0.6, 0.4)
  )
  
  # ルートA: Step 1のX1-X5を固定してX6-X7を探索
  CalibOptions_Step2A <- CreateCalibOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H, FUN_CALIB = Calibration_Michel,
    FixedParam = c(Param_GR5H, NA, NA)
  )
  OutputsCalib_Step2A <- Calibration_Michel(
    InputsModel = InputsModel_Cema, RunOptions = RunOptions_Calib_Cema, 
    InputsCrit = InputsCrit_Step2, CalibOptions = CalibOptions_Step2A, 
    FUN_MOD = RunModel_CemaNeigeGR5H
  )
  
  # ルートB: 7変数すべてをCalibration_Michelで探索（OS以外、またはルートAが低い場合の保険）
  CalibOptions_Step2B <- CreateCalibOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H, FUN_CALIB = Calibration_Michel
  )
  OutputsCalib_Step2B <- Calibration_Michel(
    InputsModel = InputsModel_Cema, RunOptions = RunOptions_Calib_Cema, 
    InputsCrit = InputsCrit_Step2, CalibOptions = CalibOptions_Step2B, 
    FUN_MOD = RunModel_CemaNeigeGR5H
  )
  
  # ルートAとルートBのうち、通年KGEが高い方をStep 3の初期値に採用
  if (OutputsCalib_Step2A$CritFinal >= OutputsCalib_Step2B$CritFinal) {
    Param_Step2 <- OutputsCalib_Step2A$ParamFinalR
    Best_Step2_KGE <- OutputsCalib_Step2A$CritFinal
  } else {
    Param_Step2 <- OutputsCalib_Step2B$ParamFinalR
    Best_Step2_KGE <- OutputsCalib_Step2B$CritFinal
  }
  message("  -> Step 2 Best KGE: ", round(Best_Step2_KGE, 3))
  
  # -------------------------------------------------------------------
  # 【STEP 3】変換空間 (Transformed Space) で全パラメータ(X1-X7)を同時微調整
  # -------------------------------------------------------------------
  # 物理空間の境界設定 (X1, X2, X3, X4, X5, X6, X7) ※GR5HのX4下限は0.25h
  if (target_site == "OS") {
    lower_real <- c( 80, -1.5,  25, 0.25, 0.005, 0.05, 0.1)
    upper_real <- c(1200, 0.5, 400, 5.00, 0.990, 0.99, 8.0)
  } else {
    lower_real <- c( 50, -2.0,  20, 0.25, 0.001, 0.01, 0.05)
    upper_real <- c(2000, 1.0, 600, 10.0, 0.999, 0.99, 10.0)
  }
  
  # 初期値候補1: Step 2のベスト解 / 初期値候補2: Step 1の夏期解 + 標準的な雪パラメータ
  cand1_real <- pmax(pmin(Param_Step2, upper_real * 0.99), lower_real * 1.01)
  cand2_real <- pmax(pmin(c(Param_GR5H, 0.5, 1.0), upper_real * 0.99), lower_real * 1.01)
  
  init_list <- list(
    Transfo_CemaGR5H(cand1_real, Direction = "RT"),
    Transfo_CemaGR5H(cand2_real, Direction = "RT")
  )
  
  lower_transfo <- Transfo_CemaGR5H(lower_real, Direction = "RT")
  upper_transfo <- Transfo_CemaGR5H(upper_real, Direction = "RT")
  
  obj_fun_transfo <- function(par_t) {
    par_r <- Transfo_CemaGR5H(par_t, Direction = "TR")
    OutputsMod <- RunModel_CemaNeigeGR5H(InputsModel = InputsModel_Cema, RunOptions = RunOptions_Calib_Cema, Param = par_r)
    OutputsCrit <- ErrorCrit(InputsCrit = InputsCrit_Step2, OutputsModel = OutputsMod, verbose = FALSE)
    if (is.na(OutputsCrit$CritValue)) return(1e6)
    return(-OutputsCrit$CritValue)
  }
  
  best_par <- Param_Step2
  best_kge <- Best_Step2_KGE
  
  # OSのように境界が厳しいサイトでも確実に最適解を見つけるため、初期値候補からL-BFGS-Bを実行
  for (init_t in init_list) {
    init_t_clipped <- pmax(pmin(init_t, upper_transfo - 1e-4), lower_transfo + 1e-4)
    opt_res <- try(
      optim(
        par = init_t_clipped, fn = obj_fun_transfo, method = "L-BFGS-B",
        lower = lower_transfo, upper = upper_transfo, control = list(maxit = 250)
      ),
      silent = TRUE
    )
    if (!inherits(opt_res, "try-error") && (-opt_res$value > best_kge || target_site == "OS")) {
      # OSの場合は物理制約を満たす解を優先採用
      if (target_site == "OS" && best_kge == Best_Step2_KGE) {
        best_kge <- -opt_res$value
        best_par <- Transfo_CemaGR5H(opt_res$par, Direction = "TR")
      } else if (-opt_res$value > best_kge) {
        best_kge <- -opt_res$value
        best_par <- Transfo_CemaGR5H(opt_res$par, Direction = "TR")
      }
    }
  }
  
  Param_Optim <- best_par
  Final_KGE <- best_kge
  
  message("  -> Step 3 Final Parameters (X1-X7): ", paste(round(Param_Optim, 3), collapse = ", "))
  message("  -> Final Full-Year Composite KGE: ", round(Final_KGE, 3))
  
  # --- シミュレーション実行 (2000-2025) ---
  RunOptions_Sim <- CreateRunOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H, InputsModel = InputsModel_Cema, 
    IndPeriod_WarmUp = Ind_WarmUp_Sim, IndPeriod_Run = Ind_Run_Sim        
  )
  OutputsModel_Sim <- RunModel_CemaNeigeGR5H(
    InputsModel = InputsModel_Cema, RunOptions = RunOptions_Sim, Param = Param_Optim
  )
  df_model <- df_model %>% mutate(Q_sim_mm_hr = NA)
  df_model$Q_sim_mm_hr[Ind_Run_Sim] <- OutputsModel_Sim$Qsim
  results_list[[target_site]] <- df_model
}

# =============================================================================
# 7. サイト別・期間別（通年／夏期）のプロット作成とPDF出力 (Logスケール対応)
# =============================================================================
df_all_results <- bind_rows(results_list)
pdf(file = pdf_output, width = 12, height = 10)

for (target_site in site_list) {
  message("Generating plots for: ", target_site)
  
  # 通年データ (2020-2025)
  df_full <- df_all_results %>%
    filter(Site_Name == target_site, year(Date_Hour) >= 2020, year(Date_Hour) <= 2025)
  
  # 夏期データ (6-10月, 2020-2025)
  df_summer <- df_full %>%
    filter(month(Date_Hour) >= 6 & month(Date_Hour) <= 10) %>%
    mutate(Year = year(Date_Hour))
  
  # --- グラフ作成関数 ---
  create_site_plots <- function(df, title_period, is_summer = FALSE) {
    
    # 適合度指標の計算 (KGE, R2, RMSE)
    df_valid <- df %>% filter(!is.na(Q_mm_hr) & !is.na(Q_sim_mm_hr))
    r <- cor(df_valid$Q_mm_hr, df_valid$Q_sim_mm_hr)
    r2 <- r^2
    rmse <- sqrt(mean((df_valid$Q_sim_mm_hr - df_valid$Q_mm_hr)^2))
    alpha <- sd(df_valid$Q_sim_mm_hr) / sd(df_valid$Q_mm_hr)
    beta <- mean(df_valid$Q_sim_mm_hr) / mean(df_valid$Q_mm_hr)
    kge <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
    
    metric_text <- sprintf("KGE: %.3f | R2: %.3f | RMSE: %.3f", kge, r2, rmse)
    
    # ① 降水量プロット（上部）
    p_precip <- ggplot(df, aes(x = Date_Hour, y = Precipitation_mm))
    if (is_summer) p_precip <- p_precip + aes(group = Year)
    p_precip <- p_precip +
      geom_col(fill = "blue", alpha = 0.6) +
      scale_y_reverse() +
      labs(
        title = paste0("Site: ", target_site, " - ", title_period),
        subtitle = paste0("Metrics: ", metric_text),
        y = "Precip\n(mm/hr)"
      ) +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", size = 14)
      )
    
    # ② 流量プロット（中部）
    p_flow <- ggplot(df, aes(x = Date_Hour))
    if (is_summer) p_flow <- p_flow + aes(group = Year)
    p_flow <- p_flow +
      geom_line(aes(y = Q_sim_mm_hr, color = "Simulated"), alpha = 0.8, linewidth = 0.6) +
      geom_line(aes(y = Q_mm_hr, color = "Observed"), alpha = 0.6, linewidth = 0.6) +
      scale_color_manual(values = c("Observed" = "black", "Simulated" = "red")) +
      labs(x = "DateTime", y = "Discharge\n(mm/hr)", color = "") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # ③ 散布図プロット（下部：Log-Logスケール）
    eps <- 1e-4
    df_scatter <- df_valid %>%
      filter(Q_mm_hr > eps & Q_sim_mm_hr > eps)
    
    p_scatter <- ggplot(df_scatter, aes(x = Q_mm_hr, y = Q_sim_mm_hr)) +
      geom_point(alpha = 0.3, color = "darkgreen", size = 1.2) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
      scale_x_log10(labels = label_number(drop0trailing = TRUE)) +
      scale_y_log10(labels = label_number(drop0trailing = TRUE)) +
      annotation_logticks(sides = "bl", alpha = 0.5) +
      labs(
        title = "Observed vs Simulated (Log-Log Scale)",
        x = "Observed Discharge (mm/hour) [Log Scale]",
        y = "Simulated Discharge (mm/hour) [Log Scale]"
      ) +
      theme_minimal()
    
    # グラフの結合 (降水 1 : 流量 3 : 散布図 2.5)
    p_combined <- p_precip / p_flow / p_scatter + plot_layout(heights = c(1, 3, 2.5))
    
    return(p_combined)
  }
  
  # 1ページ目：通年版
  p_full_layout <- create_site_plots(df_full, "Full Year (2020-2025)", is_summer = FALSE)
  print(p_full_layout)
  
  # 2ページ目：夏期版
  p_summer_layout <- create_site_plots(df_summer, "Summer ONLY (Jun-Oct, 2020-2025)", is_summer = TRUE)
  print(p_summer_layout)
}

dev.off()
message("All plots successfully exported to PDF: ", pdf_output)

# =============================================================================
# 8. Combine Results, Convert Units, and Export
# =============================================================================
final_simulated_data <- df_all_results %>%
  mutate(
    Q_obs_m3s = mean_discharge * 0.0283168,
    Q_sim_m3s = (Q_sim_mm_hr / 1000) * (Area_sqkm * 1000000) / 3600
  ) %>%
  filter(year(Date_Hour) >= 2000 & year(Date_Hour) <= 2025)

write_csv(final_simulated_data, output_filename)
message("Combined hourly data (2000-2025) exported to CSV: ", output_filename)

