# =============================================================================
# Data Preparation and GR5H + CemaNeige Hydrological Modeling (Hourly)
# Focus: 2020-2025 Calibration (NSE) / 2000-2025 Simulation (Hindcast)
# Improvements: Diurnal ET distribution, GR5H upgrade, NSE objective function
# =============================================================================

# 1. Initialize environment, load packages, and set working directory
rm(list = ls(all.names = TRUE))
library(tidyverse)
library(lubridate)
library(rstudioapi)
library(airGR)
library(zoo)

if (rstudioapi::isAvailable()) {
  try({
    current_path <- rstudioapi::getSourceEditorContext()$path
    if (current_path != "") {
      setwd(dirname(current_path))
      message("Working directory set to: ", getwd())
    }
  }, silent = TRUE)
}

output_filename <- "GR5H_Simulated_Discharge_WB_2000_2025.csv"

# 2. Define site dictionaries, catchment areas (sq km), and site-specific parameters
dict_eco <- tibble(
  raw_name = c("Jimmy Brook", "Obear Brook  Low", "Mitchell Brook", "West Brook 1 (Lower)"),
  unified_name = c("OL", "IL", "OS", "WB") 
)

# 各サイト固有のパラメータをまとめた辞書を作成
site_params <- tibble(
  Site_Name = c("OL", "IL", "OS", "WB"),
  Area_sqkm = c(2.46, 1.33, 1.06, 21.78), # Letcher et al. 2022 cohort
  
  # 河川ごとのスノー・アンダーキャッチ補正係数 (例: 小流域は風の影響が少なく1.05、WBは1.2など)
  snow_catch_factor = c(1.05, 1.05, 1.05, 1.20),
  
  # 河川ごとのETマルチプライヤ (例: 樹冠被覆率や日当たりの違いを反映)
  et_multiplier = c(1.1, 1.1, 1.1, 1.00)
)

# =============================================================================
# 3. Load Data & 4. Format/Standardize
# =============================================================================
# 流量データの読み込み (15分間隔)
eco_draught_data <- read_csv("EcoDrought_Continuous_MA.csv")

# 15分間隔の流量を1時間単位 (Hourly) に集約
eco_hourly <- eco_draught_data %>%
  left_join(dict_eco, by = c("Site_Name" = "raw_name")) %>%
  mutate(
    Site_Name = unified_name,
    Date_Hour = floor_date(as.POSIXct(DateTime_EST, tz = "UTC"), "hour") 
  ) %>%
  group_by(Site_Name, Date_Hour) %>%
  summarize(
    mean_discharge = mean(Discharge_Hobo_cfs, na.rm = TRUE),
    .groups = "drop"
  )

# 新規の気象データの読み込み (サイトごとのHourlyデータ)
site_list <- c("OL", "IL", "OS", "WB")
weather_list <- list()

for (site in site_list) {
  file_name <- paste0(site, "_hourly_weather.csv")
  
  if (file.exists(file_name)) {
    temp_df <- read_csv(file_name, show_col_types = FALSE) %>%
      mutate(
        Site_Name = site,
        Date_Hour = as.POSIXct(Time, tz = "UTC") # Time列を使用
      ) %>%
      rename(Temp_C = Temperature_C) %>%         # Temperature_C を Temp_C にリネーム
      select(Site_Name, Date_Hour, Temp_C, Precipitation_mm)
    
    weather_list[[site]] <- temp_df
  } else {
    warning("File not found: ", file_name)
  }
}
weather_hourly <- bind_rows(weather_list)

# =============================================================================
# 5. Merge datasets, calculate Diurnal PET, and Apply Snow Undercatch
# =============================================================================
# airGR用に欠損のない完全な時間配列を作成
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

# 日平均気温と日積算PET(Oudin法)の計算
daily_temp_summary <- master_data %>%
  mutate(Date = as.Date(Date_Hour)) %>%
  group_by(Site_Name, Date) %>%
  summarize(mean_daily_temp = mean(Temp_C, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Julian_Day = yday(Date),
    PET_daily = PE_Oudin(JD = Julian_Day, Temp = mean_daily_temp, Lat = site_latitude, LatUnit = "deg")
  )

# ETの動的配分（日変化）と各種単位変換
master_data <- master_data %>%
  mutate(
    Date = as.Date(Date_Hour),
    Hour = hour(Date_Hour)
  ) %>%
  left_join(daily_temp_summary %>% select(Site_Name, Date, PET_daily), by = c("Site_Name", "Date")) %>%
  left_join(site_params, by = "Site_Name") %>%  # <--- ここをsite_paramsに変更
  mutate(
    # サイトごとのsnow_catch_factorを参照して補正
    Precipitation_mm = if_else(Temp_C <= 0, Precipitation_mm * snow_catch_factor, Precipitation_mm),
    
    # 疑似日変化（Diurnal）ファクター: 夜間もベースライン(0.1)を持たせる
    Diurnal_Weight = if_else(Hour >= 6 & Hour <= 18, sin((Hour - 6) * pi / 12) + 0.1, 0.1),
    
    # サイトごとのet_multiplierを参照して最終的なPETを計算
    PET_hourly = (PET_daily * (Diurnal_Weight / 8.839437)) * et_multiplier,
    
    # 流量変換 (cfs -> mm/hour)
    Volume_m3_hour = mean_discharge * 0.0283168 * 3600,
    Area_m2 = Area_sqkm * 1000000,
    Q_mm_hr = (Volume_m3_hour / Area_m2) * 1000
  ) %>%
  select(-Date, -Hour, -PET_daily, -Diurnal_Weight, -snow_catch_factor, -et_multiplier) # 最後に不要な列を削除

# =============================================================================
# 6. GR5H + CemaNeige Modeling Loop (Site-by-Site)
# =============================================================================
results_list <- list()

for (target_site in site_list) {
  
  message("Processing site: ", target_site)
  
  df_model <- master_data %>%
    filter(Site_Name == target_site) %>%
    arrange(Date_Hour) %>%
    mutate(
      PET_hourly = na.approx(PET_hourly, na.rm = FALSE, rule = 2)
    )
  
  # 期間設定
  Ind_WarmUp_Calib <- which(df_model$Date_Hour >= as.POSIXct("2018-01-01 00:00:00", tz = "UTC") & 
                              df_model$Date_Hour <= as.POSIXct("2019-12-31 23:00:00", tz = "UTC"))
  Ind_Run_Calib <- which(df_model$Date_Hour >= as.POSIXct("2020-01-01 00:00:00", tz = "UTC") & 
                           df_model$Date_Hour <= as.POSIXct("2025-12-31 23:00:00", tz = "UTC"))
  
  Ind_WarmUp_Sim <- which(df_model$Date_Hour >= as.POSIXct("1999-01-01 00:00:00", tz = "UTC") & 
                            df_model$Date_Hour <= as.POSIXct("1999-12-31 23:00:00", tz = "UTC"))
  Ind_Run_Sim <- which(df_model$Date_Hour >= as.POSIXct("2000-01-01 00:00:00", tz = "UTC") & 
                         df_model$Date_Hour <= as.POSIXct("2025-12-31 23:00:00", tz = "UTC"))
  
  # サイトに応じたCemaNeige標高設定
  if (target_site == "WB") {
    catchment_elevation <- 200
    hypso_curve <- seq(from = 150, to = 250, length.out = 101)
  } else {
    catchment_elevation <- 150 
    hypso_curve <- seq(from = 130, to = 170, length.out = 101) 
  }
  
  # GR5H入力モデルの作成
  InputsModel <- CreateInputsModel(
    FUN_MOD = RunModel_CemaNeigeGR5H,
    DatesR = df_model$Date_Hour,
    Precip = df_model$Precipitation_mm,
    PotEvap = df_model$PET_hourly,
    TempMean = df_model$Temp_C,         
    ZInputs = catchment_elevation,      
    HypsoData = hypso_curve,            
    NLayers = 5                         
  )
  
  RunOptions_Calib <- CreateRunOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H,
    InputsModel = InputsModel,
    IndPeriod_WarmUp = Ind_WarmUp_Calib, 
    IndPeriod_Run = Ind_Run_Calib        
  )
  
  # ★ 指標をNSEに変更
  InputsCrit_Calib <- CreateInputsCrit(
    FUN_CRIT = ErrorCrit_NSE, 
    InputsModel = InputsModel,
    RunOptions = RunOptions_Calib,
    Obs = df_model$Q_mm_hr[Ind_Run_Calib]
  )
  
  CalibOptions <- CreateCalibOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H,
    FUN_CALIB = Calibration_Michel
  )
  
  OutputsCalib <- Calibration_Michel(
    InputsModel = InputsModel,
    RunOptions = RunOptions_Calib,
    InputsCrit = InputsCrit_Calib,
    CalibOptions = CalibOptions,
    FUN_MOD = RunModel_CemaNeigeGR5H
  )
  
  Param_Optim <- OutputsCalib$ParamFinalR
  message(target_site, " optimized parameters (X1-X7): ", 
          paste(round(Param_Optim, 3), collapse = ", "))
  message(target_site, " calibration NSE: ", round(OutputsCalib$CritFinal, 3))
  
  # 最適化パラメータを用いた全期間のシミュレーション
  RunOptions_Sim <- CreateRunOptions(
    FUN_MOD = RunModel_CemaNeigeGR5H,
    InputsModel = InputsModel,
    IndPeriod_WarmUp = Ind_WarmUp_Sim, 
    IndPeriod_Run = Ind_Run_Sim        
  )
  
  OutputsModel_Sim <- RunModel_CemaNeigeGR5H(
    InputsModel = InputsModel,
    RunOptions = RunOptions_Sim,
    Param = Param_Optim
  )
  
  # 結果の格納
  df_model <- df_model %>% mutate(Q_sim_mm_hr = NA)
  df_model$Q_sim_mm_hr[Ind_Run_Sim] <- OutputsModel_Sim$Qsim
  results_list[[target_site]] <- df_model
  
  message(target_site, " processing complete.\n")
}

# =============================================================================
# 7. プロット (Hourly)
# =============================================================================

df_all_results <- bind_rows(results_list) %>%
  filter(!is.na(Q_sim_mm_hr))

df_plot1 <- df_all_results %>%
  filter(year(Date_Hour) >= 2020)

p1 <- ggplot(df_plot1, aes(x = Date_Hour)) +
  geom_line(aes(y = Q_sim_mm_hr, color = "Simulated (2020-2025)"), alpha = 0.7, linewidth = 0.5) +
  geom_line(aes(y = Q_mm_hr, color = "Observed (2020-2025)"), alpha = 0.7, linewidth = 0.5) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Observed (2020-2025)" = "black", "Simulated (2020-2025)" = "red")) +
  labs(
    title = "Hourly Hydrograph (GR5H): Simulation vs Observation (2020-2025)",
    x = "DateTime",
    y = "Discharge (mm/hour)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

df_calib_scatter <- df_all_results %>%
  filter(!is.na(Q_mm_hr) & !is.na(Q_sim_mm_hr))

p2 <- ggplot(df_calib_scatter, aes(x = Q_mm_hr, y = Q_sim_mm_hr)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Site_Name, scales = "free", ncol = 2) +
  labs(
    title = "Scatter Plot: Observed vs Simulated Discharge (Hourly Calibration)",
    x = "Observed Discharge (Qob, mm/hour)",
    y = "Simulated Discharge (Qsim, mm/hour)"
  ) +
  theme_minimal()

print(p2)

# =============================================================================
# 8. Combine Results, Convert Units, and Export
# =============================================================================

final_simulated_data <- df_all_results %>%
  mutate(
    # hourly mm から m3/s へ逆算 (3600秒で割る)
    Q_obs_m3s = mean_discharge * 0.0283168,
    Q_sim_m3s = (Q_sim_mm_hr / 1000) * (Area_sqkm * 1000000) / 3600
  ) %>%
  # 出力対象を2000〜2025年に絞る
  filter(year(Date_Hour) >= 2000 & year(Date_Hour) <= 2025)

write_csv(final_simulated_data, output_filename)

message("Combined hourly data (2000-2025) exported to: ", output_filename)

