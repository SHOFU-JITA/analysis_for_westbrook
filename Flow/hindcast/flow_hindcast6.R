# =============================================================================
# Data Preparation and GR4J + CemaNeige Hydrological Modeling
# Focus: Summer 2020 (June - Oct) Calibration
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

output_filename <- "GR4J_Simulated_Discharge_WB_Summer2020.csv"

# 2. Define site dictionaries and catchment areas (sq km)
dict_wb <- tibble(
  raw_name = c("wb jimmy", "wb obear", "wb mitchell", "west brook"),
  unified_name = c("OL", "IL", "OS", "WB")
)

dict_eco <- tibble(
  raw_name = c("Jimmy Brook", "Obear Brook  Low", "Mitchell Brook", "West Brook 1 (Lower)"),
  unified_name = c("OL", "IL", "OS", "WB") 
)

area_dictionary <- tibble(
  Site_Name = c("OL", "IL", "OS", "WB"),
  Area_sqkm = c(2.46, 1.33, 1.06, 21.78) # Letcher et al. 2022 cohort
)

# 3. Load data
wb_env_data <- read_csv("../../rawdata/envDataWB.csv")
eco_draught_data <- read_csv("EcoDrought_Continuous_MA.csv")
precipitation_data <- read_csv("Westbrook_MA_precipitation.csv")

# 4. Format and standardize site names and units
# Past environmental data
wb_daily <- wb_env_data %>%
  left_join(dict_wb, by = c("river" = "raw_name")) %>%
  mutate(
    Site_Name = unified_name,
    Date = as.Date(date) 
  ) %>%
  group_by(Site_Name, Date) %>%
  summarize(
    mean_temp_C = mean(temperature, na.rm = TRUE),
    max_temp_C = mean(daily_max_temp, na.rm = TRUE),
    min_temp_C = mean(daily_min_temp, na.rm = TRUE),
    .groups = "drop"
  )

# Recent discharge and environmental data
eco_daily <- eco_draught_data %>%
  left_join(dict_eco, by = c("Site_Name" = "raw_name")) %>%
  mutate(
    Site_Name = unified_name,
    Date = as.Date(DateTime_EST) 
  ) %>%
  group_by(Site_Name, Date) %>%
  summarize(
    mean_discharge = mean(Discharge_Hobo_cfs, na.rm = TRUE),
    mean_temp_C = (mean(AirTemperature_HOBO_degF, na.rm = TRUE) - 32) * 5 / 9,
    .groups = "drop"
  )

# Precipitation data
precip_daily <- precipitation_data %>%
  mutate(Date = ymd(Date))

# =============================================================================
# 5. Merge datasets and calculate PET
# =============================================================================
master_data <- eco_daily %>%
  full_join(wb_daily, by = c("Site_Name", "Date")) %>%
  full_join(precip_daily, by = "Date") %>%
  mutate(Temp_C = coalesce(mean_temp_C.x, mean_temp_C.y)) %>%
  select(-mean_temp_C.x, -mean_temp_C.y) %>%
  arrange(Site_Name, Date) %>%
  # サイトごとにグループ化して、気温の欠損値を線形補間で埋める
  group_by(Site_Name) %>%
  mutate(Temp_C = na.approx(Temp_C, na.rm = FALSE, rule = 2)) %>%
  ungroup()

site_latitude <- 42.431868 

# ETと降水量のサイト別・最適化設定
master_data <- master_data %>%
  mutate(
    Julian_Day = yday(Date),
    
    # 【1】ETの最終設定（最も成績が良かった値に固定）
    site_pet_mult = case_when(
      Site_Name == "IL" ~ 1.3,   # 0.667を出したベスト値
      Site_Name == "OL" ~ 1.3,   # 0.621を出したベスト値
      Site_Name == "OS" ~ 1.4,   # 1.3(0.699)と1.5(0.704)の中間を狙い0.7超えを確実にする
      Site_Name == "WB" ~ 1.5,   # X1が劇的に改善したベスト値
      TRUE ~ 1.0
    ),
    
    # 【2】降水量の最終設定（最も成績が良かった値に固定）
    site_precip_mult = case_when(
      Site_Name == "IL" ~ 1.1,   # 0.667を出したベスト値
      Site_Name == "OL" ~ 1.05,  # 水収支X2が0.000になった奇跡の値
      Site_Name == "OS" ~ 1.15,  # 0.7超えを達成した不動のベスト値
      Site_Name == "WB" ~ 0.95,  # 水没状態が解消されたベスト値
      TRUE ~ 1.0
    ),
    
    PET_mm = PE_Oudin(JD = Julian_Day, Temp = Temp_C, Lat = site_latitude, LatUnit = "deg") * site_pet_mult,
    
    # 降水量に補正係数を掛ける
    Precipitation_mm = Precipitation_mm * site_precip_mult
  )

# Convert discharge (cfs) to runoff depth (mm/day)
master_data <- master_data %>%
  left_join(area_dictionary, by = "Site_Name") %>%
  mutate(
    Volume_m3_day = mean_discharge * 0.0283168 * 86400,
    Area_m2 = Area_sqkm * 1000000,
    Q_mm = (Volume_m3_day / Area_m2) * 1000
  )

# =============================================================================
# 6. GR4J + CemaNeige Modeling Loop (Site-by-Site)
# =============================================================================
results_list <- list()
site_list <- unique(master_data$Site_Name)
site_list <- site_list[!is.na(site_list)]

for (target_site in site_list) {
  
  message("Processing site: ", target_site)
  
  # Prepare site-specific data
  df_model <- master_data %>%
    filter(Site_Name == target_site) %>%
    arrange(Date) %>%
    mutate(
      Precip_mm = replace_na(Precipitation_mm, 0),
      PET_mm = na.approx(PET_mm, na.rm = FALSE, rule = 2),
      Date = as.POSIXct(Date, tz = "UTC")
    )
  
  # =======================================================
  # ★【期間設定】 2020年6〜10月にフォーカス
  # 助走期間(WarmUp): 2018年1月〜2020年5月（初期状態の安定化）
  # =======================================================
  Ind_WarmUp <- which(df_model$Date >= as.POSIXct("2018-01-01", tz = "UTC") & 
                        df_model$Date <= as.POSIXct("2020-05-31", tz = "UTC"))
  
  # 評価・キャリブレーション期間(Run): 2020年6月1日〜2020年10月31日
  Ind_Run <- which(df_model$Date >= as.POSIXct("2020-06-01", tz = "UTC") & 
                     df_model$Date <= as.POSIXct("2020-10-31", tz = "UTC"))
  
  # =======================================================
  # ★【CemaNeige用の追加設定】標高データ (130m〜170m)
  # =======================================================
  catchment_elevation <- 150 # 流域の中央値（代表標高）を150mに設定
  hypso_curve <- seq(from = 130, to = 170, length.out = 101) # 130m〜170mの勾配分布
  
  # Set up inputs and options
  InputsModel <- CreateInputsModel(
    FUN_MOD = RunModel_CemaNeigeGR4J,
    DatesR = df_model$Date,
    Precip = df_model$Precip_mm,
    PotEvap = df_model$PET_mm,
    TempMean = df_model$Temp_C,         
    ZInputs = catchment_elevation,      
    HypsoData = hypso_curve,            
    NLayers = 5                         
  )
  
  RunOptions <- CreateRunOptions(
    FUN_MOD = RunModel_CemaNeigeGR4J,
    InputsModel = InputsModel,
    IndPeriod_WarmUp = Ind_WarmUp, 
    IndPeriod_Run = Ind_Run        
  )
  
  # ★ 低水域のフィッティングを重視するため transfo = "log" に変更
  InputsCrit <- CreateInputsCrit(
    FUN_CRIT = ErrorCrit_NSE, 
    InputsModel = InputsModel,
    RunOptions = RunOptions,
    Obs = df_model$Q_mm[Ind_Run],
    transfo = "log"
  )
  
  CalibOptions <- CreateCalibOptions(
    FUN_MOD = RunModel_CemaNeigeGR4J,
    FUN_CALIB = Calibration_Michel
  )
  
  # Run calibration
  OutputsCalib <- Calibration_Michel(
    InputsModel = InputsModel,
    RunOptions = RunOptions,
    InputsCrit = InputsCrit,
    CalibOptions = CalibOptions,
    FUN_MOD = RunModel_CemaNeigeGR4J
  )
  
  Param_Optim <- OutputsCalib$ParamFinalR
  message(target_site, " optimized parameters (X1-X6): ", 
          paste(round(Param_Optim, 3), collapse = ", "))
  
  # Run simulation over the full period
  OutputsModel <- RunModel_CemaNeigeGR4J(
    InputsModel = InputsModel,
    RunOptions = RunOptions,
    Param = Param_Optim
  )
  
  # Store results
  df_model <- df_model %>% mutate(Q_sim_mm = NA)
  df_model$Q_sim_mm[Ind_Run] <- OutputsModel$Qsim
  results_list[[target_site]] <- df_model
  
  message(target_site, " processing complete.\n")
}

# =============================================================================
# 7. 観測流量(Qob)と模擬流量(Qsim)のプロット
# =============================================================================

# results_listを1つのデータフレームに結合
df_all_results <- bind_rows(results_list) %>%
  filter(!is.na(Q_mm) & !is.na(Q_sim_mm))

# ---------------------------------------------------------
# プロット1: 時系列プロット (2020年6月〜10月にズームイン)
# ---------------------------------------------------------
p1 <- ggplot(df_all_results, aes(x = Date)) +
  geom_line(aes(y = Q_mm, color = "Observed (Qob)"), alpha = 0.6, linewidth = 0.8) +
  geom_line(aes(y = Q_sim_mm, color = "Simulated (Qsim)"), alpha = 0.8, linewidth = 0.6) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Observed (Qob)" = "black", "Simulated (Qsim)" = "red")) +
  labs(
    title = "Hydrograph: Summer 2020 Focus (June - Oct)",
    x = "Date",
    y = "Discharge (mm/day)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  # X軸の範囲を2020年6月1日〜10月31日に限定
  coord_cartesian(xlim = c(as.POSIXct("2020-06-01", tz="UTC"), 
                           as.POSIXct("2020-10-31", tz="UTC")))

print(p1)

# ---------------------------------------------------------
# プロット2: 散布図 (1:1 ライン付き)
# ---------------------------------------------------------
p2 <- ggplot(df_all_results, aes(x = Q_mm, y = Q_sim_mm)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Site_Name, scales = "free", ncol = 2) +
  labs(
    title = "Scatter Plot: Observed vs Simulated Discharge (Summer 2020)",
    x = "Observed Discharge (Qob, mm/day)",
    y = "Simulated Discharge (Qsim, mm/day)"
  ) +
  theme_minimal()

print(p2)

# =============================================================================
# 8. Combine Results, Convert Units, and Export
# =============================================================================

final_simulated_data <- bind_rows(results_list) %>%
  mutate(
    Date = as.Date(Date),
    Q_obs_m3s = mean_discharge * 0.0283168,
    Q_sim_m3s = (Q_sim_mm / 1000) * (Area_sqkm * 1000000) / 86400
  )

write_csv(final_simulated_data, output_filename)

message("Combined data exported to: ", output_filename)