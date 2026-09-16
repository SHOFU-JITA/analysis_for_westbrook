# =============================================================================
# Flow Data Analysis: Trends, Seasonal Analysis, and Extreme Events
# =============================================================================

# 1. Initialize environment and load packages
rm(list = ls(all.names = TRUE))
library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)

# =============================================================================
# Configuration (設定)
# =============================================================================
# ファイルパスと出力PDF名の設定
csv_file_path   <- "../hindcast/GR4J_Simulated_Discharge_WB.csv"
output_pdf_path <- "Flow_Analysis_Temporal_and_Extremes.pdf" # 出力するPDF名
csv_output_path <- "Seasonal_Extreme_Flow_Days.csv"          # 出力するCSV名

# 2. Set Working Directory
if (rstudioapi::isAvailable()) {
  try({
    current_path <- rstudioapi::getSourceEditorContext()$path
    if (current_path != "") {
      setwd(dirname(current_path))
      message("Working directory set to: ", getwd())
    }
  }, silent = TRUE)
}
SCRIPT_DIR <- getwd()

if (!file.exists(csv_file_path)) stop("CSV file not found.")

# =============================================================================
# Data Loading & Preprocessing (Seasonの定義を含む)
# =============================================================================
# データの読み込み、フォーマット変換、対象期間のフィルタリング、およびSeasonの定義
df <- read_csv(csv_file_path) %>%
  mutate(
    Date = as.Date(Date),
    year = year(Date),
    # ここでSeasonを定義 (1:Spring, 2:Summer, 3:Autumn, 4:Winter)
    season_num = case_when(
      (month(Date) == 2 & day(Date) >= 6) | (month(Date) > 2 & month(Date) < 5) | (month(Date) == 5 & day(Date) <= 5) ~ 1,
      (month(Date) == 5 & day(Date) >= 6) | (month(Date) > 5 & month(Date) < 8) | (month(Date) == 8 & day(Date) <= 5) ~ 2,
      (month(Date) == 8 & day(Date) >= 6) | (month(Date) > 8 & month(Date) < 11) | (month(Date) == 11 & day(Date) <= 5) ~ 3,
      TRUE ~ 4
    ),
    # season列を 1, 2, 3, 4 の数値として設定
    season = season_num
  ) %>%
  filter(year >= 2000 & year <= 2015, !is.na(Q_sim_m3s))


# =============================================================================
# Objective 1: Weekly flow visualization
# =============================================================================
message("Calculating and plotting weekly flow...")
weekly_flow <- df %>%
  mutate(year_week = floor_date(Date, "week")) %>%
  group_by(Site_Name, year_week) %>%
  summarise(mean_flow = mean(Q_sim_m3s, na.rm = TRUE), .groups = 'drop')

plot_weekly <- ggplot(weekly_flow, aes(x = year_week, y = mean_flow)) +
  geom_line(aes(color = Site_Name), alpha = 0.8) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 2) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Weekly Mean Discharge",
    x = "Date (Weekly)", y = "Mean Flow (m3/s)"
  ) +
  theme_bw() + theme(legend.position = "none")


# =============================================================================
# Objective 2: Monthly flow visualization with Standard Deviation
# =============================================================================
message("Calculating and plotting monthly flow with SD...")
monthly_flow <- df %>%
  mutate(year_month = floor_date(Date, "month")) %>%
  group_by(Site_Name, year_month) %>%
  summarise(
    mean_flow = mean(Q_sim_m3s, na.rm = TRUE),
    sd_flow = sd(Q_sim_m3s, na.rm = TRUE),
    .groups = 'drop'
  )

plot_monthly <- ggplot(monthly_flow, aes(x = year_month, y = mean_flow, color = Site_Name)) +
  geom_line(alpha = 0.8) +
  geom_ribbon(aes(ymin = pmax(0, mean_flow - sd_flow), ymax = mean_flow + sd_flow, fill = Site_Name), alpha = 0.2, color = NA) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 2) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Monthly Mean Discharge with Standard Deviation",
    subtitle = "Ribbon represents ±1 Standard Deviation (lower bound capped at 0)",
    x = "Date (Monthly)", y = "Mean Flow (m3/s)"
  ) +
  theme_bw() + theme(legend.position = "none")


# =============================================================================
# Objective 3: Annual flow visualization with Standard Deviation
# =============================================================================
message("Calculating and plotting annual flow with SD...")
annual_flow <- df %>%
  group_by(Site_Name, year) %>%
  summarise(
    mean_flow = mean(Q_sim_m3s, na.rm = TRUE),
    sd_flow = sd(Q_sim_m3s, na.rm = TRUE),
    .groups = 'drop'
  )

plot_annual <- ggplot(annual_flow, aes(x = year, y = mean_flow, color = Site_Name)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, mean_flow - sd_flow), ymax = mean_flow + sd_flow), width = 0.3, alpha = 0.6) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 2) +
  labs(
    title = "Annual Mean Discharge with Standard Deviation",
    subtitle = "Error bars represent ±1 Standard Deviation (lower bound capped at 0)",
    x = "Year", y = "Annual Mean Flow (m3/s)"
  ) +
  theme_bw() + theme(legend.position = "none")


# =============================================================================
# Objective 4: Seasonal Extreme Events (Droughts & Floods)
# =============================================================================
message("Calculating seasonal extreme events...")

# 全期間に基づく各サイトの5%パーセンタイル（渇水）、95%パーセンタイル（増水）の閾値を算出
thresholds <- df %>%
  filter(Q_sim_m3s > 0) %>%
  group_by(Site_Name) %>%
  summarise(
    high_threshold = quantile(Q_sim_m3s, 0.95, na.rm = TRUE),
    low_threshold  = quantile(Q_sim_m3s, 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

# 各年・季節・河川ごとの極端な流量の日数と平均流量をカウント・計算
seasonal_extreme_days <- df %>%
  left_join(thresholds, by = "Site_Name") %>%
  group_by(Site_Name, year, season) %>%
  summarise(
    Mean_Flow      = mean(Q_sim_m3s, na.rm = TRUE), # 季節ごとの平均流量を追加
    High_Flow_Days = sum(Q_sim_m3s >= high_threshold, na.rm = TRUE),
    Low_Flow_Days  = sum(Q_sim_m3s <= low_threshold, na.rm = TRUE),
    .groups = 'drop'
  )

# CSVとして出力 (Mean_Flow 列が含まれます)
write_csv(seasonal_extreme_days, file.path(SCRIPT_DIR, csv_output_path))
message(paste("Extreme flow data saved to:", csv_output_path))

# 可視化のためにLongフォーマットへ変換（グラフに不要な Mean_Flow は除外）
seasonal_extreme_long <- seasonal_extreme_days %>%
  select(-Mean_Flow) %>% 
  pivot_longer(
    cols = c(High_Flow_Days, Low_Flow_Days), 
    names_to = "Event_Type", 
    values_to = "Days"
  )

# プロット作成: 行にサイト、列に季節を配置して一覧表示
plot_extreme_seasonal <- ggplot(seasonal_extreme_long, aes(x = year, y = Days, fill = Event_Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  facet_grid(Site_Name ~ season, scales = "free_y") +
  scale_fill_manual(
    values = c("High_Flow_Days" = "blue", "Low_Flow_Days" = "red"),
    labels = c("High Flow (>= 95th)", "Low Flow (<= 5th)")
  ) +
  labs(
    title = "Annual Extreme Flow Days by Season (1:Spring, 2:Summer, 3:Autumn, 4:Winter)",
    subtitle = paste("Thresholds based on overall 5th and 95th percentiles per site"),
    x = "Year", y = "Number of Days",
    fill = "Extreme Event"
  ) +
  theme_bw() + theme(legend.position = "bottom")


# =============================================================================
# Output All Plots to a Single PDF
# =============================================================================
message("Saving all plots to a single PDF...")

# PDFデバイスを開く (横長の方が見やすいため width=14 に設定)
pdf(file = file.path(SCRIPT_DIR, output_pdf_path), width = 14, height = 9)

print(plot_weekly)
print(plot_monthly)
print(plot_annual)
print(plot_extreme_seasonal)

# デバイスを閉じて保存完了
dev.off()

message(paste("All analyses complete. PDF saved to:", output_pdf_path))

