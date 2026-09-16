# =============================================================================
# Flow Data Analysis: Trends, Weekly/Monthly/Seasonal Visualization, Extreme Events
# (Incorporating Detrending for Spatial Synchrony Analysis)
# =============================================================================

# 1. Initialize environment and load packages
rm(list = ls(all.names = TRUE))
library(dplyr)
library(lubridate)
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)

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

csv_file_path <- "../hindcast/GR5H_Simulated_Discharge_WB_2000_2025.csv"
if (!file.exists(csv_file_path)) stop("CSV file not found.")

# 共通のカスタムテーマ（モノクロ・文字サイズ調整版）
custom_theme <- theme_bw(base_size = 18) + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text = element_text(size = 15, color = "black"),
    strip.text = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

# =============================================================================
# Format columns: Convert Hourly to Daily, filter data (2000-2015), add Season
# =============================================================================
message("時間単位(Hourly)データを日平均(Daily)に集約しています...")

df <- read_csv(csv_file_path) %>%
  # 1. Date_Hour (時間) から Date (日) を抽出
  mutate(Date = as.Date(Date_Hour)) %>%
  # 2. 日ごとにグループ化し、日平均流量を計算
  group_by(Site_Name, Date) %>%
  summarise(
    Q_sim_m3s = mean(Q_sim_m3s, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # 3. 年・月・季節などの情報を追加
  mutate(
    year = year(Date),
    month = month(Date),
    # 定義: 春(3-5), 夏(6-8), 秋(9-11), 冬(12-2)
    Season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn"
    ),
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    bio_year = if_else(month == 12, year + 1, year)
  ) %>%
  filter(year >= 2000 & year <= 2015, !is.na(Q_sim_m3s))


# =============================================================================
# Plot 1-3: Weekly, Monthly, Seasonal flow visualization
# =============================================================================
message("Weekly, Monthly, Seasonal プロットを作成中...")

# Weekly Plot
weekly_flow <- df %>%
  mutate(year_week = floor_date(Date, "week")) %>%
  group_by(Site_Name, year_week) %>%
  summarise(
    mean_flow = mean(Q_sim_m3s, na.rm = TRUE),
    sd_flow   = sd(Q_sim_m3s, na.rm = TRUE),
    n_obs     = sum(!is.na(Q_sim_m3s)),
    .groups = 'drop'
  ) %>%
  mutate(se_flow = sd_flow / sqrt(n_obs))

plot_weekly <- ggplot(weekly_flow, aes(x = year_week, y = mean_flow)) +
  geom_errorbar(aes(ymin = mean_flow - se_flow, ymax = mean_flow + se_flow), 
                color = "grey50", width = 3, alpha = 0.7) +
  geom_line(color = "black", linewidth = 0.6) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Weekly Mean Discharge Time Series (with SE)",
    subtitle = "Simulated Discharge Data (2000-2015)",
    x = "Year", y = "Weekly Mean Flow (m3/s)"
  ) + custom_theme

# Monthly Plot
monthly_flow <- df %>%
  mutate(year_month = floor_date(Date, "month")) %>%
  group_by(Site_Name, year_month) %>%
  summarise(
    mean_flow = mean(Q_sim_m3s, na.rm = TRUE),
    sd_flow   = sd(Q_sim_m3s, na.rm = TRUE),
    n_obs     = sum(!is.na(Q_sim_m3s)),
    .groups = 'drop'
  ) %>%
  mutate(se_flow = sd_flow / sqrt(n_obs))

plot_monthly <- ggplot(monthly_flow, aes(x = year_month, y = mean_flow)) +
  geom_errorbar(aes(ymin = mean_flow - se_flow, ymax = mean_flow + se_flow), 
                color = "grey50", width = 15, alpha = 0.7) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(color = "black", size = 1.5) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Monthly Mean Discharge Time Series (with SE)",
    subtitle = "Simulated Discharge Data (2000-2015)",
    x = "Year", y = "Monthly Mean Flow (m3/s)"
  ) + custom_theme

# Seasonal Plot
seasonal_flow <- df %>%
  group_by(Site_Name, bio_year, Season) %>%
  summarise(
    mean_flow = mean(Q_sim_m3s, na.rm = TRUE),
    sd_flow   = sd(Q_sim_m3s, na.rm = TRUE),
    n_obs     = sum(!is.na(Q_sim_m3s)),
    season_date = mean(Date, na.rm = TRUE),
    .groups   = 'drop'
  ) %>%
  mutate(se_flow = sd_flow / sqrt(n_obs))

plot_seasonal <- ggplot(seasonal_flow, aes(x = season_date, y = mean_flow)) +
  geom_errorbar(aes(ymin = mean_flow - se_flow, ymax = mean_flow + se_flow), 
                color = "black", width = 45, alpha = 0.7) +
  geom_line(color = "black", linewidth = 0.8) + 
  geom_point(color = "black", size = 2) +      
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    title = "Seasonal Mean Discharge Time Series (with SE)",
    subtitle = "Simulated Discharge Data (2000-2015)",
    x = "Year", y = "Seasonal Mean Flow (m3/s)" 
  ) + custom_theme


# =============================================================================
# Plot 4-5: Extreme Flow Events 
# =============================================================================
message("Extreme Events プロットを作成中...")

thresholds <- df %>%
  filter(Q_sim_m3s > 0) %>%
  group_by(Site_Name) %>%
  summarise(
    high_threshold = quantile(Q_sim_m3s, 0.95, na.rm = TRUE),
    low_threshold = quantile(Q_sim_m3s, 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

extreme_events <- df %>%
  left_join(thresholds, by = "Site_Name") %>%
  group_by(Site_Name, year) %>%
  summarise(
    High_Flow_Days = sum(Q_sim_m3s >= high_threshold, na.rm = TRUE),
    Low_Flow_Days = sum(Q_sim_m3s <= low_threshold, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(High_Flow_Days, Low_Flow_Days), names_to = "Event_Type", values_to = "Days")

plot_extreme <- ggplot(extreme_events, aes(x = year, y = Days, group = Event_Type, shape = Event_Type, linetype = Event_Type)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(color = "black", size = 3) +
  facet_wrap(~ Site_Name, scales = "free_y", ncol = 1) +
  labs(
    title = "Annual Frequency of Extreme Flow Events",
    subtitle = "Solid: High Flow (>95th) | Dashed: Low Flow (<5th)",
    x = "Year", y = "Number of Days per Year"
  ) + custom_theme + theme(legend.position = "top", legend.title = element_blank())

seasonal_extremes <- df %>%
  left_join(thresholds, by = "Site_Name") %>%
  group_by(Site_Name, bio_year, Season) %>%
  summarise(
    High_Flow_Days = sum(Q_sim_m3s >= high_threshold, na.rm = TRUE),
    Low_Flow_Days = sum(Q_sim_m3s <= low_threshold, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(High_Flow_Days, Low_Flow_Days), names_to = "Event_Type", values_to = "Days")

plot_seasonal_extreme <- ggplot(seasonal_extremes, aes(x = bio_year, y = Days, group = Event_Type, shape = Event_Type, linetype = Event_Type)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(color = "black", size = 2) +
  facet_grid(Season ~ Site_Name, scales = "free_y") +
  labs(
    title = "Seasonal Frequency of Extreme Flow Events",
    subtitle = "High and Low flow days categorized by season",
    x = "Year", y = "Days per Season"
  ) + custom_theme + theme(legend.position = "top", legend.title = element_blank())

# Save Plots
pdf_filename <- file.path(SCRIPT_DIR, "WB_flow.pdf")
pdf(file = pdf_filename, width = 10, height = 8)
print(plot_weekly)
print(plot_monthly)
print(plot_seasonal)
print(plot_extreme)
print(plot_seasonal_extreme)
dev.off()
message(paste("All 5 graphs successfully saved to:", pdf_filename))

# =============================================================================
# Section 6: Detrended Annual Mean Flow Cross-site Correlation (全期間・年間ベース)
# =============================================================================
message("Detrending annual mean flows and calculating correlations...")

# 1. 年平均流量の算出とDetrending (線形回帰によるトレンド除去・残差の抽出)
annual_flow <- df %>%
  group_by(Site_Name, year) %>%
  summarise(
    annual_mean_flow = mean(Q_sim_m3s, na.rm = TRUE),
    .groups = 'drop'
  )

detrended_annual <- annual_flow %>%
  group_by(Site_Name) %>%
  mutate(
    # nが3以上かつ分散が0より大きい場合のみ回帰を実行
    mean_flow_res = if (n() >= 3 && sd(annual_mean_flow, na.rm = TRUE) > 0) {
      residuals(lm(annual_mean_flow ~ year, na.action = na.exclude))
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

# 2. 全期間（年間ベース）のペアワイズ相関を計算
wide_res <- detrended_annual %>%
  select(year, Site_Name, mean_flow_res) %>%
  pivot_wider(names_from = Site_Name, values_from = mean_flow_res) %>%
  select(-year)

site_cols <- colnames(wide_res)
n_s <- length(site_cols)
cor_results_mean <- list()

for (i in 1:(n_s - 1)) {
  for (j in (i + 1):n_s) {
    s1 <- site_cols[i]
    s2 <- site_cols[j]
    v1 <- wide_res[[s1]]
    v2 <- wide_res[[s2]]
    
    valid_idx <- complete.cases(v1, v2)
    if (sum(valid_idx) >= 3 && sd(v1[valid_idx]) > 0 && sd(v2[valid_idx]) > 0) {
      test <- cor.test(v1, v2, method = "pearson")
      cor_results_mean[[length(cor_results_mean) + 1]] <- data.frame(
        Period = "Annual", site_1 = s1, site_2 = s2,
        r_value = test$estimate, p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

cor_export_data <- do.call(rbind, cor_results_mean) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_value_adj < 0.001 ~ "***",
      p_value_adj < 0.01  ~ "**",
      p_value_adj < 0.05  ~ "*",
      p_value_adj < 0.1   ~ ".",
      TRUE                ~ ""
    )
  ) %>%
  arrange(desc(r_value))

csv_out_path <- file.path(SCRIPT_DIR, "Detrended_Annual_Flow_Correlation.csv")
write_csv(cor_export_data, csv_out_path)
message(paste("Detrended annual mean correlation saved to:", csv_out_path))

# =============================================================================
# Section 6: Detrended Seasonal Mean Flow Cross-site Correlation
# =============================================================================
message("Detrending seasonal mean flows and calculating correlations...")

# 1. Detrending (線形回帰によるトレンド除去・残差の抽出)
detrended_seasonal <- seasonal_flow %>%
  group_by(Site_Name, Season) %>%
  mutate(
    # nが3以上かつ分散が0より大きい場合のみ回帰を実行
    mean_flow_res = if (n() >= 3 && sd(mean_flow, na.rm = TRUE) > 0) {
      residuals(lm(mean_flow ~ bio_year, na.action = na.exclude))
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

cor_results_mean <- list()
seasons_list <- levels(detrended_seasonal$Season)

# 2. 季節ごとにペアワイズ相関を計算
for (sea in seasons_list) {
  wide_res <- detrended_seasonal %>%
    filter(Season == sea) %>%
    select(bio_year, Site_Name, mean_flow_res) %>%
    pivot_wider(names_from = Site_Name, values_from = mean_flow_res) %>%
    select(-bio_year)
  
  site_cols <- colnames(wide_res)
  n_s <- length(site_cols)
  
  for (i in 1:(n_s - 1)) {
    for (j in (i + 1):n_s) {
      s1 <- site_cols[i]
      s2 <- site_cols[j]
      v1 <- wide_res[[s1]]
      v2 <- wide_res[[s2]]
      
      valid_idx <- complete.cases(v1, v2)
      if (sum(valid_idx) >= 3 && sd(v1[valid_idx]) > 0 && sd(v2[valid_idx]) > 0) {
        test <- cor.test(v1, v2, method = "pearson")
        cor_results_mean[[length(cor_results_mean) + 1]] <- data.frame(
          Season = sea, site_1 = s1, site_2 = s2,
          r_value = test$estimate, p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

cor_export_data <- do.call(rbind, cor_results_mean) %>%
  group_by(Season) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_value_adj < 0.001 ~ "***",
      p_value_adj < 0.01  ~ "**",
      p_value_adj < 0.05  ~ "*",
      p_value_adj < 0.1   ~ ".",
      TRUE                ~ ""
    )
  ) %>%
  ungroup() %>%
  arrange(Season, desc(r_value))

csv_out_path <- file.path(SCRIPT_DIR, "Detrended_Seasonal_Flow_Correlation.csv")
write_csv(cor_export_data, csv_out_path)
message(paste("Detrended seasonal mean correlation saved to:", csv_out_path))


# =============================================================================
# Section 7: Detrended Extreme Events Cross-site Correlation (Annual & Seasonal)
# =============================================================================
message("Detrending extreme events and calculating correlations...")

# 相関計算用関数（内部でトレンド除去を行う）
calculate_event_correlations <- function(df_long, context_label, time_col) {
  
  # トレンド除去: Event日数 ~ 年
  df_detrended <- df_long %>%
    group_by(Site_Name) %>%
    mutate(
      Days_res = if(n() >= 3 && sd(Days, na.rm = TRUE) > 0) {
        residuals(lm(as.formula(paste("Days ~", time_col)), na.action = na.exclude))
      } else {
        NA_real_
      }
    ) %>%
    ungroup()
  
  # 残差を使ってWide型へ
  wide_df <- df_detrended %>%
    select(all_of(time_col), Site_Name, Days_res) %>%
    pivot_wider(names_from = Site_Name, values_from = Days_res)
  
  site_names <- setdiff(names(wide_df), time_col)
  n_col <- length(site_names)
  if (n_col < 2) return(NULL)
  
  results <- list()
  for (i in 1:(n_col - 1)) {
    for (j in (i + 1):n_col) {
      site1 <- site_names[i]
      site2 <- site_names[j]
      v1 <- wide_df[[site1]]
      v2 <- wide_df[[site2]]
      
      valid_idx <- complete.cases(v1, v2)
      if (sum(valid_idx) >= 3 && sd(v1[valid_idx]) > 0 && sd(v2[valid_idx]) > 0) {
        test <- cor.test(v1, v2, method = "pearson")
        results[[length(results) + 1]] <- data.frame(
          context = context_label, site_1 = site1, site_2 = site2,
          r_value = test$estimate, p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(results) == 0) return(NULL)
  return(do.call(rbind, results))
}

extreme_cor_list <- list()

extreme_cor_list[[1]] <- calculate_event_correlations(
  extreme_events %>% filter(Event_Type == "High_Flow_Days"), 
  "Annual_High_Flow", "year"
)
extreme_cor_list[[2]] <- calculate_event_correlations(
  extreme_events %>% filter(Event_Type == "Low_Flow_Days"), 
  "Annual_Low_Flow", "year"
)

event_types <- c("High_Flow_Days", "Low_Flow_Days")
for (sea in seasons_list) {
  for (ev in event_types) {
    context_name <- paste0("Seasonal_", sea, "_", ev)
    res <- calculate_event_correlations(
      seasonal_extremes %>% filter(Season == sea, Event_Type == ev),
      context_name, "bio_year"
    )
    if (!is.null(res)) extreme_cor_list[[length(extreme_cor_list) + 1]] <- res
  }
}

extreme_cor_df <- do.call(rbind, extreme_cor_list) %>%
  filter(!is.na(p_value)) %>% 
  group_by(context) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_value_adj < 0.001 ~ "***",
      p_value_adj < 0.01  ~ "**",
      p_value_adj < 0.05  ~ "*",
      p_value_adj < 0.1   ~ ".",
      TRUE                ~ ""
    )
  ) %>%
  ungroup() %>%
  arrange(context, desc(r_value))

csv_extreme_out_path <- file.path(SCRIPT_DIR, "Detrended_Extreme_Events_Correlation.csv")
write_csv(extreme_cor_df, csv_extreme_out_path)

message(paste("Detrended extreme events correlation saved to:", csv_extreme_out_path))
message("Analysis Complete.")

