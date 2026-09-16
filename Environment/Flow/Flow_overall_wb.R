# ===== 0. 環境リセット=====
rm(list = ls(all.names = TRUE))

# ===== 1. パッケージ =====
library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(purrr)
library(corrplot)
library(ggplot2)

# ===== 2. データ読み込み =====
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}
csv_file_path <- "~/Westbrook_カワマスP/code_SF/rawdata/envDataWB.csv"
if (!file.exists(csv_file_path)) stop(paste("指定ファイルが見つかりません:", csv_file_path))
message("データを読み込んでいます...")
df <- read_csv(csv_file_path)
pdf_filename <- "Flow_index_WB.pdf"

# ===== 3. データ前処理 =====
df <- df %>%
  mutate(
    date = as.Date(date),
    year = year(date)
  ) %>%
  filter(year >= 2000 & year <= 2015)

# ===== 4. 流況分析 (FDC) データ作成 (川ごと) =====
message("川ごとに流況分析データを計算しています...")

# river列でグループ化し、NAと0を除外
fdc_data_grouped <- df %>%
  filter(!is.na(flowByRiver) & flowByRiver > 0) %>%
  group_by(river)

# グループ内でソートし、順位(rank)と超過確率(exceedance_prob_pct)を計算
# P = (Rank / (N + 1)) * 100
fdc_data <- fdc_data_grouped %>%
  # .by_group = TRUE を指定し、グループ(river)内で降順ソート
  arrange(desc(flowByRiver), .by_group = TRUE) %>% 
  mutate(
    rank = 1:n(),                  # n() はグループ内のデータ数を返す
    n_obs = n(),                   
    exceedance_prob_pct = (rank / (n_obs + 1)) * 100
  ) %>%
  ungroup() # この後の計算のため、グループ化を一旦解除

message("FDCデータ計算完了。")

# ===== 5. 主要な流況指標の計算 (川ごと) =====
# 【変更点】指標のラベルを英語に変更
prob_targets_pct <- c(
  "Q95 (High)"   = (95 / 365.25) * 100,
  "Q185 (Normal)" = (185 / 365.25) * 100,
  "Q275 (Low)"    = (275 / 365.25) * 100,
  "Q355 (Drought)" = (355 / 365.25) * 100
)

message("川ごとに主要な流況指標を計算しています...")

# 川ごとに指標を計算
key_metrics_table <- fdc_data %>%
  group_by(river) %>%
  group_modify(~ {
    # .x は現在の川のFDCデータ
    key_flows <- approx(
      x = .x$exceedance_prob_pct,
      y = .x$flowByRiver,
      xout = prob_targets_pct
    )$y
    
    tibble(
      Metric = names(prob_targets_pct),
      Exceedance_Prob_Pct = prob_targets_pct,
      Calculated_Flow = key_flows
    )
  }) %>%
  ungroup()

message("主要な流況指標 (川ごと):")
print(key_metrics_table)

# ===== 6. 流況曲線 (FDC) のプロット (川ごと) =====
message("川ごとに流況曲線をプロットしています...")

# ベースとなるプロットを作成
fdc_plot_base <- ggplot(fdc_data, aes(x = exceedance_prob_pct, y = flowByRiver)) +
  geom_line(aes(color = river), linewidth = 1.1) +
  scale_y_log10(
    labels = scales::comma
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    labels = scales::percent_format(scale = 1)
  ) +
  geom_point(
    data = key_metrics_table, 
    aes(x = Exceedance_Prob_Pct, y = Calculated_Flow, color = river), 
    size = 2.5
  ) +
  geom_text(
    data = key_metrics_table,
    aes(x = Exceedance_Prob_Pct, y = Calculated_Flow, label = Metric),
    color = "black",
    vjust = -1.1,
    size = 3.0,
    check_overlap = TRUE
  ) +
  
  # 【変更点】labs()内のラベルを英語に変更
  labs(
    title = "Flow Duration Curve (FDC) by River",
    subtitle = "Westbrook Stream Data (2000-2015)",
    x = "Exceedance Probability (%)",
    y = "Daily Mean Discharge (flowByRiver) [Log Scale]",
    color = "River" # 凡例のタイトル
  ) +
  
  theme_bw() +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )

# --- プロットを表示 ---
fdc_plot_faceted <- fdc_plot_base +
  facet_wrap(~ river, scales = "free_y") 

print(fdc_plot_faceted)

# ===== 7. 流況指数 (FCI) の計算 =====
message("流況指数 (FCI) を計算しています...")

# ステップ 7.1: 平年値（日別平均流量）の計算
df_with_md <- df %>%
  filter(!is.na(flowByRiver)) %>%
  mutate(
    month_day = format(date, "%m-%d")
  )

message("平年値（日別平均）を計算中...")
normal_flows <- df_with_md %>%
  group_by(river, month_day) %>%
  summarise(
    normal_flow = mean(flowByRiver, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(normal_flow > 0)

# ステップ 7.2: FCIの計算
fci_data <- df_with_md %>%
  left_join(normal_flows, by = c("river", "month_day")) %>%
  filter(!is.na(normal_flow)) %>%
  mutate(
    fci = (flowByRiver / normal_flow) * 100
  )

message("FCI計算完了。")

# ===== 8. 流況指数 (FCI) のプロット (時系列) =====
message("流況指数 (FCI) を時系列でプロットしています...")

fci_plot <- ggplot(fci_data, aes(x = date, y = fci)) +
  geom_line(aes(color = river), alpha = 0.8) +
  geom_hline(
    yintercept = 100, 
    linetype = "dashed", 
    color = "red", 
    linewidth = 1
  ) +
  facet_wrap(~ river, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  
  # 【変更点】labs()内のラベルを英語に変更
  labs(
    title = "Flow Condition Index (FCI) Time Series",
    subtitle = "Westbrook Stream Data (2000-2015) | Red dashed line = Normal (100%)",
    x = "Date",
    y = "FCI (Percent of Normal Flow, %)",
    color = "River"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "none"
  )

print(fci_plot)



# ===== 8.5. 時系列流量プロット (日ごと・月ごと・季節ごと) =====
message("時系列流量プロットを作成しています...")

# データの準備 (NAを除外)
flow_data <- df %>%
  filter(!is.na(flowByRiver))

# --- 1. 日ごとのプロット (Daily) ---
plot_daily <- ggplot(flow_data, aes(x = date, y = flowByRiver)) +
  geom_line(aes(color = river), alpha = 0.7) +
  facet_wrap(~ river, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Daily Mean Discharge Time Series",
    subtitle = "Westbrook Stream Data (2000-2015)",
    x = "Date",
    y = "Daily Discharge (flowByRiver)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

message("日ごとプロット作成完了。")

# --- 2. 月ごとの集計とプロット (Monthly) ---
monthly_flow <- flow_data %>%
  # floor_date() で日付をその月の初日に丸める
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(river, year_month) %>%
  summarise(mean_flow = mean(flowByRiver, na.rm = TRUE), .groups = 'drop')

plot_monthly <- ggplot(monthly_flow, aes(x = year_month, y = mean_flow)) +
  geom_line(aes(color = river)) +
  facet_wrap(~ river, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "Monthly Mean Discharge Time Series",
    subtitle = "Westbrook Stream Data (2000-2015)",
    x = "Month",
    y = "Monthly Mean Discharge (flowByRiver)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

message("月ごとプロット作成完了。")

# --- 3. 季節ごとの集計とプロット (Seasonal) [モノクロ・文字サイズ調整版] ---
message("季節ごとプロット（モノクロ・SE版）を作成しています...")

# 季節の定義と集計（データ処理部分は変更なし）
seasonal_flow <- flow_data %>%
  mutate(
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn")),
    bio_year = if_else(month == 12, year + 1, year)
  ) %>%
  group_by(river, bio_year, season) %>%
  summarise(
    mean_flow = mean(flowByRiver, na.rm = TRUE),
    sd_flow   = sd(flowByRiver, na.rm = TRUE),
    n_obs     = sum(!is.na(flowByRiver)),
    season_date = mean(date, na.rm = TRUE),
    .groups   = 'drop'
  ) %>%
  mutate(
    se_flow = sd_flow / sqrt(n_obs)
  )

# プロット作成（色分けを削除）
plot_seasonal <- ggplot(seasonal_flow, aes(x = season_date, y = mean_flow)) +
  # エラーバー: 黒色に統一
  geom_errorbar(aes(ymin = mean_flow - se_flow, 
                    ymax = mean_flow + se_flow), 
                color = "black",  # ここで色を指定
                width = 45, alpha = 0.7) +
  
  # 折れ線: 黒色に統一
  geom_line(color = "black", linewidth = 0.8) + 
  
  # 点: 黒色に統一
  geom_point(color = "black", size = 2) +      
  
  facet_wrap(~ river, scales = "free_y", ncol = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  
  coord_cartesian(ylim = c(0, NA)) +
  
  labs(
    title = "Seasonal Mean Discharge Time Series (with SE)",
    subtitle = "Westbrook Stream Data (2000-2015)",
    x = "Year",
    y = "Stream flow by river (m^3/s)" # Superscript formatting depends on output device, simplified here
  ) +
  
  theme_bw(base_size = 18) + 
  theme(
    legend.position = "none", # 凡例は不要なので非表示のまま
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.text = element_text(size = 15, color = "black"),
    strip.text = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "grey90")
  )

print(plot_seasonal)
message("季節ごとプロット（モノクロ）作成完了。")

# ===== 9. PDFへの保存 (修正版) =====
message("全グラフをPDFファイルに保存しています...")

# 【変更点】SCRIPT_DIRを使って保存先をスクリプトと同じフォルダに指定
pdf_filename <- file.path(SCRIPT_DIR, "WB_flow.pdf")

# 1. 保存するPDFファイルを開く
pdf(file = pdf_filename, width = 10, height = 8)

# 2. 1ページ目: FDC
message("1ページ目: 流況曲線 (FDC) を描画中...")
print(fdc_plot_faceted)

# 3. 2ページ目: FCI
message("2ページ目: 流況指数 (FCI) を描画中...")
print(fci_plot)

# 4. 3ページ目: 日ごと流量
message("3ページ目: 日ごと流量 (Daily) を描画中...")
print(plot_daily)

# 5. 4ページ目: 月ごと流量
message("4ページ目: 月ごと流量 (Monthly) を描画中...")
print(plot_monthly)

# 6. 5ページ目: 季節ごと流量
message("5ページ目: 季節ごと流量 (Seasonal) を描画中...")
print(plot_seasonal)

# 7. ファイルを閉じて保存を完了
dev.off()

message(paste("All 5 graphs saved to:", pdf_filename))

