# ===== 0. 環境リセット =====
rm(list = ls(all.names = TRUE))

# ===== 1. パッケージ =====
library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(purrr)
library(corrplot)
library(ggplot2)
library(tibble) # 行名操作用

# ===== 2. ディレクトリとデータ準備 =====
# スクリプトの保存場所を特定
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}
message(paste("保存先ディレクトリ:", SCRIPT_DIR))

# データ読み込み
csv_file_path <- "~/Westbrook_カワマスP/code_SF/rawdata/envDataWB.csv"
if (!file.exists(csv_file_path)) stop(paste("指定ファイルが見つかりません:", csv_file_path))

message("データを読み込んでいます...")
df <- read_csv(csv_file_path)

# ===== 3. データ前処理 =====
df <- df %>%
  mutate(
    date = as.Date(date),
    year = year(date)
  ) %>%
  filter(year >= 2000 & year <= 2015)

# 季節列の作成
df <- df %>%
  mutate(
    month_num = month(date),
    season = case_when(
      month_num %in% 3:5 ~ 1,   # 春
      month_num %in% 6:8 ~ 2,   # 夏
      month_num %in% 9:11 ~ 3,  # 秋
      TRUE ~ 4                  # 冬
    ),
    season_label = factor(season,
                          levels = c(1, 2, 3, 4),
                          labels = c("Spring", "Summer", "Autumn", "Winter"))
  )

# ===== 4. flowByRiver の存在確認 =====
if (!"flowByRiver" %in% names(df)) stop("The 'flowByRiver' column is missing from the data.")

# ===== A. 相関解析 (flowByRiver) =====
message("===== A. 相関解析開始 =====")

# 基本データの集計
seasonal_summary_stats <- df %>%
  filter(!is.na(flowByRiver)) %>%
  group_by(river, year, season_label) %>%
  summarise(
    avg_flow = mean(flowByRiver, na.rm = TRUE),
    n_obs = n(),
    sd_flow = sd(flowByRiver, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_flow = ifelse(n_obs > 0, sd_flow / sqrt(n_obs), NA_real_),
    ci_lower = avg_flow - (qt(0.975, df = n_obs - 1) * se_flow),
    ci_upper = avg_flow + (qt(0.975, df = n_obs - 1) * se_flow)
  )

# ---- CSV出力 1: 季節別・河川別集計データ ----
file_stats <- file.path(SCRIPT_DIR, "seasonal_flow_stats.csv")
write_csv(seasonal_summary_stats, file_stats)
message(paste("集計データを保存:", file_stats))


# ---- データ整形（ワイド形式） ----
# アプローチA: 直接相関用
wide_flow_data <- seasonal_summary_stats %>%
  select(year, season_label, river, avg_flow) %>%
  pivot_wider(names_from = river, values_from = avg_flow) %>%
  arrange(season_label, year)

# アプローチB: トレンド除去後用
detrended_flow <- seasonal_summary_stats %>%
  group_by(river, season_label) %>%
  filter(n() > 2) %>%
  mutate(flow_residual = residuals(lm(avg_flow ~ year, data = cur_data()))) %>%
  ungroup() %>%
  select(year, season_label, river, flow_residual)

wide_detrended_data <- detrended_flow %>%
  pivot_wider(names_from = river, values_from = flow_residual)

# ---- 相関行列の計算 ----
# 関数: データフレームを受け取り相関行列を計算してリストに格納
calc_cor_matrix <- function(data) {
  cor(data %>% select(-year), use = "pairwise.complete.obs")
}

correlation_results_A <- wide_flow_data %>%
  group_by(season_label) %>%
  nest() %>%
  mutate(correlation_matrix = map(data, calc_cor_matrix))

correlation_results_B <- wide_detrended_data %>%
  group_by(season_label) %>%
  nest() %>%
  mutate(correlation_matrix = map(data, calc_cor_matrix))

# ---- CSV出力 2: 相関行列（ロング形式） ----
# 行列を行リスト形式（River1, River2, Correlation）に変換する関数
matrix_to_long <- function(cor_mat) {
  cor_mat %>%
    as.data.frame() %>%
    rownames_to_column(var = "river_1") %>%
    pivot_longer(cols = -river_1, names_to = "river_2", values_to = "correlation")
}

# A: 直接相関のCSV化
cor_long_A <- correlation_results_A %>%
  mutate(long_data = map(correlation_matrix, matrix_to_long)) %>%
  select(season_label, long_data) %>%
  unnest(long_data)

file_cor_A <- file.path(SCRIPT_DIR, "flow_correlation_direct.csv")
write_csv(cor_long_A, file_cor_A)
message(paste("相関データ(A)を保存:", file_cor_A))

# B: トレンド除去相関のCSV化
cor_long_B <- correlation_results_B %>%
  mutate(long_data = map(correlation_matrix, matrix_to_long)) %>%
  select(season_label, long_data) %>%
  unnest(long_data)

file_cor_B <- file.path(SCRIPT_DIR, "flow_correlation_detrended.csv")
write_csv(cor_long_B, file_cor_B)
message(paste("相関データ(B)を保存:", file_cor_B))


# ---- 相関行列プロット (PDF保存) ----
pdf_cor_file <- file.path(SCRIPT_DIR, "flowByRiver_synchrony_matrices.pdf") # パス修正
pdf(pdf_cor_file, width = 12, height = 12)
par(mfrow = c(2, 2))

# Aのプロット
correlation_results_A %>%
  pwalk(~ corrplot(..3, method = "number", type = "upper",
                   title = paste("A: Direct Correlation -", ..1),
                   mar = c(0, 0, 1, 0),
                   tl.cex = 0.8, number.cex = 0.8))

# Bのプロット
correlation_results_B %>%
  pwalk(~ corrplot(..3, method = "number", type = "upper",
                   title = paste("B: Detrended Correlation -", ..1),
                   mar = c(0, 0, 1, 0),
                   tl.cex = 0.8, number.cex = 0.8))
dev.off()
par(mfrow = c(1, 1))
message(paste("相関プロットを保存:", pdf_cor_file))


# ===== B. 季節別平均流量の時系列プロット =====
message("===== B. 季節別平均流量の可視化 =====")

# ---- プロット1: 各季節×河川別（詳細） ----
flow_plot_detailed <- ggplot(seasonal_summary_stats, aes(x = year, y = avg_flow, color = season_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  facet_grid(rows = vars(season_label), cols = vars(river), scales = "free_y") +
  labs(
    title = "Seasonal Average Flow (flowByRiver) by River and Season (2000–2015)",
    subtitle = "Sections 1–10. Error bars: 95% CI. Lines: linear regression.",
    x = "Year", y = "Average Flow (m³/s)", color = "Season"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "top"
  )

file_plot1 <- file.path(SCRIPT_DIR, "seasonal_flowByRiver_change_detailed.pdf") # パス修正
ggsave(file_plot1, plot = flow_plot_detailed, width = 14, height = 10, units = "in", device = "pdf")
message(paste("プロット1を保存:", file_plot1))

# ---- プロット2: 河川別・折れ線なし ----
flow_plot_by_river_no_lines <- ggplot(seasonal_summary_stats, aes(x = year, y = avg_flow, color = season_label)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.25, linewidth = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~river, scales = "free_y", ncol = 2) +
  labs(
    title = "Annual Average Flow Trends by River (flowByRiver, 2000–2015)",
    subtitle = "Sections 1–10. Points: seasonal averages. Error bars: 95% CI.",
    x = "Year", y = "Average Flow (m³/s)", color = "Season"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "top"
  )

file_plot2 <- file.path(SCRIPT_DIR, "seasonal_flowByRiver_by_river_no_lines.pdf") # パス修正
ggsave(file_plot2, plot = flow_plot_by_river_no_lines, width = 12, height = 8, units = "in", device = "pdf")
message(paste("プロット2を保存:", file_plot2))

message("===== 全ての解析と保存が完了しました =====")

# ---- P値と星付き相関を計算する関数の定義 ----
calc_cor_and_p <- function(data) {
  df_numeric <- data %>% select(-year)
  n_col <- ncol(df_numeric)
  
  # 相関係数行列
  r_mat <- cor(df_numeric, use = "pairwise.complete.obs")
  
  # P値行列の初期化
  p_mat <- matrix(NA, n_col, n_col)
  colnames(p_mat) <- colnames(r_mat)
  rownames(p_mat) <- colnames(r_mat)
  
  # 全ペアに対してcor.testを実行
  for (i in 1:n_col) {
    for (j in 1:n_col) {
      if (i == j) {
        p_mat[i, j] <- 0
      } else {
        try({
          test <- cor.test(df_numeric[[i]], df_numeric[[j]], method = "pearson")
          p_mat[i, j] <- test$p.value
        }, silent = TRUE)
      }
    }
  }
  
  # 星付きフォーマットに変換する内部関数
  format_sig <- function(r, p) {
    if (is.na(r) || is.na(p)) return(NA_character_)
    sig <- case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.1   ~ ".",
      TRUE      ~ ""
    )
    return(paste0(sprintf("%.2f", r), sig))
  }
  
  # フォーマット済み行列の作成
  formatted_mat <- matrix(
    mapply(format_sig, as.vector(r_mat), as.vector(p_mat)),
    nrow = n_col, dimnames = list(rownames(r_mat), colnames(r_mat))
  )
  
  return(list(r = r_mat, p = p_mat, formatted = formatted_mat))
}

# ---- 実行: 季節ごとに相関とP値を計算 ----
correlation_results_A <- wide_flow_data %>%
  group_by(season_label) %>%
  nest() %>%
  mutate(res = map(data, calc_cor_and_p))

correlation_results_B <- wide_detrended_data %>%
  group_by(season_label) %>%
  nest() %>%
  mutate(res = map(data, calc_cor_and_p))

# ---- CSV出力用: リスト構造をロング形式に変換する関数 ----
extract_to_long <- function(res_list, type = "formatted") {
  mat <- res_list[[type]]
  as.data.frame(mat) %>%
    rownames_to_column(var = "river_1") %>%
    pivot_longer(cols = -river_1, names_to = "river_2", values_to = type)
}

# ---- CSV保存処理 (A: 直接相関) ----
# 1. 星付き相関
cor_sig_A <- correlation_results_A %>%
  mutate(long_data = map(res, ~extract_to_long(.x, "formatted"))) %>%
  select(season_label, long_data) %>% unnest(long_data)
write_csv(cor_sig_A, file.path(SCRIPT_DIR, "flow_correlation_direct_with_sig.csv"))

# 2. P値のみ
cor_p_A <- correlation_results_A %>%
  mutate(long_data = map(res, ~extract_to_long(.x, "p"))) %>%
  select(season_label, long_data) %>% unnest(long_data)
write_csv(cor_p_A, file.path(SCRIPT_DIR, "flow_correlation_direct_p_values.csv"))

# ---- CSV保存処理 (B: トレンド除去後) ----
# 1. 星付き相関
cor_sig_B <- correlation_results_B %>%
  mutate(long_data = map(res, ~extract_to_long(.x, "formatted"))) %>%
  select(season_label, long_data) %>% unnest(long_data)
write_csv(cor_sig_B, file.path(SCRIPT_DIR, "flow_correlation_detrended_with_sig.csv"))

# 2. P値のみ
cor_p_B <- correlation_results_B %>%
  mutate(long_data = map(res, ~extract_to_long(.x, "p"))) %>%
  select(season_label, long_data) %>% unnest(long_data)
write_csv(cor_p_B, file.path(SCRIPT_DIR, "flow_correlation_detrended_p_values.csv"))

message("P値を含むCSVの保存が完了しました。")

