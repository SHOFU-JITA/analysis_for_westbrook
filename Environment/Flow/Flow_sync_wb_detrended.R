# =============================================================================
# 流域流量データの相関分析スクリプト (トレンド除去・P値付き出力版)
# =============================================================================

# 0. 環境リセット
rm(list = ls(all.names = TRUE))

# 1. パッケージ読み込み
library(dplyr)
library(lubridate)
library(readr)
library(tidyr)
library(purrr)
library(corrplot)
library(ggplot2)
library(tibble)
library(scales)

# =============================================================================
# 2. 基本設定と関数定義
# =============================================================================

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}

# --- 関数: 相関係数と補正済みP値の計算 ---
calc_cor_and_p <- function(data, adjust_method = "BH") {
  df_numeric <- data %>% select(-year)
  n_col <- ncol(df_numeric)
  r_mat <- cor(df_numeric, use = "pairwise.complete.obs")
  
  # 生のP値行列を作成
  p_mat_raw <- matrix(NA, n_col, n_col, dimnames = list(colnames(r_mat), colnames(r_mat)))
  
  for (i in 1:n_col) {
    for (j in 1:n_col) {
      if (i == j) {
        p_mat_raw[i, j] <- NA 
      } else {
        try({
          test <- cor.test(df_numeric[[i]], df_numeric[[j]], method = "pearson")
          p_mat_raw[i, j] <- test$p.value
        }, silent = TRUE)
      }
    }
  }
  
  # 多重比較の補正 (BH法)
  p_vec_adj <- p.adjust(as.vector(p_mat_raw), method = adjust_method)
  p_mat_adj <- matrix(p_vec_adj, nrow = n_col, dimnames = dimnames(p_mat_raw))
  diag(p_mat_adj) <- 0 
  
  # 星付きフォーマット用関数
  format_sig <- function(r, p) {
    if (is.na(r) || is.na(p)) return(NA_character_)
    sig <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
    paste0(sprintf("%.2f", r), sig)
  }
  
  formatted_mat <- matrix(
    mapply(format_sig, as.vector(r_mat), as.vector(p_mat_adj)),
    nrow = n_col, dimnames = dimnames(r_mat)
  )
  
  # 全ての結果をリストで返す
  return(list(r = r_mat, p_raw = p_mat_raw, p_adj = p_mat_adj, formatted = formatted_mat))
}

# --- 関数: 行列をロング形式に変換 (内部処理用) ---
matrix_to_long_internal <- function(mat, value_name) {
  as.data.frame(mat) %>%
    rownames_to_column(var = "river_1") %>%
    pivot_longer(cols = -river_1, names_to = "river_2", values_to = value_name)
}

# =============================================================================
# 3. データ処理と解析
# =============================================================================
csv_file_path <- "~/Westbrook_カワマスP/code_SF/rawdata/envDataWB.csv"
if (!file.exists(csv_file_path)) stop("CSVファイルが見つかりません。")

df <- read_csv(csv_file_path) %>%
  mutate(date = as.Date(date), year = year(date)) %>%
  filter(year >= 2000 & year <= 2015) %>%
  mutate(
    month_num = month(date),
    season_label = factor(case_when(
      month_num %in% 3:5 ~ "Spring", month_num %in% 6:8 ~ "Summer",
      month_num %in% 9:11 ~ "Autumn", TRUE ~ "Winter"
    ), levels = c("Spring", "Summer", "Autumn", "Winter"))
  )

seasonal_summary <- df %>%
  filter(!is.na(flowByRiver)) %>%
  group_by(river, year, season_label) %>%
  summarise(avg_flow = mean(flowByRiver, na.rm = TRUE), n = n(), sd = sd(flowByRiver, na.rm = TRUE), .groups = "drop") %>%
  mutate(se = sd / sqrt(n), ci_low = avg_flow - (qt(0.975, df = n-1) * se), ci_up = avg_flow + (qt(0.975, df = n-1) * se))

# トレンド除去データの作成
wide_detrended <- seasonal_summary %>%
  group_by(river, season_label) %>%
  filter(n() > 2) %>%
  mutate(resid = residuals(lm(avg_flow ~ year))) %>%
  ungroup() %>%
  select(year, season_label, river, resid) %>%
  pivot_wider(names_from = river, values_from = resid)

# 解析実行
analysis_results <- wide_detrended %>%
  group_by(season_label) %>%
  nest() %>%
  mutate(res = map(data, ~calc_cor_and_p(.x, adjust_method = "BH")))

# =============================================================================
# 4. CSV出力 (重複ペアを除去してP値を含める)
# =============================================================================

csv_export_data <- analysis_results %>%
  mutate(long_combined = map(res, function(res_obj) {
    # 各要素をロング形式にする
    df_r   <- matrix_to_long_internal(res_obj$r, "r_value")
    df_fmt <- matrix_to_long_internal(res_obj$formatted, "r_with_sig")
    df_raw <- matrix_to_long_internal(res_obj$p_raw, "p_raw")
    df_adj <- matrix_to_long_internal(res_obj$p_adj, "p_adjusted_BH")
    
    df_r %>%
      left_join(df_fmt, by = c("river_1", "river_2")) %>%
      left_join(df_raw, by = c("river_1", "river_2")) %>%
      left_join(df_adj, by = c("river_1", "river_2")) %>%
      # --- ここで重複と自己相関を排除 ---
      # アルファベット順で river_1 < river_2 となるものだけを残すことで、
      # 「自己相関」と「入れ替わりの重複」を同時に消去できます
      filter(as.character(river_1) < as.character(river_2))
  })) %>%
  select(season_label, long_combined) %>%
  unnest(long_combined)

# ファイルに書き出し
write_csv(csv_export_data, file.path(SCRIPT_DIR, "flow_correlation_detrended_unique_pairs.csv"))

message("CSV出力完了: 自己相関および重複ペアを排除しました。")