# =============================================================================
# 0. Environment Reset & Package Loading
# =============================================================================
rm(list = ls(all.names = TRUE))

packages <- c("readr", "dplyr", "lubridate", "ggplot2", "tidyr", "scales", 
              "broom", "purrr", "nnet", "emmeans", "patchwork", "ggh4x", "RColorBrewer")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(SCRIPT_DIR)
} else {
  SCRIPT_DIR <- getwd()
}

# =============================================================================
# 1. Basic Settings & Data Loading
# =============================================================================
# File path
file_path <- "../Clustering/WB_bkt_clustered_rev2.csv"
OUTPUT_FILENAME <- "Combined_Cluster_Dynamics_k4_vs_k9.pdf"

years_wb <- 2000:2015     
years_others <- 2002:2015 
min_year <- min(years_wb) 

river_levels <- c("Metapopulation", "WB", "OL", "OS", "IL")

# -----------------------------------------------------------------------------
# カラーパレットのカスタマイズ (クラスターと色の対応を入れ替え)
# -----------------------------------------------------------------------------
# Set1の全9色を取得
# (1:赤, 2:青, 3:緑, 4:紫, 5:オレンジ, 6:黄, 7:茶, 8:ピンク, 9:グレー)
set1_colors <- brewer.pal(9, "Set1")

# k=4 の色の割り当て (例: cluster 1 = 4番目[紫], cluster 2 = 3番目[緑], 3=1番目[赤], 4=2番目[青])
order_k4 <- c(4, 3, 1, 2)
COLORS_K4 <- set1_colors[order_k4]

# k=9 の色の割り当て (例: cluster 1 = 4番目, cluster 2 = 3番目, 残りは適宜割り当て)
order_k9 <- c(4, 3, 1, 2, 5, 6, 7, 8, 9)
COLORS_K9 <- set1_colors[order_k9]
# -----------------------------------------------------------------------------

# 共通テーマ
shared_theme <- theme_bw(base_size = 20) +
  theme(
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold", color = "black"),
    strip.background = element_blank(),  
    panel.border = element_blank(),                        
    axis.line = element_line(color = "black", linewidth = 1), 
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = "bottom"
  )

# 河川名の統一関数
rename_river_func <- function(df, col_name) {
  df %>%
    mutate(!!sym(col_name) := case_when(
      !!sym(col_name) == "west brook"  ~ "WB",
      !!sym(col_name) == "wb jimmy"    ~ "OL",
      !!sym(col_name) == "wb mitchell" ~ "OS",
      !!sym(col_name) == "wb obear"    ~ "IL",
      !!sym(col_name) == "All Rivers"  ~ "Metapopulation", 
      TRUE                             ~ as.character(!!sym(col_name))
    )) %>%
    filter(!!sym(col_name) %in% river_levels) %>%
    mutate(!!sym(col_name) := factor(!!sym(col_name), levels = river_levels))
}

# データ読み込み
df_raw <- read_csv(file_path, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

# =============================================================================
# 2. グラフ作成関数の定義 (データ処理・モデリング・作図を一括実行)
# =============================================================================
create_cluster_plot <- function(data_raw, k_val, cluster_colors, is_left_plot = TRUE) {
  
  TARGET_CLUSTER_COL <- paste0("cluster_k", k_val)
  ref_cluster <- as.character(k_val)
  
  # ターゲットとなるクラスター列の抽出とファクタ化
  df_k <- data_raw %>%
    rename(ocluster_rename_dynamic = !!sym(TARGET_CLUSTER_COL)) %>%
    mutate(ocluster_rename_dynamic = relevel(factor(ocluster_rename_dynamic, levels = as.character(1:k_val)), ref = ref_cluster))
  
  # Multiタグの特定
  river_tags <- df_k %>% filter(cohort %in% years_wb) %>% distinct(tag, river)
  multi_tags <- river_tags %>% count(tag) %>% filter(n > 1) %>% pull(tag)
  
  df_model_base <- df_k %>%
    mutate(
      river_grp = if_else(tag %in% multi_tags, "Multi", river),
      cohort_numeric = as.numeric(as.character(cohort))
    )
  
  df_model_all <- df_model_base %>% mutate(river_grp = "All Rivers")
  df_model <- bind_rows(df_model_base, df_model_all) %>%
    filter(
      (river_grp == "west brook" & cohort_numeric %in% years_wb) |
        (river_grp != "west brook" & cohort_numeric %in% years_others)
    ) %>%
    mutate(year_centered = cohort_numeric - min_year)
  
  # 実測割合の計算
  proportion_data <- df_model %>%
    group_by(river_grp, cohort, cohort_numeric, year_centered, ocluster_rename_dynamic) %>%
    summarise(n_cluster_tags = n_distinct(tag), .groups = "drop") %>%
    group_by(river_grp, cohort) %>%
    mutate(n_total_tags = sum(n_cluster_tags)) %>%
    ungroup() %>%
    mutate(proportion = n_cluster_tags / n_total_tags)
  
  plot_data_full <- proportion_data %>%
    complete(river_grp, cohort_numeric = years_wb, ocluster_rename_dynamic, fill = list(proportion = 0)) %>%
    mutate(
      year_centered = cohort_numeric - min_year,
      proportion = if_else(river_grp != "west brook" & cohort_numeric < 2002, NA_real_, proportion)
    )
  
  # 多様度指数 (Metapopulationのみ計算)
  simpson_data <- plot_data_full %>%
    filter(!is.na(proportion)) %>% 
    group_by(river_grp, cohort_numeric) %>%
    summarise(simpson_index = 1 - sum(proportion^2, na.rm = TRUE), .groups = "drop") %>%
    filter(river_grp == "All Rivers")
  
  # モデリング1: All Rivers (独立)
  multinom_results_all <- df_model %>%
    filter(river_grp == "All Rivers") %>%
    group_by(river_grp) %>% nest() %>%
    mutate(
      model = map(data, ~multinom(ocluster_rename_dynamic ~ year_centered, data = ., trace = FALSE)),
      predicted = map(model, ~{
        new_data <- tibble(cohort_numeric = years_others, year_centered = years_others - min_year)
        predict(., newdata = new_data, type = "probs") %>%
          as.data.frame() %>% mutate(cohort_numeric = new_data$cohort_numeric) %>%
          pivot_longer(cols = -cohort_numeric, names_to = "ocluster_rename_dynamic", values_to = "predicted_prop")
      })
    )
  
  predicted_data_all <- multinom_results_all %>% select(river_grp, predicted) %>% unnest(predicted)
  plot_data_all <- plot_data_full %>%
    filter(river_grp == "All Rivers") %>%
    left_join(predicted_data_all, by = c("river_grp", "cohort_numeric", "ocluster_rename_dynamic"))
  
  # モデリング2: 支流 (交互作用)
  df_interaction <- df_model %>% filter(river_grp != "All Rivers")
  
  full_model <- multinom(
    ocluster_rename_dynamic ~ year_centered * river_grp, 
    data = df_interaction, 
    trace = FALSE, 
    maxit = 1000, 
    Hess = TRUE
  )
  
  new_data_interaction <- plot_data_full %>% 
    filter(river_grp != "All Rivers") %>% 
    filter((river_grp == "west brook" & cohort_numeric %in% years_wb) |
             (river_grp != "west brook" & cohort_numeric %in% years_others)) %>%
    distinct(cohort_numeric, year_centered, river_grp)
  
  predicted_interaction <- as.data.frame(predict(full_model, newdata = new_data_interaction, type = "probs")) %>%
    mutate(cohort_numeric = new_data_interaction$cohort_numeric, river_grp = new_data_interaction$river_grp) %>%
    pivot_longer(cols = -c(cohort_numeric, river_grp), names_to = "ocluster_rename_dynamic", values_to = "predicted_prop")
  
  plot_data_interaction <- plot_data_full %>%
    filter(river_grp != "All Rivers") %>%
    left_join(predicted_interaction, by = c("river_grp", "cohort_numeric", "ocluster_rename_dynamic"))
  
  # データの結合と名称統一
  data_vertical <- bind_rows(
    plot_data_all %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop),
    plot_data_interaction %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop)
  )
  
  data_vertical <- rename_river_func(data_vertical, "river_grp")
  simpson_data  <- rename_river_func(simpson_data, "river_grp")
  
  n_labels <- df_model %>%
    rename_river_func("river_grp") %>%
    group_by(river_grp) %>%
    summarise(total_n = n_distinct(tag), .groups = "drop") %>%
    mutate(label = paste0("n = ", total_n))
  
  # X軸のスケール設定 (パッチワーク用)
  x_scale_show <- scale_x_continuous(breaks = seq(2000, 2015, by = 2), limits = c(min(years_wb), max(years_wb)))
  x_scale_hide <- scale_x_continuous(breaks = seq(2000, 2015, by = 2), limits = c(min(years_wb), max(years_wb)), labels = NULL)
  x_scales_list <- list(x_scale_show, x_scale_hide, x_scale_hide, x_scale_hide, x_scale_show)
  
  # 作図
  p <- ggplot(data_vertical, aes(x = cohort_numeric)) +
    geom_point(aes(y = proportion, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), size = 3, alpha = 0.6) +
    geom_line(aes(y = predicted_prop, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), linewidth = 2) +
    
    # Diversity Index (Metapopulationのみに描画される)
    geom_line(data = simpson_data, aes(y = simpson_index), color = "black", linetype = "dashed", linewidth = 1.2) +
    geom_point(data = simpson_data, aes(y = simpson_index), color = "black", shape = 17, size = 4) +
    
    facet_wrap(~ river_grp, ncol = 1, scales = "free_x") + 
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    ggh4x::facetted_pos_scales(x = x_scales_list) + 
    geom_text(data = n_labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.2, vjust = 1.5, size = 6, fontface = "bold", color = "black") +
    
    labs(
      x = "Year (Cohort)", 
      y = if(is_left_plot) "Proportion / Diversity Index" else "", # 右側のグラフはY軸ラベルを省略してスッキリさせる
      color = paste0("Cluster (k=", k_val, ")")
    ) +
    scale_color_manual(
      values = cluster_colors[1:k_val],
      # ★ 5個で折り返して2行にする設定 (byrow = TRUE で左から右へ並べる)
      guide = guide_legend(ncol = 5, byrow = TRUE) 
    ) +
    shared_theme +
    ggh4x::force_panelsizes(rows = c(1, 1, 1, 1, 1))
  
  return(p)
}

# =============================================================================
# 3. グラフ生成と結合
# =============================================================================
cat("\nProcessing k=4...\n")
p_left <- create_cluster_plot(df_raw, k_val = 4, cluster_colors = COLORS_K4, is_left_plot = TRUE)

cat("Processing k=9...\n")
p_right <- create_cluster_plot(df_raw, k_val = 9, cluster_colors = COLORS_K9, is_left_plot = FALSE)

# 左右のグラフを結合
# guides="keep"にすることで、それぞれのクラスター凡例（k=4とk=9）を並べて表示させます
combined_plot <- p_left + p_right + 
  plot_layout(ncol = 2, guides = "keep") & 
  theme(legend.position = "bottom")

# =============================================================================
# 4. Output to PDF
# =============================================================================
n_panels <- length(river_levels)
PANEL_HEIGHT <- 3.5 
total_pdf_height <- (PANEL_HEIGHT * n_panels) + 2

output_path <- file.path(SCRIPT_DIR, OUTPUT_FILENAME)
ggsave(filename = output_path, plot = combined_plot, width = 14, height = total_pdf_height)

cat("\n✅ Saved successfully: ", output_path, "\n")

