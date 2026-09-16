# ===================================================================
# Westbrook Portfolio Effect Analysis (Combined: Static & Dynamic)
  # ===================================================================
# ===== 0. Reset environment & User Configuration =====
rm(list = ls(all.names = TRUE))

CLUSTER_NUM <- 3
WINDOW_SIZE <- 10
TARGET_CLUSTER_COL <- paste0("cluster_k", CLUSTER_NUM)
CLUSTER_IDS <- c("1", "2", "3") # Adjust according to CLUSTER_NUM

# 色指定を名前付きベクトルに変更（一部のクラスターのみを描画する際に色がずれないようにするため）
CLUSTER_COLORS <- c("1" = "steelblue", "2" = "aquamarine3", "3" = "salmon3")

INPUT_FILE <- "../csv/WB_combined_bkt_rev2.csv"
OUTPUT_PDF <- paste0("Westbrook_Combined_All_Plots_", WINDOW_SIZE, "yr_window_k", CLUSTER_NUM, ".pdf")

# ===== 1. Package Loading =====
packages <- c("readr", "dplyr", "lubridate", "tidyr", "ggplot2", "ggrepel", "patchwork")

if (!require("rstudioapi")) install.packages("rstudioapi") 
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ===== 2. Set Working Directory to Script Location =====
if (rstudioapi::isAvailable()) {
  try({
    current_path <- rstudioapi::getSourceEditorContext()$path
    if (current_path != "") {
      setwd(dirname(current_path))
      message("Working directory set to: ", getwd())
    }
  }, silent = TRUE)
}

# ===== 3. Data Loading and Preprocessing (Common) =====
data <- read_csv(INPUT_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

if (!TARGET_CLUSTER_COL %in% names(data)) {
  stop(paste("Error: Column", TARGET_CLUSTER_COL, "not found in the data."))
}

# Rename cluster column and recode river names
data <- data %>%
  rename(cluster_id = all_of(TARGET_CLUSTER_COL)) %>%
  mutate(
    year = year(date),
    river = case_when(
      river == "west brook"  ~ "WB",
      river == "wb jimmy"    ~ "OL",
      river == "wb mitchell" ~ "OS",
      river == "wb obear"    ~ "IL",
      TRUE                   ~ river
    )
  ) %>%
  filter(!is.na(cluster_id))

# Expand residence period to seasonal units and filter S3 (2002-2015)
alive_by_season_filtered <- data %>%
  arrange(tag, date) %>%
  group_by(tag) %>%
  mutate(
    next_date = lead(date),
    end_date = if_else(!is.na(next_date), next_date, date)
  ) %>%
  rowwise() %>%
  mutate(season_seq = list({
    dseq <- seq(date, end_date, by = "day")
    unique(paste0(year(dseq), "-S", quarter(dseq)))
  })) %>%
  unnest(season_seq) %>%
  ungroup() %>%
  select(tag, river, cluster_id, season_seq) %>%
  distinct() %>%
  group_by(river, cluster_id, season_seq) %>%
  summarise(n_alive = n_distinct(tag), .groups = "drop") %>%
  mutate(year = as.numeric(substr(season_seq, 1, 4))) %>%
  filter(year >= 2002 & year <= 2015) %>%
  filter(grepl("-S3$", season_seq))

# ===================================================================
# Analysis Part 1: Static Portfolio Effect Decomposition
# ===================================================================
cv_by_river_cluster <- alive_by_season_filtered %>%
  group_by(river, cluster_id) %>%
  summarise(mean_alive = mean(n_alive), cv = sd(n_alive) / mean(n_alive), .groups = 'drop')

cv_by_river <- alive_by_season_filtered %>%
  group_by(river, season_seq) %>%
  summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
  group_by(river) %>%
  summarise(mean_alive = mean(n_alive_sum), cv = sd(n_alive_sum) / mean(n_alive_sum), .groups = 'drop')

cv_by_cluster <- alive_by_season_filtered %>%
  group_by(cluster_id, season_seq) %>%
  summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
  group_by(cluster_id) %>%
  summarise(mean_alive = mean(n_alive_sum), cv = sd(n_alive_sum) / mean(n_alive_sum), .groups = 'drop')

cv_total <- alive_by_season_filtered %>%
  group_by(season_seq) %>%
  summarise(total_n_alive = sum(n_alive), .groups = "drop") %>%
  summarise(cv = sd(total_n_alive) / mean(total_n_alive)) %>%
  pull(cv)

cv_w_mean_river_cluster <- weighted.mean(cv_by_river_cluster$cv, w = cv_by_river_cluster$mean_alive, na.rm = TRUE)
cv_w_mean_river <- weighted.mean(cv_by_river$cv, w = cv_by_river$mean_alive, na.rm = TRUE)
cv_w_mean_cluster <- weighted.mean(cv_by_cluster$cv, w = cv_by_cluster$mean_alive, na.rm = TRUE)

pe_cluster_effect <- cv_w_mean_river_cluster - cv_w_mean_river
pe_river_effect <- cv_w_mean_river_cluster - cv_w_mean_cluster

# ===================================================================
# Helper Function: Partial Spearman Correlation
# ===================================================================
# 裾依存的同期性（Tail-dependent synchrony）を計算するためのヘルパー関数
calc_partial_spearman <- function(x, y, b_l, b_u) {
  # 経験的累積分布関数（ECDF）に基づくランク（分位数）の計算
  u <- rank(x, ties.method = "average") / (length(x) + 1)
  v <- rank(y, ties.method = "average") / (length(y) + 1)
  
  # 指定された分位数範囲に両方の変数が該当するインデックスを抽出
  idx <- which(u >= b_l & u <= b_u & v >= b_l & v <= b_u)
  
  # サンプルサイズが小さすぎる場合はNAを返す
  if (length(idx) < 3) return(NA_real_)
  
  return(cor(x[idx], y[idx], method = "spearman"))
}

# ===================================================================
# Analysis Part 2: Dynamic Portfolio Effect & Synchrony (Sliding Window)
# ===================================================================
min_year <- min(alive_by_season_filtered$year)
max_year <- max(alive_by_season_filtered$year)
start_years <- min_year:(max_year - WINDOW_SIZE + 1)
windows <- lapply(start_years, function(y) y:(y + WINDOW_SIZE - 1))
names(windows) <- paste0(start_years, "-", start_years + WINDOW_SIZE - 1)

sync_base <- alive_by_season_filtered %>%
  group_by(cluster_id, year) %>%
  summarise(n_alive = sum(n_alive), .groups = "drop")

r_one_vs_others_all <- bind_rows(lapply(CLUSTER_IDS, function(cid) {
  sync_base %>%
    mutate(is_target = if_else(cluster_id == cid, "Target", "Others")) %>%
    group_by(is_target, year) %>%
    summarise(n_alive = sum(n_alive), .groups = "drop") %>%
    arrange(is_target, year) %>%
    group_by(is_target) %>%
    mutate(
      n_prev = lag(n_alive, order_by = year),
      year_diff = year - lag(year, order_by = year),
      r_rate = if_else(year_diff == 1 & n_prev > 0 & n_alive > 0, log(n_alive / n_prev), NA_real_)
    ) %>%
    filter(!is.na(r_rate)) %>%
    mutate(target_cluster = cid) %>%
    ungroup()
}))

contribution_results <- list(); sync_results <- list(); pe_results <- list()

for (w_name in names(windows)) {
  target_years <- windows[[w_name]]
  window_data <- alive_by_season_filtered %>% filter(year %in% target_years)
  if (n_distinct(window_data$season_seq) < 2) next
  
  # Analysis A: Portfolio Effect
  cv_rc <- window_data %>% group_by(river, cluster_id) %>% summarise(m = mean(n_alive), cv = sd(n_alive) / mean(n_alive), .groups = 'drop')
  cv_r <- window_data %>% group_by(river, season_seq) %>% summarise(n = sum(n_alive), .groups = "drop") %>% group_by(river) %>% summarise(m = mean(n), cv = sd(n) / mean(n), .groups = 'drop')
  cv_c <- window_data %>% group_by(cluster_id, season_seq) %>% summarise(n = sum(n_alive), .groups = "drop") %>% group_by(cluster_id) %>% summarise(m = mean(n), cv = sd(n) / mean(n), .groups = 'drop')
  cv_t <- window_data %>% group_by(season_seq) %>% summarise(tot = sum(n_alive), .groups = 'drop') %>% summarise(cv = sd(tot) / mean(tot)) %>% pull(cv)
  
  wm_rc <- weighted.mean(cv_rc$cv, w = cv_rc$m, na.rm = TRUE)
  wm_r  <- weighted.mean(cv_r$cv, w = cv_r$m, na.rm = TRUE)
  wm_c  <- weighted.mean(cv_c$cv, w = cv_c$m, na.rm = TRUE)
  
  pe_results[[w_name]] <- data.frame(Window = w_name, `Cluster Effect` = wm_rc - wm_r, `River Effect` = wm_rc - wm_c, `Overall Effect` = wm_rc - cv_t, check.names = FALSE)
  
  # Analysis B: Cluster Contribution Effect
  for (cid in CLUSTER_IDS) {
    temp_cv <- window_data %>%
      mutate(merged_cluster = if_else(cluster_id == cid, as.character(cid), "Others")) %>%
      group_by(river, merged_cluster, season_seq) %>% summarise(n = sum(n_alive), .groups = "drop") %>%
      group_by(river, merged_cluster) %>% summarise(m = mean(n), cv = sd(n) / mean(n), .groups = "drop")
    
    wm_temp <- weighted.mean(temp_cv$cv, w = temp_cv$m, na.rm = TRUE)
    contribution_results[[paste(w_name, cid)]] <- data.frame(Window = w_name, cluster_id = cid, contribution_effect = wm_temp - wm_r)
  }
  
  # Analysis C: Tail-dependent Synchrony (Cluster 3 vs Othersのみに限定)
  window_r_data <- r_one_vs_others_all %>% filter(year %in% target_years)
  for (cid in c("3")) { 
    df_wide <- window_r_data %>% filter(target_cluster == cid) %>% pivot_wider(id_cols = year, names_from = is_target, values_from = r_rate)
    
    r_val <- NA_real_; p_val <- NA_real_
    r_upper <- NA_real_; r_lower <- NA_real_; tail_dep_strength <- NA_real_
    
    if (all(c("Target", "Others") %in% names(df_wide))) {
      valid_data <- df_wide %>% filter(is.finite(Target) & is.finite(Others))
      if (nrow(valid_data) >= 3) {
        # 従来の全体のピアソン相関
        test_res <- cor.test(valid_data$Target, valid_data$Others, method = "pearson")
        r_val <- test_res$estimate
        p_val <- test_res$p.value
        
        # 裾依存的な部分スピアマン相関の計算
        r_lower <- calc_partial_spearman(valid_data$Target, valid_data$Others, 0.0, 0.5)
        r_upper <- calc_partial_spearman(valid_data$Target, valid_data$Others, 0.5, 1.0)
        
        # 裾依存性の強度 (Upper - Lower)
        if(!is.na(r_upper) & !is.na(r_lower)) {
          tail_dep_strength <- r_upper - r_lower
        }
      }
    }
    sync_results[[paste(w_name, cid)]] <- data.frame(
      Window = w_name, 
      cluster_id = cid, 
      correlation = r_val, 
      p_value = p_val,
      r_lower = r_lower,
      r_upper = r_upper,
      tail_dependence = tail_dep_strength
    )
  }
}

df_contrib <- bind_rows(contribution_results) %>% mutate(Window = factor(Window, levels = names(windows)))
df_sync <- bind_rows(sync_results) %>% 
  mutate(Window = factor(Window, levels = names(windows)),
         p_adj = p.adjust(p_value, method = "BH"),
         sig_adj = case_when(is.na(p_adj) ~ "", p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", p_adj < 0.1 ~ ".", TRUE ~ ""))
df_pe <- bind_rows(pe_results) %>% mutate(Window = factor(Window, levels = names(windows))) %>%
  pivot_longer(cols = c("Cluster Effect", "River Effect", "Overall Effect"), names_to = "effect_type", values_to = "pe_value")

# ===================================================================
# Visualization Part
# ===================================================================
total_mean_abundance <- sum(cv_by_river_cluster$mean_alive, na.rm = TRUE)
river_cluster_plot_data <- cv_by_river_cluster %>% mutate(proportion = mean_alive / total_mean_abundance, cluster_label = paste0("c", cluster_id))
max_prop <- max(river_cluster_plot_data$proportion, na.rm = TRUE)
pos_vline <- max_prop * 1.1; pos_summary_pts <- max_prop * 1.2; pos_arrow_clus <- max_prop * 1.15; pos_arrow_river <- max_prop * 1.6; x_limit_max <- max_prop * 1.8

summary_metrics <- data.frame(
  metric = factor(c("Expected CV (No Diversity)", "Expected CV (Cluster Diversity Only)", "Expected CV (River Diversity Only)", "Observed Total CV"),
                  levels = c("Expected CV (No Diversity)", "Expected CV (Cluster Diversity Only)", "Expected CV (River Diversity Only)", "Observed Total CV")),
  cv = c(cv_w_mean_river_cluster, cv_w_mean_river, cv_w_mean_cluster, cv_total),
  x_pos = pos_summary_pts
)

ref_lines <- data.frame(y = c(cv_w_mean_river_cluster, cv_w_mean_river, cv_w_mean_cluster), x_start = pos_vline, x_end = x_limit_max)

# Dynamic Theme Setting (共通のL字テーマ・文字を大きく・斜めを解除)
custom_theme <- theme_classic(base_size = 18) +
  theme(axis.text = element_text(color = "black", face = "bold", size = 16),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", face = "bold", size = 16), 
        axis.title = element_text(color = "black", face = "bold", size = 20),
        axis.line = element_line(linewidth = 1, color = "black"), 
        legend.position = "bottom", 
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 22))

# Plot a
p_static <- ggplot() +
  geom_segment(data = ref_lines, aes(x = x_start, xend = x_end, y = y, yend = y), linetype = "dashed", color = "gray60", size = 0.6) +
  geom_vline(xintercept = pos_vline, linetype = "dotted", color = "grey50", size = 0.8) +
  geom_point(data = river_cluster_plot_data, aes(x = proportion, y = cv, color = river), size = 10, alpha = 0.8) +
  geom_text(data = river_cluster_plot_data, aes(x = proportion, y = cv, label = cluster_label), color = "white", size = 7, fontface = "bold") +
  geom_point(data = summary_metrics, aes(x = x_pos, y = cv), size = 8, shape = 18, color = "black") +
  geom_text(data = summary_metrics, aes(x = x_pos + 0.008, y = cv, label = metric), hjust = 0, size = 6, color = "black") +
  annotate("segment", x = pos_arrow_clus, xend = pos_arrow_clus, y = cv_w_mean_river_cluster, yend = cv_w_mean_river, arrow = arrow(length = unit(0.2, "cm"), ends = "both"), colour = "blue", size = 1) +
  annotate("text", x = pos_arrow_clus + 0.008, y = (cv_w_mean_river_cluster + cv_w_mean_river) / 2, label = paste0("Cluster Effect\n(ΔCV = ", round(pe_cluster_effect, 3), ")"), color = "blue", hjust = 0, size = 5, fontface = "bold") +
  annotate("segment", x = pos_arrow_river, xend = pos_arrow_river, y = cv_w_mean_river_cluster, yend = cv_w_mean_cluster, arrow = arrow(length = unit(0.2, "cm"), ends = "both"), colour = "red", size = 1) +
  annotate("text", x = pos_arrow_river + 0.008, y = (cv_w_mean_river_cluster + cv_w_mean_cluster) / 2, label = paste0("River Effect\n(ΔCV = ", round(pe_river_effect, 3), ")"), color = "red", hjust = 0, size = 5, fontface = "bold") +
  labs(title = "a", x = "Proportion of Total Mean Abundance", y = "Coefficient of Variation (CV)", color = "River") +
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(color = "black", face = "bold", size = 16),
        axis.title = element_text(color = "black", face = "bold", size = 20),
        axis.line = element_line(linewidth = 1, color = "black"),
        plot.title = element_text(face = "bold", size = 22),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 16)) +
  coord_cartesian(xlim = c(0, x_limit_max))

# Plot c
p_contrib <- ggplot(df_contrib, aes(x = Window, y = contribution_effect, group = cluster_id, color = as.factor(cluster_id))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = CLUSTER_COLORS) + custom_theme +
  labs(title = "c", x = paste0(WINDOW_SIZE, "-Year Window"), y = "CV Reduction", color = "Cluster")

# Plot d (Tail Dependence Plotに変更)
p_sync <- ggplot(df_sync, aes(x = Window, y = tail_dependence, group = cluster_id, color = as.factor(cluster_id))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_line(linewidth = 1.2) + 
  geom_point(size = 4) +
  scale_color_manual(values = CLUSTER_COLORS) + 
  custom_theme +
  coord_cartesian(ylim = c(-1.0, 1.0)) + 
  labs(
    title = "d: Tail Dependence in Synchrony", 
    x = paste0(WINDOW_SIZE, "-Year Window"), 
    y = "Tail Dependence Strength (Upper - Lower)", 
    color = "Cluster"
  )

# Plot b
p_pe <- ggplot(df_pe, aes(x = Window, y = pe_value, group = effect_type, shape = effect_type, linetype = effect_type)) +
  geom_line(linewidth = 1.2, color = "black") + geom_point(size = 4, color = "black", fill = "black", stroke = 1.2) +
  scale_shape_manual(values = c("Cluster Effect" = 22, "River Effect" = 21, "Overall Effect" = 24)) +
  scale_linetype_manual(values = c("Cluster Effect" = "solid", "River Effect" = "solid", "Overall Effect" = "solid")) +
  custom_theme +
  labs(title = "b", x = paste0(WINDOW_SIZE, "-Year Window"), y = "CV Reduction", shape = "Effect Type", linetype = "Effect Type")

# ===================================================================
# Combine and Save Plot (Using patchwork)
# ===================================================================
cat("\nCombining plots and saving to:", OUTPUT_PDF, "\n")

# Combine 4 plots into a 2x2 grid
combined_plot <- (p_static + p_pe) / (p_contrib + p_sync)

ggsave(
  filename = OUTPUT_PDF,
  plot = combined_plot,
  device = cairo_pdf,
  width = 20,
  height = 14,
  units = "in",
  dpi = 300
)

message("All Analysis Completed Successfully.")

