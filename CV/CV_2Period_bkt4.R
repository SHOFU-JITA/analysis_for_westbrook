# ==============================================================================
# Brook Trout: Portfolio Effect Analysis (2-period comparison, lightweight)
# ==============================================================================

# ===== 0. Reset environment =====
rm(list = ls(all.names = TRUE))

# 1. Input file path
input_file_path <- "../csv/WB_combined_bkt_rev2.csv"

# 2. Number of clusters (k)
cluster_n <- 4

# 3. Set output directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  message("Output directory set to script location: ", getwd())
} else {
  message("Output directory is current working directory: ", getwd())
}

# ===== 1. Load libraries ======================================================
library(tidyverse)
library(lubridate)

# ===== 2. Data loading and preprocessing ======================================
data <- read_csv(input_file_path)
data$date <- as.Date(data$date)
cluster_col <- paste0("cluster_k", cluster_n)

# Extract year and exclude NAs based on the target cluster column
data <- data %>%
  mutate(year = year(date)) %>%
  filter(!is.na(!!sym(cluster_col)))

# ===== 3. Expand data (residence period to seasonal units) ====================
stay_data <- data %>%
  group_by(tag, river) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    .groups = "drop"
  )

expand_season <- stay_data %>%
  rowwise() %>%
  mutate(
    season_seq = list({
      dseq <- seq(start_date, end_date, by = "day")
      unique(paste0(year(dseq), "-S", quarter(dseq)))
    })
  ) %>%
  unnest(season_seq) %>%
  ungroup() %>%
  left_join(data %>% select(tag, !!sym(cluster_col)) %>% distinct(), by = "tag") %>%
  select(tag, river, !!sym(cluster_col), season_seq) %>%
  distinct()

alive_by_season <- expand_season %>%
  group_by(river, !!sym(cluster_col), season_seq) %>%
  summarise(n_alive = n_distinct(tag), .groups = "drop")

base_data_s3_all <- alive_by_season %>%
  mutate(year = as.numeric(substr(season_seq, 1, 4))) %>%
  filter(grepl("-S3$", season_seq))

# ===== 4. Portfolio effect calculation function (Delta CV) ====================
calculate_pe_for_window <- function(window_data, label) {
  
  cluster_col_name <- names(window_data)[grepl("^cluster_k", names(window_data))][1]
  
  if (n_distinct(window_data$season_seq) < 2) {
    return(data.frame(
      period = label,
      pe_cluster_effect = NA_real_,
      pe_river_effect = NA_real_,
      pe_overall_effect = NA_real_
    ))
  }
  
  # 1. River x Cluster CV
  cv_by_river_cluster <- window_data %>%
    group_by(river, !!sym(cluster_col_name)) %>%
    summarise(mean_alive = mean(n_alive),
              cv = sd(n_alive) / mean(n_alive), .groups = 'drop')
  
  # 2. River CV
  cv_by_river <- window_data %>%
    group_by(river, season_seq) %>%
    summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
    group_by(river) %>%
    summarise(mean_alive = mean(n_alive_sum),
              cv = sd(n_alive_sum) / mean(n_alive_sum), .groups = 'drop')
  
  # 3. Cluster CV
  cv_by_cluster <- window_data %>%
    group_by(!!sym(cluster_col_name), season_seq) %>%
    summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
    group_by(!!sym(cluster_col_name)) %>%
    summarise(mean_alive = mean(n_alive_sum),
              cv = sd(n_alive_sum) / mean(n_alive_sum), .groups = 'drop')
  
  # 4. Total CV
  cv_total <- window_data %>%
    group_by(season_seq) %>%
    summarise(total_n_alive = sum(n_alive), .groups = 'drop') %>%
    summarise(cv = sd(total_n_alive) / mean(total_n_alive)) %>%
    pull(cv)
  
  # Weighted mean
  cv_w_mean_river_cluster <- weighted.mean(cv_by_river_cluster$cv, w = cv_by_river_cluster$mean_alive, na.rm = TRUE)
  cv_w_mean_river <- weighted.mean(cv_by_river$cv, w = cv_by_river$mean_alive, na.rm = TRUE)
  cv_w_mean_cluster <- weighted.mean(cv_by_cluster$cv, w = cv_by_cluster$mean_alive, na.rm = TRUE)
  
  # PE (Delta CV)
  pe_cluster_effect <- cv_w_mean_river_cluster - cv_w_mean_river
  pe_river_effect <- cv_w_mean_river_cluster - cv_w_mean_cluster
  pe_overall_effect <- cv_w_mean_river_cluster - cv_total
  
  return(data.frame(
    period = label,
    pe_cluster_effect = pe_cluster_effect,
    pe_river_effect = pe_river_effect,
    pe_overall_effect = pe_overall_effect
  ))
}

# ===== 5. Execute PE analysis for 2 periods ===================================
data_p1 <- base_data_s3_all %>% filter(year >= 2002, year <= 2008)
pe_p1 <- calculate_pe_for_window(data_p1, "2002–2008")

data_p2 <- base_data_s3_all %>% filter(year >= 2009, year <= 2015)
pe_p2 <- calculate_pe_for_window(data_p2, "2009–2015")

pe_summary <- bind_rows(pe_p1, pe_p2)

# ===== 6. Save results (CSV & PDF) ============================================

# --- 6a. CSV Output ---
csv_filename <- paste0("CV_bkt", cluster_n, ".csv")
write_csv(pe_summary, csv_filename)
message("✅ CSV saved successfully: ", csv_filename)

# --- 6b. Create plot (Cluster Effect only, black, large points, fixed Y-axis) ---
plot_data <- pe_summary %>%
  select(period, pe_cluster_effect) %>%
  pivot_longer(cols = pe_cluster_effect, names_to = "effect_type", values_to = "value")

p <- ggplot(plot_data, aes(x = period, y = value, group = 1)) +
  geom_line(linewidth = 1.2, color = "black") + 
  geom_point(size = 5, color = "black") +
  # Specify Y-axis range (0 to 0.25)
  scale_y_continuous(limits = c(0, 0.25)) +
  theme_bw(base_size = 18) + 
  labs(
    x = "Period",
    y = "Cluster Portfolio Effect"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(color = "black", size = 20, face = "bold"),
    axis.text = element_text(color = "black", size = 18, face = "bold"),
    plot.title = element_text(color = "black", size = 22, face = "bold", hjust = 0.5)
  )

# --- 6c. Save PDF ---
pdf_filename <- paste0("CV_bkt", cluster_n, "_2periods_cluster_only.pdf")
ggsave(
  filename = pdf_filename,
  plot = p,
  device = "pdf",
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)
message("✅ PDF saved successfully: ", pdf_filename)

