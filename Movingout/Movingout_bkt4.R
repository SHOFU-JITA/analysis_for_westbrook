# ===== 0. Setup & Load Packages =====
rm(list = ls(all.names = TRUE))

library(dplyr)
library(lubridate)
library(tidyr)

# Output directory 
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  output_dir <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/")
} else {
  output_dir <- "./"
}

# Input file path
file_path <- "../csv/WB_combined_bkt_rev2.csv"

# ===== 1. Analysis Settings =====
# Target cluster column name
target_cluster_col <- "cluster_k4"

# ===== 2. Data Loading & Preprocessing =====
df <- read.csv(file_path)
df$date <- as.Date(df$date)

# Define the representative cluster for each individual (mode)
df$cluster_id <- df[[target_cluster_col]]

cluster_def <- df %>%
  filter(!is.na(cluster_id)) %>%
  group_by(tag) %>%
  summarise(cluster_id = names(which.max(table(cluster_id)))) %>%
  ungroup()

# Merge representative cluster_id back to the main dataframe
df <- df %>%
  select(-cluster_id) %>% 
  left_join(cluster_def, by = "tag")

# ==============================================================================
# Common Setup: Identify individuals using specific sections
# ==============================================================================

# Tags in the upper section (5400-6000)
tags_in_A <- df %>%
  filter(pittag == 1, riverMeter >= 5400, riverMeter <= 6000) %>%
  distinct(tag) %>%
  pull(tag)

# Tags in the lower section (4150-4393)
tags_in_B <- df %>%
  filter(pittag == 1, riverMeter >= 4150, riverMeter <= 4393) %>%
  distinct(tag) %>%
  pull(tag)

# Base seasons
unique_seasons <- sort(unique(df$season))

# ==============================================================================
# Part 1: Aggregation by Cluster
# ==============================================================================

# Denominator: Total individuals per cluster
denom_cluster <- cluster_def %>%
  group_by(cluster_id) %>%
  summarise(n_total = n(), .groups = "drop")

# Create a base grid
base_grid_cluster <- denom_cluster %>%
  crossing(season = unique_seasons)

# --- Above ---
movers_A_cluster <- df %>%
  filter(tag %in% tags_in_A) %>%
  group_by(tag) %>%
  filter(date == max(date, na.rm = TRUE)) %>% # Last observation
  ungroup() %>%
  filter(area == "above") %>%
  group_by(cluster_id, season) %>%
  summarise(n_above = n_distinct(tag), .groups = "drop")

# --- Below ---
movers_B_cluster <- df %>%
  filter(tag %in% tags_in_B) %>%
  group_by(tag) %>%
  filter(date == max(date, na.rm = TRUE)) %>% # Last observation
  ungroup() %>%
  filter(area == "below") %>%
  group_by(cluster_id, season) %>%
  summarise(n_below = n_distinct(tag), .groups = "drop")

# --- Merge & Calculate Totals ---
result_cluster <- base_grid_cluster %>%
  left_join(movers_A_cluster, by = c("cluster_id", "season")) %>%
  left_join(movers_B_cluster, by = c("cluster_id", "season")) %>%
  mutate(
    n_above = replace_na(n_above, 0),
    n_below = replace_na(n_below, 0),
    n_move_sum = n_above + n_below,
    ratio_above = n_above / n_total,
    ratio_below = n_below / n_total,
    ratio_move_sum = n_move_sum / n_total
  ) %>%
  rename(!!target_cluster_col := cluster_id) %>%
  arrange(!!sym(target_cluster_col), season)

# ==============================================================================
# Part 2: Overall Aggregation (Ignoring Clusters)
# ==============================================================================

# Denominator: Total unique individuals
n_grand_total <- n_distinct(df$tag)

# Base: Seasons only
base_grid_overall <- data.frame(season = unique_seasons) %>%
  mutate(n_total = n_grand_total)

# --- Above ---
movers_A_overall <- df %>%
  filter(tag %in% tags_in_A) %>%
  group_by(tag) %>%
  filter(date == max(date, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(area == "above") %>%
  group_by(season) %>%
  summarise(n_above = n_distinct(tag), .groups = "drop")

# --- Below ---
movers_B_overall <- df %>%
  filter(tag %in% tags_in_B) %>%
  group_by(tag) %>%
  filter(date == max(date, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(area == "below") %>%
  group_by(season) %>%
  summarise(n_below = n_distinct(tag), .groups = "drop")

# --- Merge & Calculate Totals ---
result_overall <- base_grid_overall %>%
  left_join(movers_A_overall, by = "season") %>%
  left_join(movers_B_overall, by = "season") %>%
  mutate(
    n_above = replace_na(n_above, 0),
    n_below = replace_na(n_below, 0),
    n_move_sum = n_above + n_below,
    ratio_above = n_above / n_total,
    ratio_below = n_below / n_total,
    ratio_move_sum = n_move_sum / n_total
  ) %>%
  arrange(season)

# ==============================================================================
# Part 3: Output Results
# ==============================================================================

cat("\n=== 1. Move-out ratios by cluster and season ===\n")
cat("Target cluster column:", target_cluster_col, "\n")
print(result_cluster)

cat("\n=== 2. Overall move-out ratios by season ===\n")
print(result_overall)

# Export to CSV
file_cluster <- paste0(output_dir, "moveout_by_", target_cluster_col, "_bkt.csv")
file_overall <- paste0(output_dir, "moveout_overall_bkt.csv")

write.csv(result_cluster, file = file_cluster, row.names = FALSE)
write.csv(result_overall, file = file_overall, row.names = FALSE)

cat("\n✅ CSV output completed:\n")
cat("  1. ", file_cluster, "\n")
cat("  2. ", file_overall, "\n")

