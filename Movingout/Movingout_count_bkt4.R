# ===== 0. Load Packages & Environment Reset =====
rm(list = ls(all.names = TRUE))
library(dplyr)
library(lubridate)
library(tidyr)

## ===== 1. Automatically Get Script Directory =====
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  current_dir <- dirname(script_path)
} else {
  current_dir <- getwd() 
  message("Warning: rstudioapi is not available. Using working directory (getwd).")
}

# ===== 2. Set Output Directory =====
output_dir <- current_dir

# File Path
file_path <- "../csv/WB_combined_bkt_rev2.csv"

# ===== 3. Analysis Configuration =====
target_cluster_col <- "cluster_k4"

# ===== 4. Load Data & Preprocessing =====
df <- read.csv(file_path)
df$date <- as.Date(df$date)

# Standardize river names to short abbreviations (if river column exists)
if ("river" %in% colnames(df)) {
  df <- df %>%
    mutate(river = case_when(
      tolower(river) == "west brook"  ~ "WB",
      tolower(river) == "wb jimmy"     ~ "OL",
      tolower(river) == "wb mitchell"  ~ "OS",
      tolower(river) == "wb obear"     ~ "IL",
      TRUE ~ river
    ))
}

# Determine representative cluster (Fix one cluster per individual fish)
df$cluster_id_tmp <- df[[target_cluster_col]]
cluster_def <- df %>%
  filter(!is.na(cluster_id_tmp)) %>%
  group_by(tag) %>%
  summarise(cluster_fixed = names(which.max(table(cluster_id_tmp)))) %>%
  ungroup()

# Combine representative cluster values
df <- df %>%
  left_join(cluster_def, by = "tag")

# ==============================================================================
# 5. Identify Emigrating Individuals (Create Tag List)
# ==============================================================================

# A: Individuals that used the upstream section (5400-6000) at least once
tags_in_A <- df %>%
  filter(pittag == 1, riverMeter >= 5400, riverMeter <= 6000) %>%
  distinct(tag) %>%
  pull(tag)

# B: Individuals that used the downstream section (4150-4393) at least once
tags_in_B <- df %>%
  filter(pittag == 1, riverMeter >= 4150, riverMeter <= 4393) %>%
  distinct(tag) %>%
  pull(tag)

# --- Emigrant list to upstream (above) ---
list_above <- df %>%
  filter(tag %in% tags_in_A) %>%
  group_by(tag) %>%
  filter(date == max(date, na.rm = TRUE)) %>% # Extract final observation date
  ungroup() %>%
  filter(area == "above") %>%
  select(tag, cluster_fixed, season, date, area, riverMeter) %>%
  mutate(move_direction = "to_above")

# --- Emigrant list to downstream (below) ---
list_below <- df %>%
  filter(tag %in% tags_in_B) %>%
  group_by(tag) %>%
  filter(date == max(date, na.rm = TRUE)) %>% # Extract final observation date
  ungroup() %>%
  filter(area == "below") %>%
  select(tag, cluster_fixed, season, date, area, riverMeter) %>%
  mutate(move_direction = "to_below")

# --- Combine lists ---
moveout_tag_list <- bind_rows(list_above, list_below) %>%
  rename(!!target_cluster_col := cluster_fixed) %>%
  arrange(season, !!sym(target_cluster_col), tag)

# ==============================================================================
# 6. Export Results
# ==============================================================================

cat("\n=== Emigrant Tag List (Preview) ===\n")
print(head(moveout_tag_list))

# Export to CSV
filename <- paste0("moveout_tag_list_", target_cluster_col, "_bkt.csv")
file_taglist <- file.path(output_dir, filename)

write.csv(moveout_tag_list, file = file_taglist, row.names = FALSE)

cat("\n✅ CSV Export Completed:\n")
cat("  - ", file_taglist, "\n")
cat("  - Total Emigrants: ", nrow(moveout_tag_list), " tags\n")

