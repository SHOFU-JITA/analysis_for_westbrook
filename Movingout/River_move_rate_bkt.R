# ==============================================================================
# ===== 0. Environment Reset and Setup =====
# ==============================================================================
rm(list = ls(all.names = TRUE)) # Clear workspace

# --- Libraries ---
library(tidyverse)
library(lubridate)

# --- Constants and File Paths ---
INPUT_CSV_PATH <- "../csv/WB_combined_bkt_rev2.csv"

# --- Output Files ---
# 1. Movement counts summary (Annual)
OUTPUT_CSV_OVERALL <- "WB_movements_overall_annual_bkt.csv"
OUTPUT_CSV_CLUSTER <- "WB_movements_by_cluster_annual_bkt.csv"

# 2. Population stats by river (River-based)
OUTPUT_CSV_POP_OVERALL <- "WB_population_stats_by_river.csv_bkt"
OUTPUT_CSV_POP_CLUSTER <- "WB_population_stats_by_river_and_cluster.csv_bkt"

# 3. Unique individual movement type summary
OUTPUT_CSV_UNIQUE_SUMMARY <- "WB_population_unique_summary.csv_bkt"

# Rate output files
OUTPUT_CSV_RATES_OVERALL <- "WB_rates_overall.csv_bkt"
OUTPUT_CSV_RATES_CLUSTER <- "WB_rates_by_cluster.csv_bkt"

CLUSTER_COLUMN <- "cluster_k4"

# --- River Network Settings ---
MAINSTEM <- "WB"
TRIBUTARIES <- c("OS", "OL", "IL")

# ==============================================================================
# ===== 1. Data Loading and Preprocessing =====
# ==============================================================================
cat("--- 1. Data Loading and Preprocessing ---\n")

if (!file.exists(INPUT_CSV_PATH)) {
  stop(paste("Error: Input file not found:", INPUT_CSV_PATH))
}
df <- read_csv(INPUT_CSV_PATH, show_col_types = FALSE)

if (!CLUSTER_COLUMN %in% names(df)) {
  stop(paste("Error: Specified cluster column '", CLUSTER_COLUMN, "' does not exist.", sep=""))
}

# Data cleaning and renaming river codes
df <- df %>%
  mutate(
    river = case_when(
      river %in% c("west brook", "westbrook")   ~ "WB",
      river %in% c("wb jimmy", "wb_jimmy")       ~ "OL",
      river %in% c("wb mitchell", "wb_mitchell") ~ "OS",
      river %in% c("wb obear", "wb_obear")       ~ "IL",
      TRUE ~ river
    ),
    date = as.Date(date)
  ) %>%
  filter(!is.na(.data[[CLUSTER_COLUMN]]))

cat("Data loading complete.\n")

# ==============================================================================
# ===== 2. Create Movement Event Data =====
# ==============================================================================
cat("--- 2. Creating Base Movement Dataframe ---\n")

movement_events <- df %>%
  arrange(tag, date) %>%
  group_by(tag) %>%
  mutate(
    from_river = lag(river),
    to_river = river
  ) %>%
  ungroup() %>%
  filter(!is.na(from_river))

cat("Movement events extracted.\n")

# ==============================================================================
# ===== 3. Logic for Unifying Movements (Annual) =====
# ==============================================================================
process_movement_data_annual <- function(movement_data) {
  actual_movements <- movement_data %>%
    filter(from_river != to_river)
  
  tributary_to_tributary <- actual_movements %>%
    filter(from_river %in% TRIBUTARIES & to_river %in% TRIBUTARIES)
  
  all_movement_segments <- list()
  
  if(nrow(tributary_to_tributary) > 0) {
    all_movement_segments[[1]] <- tributary_to_tributary %>%
      mutate(to_river = MAINSTEM)
    all_movement_segments[[2]] <- tributary_to_tributary %>%
      mutate(from_river = MAINSTEM)
    all_movement_segments[[3]] <- actual_movements %>%
      filter(!(from_river %in% TRIBUTARIES & to_river %in% TRIBUTARIES))
    
    unified_movements <- bind_rows(all_movement_segments) %>%
      count(from_river, to_river, name = "weight")
  } else {
    unified_movements <- actual_movements %>%
      count(from_river, to_river, name = "weight")
  }
  return(unified_movements)
}

# ==============================================================================
# ===== 4. Aggregation and Export: Movement Overall =====
# ==============================================================================
cat("--- 3. Aggregating Movement Data (Overall) ---\n")
overall_aggregated <- process_movement_data_annual(movement_events) %>%
  mutate(category = "Overall") %>%
  select(category, from_river, to_river, count = weight)
write_csv(overall_aggregated, OUTPUT_CSV_OVERALL)

# ==============================================================================
# ===== 5. Aggregation and Export: Movement By Cluster =====
# ==============================================================================
cat("--- 4. Aggregating Movement Data (By Cluster) ---\n")
unique_clusters <- unique(df[[CLUSTER_COLUMN]]) %>% sort()
cluster_results_list <- list()

for (current_cluster in unique_clusters) {
  tags_in_cluster <- df %>% 
    filter(.data[[CLUSTER_COLUMN]] == current_cluster) %>% 
    pull(tag) %>% unique()
  cluster_sub_data <- movement_events %>% filter(tag %in% tags_in_cluster)
  processed_data <- process_movement_data_annual(cluster_sub_data)
  if(nrow(processed_data) > 0) {
    processed_data <- processed_data %>% mutate(!!CLUSTER_COLUMN := current_cluster)
    cluster_results_list[[as.character(current_cluster)]] <- processed_data
  }
}

if (length(cluster_results_list) > 0) {
  final_cluster_df <- bind_rows(cluster_results_list) %>%
    select(!!sym(CLUSTER_COLUMN), from_river, to_river, count = weight)
  write_csv(final_cluster_df, OUTPUT_CSV_CLUSTER)
}

# ==============================================================================
# ===== 6. Population Analysis Setup =====
# ==============================================================================
cat("\n--- 5. Analyzing Population Stats (Residents vs Mobile) ---\n")

# 1. Classify movement type for each individual
tag_classification <- df %>%
  group_by(tag) %>%
  summarise(
    n_rivers_visited = n_distinct(river),
    .groups = "drop"
  ) %>%
  mutate(
    move_type = ifelse(n_rivers_visited == 1, "Resident", "Mobile")
  )

# 2. Join movement type to the main dataframe
df_with_type <- df %>%
  left_join(tag_classification, by = "tag")

# ==============================================================================
# ===== 7. Export Population Stats: By River (Overall) =====
# ==============================================================================
cat("--- 6. Exporting River-based Stats (Overall) ---\n")
pop_stats_overall <- df_with_type %>%
  group_by(river) %>%
  summarise(
    total_observed = n_distinct(tag),
    resident_count = n_distinct(tag[move_type == "Resident"]),
    mobile_count   = n_distinct(tag[move_type == "Mobile"]),
    .groups = "drop"
  ) %>%
  mutate(category = "Overall") %>%
  select(category, river, total_observed, resident_count, mobile_count)
write_csv(pop_stats_overall, OUTPUT_CSV_POP_OVERALL)

# ==============================================================================
# ===== 8. Export Population Stats: By River (By Cluster) =====
# ==============================================================================
cat("--- 7. Exporting River-based Stats (By Cluster) ---\n")
pop_stats_cluster <- df_with_type %>%
  group_by(!!sym(CLUSTER_COLUMN), river) %>%
  summarise(
    total_observed = n_distinct(tag),
    resident_count = n_distinct(tag[move_type == "Resident"]),
    mobile_count   = n_distinct(tag[move_type == "Mobile"]),
    .groups = "drop"
  ) %>%
  arrange(!!sym(CLUSTER_COLUMN), river)
write_csv(pop_stats_cluster, OUTPUT_CSV_POP_CLUSTER)


# ==============================================================================
# ===== 9. Total Unique Population Summary =====
# ==============================================================================
cat("\n--- 8. Exporting Total Unique Population Summary ---\n")

# Combine cluster info with tag classification data (Assumes 1 tag belongs to 1 cluster)
tag_cluster_map <- df %>% distinct(tag, !!sym(CLUSTER_COLUMN))
tag_classification_full <- tag_classification %>%
  left_join(tag_cluster_map, by = "tag")

# 9.1 Overall summary
summary_overall <- tag_classification_full %>%
  count(move_type) %>%
  pivot_wider(names_from = move_type, values_from = n, values_fill = 0) %>%
  mutate(
    Group = "Overall",
    Total_Tags = Resident + Mobile,
    Mobile_Rate = Mobile / Total_Tags
  ) %>%
  select(Group, Total_Tags, Resident, Mobile, Mobile_Rate)

# 9.2 Summary by cluster
summary_cluster <- tag_classification_full %>%
  group_by(!!sym(CLUSTER_COLUMN), move_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = move_type, values_from = n, values_fill = 0) %>%
  mutate(
    Total_Tags = Resident + Mobile,
    Mobile_Rate = Mobile / Total_Tags
  ) %>%
  rename(Group = !!sym(CLUSTER_COLUMN)) %>% 
  mutate(Group = as.character(Group)) %>%
  select(Group, Total_Tags, Resident, Mobile, Mobile_Rate) %>%
  arrange(Group)

# 9.3 Combine and export
final_unique_summary <- bind_rows(summary_overall, summary_cluster)

write_csv(final_unique_summary, OUTPUT_CSV_UNIQUE_SUMMARY)
cat("Exported:", OUTPUT_CSV_UNIQUE_SUMMARY, "\n")
print(final_unique_summary)

cat("\nAll processing complete.\n")

# ==============================================================================
# ===== 10. Calculate Emigration and Immigration Rates =====
# ==============================================================================
cat("\n--- 9. Calculating Emigration and Immigration Rates ---\n")

# --- Function: Rate calculation logic ---
calculate_rates <- function(movement_df, pop_df, group_cols = NULL) {
  
  # 1. Aggregate Emigration Count
  group_keys_out <- c(group_cols, "from_river")
  outflow <- movement_df %>%
    group_by(across(all_of(group_keys_out))) %>%
    summarise(emigration_count = sum(count), .groups = "drop") %>%
    rename(river = from_river)
  
  # 2. Aggregate Immigration Count
  group_keys_in <- c(group_cols, "to_river")
  inflow <- movement_df %>%
    group_by(across(all_of(group_keys_in))) %>%
    summarise(immigration_count = sum(count), .groups = "drop") %>%
    rename(river = to_river)
  
  # 3. Join with population data and calculate rates
  join_keys <- c("river")
  if (!is.null(group_cols)) {
    join_keys <- c(group_cols, "river")
  }
  
  rates_df <- pop_df %>%
    select(all_of(c(join_keys, "total_observed"))) %>%
    left_join(outflow, by = join_keys) %>%
    left_join(inflow, by = join_keys) %>%
    mutate(
      emigration_count = replace_na(emigration_count, 0),
      immigration_count = replace_na(immigration_count, 0)
    ) %>%
    mutate(
      emigration_rate = emigration_count / total_observed,
      immigration_rate = immigration_count / total_observed
    )
  
  return(rates_df)
}

# ----- 10.1 Overall Rates -----
if (exists("overall_aggregated") && exists("pop_stats_overall")) {
  
  rates_overall <- calculate_rates(overall_aggregated, pop_stats_overall, group_cols = NULL) %>%
    mutate(category = "Overall") %>%
    select(category, river, total_observed, 
           emigration_count, emigration_rate, 
           immigration_count, immigration_rate)
  
  write_csv(rates_overall, OUTPUT_CSV_RATES_OVERALL)
  cat("Exported:", OUTPUT_CSV_RATES_OVERALL, "\n")
  print(rates_overall)
}

# ----- 10.2 Rates by Cluster -----
if (exists("final_cluster_df") && exists("pop_stats_cluster")) {
  
  rates_cluster <- calculate_rates(final_cluster_df, pop_stats_cluster, group_cols = CLUSTER_COLUMN) %>%
    select(!!sym(CLUSTER_COLUMN), river, total_observed, 
           emigration_count, emigration_rate, 
           immigration_count, immigration_rate) %>%
    arrange(!!sym(CLUSTER_COLUMN), river)
  
  write_csv(rates_cluster, OUTPUT_CSV_RATES_CLUSTER)
  cat("Exported:", OUTPUT_CSV_RATES_CLUSTER, "\n")
}

cat("\nRate calculation complete.\n")

