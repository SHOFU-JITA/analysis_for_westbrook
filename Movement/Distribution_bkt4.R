# ==============================================================================
# Fish Distribution and Behavior Analysis by Cluster
# ==============================================================================

# ===== 0. Environment Reset =====
rm(list = ls(all.names = TRUE))

# --- Load Packages ---
library(readr)
library(dplyr)
library(tidyr) 

# --- Constants and File Paths ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}

INPUT_CSV_PATH <- "../csv/WB_combined_bkt_rev2.csv"
OUTPUT_CSV_PATH <- file.path(SCRIPT_DIR, "WB_Distribution_bkt4.csv") 
CLUSTER_COLUMN <- "cluster_k4"

# ===== 1. Load Data and Preprocess =====
cat("--- 1. Load Data ---\n")
if (!file.exists(INPUT_CSV_PATH)) {
  stop(paste("Error: Input file not found:", INPUT_CSV_PATH))
}

df <- read_csv(INPUT_CSV_PATH, show_col_types = FALSE)

# Standardize river names to short abbreviations
df <- df %>%
  mutate(
    river = case_when(
      tolower(river) == "west brook"   ~ "WB",
      tolower(river) == "wb jimmy"     ~ "OL",
      tolower(river) == "wb mitchell"  ~ "OS",
      tolower(river) == "wb obear"     ~ "IL",
      TRUE ~ river
    )
  )

# Verify cluster column existence
if (!CLUSTER_COLUMN %in% names(df)) {
  stop(paste("Error: Column '", CLUSTER_COLUMN, "' does not exist.", sep=""))
}

# ===== 2. Determine Behavior Types and Aggregate =====
cat("--- 2. Determine Behavior Types and Aggregate ---\n")

# Step 1: Determine behavior type per individual fish
fish_behavior <- df %>%
  group_by(tag, !!sym(CLUSTER_COLUMN)) %>%
  summarize(
    n_rivers = n_distinct(river),
    # Categorize as "Multi-River" if detected in 2+ rivers; otherwise, use the river name
    behavior_type = if_else(n_rivers > 1, "Multi-River", first(river)),
    .groups = "drop"
  )

# Step 2: Count per cluster (Pivot to wide format)
cluster_counts <- fish_behavior %>%
  count(!!sym(CLUSTER_COLUMN), behavior_type) %>%
  pivot_wider(
    names_from = behavior_type, 
    values_from = n, 
    values_fill = 0
  )

# Step 3: Create and combine Total row
df_main <- cluster_counts %>%
  mutate(!!sym(CLUSTER_COLUMN) := as.character(!!sym(CLUSTER_COLUMN)))

# Calculate Total row by summing all numeric columns
df_total <- df_main %>%
  summarise(across(where(is.numeric), sum)) %>%
  mutate(!!sym(CLUSTER_COLUMN) := "total")

# Combine main data and Total row
df_combined <- bind_rows(df_main, df_total)

# Step 4: Calculate total individuals and proportions
final_stats <- df_combined %>%
  
  mutate(Total_Ind = rowSums(across(where(is.numeric)))) %>%
  
  mutate(
    across(
      .cols = where(is.numeric) & !any_of("Total_Ind"), 
      .fns = ~ .x / Total_Ind,
      .names = "{.col}_Prop"
    )
  ) %>%
  
  mutate(
    across(
      .cols = where(is.numeric) & !ends_with("_Prop"), 
      .fns = ~ .x / last(.x), 
      .names = "{.col}_ClusterProp"
    )
  )

# Step 5: Reorder columns
# Desired order: Cluster, Total_Ind metrics, Multi-River metrics, followed by alphabetical river metrics
final_stats <- final_stats %>%
  select(
    !!sym(CLUSTER_COLUMN),
    starts_with("Total_Ind"),
    starts_with("Multi-River"),
    sort(names(select(., -c(!!sym(CLUSTER_COLUMN), starts_with("Total_Ind"), starts_with("Multi-River")))))
  )

# ===== 3. Export to CSV =====
cat("--- 3. Export CSV ---\n")
write_csv(final_stats, OUTPUT_CSV_PATH)
cat(paste("Aggregation completed. File saved to:", OUTPUT_CSV_PATH, "\n"))

# Print results to console
print(final_stats)

