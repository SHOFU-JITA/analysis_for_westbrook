# ===== 0. Reset environment =====
rm(list = ls(all.names = TRUE))

# ===== 1. Load packages =====
library(tidyverse)
library(lubridate)
library(ggplot2)
library(chron)
library(readr)
library(dplyr)
library(stringr)
  
# ===== 2. Read CSV files =====
path_df1 <- "../csv/WB_stationaryantenna.csv"
path_df2 <- "../Clustering/WB_bkt_clustered_rev2.csv"
output_dir <- "../csv/"

tryCatch({
  df1_raw <- read.csv(path_df1, stringsAsFactors = FALSE)
  df2_raw <- read.csv(path_df2, stringsAsFactors = FALSE)
  cat("Successfully loaded files.\n")
}, error = function(e) {
  stop("Failed to load CSV files. Please check the file paths.\nError: ", e$message)
})

# ===== 3. Define common column names =====
common_cols <- c(
  "species", "river", "tag", "area", "date", "year", "season", "riverMeter",
  "cohort", "sampleName", "sampleNumber",
  "cluster_k1","cluster_k2","cluster_k3","cluster_k4","cluster_k5",
  "cluster_k6","cluster_k7","cluster_k8","cluster_k9","cluster_k10",
  "observedLength", "observedWeight",
  "pittag"
)

# ===== 4. Preprocess data frames =====

# --- Process df1 (Pittag / Antenna) ---
df1 <- df1_raw %>%
  filter(species == "bkt") %>%
  mutate(
    date = as.Date(parse_date_time(date, orders = c("mdy", "ymd", "dmy", "mdy HMS", "ymd HMS"))),
    tag = if("tag" %in% names(.)) str_trim(tolower(tag)) else NA_character_,
    cohort = NA_character_,
    sampleName = NA_integer_,
    sampleNumber = NA_integer_,
    observedLength = NA_real_,
    observedWeight = NA_real_,
    cluster_k1 = NA_real_, cluster_k2 = NA_real_, cluster_k3 = NA_real_, cluster_k4 = NA_real_, cluster_k5 = NA_real_,
    cluster_k6 = NA_real_, cluster_k7 = NA_real_, cluster_k8 = NA_real_, cluster_k9 = NA_real_, cluster_k10 = NA_real_,
    pittag = 1,
    area = NA_character_,
    season = NA_real_
  )

# --- Process df2 (Clustered / Capture) ---
df2_intermediate <- df2_raw %>%
  filter(species == "bkt") %>%
  mutate(
    date = as.Date(parse_date_time(date, orders = c("mdy", "ymd", "dmy", "mdy HMS", "ymd HMS"))),
    tag = if("tag" %in% names(.)) str_trim(tolower(tag)) else NA_character_,
    pittag = 0
  )

# Fill missing columns in df2 with NA
missing_cols_in_df2 <- setdiff(common_cols, names(df2_intermediate))
for (col_name in missing_cols_in_df2) {
  if (col_name != "date") {
    df2_intermediate[[col_name]] <- NA
  }
}
df2 <- df2_intermediate

# --- Select columns and convert types ---
df1 <- df1 %>% select(all_of(common_cols)) %>% mutate(cohort = as.character(cohort))
df2 <- df2 %>% select(all_of(common_cols)) %>% mutate(cohort = as.character(cohort))

# ===== 5. Combine and format data =====
df_combined <- bind_rows(df1, df2)

# Fill missing information per tag
df_filled <- df_combined %>%
  arrange(tag, pittag) %>%
  group_by(tag) %>%
  fill(
    cohort,
    cluster_k1, cluster_k2, cluster_k3, cluster_k4, cluster_k5, cluster_k6,
    cluster_k7, cluster_k8, cluster_k9, cluster_k10,
    .direction = "downup"
  ) %>%
  ungroup() %>%
  arrange(tag, date)

# ===== 6. Filter out individuals without cluster information =====
initial_count <- nrow(df_filled)

df_filled <- df_filled %>%
  filter(!is.na(cluster_k4))

removed_count <- initial_count - nrow(df_filled)
cat(paste0("\n[Info] Rows removed due to missing cluster info: ", removed_count, "\n"))

# ===== 7. Define Season and Area =====
df_filled <- df_filled %>%
  mutate(
    season = case_when(
      (month(date) == 2 & day(date) >= 6) | (month(date) > 2 & month(date) < 5) | (month(date) == 5 & day(date) <= 5) ~ 1,
      (month(date) == 5 & day(date) >= 6) | (month(date) > 5 & month(date) < 8) | (month(date) == 8 & day(date) <= 5) ~ 2,
      (month(date) == 8 & day(date) >= 6) | (month(date) > 8 & month(date) < 11) | (month(date) == 11 & day(date) <= 5) ~ 3,
      TRUE ~ 4
    ),
    area = case_when(
      !is.na(area) ~ area,
      river %in% c("wb mitchell", "wb jimmy", "wb obear") ~ "trib",
      river == "west brook" & riverMeter <= 4992 ~ "below",
      river == "west brook" & riverMeter >= 5354 ~ "above",
      river == "west brook" & riverMeter > 4992 & riverMeter < 5354 ~ "inside",
      TRUE ~ NA_character_
    )
  )

# ===== 8. Verify and save results =====
cat("\n--- Final Data Summary ---\n")
print(dim(df_filled))
cat("\n--- NA Check for cluster_k4 (Should be 0) ---\n")
print(sum(is.na(df_filled$cluster_k4)))

output_filename <- "WB_combined_bkt_rev2.csv"
output_path <- file.path(output_dir, output_filename)

write.csv(df_filled, output_path, row.names = FALSE)
cat(paste("\nSuccessfully saved to:", output_path, "\n"))

