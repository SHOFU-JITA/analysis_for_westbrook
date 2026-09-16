# ==============================================================================
# Population Estimation & Synchrony Analysis via POPAN Model (Jolly-Seber)
# Customized: S3 Only, Specific Years, Total Calculation, Simplified Synchrony
# Plus: Conditional Error Bars, Overlay Plot, Output to Script Directory
# Final Update: Export Formatted Synchrony Matrix to CSV
# ==============================================================================

# ===== 0. Environment Reset =====
rm(list = ls(all.names = TRUE))
graphics.off()

# ===== 1. Load Packages =====
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RMark)
library(lubridate)
library(stringr)
library(tibble)

# ===== Configuration: Modify according to your data =====
RIVER_COL_NAME <- "river" 
FILE_SUFFIX <- "ats" # For file identification
FILE_PATH <- "../rawdata/cdWB_all.csv" # Input file path
TARGET_SPECIES <- "ats"

# [Output Directory Setup] Output to the directory where this script is saved
# * Requires running in RStudio and the script must be saved
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  OUTPUT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  OUTPUT_DIR <- getwd() # Fallback: Current directory
}
message(paste("Output Directory:", OUTPUT_DIR))

# [Important] Specify the correct name for the main river
MAIN_RIVER_NAME <- "WB" 

# ==============================================================================
# ===== 2. Load and Preprocess Data & Impute Missing Values =====
# ==============================================================================
if (!file.exists(FILE_PATH)) stop("Specified file not found.")
df_raw <- read_csv(FILE_PATH)

# --- 2-1. Filter by specific species ---
if ("species" %in% colnames(df_raw)) {
  message(paste("Filtering data for species:", TARGET_SPECIES))
  df_raw <- df_raw %>% filter(species == TARGET_SPECIES)
} else {
  warning("Warning: 'species' column not found. Skipping filtering.")
}

# Convert to Date type
df_raw$date <- as.Date(df_raw$date)

# --- 2-2. Define and Impute Season and Area ---
message("Imputing missing season and area values...")

df_raw <- df_raw %>%
  mutate(
    # Impute season from date only if it is NA
    season = case_when(
      !is.na(season) ~ as.numeric(season), # Keep existing values
      (month(date) == 2 & day(date) >= 6) | (month(date) > 2 & month(date) < 5) | (month(date) == 5 & day(date) <= 5) ~ 1,
      (month(date) == 5 & day(date) >= 6) | (month(date) > 5 & month(date) < 8) | (month(date) == 8 & day(date) <= 5) ~ 2,
      (month(date) == 8 & day(date) >= 6) | (month(date) > 8 & month(date) < 11) | (month(date) == 11 & day(date) <= 5) ~ 3,
      TRUE ~ 4
    ),
    # Define area (if 'riverMeter' column exists)
    area = case_when(
      !is.na(area) ~ area,
      river %in% c("wb mitchell", "wb jimmy", "wb obear") ~ "trib",
      river == "west brook" & "riverMeter" %in% colnames(df_raw) & riverMeter <= 4992 ~ "below",
      river == "west brook" & "riverMeter" %in% colnames(df_raw) & riverMeter >= 5354 ~ "above",
      river == "west brook" & "riverMeter" %in% colnames(df_raw) & riverMeter > 4992 & riverMeter < 5354 ~ "inside",
      TRUE ~ NA_character_
    )
  )

# Check if season was imputed correctly (for debugging)
if (any(is.na(df_raw$season))) {
  warning("Warning: Some rows still have missing season values. Please check the 'date' column.")
}

# --- 2-3. Rename river names to short abbreviations ---
df_raw <- df_raw %>%
  mutate(!!sym(RIVER_COL_NAME) := case_when(
    tolower(!!sym(RIVER_COL_NAME)) == "west brook"  ~ "WB",
    tolower(!!sym(RIVER_COL_NAME)) == "wb jimmy"     ~ "OL",
    tolower(!!sym(RIVER_COL_NAME)) == "wb mitchell"  ~ "OS",
    tolower(!!sym(RIVER_COL_NAME)) == "wb obear"     ~ "IL",
    TRUE ~ !!sym(RIVER_COL_NAME)
  ))

message(paste("Data processed. Total rows:", nrow(df_raw)))
message("Preparing Capture History...")

# ===== 3. Create Seasonal Capture History =====
df <- df_raw %>%
  mutate(
    year = year(date),
    time_label = paste0(year, "-S", season) 
  ) %>%
  select(tag, !!sym(RIVER_COL_NAME), time_label) %>%
  distinct()

capture_history_wide <- df %>%
  mutate(capture = 1) %>%
  pivot_wider(names_from = time_label, values_from = capture, values_fill = list(capture = 0))

ch_cols <- names(capture_history_wide)[!names(capture_history_wide) %in% c("tag", RIVER_COL_NAME)]
ch_order <- data.frame(col_name = ch_cols) %>%
  mutate(
    y = as.numeric(str_extract(col_name, "^[0-9]{4}")),
    s = as.numeric(str_extract(col_name, "[0-9]+$"))
  ) %>%
  arrange(y, s)
sorted_cols <- ch_order$col_name

capture_history_sorted <- capture_history_wide %>%
  select(tag, !!sym(RIVER_COL_NAME), all_of(sorted_cols))

message("Capture History created.")

# ==============================================================================
# ===== 4. Define Function for POPAN Estimation (Updated) =====
# ==============================================================================

estimate_N_POPAN <- function(data, river_name, sorted_cols) {
  ch_data <- data
  ch_cols_local <- sorted_cols
  
  if (length(ch_cols_local) == 0) return(NULL)
  
  # Aggregate actual number of captures (n)
  n_captured_v <- colSums(ch_data[, ch_cols_local, drop = FALSE])
  
  # Create CH string for MARK
  ch_data$ch <- apply(ch_data[, ch_cols_local, drop = FALSE], 1, paste0, collapse = "")
  ch_data_mark <- ch_data %>% select(tag, ch) %>% filter(grepl("1", ch))
  
  if (nrow(ch_data_mark) < 10) return(NULL)
  
  processed <- tryCatch({
    process.data(ch_data_mark, model = "POPAN", begin.time = 1)
  }, error = function(e) return(NULL))
  
  if (is.null(processed)) return(NULL)
  
  design <- make.design.data(processed)
  
  # --- C. Fix recruitment probability (pent) to 0 from 2005 onwards ---
  year_labels <- as.numeric(str_extract(ch_cols_local, "^[0-9]{4}"))
  pent_fix_indices <- which(year_labels[-length(year_labels)] >= 2005)
  
  fixed_list <- list()
  if (length(pent_fix_indices) > 0) {
    fixed_list$pent <- list(index = pent_fix_indices, value = 0)
  }
  
  model_params <- list(
    Phi  = list(formula = ~time),
    p    = list(formula = ~time),
    pent = list(formula = ~1, fixed = fixed_list$pent),
    N    = list(formula = ~1)
  )
  
  model <- tryCatch({
    mark(processed, design, model.parameters = model_params,
         silent = TRUE, output = FALSE, delete = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(model)) return(NULL)
  
  derived <- tryCatch({ popan.derived(processed, model) }, error = function(e) return(NULL))
  if (is.null(derived)) return(NULL)
  
  # --- B. Set estimates to NA if actual captures < 10 ---
  est_df <- derived$N %>%
    mutate(
      river = river_name,
      season_label = ch_cols_local,
      n_captured = n_captured_v
    ) %>%
    mutate(
      across(c(N, se, LCL, UCL), ~ ifelse(n_captured < 10, NA, .))
    ) %>%
    rename(N_est = N, N_se = se, N_lcl = LCL, N_ucl = UCL)
  
  return(est_df)
}

# ==============================================================================
# ===== 5. Run Analysis by River =====
# ==============================================================================
river_list <- sort(unique(capture_history_sorted[[RIVER_COL_NAME]]))
all_estimates <- list()

message("\nStarting POPAN Analysis...")

for (riv in river_list) {
  cat(paste0("Analyzing River: ", riv, " ... "))
  river_data <- capture_history_sorted %>% filter(!!sym(RIVER_COL_NAME) == riv)
  
  res <- estimate_N_POPAN(river_data, riv, sorted_cols)
  
  if (!is.null(res)) {
    all_estimates[[length(all_estimates) + 1]] <- res
    cat("Done.\n")
  } else {
    cat("Skipped.\n")
  }
}

# Combine the list into a single dataframe to create base_results
if (length(all_estimates) > 0) {
  base_results <- bind_rows(all_estimates)
  message("Success: base_results created.")
} else {
  stop("Error: No valid analysis results were obtained. Please check your data.")
}

# ===== 6. Process Results (Filter S3 & Years, Add Total) =====
message("\nProcessing results for S3 and calculating Totals...")

# 6-1. Restore year and season columns, and extract S3 only
results_s3 <- base_results %>%
  separate(season_label, into = c("year_str", "season_str"), sep = "-", remove = FALSE) %>%
  mutate(
    year = as.numeric(year_str),
    season = as.numeric(str_replace(season_str, "S", ""))
  ) %>%
  filter(season == 3)

# 6-2. Period filtering per river
filtered_results <- results_s3 %>%
  filter(
    (river == MAIN_RIVER_NAME & year >= 2000 & year <= 2008) |
      (river != MAIN_RIVER_NAME & year >= 2002 & year <= 2008)
  )

# 6-3. Calculate Total (2002-2008 only, sum of all rivers)
total_df <- filtered_results %>%
  filter(year >= 2002, year <= 2008) %>%
  group_by(year, season, season_label) %>%
  summarise(
    N_est = sum(N_est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    river = "Total",
    N_se = NA, N_lcl = NA, N_ucl = NA,
    occasion_idx = NA
  )

# 6-4. Combine data
final_output_df <- bind_rows(filtered_results, total_df)
final_output_df$river <- factor(final_output_df$river, levels = c(river_list, "Total"))
final_output_df <- final_output_df %>% arrange(river, year)

# ===== 7. Export Results & Visualization =====

get_output_path <- function(filename) {
  file.path(OUTPUT_DIR, filename)
}

# CSV Output
csv_filename <- get_output_path(paste0("Abundance_Estimates_S3_Total_", FILE_SUFFIX, ".csv"))
write_csv(final_output_df, csv_filename)
message(paste0("CSV Saved: ", csv_filename))


# --- Create Plots ---

# Condition: Error bars drawn only for non-Total rows and where LCL > 0
errorbar_data <- final_output_df %>% 
  filter(river != "Total" & N_lcl > 0)

# Plot 1: Facet Wrap (Individual Display)
p_s3_facet <- ggplot(final_output_df, aes(x = year, y = N_est, color = river)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Conditional error bars
  geom_errorbar(data = errorbar_data, 
                aes(ymin = N_lcl, ymax = N_ucl), width = 0.2, alpha = 0.5) +
  facet_wrap(~river, scales = "free_y", ncol = 1) +
  labs(
    title = "Estimated Population Abundance (Season 3) - Faceted",
    subtitle = "Error bars hidden if LCL <= 0. Main River: 2000-2015, Others/Total: 2002-2015",
    x = "Year",
    y = "Estimated Abundance (N)"
  ) +
  scale_x_continuous(breaks = 2000:2008) +
  theme_bw() +
  theme(legend.position = "none")

# Plot 2: Overlay (Overlapped Display)
p_s3_overlay <- ggplot(final_output_df, 
                       aes(x = year, y = N_est, color = river, group = river)) +
  geom_line(linewidth = 2, alpha = 0.8) +
  geom_point(size = 4) +
  geom_errorbar(data = errorbar_data, 
                aes(ymin = N_lcl, ymax = N_ucl), width = 0.5, alpha = 0.3) +
  labs(
    title = "Estimated Population Abundance (Season 3) - Overlaid",
    subtitle = "Comparing trends across rivers (Including Total)",
    x = "Year",
    y = "Estimated Abundance (N)"
  ) +
  scale_x_continuous(breaks = 2000:2008) +
  scale_y_continuous(breaks = seq(0, 1400, 200)) + 
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 14), 
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# PDF Output
pdf_filename <- get_output_path(paste0("Abundance_Plots_S3_", FILE_SUFFIX, ".pdf"))
pdf(pdf_filename, width = 10, height = 12)
print(p_s3_facet)
print(p_s3_overlay)
dev.off()
message(paste0("PDF Saved: ", pdf_filename))


# ==============================================================================
# ===== 8. Synchrony Analysis (Growth Rates r) - Matrix Output Only =====
# ==============================================================================
message("\nCalculating Synchrony Matrix (Growth Rates r)...")

sync_base <- filtered_results %>%
  filter(river != "Total", year >= 2002, year <= 2008) %>%
  select(river, year, N_est)

r_df <- sync_base %>%
  arrange(river, year) %>%
  group_by(river) %>%
  mutate(
    N_prev = lag(N_est, order_by = year),
    year_diff = year - lag(year, order_by = year),
    r_rate = if_else(year_diff == 1, log(N_est / N_prev), NA_real_)
  ) %>%
  filter(!is.na(r_rate)) %>%
  select(river, year, r_rate) %>%
  ungroup()

r_wide <- r_df %>%
  pivot_wider(names_from = river, values_from = r_rate) %>%
  arrange(year) %>%
  select(-year) %>%
  as.data.frame()

cor_r <- cor(r_wide, use = "pairwise.complete.obs", method = "pearson")

get_p_value_matrix <- function(data_mat) {
  n <- ncol(data_mat)
  p_mat <- matrix(NA, n, n)
  colnames(p_mat) <- colnames(data_mat)
  rownames(p_mat) <- colnames(data_mat)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tryCatch({
        test <- cor.test(data_mat[, i], data_mat[, j], method = "pearson")
        p_mat[i, j] <- test$p.value
        p_mat[j, i] <- test$p.value
      }, error = function(e) {
        p_mat[i, j] <<- NA; p_mat[j, i] <<- NA
      })
    }
  }
  diag(p_mat) <- 0
  return(p_mat)
}

p_vals_r <- get_p_value_matrix(r_wide)

format_significance <- function(r, p) {
  if(is.na(r) || is.na(p)) return("")
  sig_mark <- case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".", 
    TRUE      ~ ""
  )
  return(paste0(sprintf("%.2f", r), sig_mark))
}

r_vec <- as.vector(as.matrix(cor_r))
p_vec <- as.vector(p_vals_r)
r_vec[is.na(r_vec)] <- 0
p_vec[is.na(p_vec)] <- 1

formatted_vals <- mapply(format_significance, r_vec, p_vec)
formatted_r_matrix <- matrix(formatted_vals, nrow = nrow(cor_r))
colnames(formatted_r_matrix) <- colnames(cor_r)
rownames(formatted_r_matrix) <- rownames(cor_r)

# Console Output
message("\n[Result] Synchrony Matrix based on Growth Rates (r):")
res_sync_df <- as.data.frame(formatted_r_matrix) %>% 
  rownames_to_column("river")

print(res_sync_df)

# --- 8-4. CSV Output (2 types) ---

# 1. Matrix with Significance Stars (Formatted)
sync_csv_formatted <- get_output_path(paste0("Synchrony_Matrix_Formatted_", FILE_SUFFIX, ".csv"))
write_csv(res_sync_df, sync_csv_formatted)
message(paste0("Formatted Matrix Saved: ", sync_csv_formatted))

# 2. P-value Matrix (Raw P-values)
p_val_output <- as.data.frame(p_vals_r) %>% 
  rownames_to_column("river")

sync_csv_pvals <- get_output_path(paste0("Synchrony_Matrix_P_Values_", FILE_SUFFIX, ".csv"))
write_csv(p_val_output, sync_csv_pvals)
message(paste0("P-value Matrix Saved:   ", sync_csv_pvals))

message("\nAll Analysis Completed.")

