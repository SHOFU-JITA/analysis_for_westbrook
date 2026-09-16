# ===================================================================
# Westbrook Portfolio Effect Analysis 
# ===================================================================
rm(list = ls(all.names = TRUE))
# ===== 0. User Configuration =====
TARGET_CLUSTER_COL <- "cluster_k3"
OUTPUT_FILENAME     <- "Population_dynamics_cluster_k3.pdf"

# Specify cluster IDs (adjust this list if the number of clusters changes)
CLUSTER_IDS <- c("1", "2", "3") 

# Output filenames for Synchrony Analysis (Prefixes and Plot)
SYNC_FORMATTED_PREFIX     <- "Synchrony_Cluster_"
SYNC_PVALS_PREFIX         <- "Synchrony_Cluster_P_"
SYNC_ONE_VS_OTHERS_PREFIX <- "Synchrony_cluster_comparison"
SYNC_OVO_PLOT_FILENAME    <- "Synchrony_cluster_comparison_k3.pdf"

# ===================================================================
# Automatic Processing Below
# ===================================================================

# ===== 1. Environment Reset and Package Loading =====


packages <- c("readr", "dplyr", "lubridate", "tidyr", "ggplot2", "gridExtra")
if (!require("rstudioapi")) install.packages("rstudioapi") 

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ===== 2. Set Working Directory =====
if (rstudioapi::isAvailable()) {
  try({
    current_path = rstudioapi::getSourceEditorContext()$path
    if (current_path != "") {
      setwd(dirname(current_path))
      message("Working directory set to: ", getwd())
    }
  }, silent = TRUE)
}

# ===== 3. Data Loading and Column Renaming =====
# Note: Please update the file path according to your environment.
data <- read_csv("../csv/WB_combined_bkt_rev2.csv")

data$date <- as.Date(data$date)

if (!TARGET_CLUSTER_COL %in% names(data)) {
  stop(paste("Error: Column", TARGET_CLUSTER_COL, "not found in the data."))
}

data <- data %>%
  rename(cluster_id = all_of(TARGET_CLUSTER_COL))

# ===== 4. Season Definition and NA Filtering =====
data <- data %>%
  mutate(
    year = year(date),
    season = quarter(date, with_year = FALSE)
  ) %>%
  filter(!is.na(cluster_id))

# ===== 5. Impute Stay Duration =====
stay_data <- data %>%
  arrange(tag, date) %>%
  group_by(tag) %>%
  mutate(
    next_date = lead(date),
    end_date = if_else(!is.na(next_date), next_date, date)
  ) %>%
  ungroup()

# ===== 6. Expand to Seasonal Scale =====
expand_season <- stay_data %>%
  rowwise() %>%
  mutate(
    season_seq = list({
      dseq <- seq(date, end_date, by = "day")
      unique(paste0(year(dseq), "-S", quarter(dseq)))
    })
  ) %>%
  unnest(season_seq) %>%
  ungroup() %>%
  select(tag, river, cluster_id, season_seq) %>%
  distinct()

# ===== 7. Aggregation and Filtering (S3 Only) =====
alive_by_season <- expand_season %>%
  group_by(river, cluster_id, season_seq) %>%
  summarise(n_alive = n_distinct(tag), .groups = "drop")

alive_by_season_filtered <- alive_by_season %>%
  mutate(year = as.numeric(substr(season_seq, 1, 4))) %>%
  filter(year >= 2002 & year <= 2015) %>%
  filter(grepl("-S3$", season_seq))

# Create "All Rivers" data and combine it
all_rivers_data <- alive_by_season_filtered %>%
  group_by(cluster_id, season_seq, year) %>% 
  summarise(n_alive = sum(n_alive), .groups = "drop") %>%
  mutate(river = "All Rivers") 

alive_by_season_combined <- bind_rows(alive_by_season_filtered, all_rivers_data)

# ===== 8. Visualization of Abundance Trends by River and Cluster (S3) =====

# 1. Rename rivers and convert to factors
alive_by_season_combined <- alive_by_season_combined %>%
  mutate(
    river = case_when(
      river == "west brook"   ~ "WB",
      river == "wb jimmy"    ~ "OL",
      river == "wb mitchell" ~ "OS",
      river == "wb obear"    ~ "IL",
      TRUE                   ~ river
    ),
    cluster_id = as.factor(cluster_id),
    river = factor(river, levels = c("All Rivers", "WB", "OL", "OS", "IL"))
  )

# 2. Fill missing periods with 0
alive_by_season_combined <- alive_by_season_combined %>%
  complete(river, year = 2002:2015, cluster_id, fill = list(n_alive = 0))

# 3. Create stacked area chart
p_trend <- ggplot(alive_by_season_combined, aes(x = year, y = n_alive, fill = cluster_id)) +
  geom_area(color = "black", linewidth = 0.4, alpha = 0.8) +   
  scale_fill_brewer(palette = "Set1", name = "Cluster") +      
  facet_wrap(~ river, scales = "free_y", ncol = 1) +           
  scale_x_continuous(breaks = seq(min(alive_by_season_combined$year, na.rm = TRUE), 
                                  max(alive_by_season_combined$year, na.rm = TRUE), by = 2)) +
  labs(
    x = "Year",
    y = "Number of individuals"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),                        
    axis.title  = element_text(color = "black", size = 20),                        
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),                
    legend.position = "bottom"
  )

print(p_trend)

# ===== 9. Save Results to PDF =====

# Set height per panel (inches)
PANEL_HEIGHT <- 2.0  

# Get the number of panels (rivers)
n_panels <- nlevels(alive_by_season_combined$river)

# Calculate total PDF height automatically
total_pdf_height <- (PANEL_HEIGHT * n_panels) + 3

ggsave(filename = OUTPUT_FILENAME, plot = p_trend, width = 6, height = total_pdf_height)
message("Saved plot to ", OUTPUT_FILENAME, " (Auto-adjusted height: ", total_pdf_height, " inches)")


# ==============================================================================
# ===== 10. Synchrony Analysis (Cluster Growth Rates r) =====
# ==============================================================================
message("\nPreparing Synchrony Analysis for Clusters...")

# --- 10-1. Define Helper Functions ---

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

# Core analysis and export function
analyze_and_save_synchrony <- function(r_df_input, period_label) {
  r_wide <- r_df_input %>%
    pivot_wider(names_from = cluster_id, values_from = r_rate) %>%
    arrange(year) %>%
    select(-year) %>%
    as.data.frame()
  
  if (ncol(r_wide) < 2 || nrow(r_wide) < 3) {
    message(paste("Skipping", period_label, ": Not enough data or clusters."))
    return(NULL)
  }
  
  cor_r <- cor(r_wide, use = "pairwise.complete.obs", method = "pearson")
  p_vals_r <- get_p_value_matrix(r_wide)
  
  r_vec <- as.vector(as.matrix(cor_r))
  p_vec <- as.vector(p_vals_r)
  r_vec[is.na(r_vec)] <- 0
  p_vec[is.na(p_vec)] <- 1
  
  formatted_vals <- mapply(format_significance, r_vec, p_vec)
  formatted_r_matrix <- matrix(formatted_vals, nrow = nrow(cor_r))
  colnames(formatted_r_matrix) <- colnames(cor_r)
  rownames(formatted_r_matrix) <- rownames(cor_r)
  
  message(paste0("\n[Result] Synchrony Matrix (", period_label, ") :"))
  res_sync_df <- as.data.frame(formatted_r_matrix) %>% 
    tibble::rownames_to_column("cluster_id")
  print(res_sync_df)
  
  # Export Formatted Matrix
  sync_csv_formatted <- paste0(SYNC_FORMATTED_PREFIX, period_label, "_", TARGET_CLUSTER_COL, ".csv")
  write_csv(res_sync_df, sync_csv_formatted)
  message(paste0("Formatted Matrix Saved: ", sync_csv_formatted))
  
  # Export Raw P-values
  p_val_output <- as.data.frame(p_vals_r) %>% 
    tibble::rownames_to_column("cluster_id")
  sync_csv_pvals <- paste0(SYNC_PVALS_PREFIX, period_label, "_", TARGET_CLUSTER_COL, ".csv")
  write_csv(p_val_output, sync_csv_pvals)
  message(paste0("P-value Matrix Saved:   ", sync_csv_pvals))
}

# --- 10-2. Calculate Base Growth Rates (r) ---

sync_base <- alive_by_season_combined %>%
  filter(river == "All Rivers", year >= 2002, year <= 2015) %>%
  select(cluster_id, year, n_alive)

r_all <- sync_base %>%
  arrange(cluster_id, year) %>%
  group_by(cluster_id) %>%
  mutate(
    n_prev = lag(n_alive, order_by = year),
    year_diff = year - lag(year, order_by = year),
    r_rate = if_else(year_diff == 1 & n_prev > 0 & n_alive > 0, log(n_alive / n_prev), NA_real_)
  ) %>%
  filter(!is.na(r_rate)) %>%
  select(cluster_id, year, r_rate) %>%
  ungroup()

# --- 10-3. Run Analysis for Each Period ---

# 1. Full Period (2002-2015)
analyze_and_save_synchrony(r_all, "Full_2002-2015")

# 2. Early Period (2002-2008)
r_early <- r_all %>% filter(year >= 2002 & year <= 2008)
analyze_and_save_synchrony(r_early, "Early_2002-2008")

# 3. Late Period (2009-2015)
r_late <- r_all %>% filter(year >= 2009 & year <= 2015)
analyze_and_save_synchrony(r_late, "Late_2009-2015")


# ==============================================================================
# ===== 11. One vs Others Synchrony Analysis =====
# ==============================================================================
message("\nPreparing One vs Others Synchrony Analysis for Clusters...")

# --- 11-1. Calculate Base Growth Rates for One vs Others ---

r_one_vs_others_list <- list()

for (cid in CLUSTER_IDS) {
  df_binary <- sync_base %>%
    mutate(is_target = if_else(cluster_id == cid, "Target", "Others")) %>%
    group_by(is_target, year) %>%
    summarise(n_alive = sum(n_alive), .groups = "drop")
  
  df_r <- df_binary %>%
    arrange(is_target, year) %>%
    group_by(is_target) %>%
    mutate(
      n_prev = lag(n_alive, order_by = year),
      year_diff = year - lag(year, order_by = year),
      r_rate = if_else(year_diff == 1 & n_prev > 0 & n_alive > 0, log(n_alive / n_prev), NA_real_)
    ) %>%
    filter(!is.na(r_rate)) %>%
    mutate(target_cluster = cid) %>%
    select(target_cluster, is_target, year, r_rate) %>%
    ungroup()
  
  r_one_vs_others_list[[cid]] <- df_r
}

r_one_vs_others_all <- bind_rows(r_one_vs_others_list)

# --- 11-2. Define Period Analysis Function ---

analyze_one_vs_others_period <- function(r_df_input, period_label, cluster_ids) {
  summary_list <- list()
  
  for (cid in cluster_ids) {
    df_wide <- r_df_input %>%
      filter(target_cluster == cid) %>%
      pivot_wider(names_from = is_target, values_from = r_rate) %>%
      arrange(year)
    
    if (!all(c("Target", "Others") %in% names(df_wide)) || nrow(df_wide) < 3) {
      summary_list[[cid]] <- tibble(
        cluster_id = cid,
        correlation = NA_real_,
        p_value = NA_real_,
        formatted = ""
      )
      next
    }
    
    tryCatch({
      test_res <- cor.test(df_wide$Target, df_wide$Others, method = "pearson")
      r_val <- test_res$estimate
      p_val <- test_res$p.value
      fmt_val <- format_significance(r_val, p_val)
    }, error = function(e) {
      r_val <- NA_real_
      p_val <- NA_real_
      fmt_val <- ""
    })
    
    summary_list[[cid]] <- tibble(
      cluster_id = cid,
      correlation = r_val,
      p_value = p_val,
      formatted = fmt_val
    )
  }
  
  summary_df <- bind_rows(summary_list) %>%
    mutate(period = period_label)
  
  message(paste0("\n[Result] One vs Others Synchrony (", period_label, ") :"))
  print(summary_df)
  
  # Export Summary Matrix
  out_csv <- paste0(SYNC_ONE_VS_OTHERS_PREFIX, period_label, "_", TARGET_CLUSTER_COL, ".csv")
  write_csv(summary_df %>% select(-period), out_csv)
  message(paste0("One vs Others Summary Saved: ", out_csv))
  
  return(summary_df)
}

# --- 11-3. Run Analysis for Each Period ---

# 1. Full Period (2002-2015)
res_full <- analyze_one_vs_others_period(r_one_vs_others_all, "Full_2002-2015", CLUSTER_IDS)

# 2. Early Period (2002-2008)
r_ovo_early <- r_one_vs_others_all %>% filter(year >= 2002 & year <= 2008)
res_early <- analyze_one_vs_others_period(r_ovo_early, "Early_2002-2008", CLUSTER_IDS)

# 3. Late Period (2009-2015)
r_ovo_late <- r_one_vs_others_all %>% filter(year >= 2009 & year <= 2015)
res_late <- analyze_one_vs_others_period(r_ovo_late, "Late_2009-2015", CLUSTER_IDS)

# Combine results for plotting (Filtering out 'Full' period and renaming for plot)
ovo_plot_data <- bind_rows(res_early, res_late) %>%
  mutate(
    Period = case_when(
      period == "Early_2002-2008" ~ "First_Half (2002-2008)",
      period == "Late_2009-2015"  ~ "Second_Half (2009-2015)"
    ),
    Period = factor(Period, levels = c("First_Half (2002-2008)", "Second_Half (2009-2015)"))
  )

# --- 11-4. Plot One vs Others Synchrony (Period Split Version) ---
message("\nGenerating One vs Others Synchrony Plot (Period Comparison)...")

p_ovo <- ggplot(ovo_plot_data, aes(x = factor(cluster_id), y = correlation, fill = Period)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_text(aes(
    label = formatted, 
    vjust = if_else(!is.na(correlation) & correlation >= 0, -0.5, 1.5)
  ), 
  position = position_dodge(width = 0.8), size = 6, fontface = "bold") +
  theme_bw(base_size = 16) +
  scale_fill_manual(values = c("First_Half (2002-2008)" = "#4E79A7", 
                               "Second_Half (2009-2015)" = "#F28E2B")) +
  labs(
    x = "Cluster",
    y = "Pearson's Correlation Coefficient (r)",
    fill = "Period"
  ) +
  scale_y_continuous(limits = c(-0.0, 1.2), breaks = seq(-0.0, 1.0, by = 0.5)) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    
    # Darken and enlarge axis titles/text and legend based on reference code
    axis.title = element_text(color = "black", size = 18, face = "bold"),
    axis.text = element_text(color = "black", size = 18, face = "bold"),
    legend.title = element_text(color = "black", size = 18, face = "bold"),
    legend.text = element_text(color = "black", size = 18, face = "bold")
  )

print(p_ovo)

ggsave(
  filename = SYNC_OVO_PLOT_FILENAME, 
  plot = p_ovo, 
  device = cairo_pdf,
  width = 11, 
  height = 6,
  units = "in",
  dpi = 300
)
message("Saved One vs Others Synchrony Plot to ", SYNC_OVO_PLOT_FILENAME)

message("\nAll Analysis Completed.")

