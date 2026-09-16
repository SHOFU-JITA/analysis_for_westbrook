# ===================================================================
# Westbrook Portfolio Effect Analysis (Dynamic Config Version)
# ===================================================================
# ===== 0. Reset environment =====
rm(list = ls(all.names = TRUE))
CLUSTER_NUM <- 3

TARGET_CLUSTER_COL <- paste0("cluster_k", CLUSTER_NUM)
OUTPUT_FILENAME    <- paste0("CV_bkt", CLUSTER_NUM, ".pdf")
OUTPUT_FILENAME_2  <- paste0("CV_2periods_cluster_bkt", CLUSTER_NUM, ".pdf")


# ===================================================================
# Automatic Processing below
# ===================================================================

# ===== 1.Package Loading =====
packages <- c("readr", "dplyr", "lubridate", "tidyr", "ggplot2", "ggrepel")

if (!require("rstudioapi")) install.packages() 

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

# ===== 3. Data Loading and Column Standardization =====
# Adjust file path as needed
data <- read_csv("../csv/WB_combined_bkt_rev2.csv")
data$date <- as.Date(data$date)

# Standardize the target cluster column to 'cluster_id' to allow switching cluster numbers without changing subsequent code
if (!TARGET_CLUSTER_COL %in% names(data)) {
  stop(paste("Error: Column", TARGET_CLUSTER_COL, "not found in the data. Please check your settings."))
}

# Rename cluster column and recode river names
data <- data %>%
  rename(cluster_id = all_of(TARGET_CLUSTER_COL)) %>%
  mutate(river = case_when(
    river == "west brook" ~ "WB",
    river == "wb jimmy" ~ "OL",
    river == "wb mitchell" ~ "OS",
    river == "wb obear" ~ "IL",
    TRUE ~ river
  ))

message("Using cluster column: ", TARGET_CLUSTER_COL)
message("Output 1 will be saved as: ", OUTPUT_FILENAME)
message("Output 2 will be saved as: ", OUTPUT_FILENAME_2)

# ===================================================================
# Calculation Part
# ===================================================================

# ===== 4. Define Seasons and Exclude NA Clusters =====
data <- data %>%
  mutate(
    year = year(date),
    season = quarter(date, with_year = FALSE)
  ) %>%
  filter(!is.na(cluster_id))

# ===== 5. Impute Residence Period =====
stay_data <- data %>%
  arrange(tag, date) %>%
  group_by(tag) %>%
  mutate(
    next_date = lead(date),
    end_date = if_else(!is.na(next_date), next_date, date)
  ) %>%
  ungroup()

# ===== 6. Expand to Seasonal Units =====
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

# ===== 7. Aggregation and Filtering (S3 only) =====
alive_by_season <- expand_season %>%
  group_by(river, cluster_id, season_seq) %>%
  summarise(n_alive = n_distinct(tag), .groups = "drop")

alive_by_season_filtered <- alive_by_season %>%
  mutate(
    year = as.numeric(substr(season_seq, 1, 4))
  ) %>%
  filter(year >= 2002 & year <= 2015) %>%
  # Extract only S3 (3rd quarter) data
  filter(grepl("-S3$", season_seq))

# ===== 8. Calculate CV at three levels =====

# 8A. River x Cluster
cv_by_river_cluster <- alive_by_season_filtered %>%
  group_by(river, cluster_id) %>%
  summarise(
    mean_alive = mean(n_alive),
    cv = sd(n_alive) / mean(n_alive),
    .groups = 'drop'
  )

# 8B. River only
cv_by_river <- alive_by_season_filtered %>%
  group_by(river, season_seq) %>%
  summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
  group_by(river) %>%
  summarise(
    mean_alive = mean(n_alive_sum),
    cv = sd(n_alive_sum) / mean(n_alive_sum),
    .groups = 'drop'
  )

# 8C. Cluster only
cv_by_cluster <- alive_by_season_filtered %>%
  group_by(cluster_id, season_seq) %>%
  summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
  group_by(cluster_id) %>%
  summarise(
    mean_alive = mean(n_alive_sum),
    cv = sd(n_alive_sum) / mean(n_alive_sum),
    .groups = 'drop'
  )

# ===== 9. Calculate Total and Weighted Mean CVs =====

# Observed total CV
cv_total <- alive_by_season_filtered %>%
  group_by(season_seq) %>%
  summarise(total_n_alive = sum(n_alive)) %>%
  summarise(cv = sd(total_n_alive) / mean(total_n_alive)) %>%
  pull(cv)

# Expected weighted mean CV
cv_w_mean_river_cluster <- weighted.mean(cv_by_river_cluster$cv, w = cv_by_river_cluster$mean_alive, na.rm = TRUE)
cv_w_mean_river <- weighted.mean(cv_by_river$cv, w = cv_by_river$mean_alive, na.rm = TRUE)
cv_w_mean_cluster <- weighted.mean(cv_by_cluster$cv, w = cv_by_cluster$mean_alive, na.rm = TRUE)

# ===== 10. Create Summary Data =====
summary_metrics <- data.frame(
  metric = c("Expected CV (No Diversity)", "Expected CV (Cluster Diversity Only)", "Expected CV (River Diversity Only)", "Observed Total CV"),
  cv = c(cv_w_mean_river_cluster, cv_w_mean_river, cv_w_mean_cluster, cv_total)
) %>%
  mutate(metric = factor(metric, levels = metric))

# Effect size (Delta CV)
pe_cluster_effect <- cv_w_mean_river_cluster - cv_w_mean_river
pe_river_effect <- cv_w_mean_river_cluster - cv_w_mean_cluster

# ===================================================================
# Visualization Part
# ===================================================================

# ===== 11. Prepare Coordinates for Plotting =====

total_mean_abundance <- sum(cv_by_river_cluster$mean_alive, na.rm = TRUE)

river_cluster_plot_data <- cv_by_river_cluster %>%
  mutate(
    proportion = mean_alive / total_mean_abundance,
    cluster_label = paste0("c", cluster_id)
  )

# Define coordinates
max_prop <- max(river_cluster_plot_data$proportion, na.rm = TRUE)
pos_vline       <- max_prop * 1.1   # Vertical separator line
pos_summary_pts <- max_prop * 1.2   # Summary points
pos_arrow_clus  <- max_prop * 1.15  # Cluster Effect arrow
pos_arrow_river <- max_prop * 1.6   # River Effect arrow
x_limit_max     <- max_prop * 1.8   # Overall max width

summary_metrics$x_pos <- pos_summary_pts

# Data for reference lines
ref_lines <- data.frame(
  y = c(cv_w_mean_river_cluster, cv_w_mean_river, cv_w_mean_cluster),
  x_start = pos_vline,
  x_end = x_limit_max
)

# ===== 12. Create Portfolio Decomposition Plot =====

portfolio_decomposition_plot <- ggplot() +
  
  # A. Background reference lines
  geom_segment(data = ref_lines,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               linetype = "dashed", color = "gray60", size = 0.6) +
  
  # B. Vertical separator line
  geom_vline(xintercept = pos_vline,
             linetype = "dotted", color = "grey50", size = 0.8) +
  
  # C. Left side: Individual plots
  geom_point(data = river_cluster_plot_data,
             aes(x = proportion, y = cv, color = river),
             size = 8, alpha = 0.8) +
  geom_text(data = river_cluster_plot_data,
            aes(x = proportion, y = cv, label = cluster_label),
            color = "white", size = 3.5, fontface = "bold") +
  
  # D. Right side: Summary points
  geom_point(data = summary_metrics, aes(x = x_pos, y = cv),
             size = 8, shape = 18, color = "black") +
  geom_text(data = summary_metrics, aes(x = x_pos + 0.008, y = cv, label = metric),
            hjust = 0, size = 5, color = "black", fontface = "plain") +
  
  # E. Annotation: Cluster Effect
  annotate("segment",
           x = pos_arrow_clus, xend = pos_arrow_clus,
           y = cv_w_mean_river_cluster, yend = cv_w_mean_river,
           arrow = arrow(length = unit(0.3, "cm"), ends = "both"),
           colour = "blue", size = 1.2) +
  annotate("text",
           x = pos_arrow_clus + 0.008,
           y = (cv_w_mean_river_cluster + cv_w_mean_river) / 2,
           label = paste0("Cluster Effect\n(ΔCV = ", round(pe_cluster_effect, 3), ")"),
           color = "blue", hjust = 0, size = 4.5, lineheight = 0.9, fontface = "bold") +
  
  # F. Annotation: River Effect
  annotate("segment",
           x = pos_arrow_river, xend = pos_arrow_river,
           y = cv_w_mean_river_cluster, yend = cv_w_mean_cluster,
           arrow = arrow(length = unit(0.3, "cm"), ends = "both"), 
           colour = "red", size = 1.2) +
  annotate("text",
           x = pos_arrow_river + 0.008,
           y = (cv_w_mean_river_cluster + cv_w_mean_cluster) / 2,
           label = paste0("River Effect\n(ΔCV = ", round(pe_river_effect, 3), ")"),
           color = "red", hjust = 0, size = 4.5, lineheight = 0.9, fontface = "bold") +
  
  # G. Theme settings
  labs(
    title = paste("Decomposition of Portfolio Effect (Brook Trout)"),
    x = "Proportion of Total Mean Abundance",
    y = "Coefficient of Variation (CV)",
    color = "River"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, color = "black") 
  ) +
  coord_cartesian(xlim = c(0, x_limit_max))

# ===== 13. Print Plot =====
print(portfolio_decomposition_plot)

# ===== 14. Save as PDF =====
# Save to current working directory
ggsave(
  filename = OUTPUT_FILENAME,
  plot = portfolio_decomposition_plot,
  device = cairo_pdf, 
  width = 15,
  height = 9,
  units = "in",
  dpi = 300
)

message("Processing complete.")
message("Saved to: ", file.path(getwd(), OUTPUT_FILENAME))

# ===================================================================
# 15. Cluster Contribution Analysis (Target vs. Rest Integrated Approach)
# ===================================================================
cat("\n\n======================================================\n")
cat("Approach 3: Cluster Effect Contribution via Target vs. Rest Integration\n")
cat("======================================================\n")

clusters <- unique(alive_by_season_filtered$cluster_id)
contribution_results <- data.frame()

# Base expected CV (all clusters integrated, river effect only)
base_w_mean_river <- cv_w_mean_river 

for (target_c in clusters) {
  
  # 1. Reclassify data into 'Target' and 'Others' and aggregate
  temp_merged_data <- alive_by_season_filtered %>%
    mutate(merged_cluster = if_else(cluster_id == target_c, as.character(target_c), "Others")) %>%
    # Sum up 'Others' clusters by river and season
    group_by(river, merged_cluster, season_seq) %>%
    summarise(n_alive = sum(n_alive), .groups = "drop")
  
  # 2. Calculate CV for 'River x Integrated Cluster'
  temp_cv_rc <- temp_merged_data %>%
    group_by(river, merged_cluster) %>%
    summarise(
      m = mean(n_alive), 
      cv = sd(n_alive) / mean(n_alive), 
      .groups = "drop"
    )
  
  # 3. Weighted mean CV in integrated state
  temp_w_mean_rc <- weighted.mean(temp_cv_rc$cv, w = temp_cv_rc$m, na.rm = TRUE)
  
  # 4. Contribution = (CV of Target vs Others) - (CV of all integrated)
  target_contribution <- temp_w_mean_rc - base_w_mean_river
  
  # Save results
  contribution_results <- rbind(contribution_results, data.frame(
    cluster_id = target_c,
    contribution_effect = target_contribution
  ))
}

print(contribution_results)
cat("* contribution_effect: CV reduction independently generated by the cluster (portfolio effect)\n")
cat("* contribution_pct: Relative percentage of each cluster's contribution\n")

# ===================================================================
# 15. Cluster Contribution Analysis (Target vs. Rest) - Period Split Version
# ===================================================================
cat("\n\n======================================================\n")
cat("Approach 3: Cluster Effect Contribution by First/Second Half Split\n")
cat("======================================================\n")

# Define periods
periods <- list(
  "First_Half (2002-2008)" = 2002:2008,
  "Second_Half (2009-2015)" = 2009:2015
)

contribution_results_periods <- data.frame()

for (period_name in names(periods)) {
  target_years <- periods[[period_name]]
  
  # 1. Extract data for the specific period
  period_data <- alive_by_season_filtered %>%
    mutate(year_num = as.numeric(substr(season_seq, 1, 4))) %>%
    filter(year_num %in% target_years)
  
  # 2. Recalculate baseline CV for the period (all clusters integrated, river only)
  cv_by_river_period <- period_data %>%
    group_by(river, season_seq) %>%
    summarise(n_alive_sum = sum(n_alive), .groups = "drop") %>%
    group_by(river) %>%
    summarise(
      m = mean(n_alive_sum),
      cv = sd(n_alive_sum) / mean(n_alive_sum),
      .groups = 'drop'
    )
  base_w_mean_river_period <- weighted.mean(cv_by_river_period$cv, w = cv_by_river_period$m, na.rm = TRUE)
  
  clusters <- unique(period_data$cluster_id)
  
  # 3. Calculate Target vs. Rest for each cluster
  for (target_c in clusters) {
    
    temp_merged_data <- period_data %>%
      mutate(merged_cluster = if_else(cluster_id == target_c, as.character(target_c), "Others")) %>%
      group_by(river, merged_cluster, season_seq) %>%
      summarise(n_alive = sum(n_alive), .groups = "drop")
    
    temp_cv_rc <- temp_merged_data %>%
      group_by(river, merged_cluster) %>%
      summarise(
        m = mean(n_alive), 
        cv = sd(n_alive) / mean(n_alive), 
        .groups = "drop"
      )
    
    temp_w_mean_rc <- weighted.mean(temp_cv_rc$cv, w = temp_cv_rc$m, na.rm = TRUE)
    
    # Calculate contribution
    target_contribution <- temp_w_mean_rc - base_w_mean_river_period
    
    contribution_results_periods <- rbind(contribution_results_periods, data.frame(
      Period = period_name,
      cluster_id = target_c,
      contribution_effect = target_contribution
    ))
  }
}

# 4. Use only Period data for plotting
combined_contribution <- contribution_results_periods[, c("Period", "cluster_id", "contribution_effect")]

# Fix display order for plots
combined_contribution$Period <- factor(
  combined_contribution$Period,
  levels = c("First_Half (2002-2008)", "Second_Half (2009-2015)")
)

print(as.data.frame(combined_contribution))

# ===================================================================
# 16. Visualization and PDF Output of Cluster Contributions (Period Comparison)
# ===================================================================
# Create plot (no title, large black axis labels)
contribution_plot <- ggplot(combined_contribution, aes(x = factor(cluster_id), y = contribution_effect, fill = Period)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Increase base text size
  theme_bw(base_size = 16) + 
  # Color settings
  scale_fill_manual(values = c("First_Half (2002-2008)" = "#4E79A7", 
                               "Second_Half (2009-2015)" = "#F28E2B")) +
  labs(
    x = "Cluster",
    y = "Contribution Effect (CV Reduction)",
    fill = "Period"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(), # Remove vertical grid lines on X-axis
    
    # Darken and enlarge axis titles
    axis.title = element_text(color = "black", size = 18, face = "bold"),
    # Darken and enlarge axis text
    axis.text = element_text(color = "black", size = 18, face = "bold"),
    # Adjust legend title and text
    legend.title = element_text(color = "black", size = 18, face = "bold"),
    legend.text = element_text(color = "black", size = 18, face = "bold")
  )

# Print plot
print(contribution_plot)

# Save as PDF
ggsave(
  filename = OUTPUT_FILENAME_2,
  plot = contribution_plot,
  device = cairo_pdf, 
  width = 11,   
  height = 6,
  units = "in",
  dpi = 300
)

cat("Processing complete.\n")
cat("Plot saved as '", file.path(getwd(), OUTPUT_FILENAME_2), "'\n", sep="")

