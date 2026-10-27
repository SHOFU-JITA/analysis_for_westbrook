# =============================================================================
# 0. Environment Reset & Package Loading
# =============================================================================
rm(list = ls(all.names = TRUE))

packages <- c("readr", "dplyr", "lubridate", "ggplot2", "tidyr", "scales", 
              "broom", "purrr", "nnet", "emmeans", "patchwork", "ggh4x", "RColorBrewer")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(SCRIPT_DIR)
} else {
  SCRIPT_DIR <- getwd()
}

# =============================================================================
# 1. Basic Settings & Data Loading
# =============================================================================
file_path <- "../Clustering/WB_bkt_clustered_rev2.csv"
OUTPUT_FILENAME <- "Combined_Cluster_Dynamics_k4_vs_k9.pdf"

years_wb <- 2000:2015     
years_others <- 2002:2015 
min_year <- min(years_wb) 

river_levels <- c("Metapopulation", "WB", "OL", "OS", "IL")

# -----------------------------------------------------------------------------
# Customize Color Palettes
# -----------------------------------------------------------------------------
# Get all 9 colors from Set1 
# (1:Red, 2:Blue, 3:Green, 4:Purple, 5:Orange, 6:Yellow, 7:Brown, 8:Pink, 9:Gray)
set1_colors <- brewer.pal(9, "Set1")

# Assign colors for k=4
order_k4 <- c(4, 3, 1, 2)
COLORS_K4 <- set1_colors[order_k4]

# Assign colors for k=9
order_k9 <- c(4, 3, 1, 2, 5, 6, 7, 8, 9)
COLORS_K9 <- set1_colors[order_k9]
# -----------------------------------------------------------------------------

# Common plotting theme
shared_theme <- theme_bw(base_size = 20) +
  theme(
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold", color = "black"),
    strip.background = element_blank(),  
    panel.border = element_blank(),                        
    axis.line = element_line(color = "black", linewidth = 1), 
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = "bottom"
  )

# Function to standardize river names
rename_river_func <- function(df, col_name) {
  df %>%
    mutate(!!sym(col_name) := case_when(
      !!sym(col_name) == "west brook"  ~ "WB",
      !!sym(col_name) == "wb jimmy"    ~ "OL",
      !!sym(col_name) == "wb mitchell" ~ "OS",
      !!sym(col_name) == "wb obear"    ~ "IL",
      !!sym(col_name) == "All Rivers"  ~ "Metapopulation", 
      TRUE                             ~ as.character(!!sym(col_name))
    )) %>%
    filter(!!sym(col_name) %in% river_levels) %>%
    mutate(!!sym(col_name) := factor(!!sym(col_name), levels = river_levels))
}

df_raw <- read_csv(file_path, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

# =============================================================================
# 2. Define Plotting Function (Data Processing, Modeling, and Plotting)
# =============================================================================
create_cluster_plot <- function(data_raw, k_val, cluster_colors, is_left_plot = TRUE) {
  
  TARGET_CLUSTER_COL <- paste0("cluster_k", k_val)
  ref_cluster <- as.character(k_val)
  
  df_k <- data_raw %>%
    rename(ocluster_rename_dynamic = !!sym(TARGET_CLUSTER_COL)) %>%
    mutate(ocluster_rename_dynamic = relevel(factor(ocluster_rename_dynamic, levels = as.character(1:k_val)), ref = ref_cluster))
  
  # Identify 'Multi' tags
  river_tags <- df_k %>% filter(cohort %in% years_wb) %>% distinct(tag, river)
  multi_tags <- river_tags %>% count(tag) %>% filter(n > 1) %>% pull(tag)
  
  df_model_base <- df_k %>%
    mutate(
      river_grp = if_else(tag %in% multi_tags, "Multi", river),
      cohort_numeric = as.numeric(as.character(cohort))
    )
  
  df_model_all <- df_model_base %>% mutate(river_grp = "All Rivers")
  df_model <- bind_rows(df_model_base, df_model_all) %>%
    filter(
      (river_grp == "west brook" & cohort_numeric %in% years_wb) |
        (river_grp != "west brook" & cohort_numeric %in% years_others)
    ) %>%
    mutate(year_centered = cohort_numeric - min_year)
  
  # Calculate observed proportions
  proportion_data <- df_model %>%
    group_by(river_grp, cohort, cohort_numeric, year_centered, ocluster_rename_dynamic) %>%
    summarise(n_cluster_tags = n_distinct(tag), .groups = "drop") %>%
    group_by(river_grp, cohort) %>%
    mutate(n_total_tags = sum(n_cluster_tags)) %>%
    ungroup() %>%
    mutate(proportion = n_cluster_tags / n_total_tags)
  
  plot_data_full <- proportion_data %>%
    complete(river_grp, cohort_numeric = years_wb, ocluster_rename_dynamic, fill = list(proportion = 0)) %>%
    mutate(
      year_centered = cohort_numeric - min_year,
      proportion = if_else(river_grp != "west brook" & cohort_numeric < 2002, NA_real_, proportion)
    )
  
  # Calculate Simpson Diversity Index (only for Metapopulation)
  simpson_data <- plot_data_full %>%
    filter(!is.na(proportion)) %>% 
    group_by(river_grp, cohort_numeric) %>%
    summarise(simpson_index = 1 - sum(proportion^2, na.rm = TRUE), .groups = "drop") %>%
    filter(river_grp == "All Rivers")
  
  # Model 1: All Rivers (Independent)
  multinom_results_all <- df_model %>%
    filter(river_grp == "All Rivers") %>%
    group_by(river_grp) %>% nest() %>%
    mutate(
      model = map(data, ~multinom(ocluster_rename_dynamic ~ year_centered, data = ., trace = FALSE)),
      predicted = map(model, ~{
        new_data <- tibble(cohort_numeric = years_others, year_centered = years_others - min_year)
        predict(., newdata = new_data, type = "probs") %>%
          as.data.frame() %>% mutate(cohort_numeric = new_data$cohort_numeric) %>%
          pivot_longer(cols = -cohort_numeric, names_to = "ocluster_rename_dynamic", values_to = "predicted_prop")
      })
    )
  
  predicted_data_all <- multinom_results_all %>% select(river_grp, predicted) %>% unnest(predicted)
  plot_data_all <- plot_data_full %>%
    filter(river_grp == "All Rivers") %>%
    left_join(predicted_data_all, by = c("river_grp", "cohort_numeric", "ocluster_rename_dynamic"))
  
  # Model 2: Tributaries (Interaction)
  df_interaction <- df_model %>% filter(river_grp != "All Rivers")
  
  full_model <- multinom(
    ocluster_rename_dynamic ~ year_centered * river_grp, 
    data = df_interaction, 
    trace = FALSE, 
    maxit = 1000, 
    Hess = TRUE
  )
  
  new_data_interaction <- plot_data_full %>% 
    filter(river_grp != "All Rivers") %>% 
    filter((river_grp == "west brook" & cohort_numeric %in% years_wb) |
             (river_grp != "west brook" & cohort_numeric %in% years_others)) %>%
    distinct(cohort_numeric, year_centered, river_grp)
  
  predicted_interaction <- as.data.frame(predict(full_model, newdata = new_data_interaction, type = "probs")) %>%
    mutate(cohort_numeric = new_data_interaction$cohort_numeric, river_grp = new_data_interaction$river_grp) %>%
    pivot_longer(cols = -c(cohort_numeric, river_grp), names_to = "ocluster_rename_dynamic", values_to = "predicted_prop")
  
  plot_data_interaction <- plot_data_full %>%
    filter(river_grp != "All Rivers") %>%
    left_join(predicted_interaction, by = c("river_grp", "cohort_numeric", "ocluster_rename_dynamic"))
  
  # Combine data and standardize names
  data_vertical <- bind_rows(
    plot_data_all %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop),
    plot_data_interaction %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop)
  )
  
  data_vertical <- rename_river_func(data_vertical, "river_grp")
  simpson_data  <- rename_river_func(simpson_data, "river_grp")
  
  n_labels <- df_model %>%
    rename_river_func("river_grp") %>%
    group_by(river_grp) %>%
    summarise(total_n = n_distinct(tag), .groups = "drop") %>%
    mutate(label = paste0("n = ", total_n))
  
  # Set X-axis scales for patchwork
  x_scale_show <- scale_x_continuous(breaks = seq(2000, 2015, by = 2), limits = c(min(years_wb), max(years_wb)))
  x_scale_hide <- scale_x_continuous(breaks = seq(2000, 2015, by = 2), limits = c(min(years_wb), max(years_wb)), labels = NULL)
  x_scales_list <- list(x_scale_show, x_scale_hide, x_scale_hide, x_scale_hide, x_scale_show)
  
  # Generate plot
  p <- ggplot(data_vertical, aes(x = cohort_numeric)) +
    geom_point(aes(y = proportion, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), size = 3, alpha = 0.6) +
    geom_line(aes(y = predicted_prop, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), linewidth = 2) +
    
    # Diversity Index (plotted only for Metapopulation)
    geom_line(data = simpson_data, aes(y = simpson_index), color = "black", linetype = "dashed", linewidth = 1.2) +
    geom_point(data = simpson_data, aes(y = simpson_index), color = "black", shape = 17, size = 4) +
    
    facet_wrap(~ river_grp, ncol = 1, scales = "free_x") + 
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    ggh4x::facetted_pos_scales(x = x_scales_list) + 
    geom_text(data = n_labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.2, vjust = 1.5, size = 6, fontface = "bold", color = "black") +
    
    labs(
      x = "Year (Cohort)", 
      y = if(is_left_plot) "Proportion / Diversity Index" else "", # Hide Y-axis label on the right plot for a cleaner look
      color = paste0("Cluster (k=", k_val, ")")
    ) +
    scale_color_manual(
      values = cluster_colors[1:k_val],
      # Wrap legend into 2 rows with 5 items per row
      guide = guide_legend(ncol = 5, byrow = TRUE) 
    ) +
    shared_theme +
    ggh4x::force_panelsizes(rows = c(1, 1, 1, 1, 1))
  
  return(p)
}

# =============================================================================
# 3. Generate and Combine Plots
# =============================================================================
cat("\nProcessing k=4...\n")
p_left <- create_cluster_plot(df_raw, k_val = 4, cluster_colors = COLORS_K4, is_left_plot = TRUE)

cat("Processing k=9...\n")
p_right <- create_cluster_plot(df_raw, k_val = 9, cluster_colors = COLORS_K9, is_left_plot = FALSE)

# Combine left and right plots
# Use guides="keep" to display legends for both k=4 and k=9
combined_plot <- p_left + p_right + 
  plot_layout(ncol = 2, guides = "keep") & 
  theme(legend.position = "bottom")

# =============================================================================
# 4. Output to PDF
# =============================================================================
n_panels <- length(river_levels)
PANEL_HEIGHT <- 3.5 
total_pdf_height <- (PANEL_HEIGHT * n_panels) + 2

output_path <- file.path(SCRIPT_DIR, OUTPUT_FILENAME)
ggsave(filename = output_path, plot = combined_plot, width = 14, height = total_pdf_height)

cat("\n✅ Saved successfully: ", output_path, "\n")
