# =============================================================================
# 0. Environment Reset & Package Loading
# =============================================================================
rm(list = ls(all.names = TRUE))

packages <- c("readr", "dplyr", "lubridate", "ggplot2", "tidyr", "scales", 
              "broom", "purrr", "nnet", "emmeans", "patchwork")

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
# File paths (Adjust these to match your environment)
file_path_1 <- "../Clustering/WB_bkt_clustered_rev2.csv"
file_path_2 <- "../csv/WB_combined_bkt_rev2.csv"

# Common settings
num_clusters <- 4
ref_cluster <- as.character(num_clusters)
TARGET_CLUSTER_COL <- paste0("cluster_k", num_clusters)

# Output filename
OUTPUT_FILENAME <- paste0("Combined_Population_Dynamics_k", num_clusters, ".pdf")

# Standardized factor levels for river groupings
river_levels <- c("All Rivers", "WB", "OL", "OS", "IL")

# Shared theme to ensure visual consistency across both plots
shared_theme <- theme_bw(base_size = 20) +
  theme(
    text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20, face = "bold", color = "black"),
    strip.background = element_rect(fill = "gray90", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    legend.position = "bottom"
  )

# =============================================================================
# 2. Left Plot Data: Cluster Proportions & Simpson Index
# =============================================================================
years_wb <- 2000:2015     
years_others <- 2002:2015 
min_year <- min(years_wb) 

df1 <- read_csv(file_path_1, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  rename(ocluster_rename_dynamic = !!sym(TARGET_CLUSTER_COL)) %>%
  mutate(ocluster_rename_dynamic = relevel(factor(ocluster_rename_dynamic, levels = as.character(1:num_clusters)), ref = ref_cluster))

river_tags <- df1 %>% filter(cohort %in% years_wb) %>% distinct(tag, river)
multi_tags <- river_tags %>% count(tag) %>% filter(n > 1) %>% pull(tag)

df_model_base <- df1 %>%
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

# Calculate Simpson Diversity Index
simpson_data <- plot_data_full %>%
  filter(!is.na(proportion)) %>% 
  group_by(river_grp, cohort_numeric) %>%
  summarise(simpson_index = 1 - sum(proportion^2, na.rm = TRUE), .groups = "drop")

# 'All Rivers' multinomial model
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

# Interaction model for individual tributaries
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

# Combine predictions and standardize river nomenclature
data_vertical <- bind_rows(
  plot_data_all %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop),
  plot_data_interaction %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop)
)

rename_river_func <- function(df, col_name) {
  df %>%
    mutate(!!sym(col_name) := case_when(
      !!sym(col_name) == "west brook"  ~ "WB",
      !!sym(col_name) == "wb jimmy"    ~ "OL",
      !!sym(col_name) == "wb mitchell" ~ "OS",
      !!sym(col_name) == "wb obear"    ~ "IL",
      TRUE                             ~ as.character(!!sym(col_name))
    )) %>%
    filter(!!sym(col_name) %in% river_levels) %>%
    mutate(!!sym(col_name) := factor(!!sym(col_name), levels = river_levels))
}

data_vertical <- rename_river_func(data_vertical, "river_grp")
simpson_data  <- rename_river_func(simpson_data, "river_grp")

n_labels <- df_model %>%
  rename_river_func("river_grp") %>%
  group_by(river_grp) %>%
  summarise(total_n = n_distinct(tag), .groups = "drop") %>%
  mutate(label = paste0("n = ", total_n))

# =============================================================================
# 3. Right Plot Data: Population Dynamics (Stacked Area)
# =============================================================================
df2 <- read_csv(file_path_2, show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  rename(cluster_id = all_of(TARGET_CLUSTER_COL)) %>%
  mutate(year = year(date), season = quarter(date, with_year = FALSE)) %>%
  filter(!is.na(cluster_id))

stay_data <- df2 %>%
  arrange(tag, date) %>% group_by(tag) %>%
  mutate(next_date = lead(date), end_date = if_else(!is.na(next_date), next_date, date)) %>%
  ungroup()

expand_season <- stay_data %>%
  rowwise() %>%
  mutate(season_seq = list(unique(paste0(year(seq(date, end_date, by = "day")), "-S", quarter(seq(date, end_date, by = "day")))))) %>%
  unnest(season_seq) %>% ungroup() %>%
  select(tag, river, cluster_id, season_seq) %>% distinct()

alive_by_season_filtered <- expand_season %>%
  group_by(river, cluster_id, season_seq) %>%
  summarise(n_alive = n_distinct(tag), .groups = "drop") %>%
  mutate(year = as.numeric(substr(season_seq, 1, 4))) %>%
  filter(
    (river == "west brook" & year >= 2000 & year <= 2015) |
      (river != "west brook" & year >= 2002 & year <= 2015)
  ) %>%
  filter(grepl("-S3$", season_seq))

all_rivers_data <- alive_by_season_filtered %>%
  filter(year >= 2002) %>% 
  group_by(cluster_id, season_seq, year) %>% 
  summarise(n_alive = sum(n_alive), .groups = "drop") %>%
  mutate(river = "All Rivers") 

alive_by_season_combined <- bind_rows(alive_by_season_filtered, all_rivers_data) %>%
  rename_river_func("river") %>%
  mutate(cluster_id = as.factor(cluster_id)) %>%
  complete(river, year = 2000:2015, cluster_id, fill = list(n_alive = 0)) %>%
  filter(
    (river == "WB" & year >= 2000) |
      (river != "WB" & year >= 2002)
  )

# =============================================================================
# 4. Plot Assembly & Patchwork Integration
# =============================================================================

# Left Plot: Cluster Proportions & Simpson Index
p_left <- ggplot(data_vertical, aes(x = cohort_numeric)) +
  geom_point(aes(y = proportion, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), size = 3, alpha = 0.6) +
  geom_line(aes(y = proportion, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), linetype = "dotted", alpha = 0.6, linewidth = 1) + 
  geom_line(aes(y = predicted_prop, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), linewidth = 2) +
  geom_line(data = simpson_data, aes(y = simpson_index), color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_point(data = simpson_data, aes(y = simpson_index), color = "black", shape = 17, size = 4) +
  facet_wrap(~ river_grp, ncol = 1) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_x_continuous(limits = c(min(years_wb), max(years_wb)), breaks = seq(2000, 2015, 2)) + 
  geom_text(data = n_labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = 1.2, vjust = 1.5, size = 6, fontface = "bold", color = "black") +
  labs(x = "Year", y = "Proportion / Simpson Index", color = "Cluster") +
  scale_color_brewer(palette = "Set1") +
  shared_theme

p_right <- ggplot(alive_by_season_combined, aes(x = year, y = n_alive, fill = cluster_id)) +
  geom_area(color = "black", linewidth = 0.4, alpha = 0.8) +   
  scale_fill_brewer(palette = "Set1", name = "Cluster") +      
  facet_wrap(~ river, scales = "free_y", ncol = 1) +           
  scale_x_continuous(breaks = seq(2000, 2015, by = 2), limits = c(min(years_wb), max(years_wb))) +
  labs(x = "Year", y = "Number of Individuals") +
  shared_theme

# Combine plots side-by-side (guides="collect" merges identical legends)
combined_plot <- p_left + p_right + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

# =============================================================================
# 5. Output to PDF
# =============================================================================
n_panels <- length(river_levels)
PANEL_HEIGHT <- 3.0  
total_pdf_height <- (PANEL_HEIGHT * n_panels) + 2

output_path <- file.path(SCRIPT_DIR, OUTPUT_FILENAME)
ggsave(filename = output_path, plot = combined_plot, width = 16, height = total_pdf_height)

cat("\n✅ Saved successfully: ", output_path, "\n")


