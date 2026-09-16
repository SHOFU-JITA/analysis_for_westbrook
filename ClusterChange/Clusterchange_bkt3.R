# =============================================================================
# 0. Load Packages
# =============================================================================
rm(list = ls(all.names = TRUE))

library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(scales)
library(broom)
library(purrr)
library(nnet)
library(emmeans)

# =============================================================================
# 1. Basic Settings and Data Loading
# =============================================================================
file_path <- "../Clustering/WB_bkt_clustered_rev2.csv"
years_wb <- 2000:2015     
years_others <- 2002:2015 
min_year <- min(years_wb) 

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}

if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}
df <- read.csv(file_path)

num_clusters <- 3
ref_cluster <- as.character(num_clusters) 

y_axis_max_proportion <- 1
base_pdf_name <- paste0("clusterchange_bktk", num_clusters, ".pdf")
pdf_file_path <- file.path(SCRIPT_DIR, base_pdf_name)

river_labeller <- as_labeller(c(
  "All Rivers"  = "All Rivers",
  "west brook"  = "WB",
  "wb jimmy"    = "OL",
  "wb mitchell" = "OS",
  "wb obear"    = "IL"
))

# =============================================================================
# 2. Data Preprocessing
# =============================================================================
cluster_col_name <- paste0("cluster_k", num_clusters)

df <- df %>%
  mutate(date = as.Date(date)) %>%
  rename(ocluster_rename_dynamic = !!sym(cluster_col_name)) %>%
  mutate(ocluster_rename_dynamic = relevel(factor(ocluster_rename_dynamic, levels = as.character(1:num_clusters)), ref = ref_cluster))

river_tags <- df %>%
  filter(cohort %in% years_wb) %>%
  distinct(tag, river)

multi_tags <- river_tags %>%
  count(tag) %>%
  filter(n > 1) %>%
  pull(tag)

df_model_base <- df %>%
  mutate(
    river_grp = if_else(tag %in% multi_tags, "Multi", river),
    cohort_numeric = as.numeric(as.character(cohort))
  )

df_model_all <- df_model_base %>% mutate(river_grp = "All Rivers")
df_model <- bind_rows(df_model_base, df_model_all)

df_model <- df_model %>%
  filter(
    (river_grp == "west brook" & cohort_numeric %in% years_wb) |
      (river_grp != "west brook" & cohort_numeric %in% years_others)
  ) %>%
  mutate(year_centered = cohort_numeric - min_year)

custom_order <- c("All Rivers", "west brook", "wb jimmy", "wb mitchell", "wb obear")
other_rivers <- sort(setdiff(unique(df_model$river_grp), c(custom_order, "Multi")))
all_river_grps <- c(custom_order, other_rivers, "Multi")
df_model$river_grp <- factor(df_model$river_grp, levels = all_river_grps)

# =============================================================================
# 3. Prepare Basic Data for Plotting (Observed Values)
# =============================================================================
proportion_data <- df_model %>%
  group_by(river_grp, cohort, cohort_numeric, year_centered, ocluster_rename_dynamic) %>%
  summarise(n_cluster_tags = n_distinct(tag), .groups = "drop") %>%
  group_by(river_grp, cohort) %>%
  mutate(n_total_tags = sum(n_cluster_tags)) %>%
  ungroup() %>%
  mutate(proportion = n_cluster_tags / n_total_tags)

plot_data_full <- proportion_data %>%
  complete(
    river_grp,
    cohort_numeric = years_wb, 
    ocluster_rename_dynamic,
    fill = list(proportion = 0)
  ) %>%
  mutate(
    year_centered = cohort_numeric - min_year,
    proportion = if_else(river_grp != "west brook" & cohort_numeric < 2002, NA_real_, proportion)
  )

# --- Calculate Simpson Diversity Index (Gini-Simpson: 1 - sum(p^2)) ---
simpson_data <- plot_data_full %>%
  filter(!is.na(proportion)) %>% 
  group_by(river_grp, cohort_numeric) %>%
  summarise(simpson_index = 1 - sum(proportion^2, na.rm = TRUE), .groups = "drop") %>%
  filter(river_grp %in% custom_order)

simpson_data$river_grp <- factor(simpson_data$river_grp, levels = custom_order)

# --- Analysis 1: Independent Model for All Rivers ---
cat("\n--- Analysis 1: 'All Rivers' Independent Model ---\n")
multinom_results_all <- df_model %>%
  filter(river_grp == "All Rivers") %>%
  group_by(river_grp) %>%
  nest() %>%
  mutate(
    model = map(data, ~multinom(ocluster_rename_dynamic ~ year_centered, data = ., trace = FALSE)),
    predicted = map(model, ~{
      new_data <- tibble(cohort_numeric = years_others, year_centered = years_others - min_year)
      predict(., newdata = new_data, type = "probs") %>%
        as.data.frame() %>%
        mutate(cohort_numeric = new_data$cohort_numeric) %>%
        pivot_longer(cols = -cohort_numeric, names_to = "ocluster_rename_dynamic", values_to = "predicted_prop")
    })
  )

predicted_data_all <- multinom_results_all %>% select(river_grp, predicted) %>% unnest(predicted)
plot_data_all <- plot_data_full %>%
  filter(river_grp == "All Rivers") %>%
  left_join(predicted_data_all, by = c("river_grp", "cohort_numeric", "ocluster_rename_dynamic"))


# --- Analysis 2: Multinomial Logit Model with River Interaction ---
cat("\n--- Analysis 2: Multinomial Logit Model with River Interaction ---\n")
df_interaction <- df_model %>% filter(river_grp != "All Rivers")

full_model <- multinom(
  ocluster_rename_dynamic ~ year_centered * river_grp, 
  data = df_interaction, trace = FALSE, maxit = 1000, Hess = TRUE
)

new_data_interaction <- plot_data_full %>% 
  filter(river_grp != "All Rivers") %>% 
  filter(
    (river_grp == "west brook" & cohort_numeric %in% years_wb) |
      (river_grp != "west brook" & cohort_numeric %in% years_others)
  ) %>%
  distinct(cohort_numeric, year_centered, river_grp)

predicted_interaction <- as.data.frame(predict(full_model, newdata = new_data_interaction, type = "probs")) %>%
  mutate(cohort_numeric = new_data_interaction$cohort_numeric, river_grp = new_data_interaction$river_grp) %>%
  pivot_longer(cols = -c(cohort_numeric, river_grp), names_to = "ocluster_rename_dynamic", values_to = "predicted_prop")

plot_data_interaction <- plot_data_full %>%
  filter(river_grp != "All Rivers") %>%
  left_join(predicted_interaction, by = c("river_grp", "cohort_numeric", "ocluster_rename_dynamic"))


# --- Analysis 3: Extract Slopes from the Integrated Model ---
cat("\n--- Analysis 3: Extracting slopes ---\n")

df_interaction_clean <- df_model %>%
  filter(river_grp != "All Rivers") %>%
  select(ocluster_rename_dynamic, year_centered, river_grp) %>%
  mutate(
    river_grp = factor(as.character(river_grp)),
    ocluster_rename_dynamic = factor(as.character(ocluster_rename_dynamic))
  ) %>%
  droplevels()


df_interaction_clean$ocluster_rename_dynamic <- relevel(df_interaction_clean$ocluster_rename_dynamic, ref = ref_cluster)

full_model_clean <- multinom(
  ocluster_rename_dynamic ~ year_centered * river_grp, 
  data = df_interaction_clean, 
  trace = FALSE, 
  Hess = TRUE, 
  maxit = 1000
)

slopes_ems <- emtrends(
  full_model_clean, 
  ~ river_grp | ocluster_rename_dynamic, 
  var = "year_centered", 
  mode = "latent",              
  data = df_interaction_clean 
)

river_slopes_global <- as.data.frame(summary(slopes_ems, infer = c(TRUE, TRUE))) %>%
  rename(
    estimate = matches("year_centered.trend"),
    y.level = matches("ocluster_rename_dynamic"),
    conf.low = matches("lower.CL|asymp.LCL|LCL"),
    conf.high = matches("upper.CL|asymp.UCL|UCL")
  ) %>%
  filter(y.level != ref_cluster) 

cat("✅ Success: Extracted slopes.\n")

cat(paste0("\nSaving analysis results to '", pdf_file_path, "'...\n"))
pdf(pdf_file_path, width = 10, height = 20)

# =============================================================================
# Analysis 4: Vertical Comparison Plot
# =============================================================================
cat("\n--- Analysis 4: Vertical Comparison Plot ---\n")

data_vertical <- bind_rows(
  plot_data_all %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop),
  plot_data_interaction %>% filter(river_grp %in% custom_order) %>% select(river_grp, cohort_numeric, ocluster_rename_dynamic, proportion, predicted_prop)
)
data_vertical$river_grp <- factor(data_vertical$river_grp, levels = custom_order)

n_labels <- df_model %>%
  filter(river_grp %in% custom_order) %>%
  group_by(river_grp) %>%
  summarise(total_n = n_distinct(tag), .groups = "drop") %>%
  mutate(label = paste0("n = ", total_n))

p_vertical_stack <- ggplot(data_vertical, aes(x = cohort_numeric)) +
  geom_point(aes(y = proportion, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), size = 3, alpha = 0.6) +
  geom_line(aes(y = proportion, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), linetype = "dotted", alpha = 0.6, linewidth = 1) + 
  geom_line(aes(y = predicted_prop, color = ocluster_rename_dynamic, group = ocluster_rename_dynamic), linewidth = 2) +
  
  geom_line(data = simpson_data, aes(y = simpson_index), 
            color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_point(data = simpson_data, aes(y = simpson_index), 
             color = "black", shape = 17, size = 4) +
  
  facet_wrap(~ river_grp, ncol = 1, labeller = river_labeller) +
  scale_y_continuous(labels = percent_format(), limits = c(0, y_axis_max_proportion)) +
  scale_x_continuous(limits = c(min(years_wb), max(years_wb))) + 
  geom_text(
    data = n_labels,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.5, vjust = 2.0, 
    size = 7,                  
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "Cohort Year", 
    y = "Proportion / Simpson Index", 
    color = "Cluster",
  ) +
  theme_minimal(base_size = 24) + 
  theme(
    text = element_text(color = "black", face = "bold"),         
    axis.text = element_text(color = "black"),    
    strip.text = element_text(color = "black", face = "bold") 
  ) +
  scale_color_brewer(palette = "Set1")

print(p_vertical_stack)

dev.off()

# =============================================================================
# Data Preparation for Slope Comparison PDF (Horizontal)
# =============================================================================
model_all_rivers <- multinom_results_all$model[[1]]
data_all_rivers <- df_model %>% 
  filter(river_grp == "All Rivers") %>% 
  mutate(ocluster_rename_dynamic = relevel(factor(as.character(ocluster_rename_dynamic)), ref = ref_cluster)) %>% 
  droplevels()

slopes_all_ems <- emtrends(
  model_all_rivers, 
  ~ 1 | ocluster_rename_dynamic, 
  var = "year_centered", 
  mode = "latent",              
  data = data_all_rivers
)

slope_all_df <- as.data.frame(summary(slopes_all_ems, infer = c(TRUE, TRUE))) %>%
  rename(
    estimate = matches("year_centered.trend"),
    y.level = matches("ocluster_rename_dynamic"),
    conf.low = matches("lower.CL|asymp.LCL|LCL"),
    conf.high = matches("upper.CL|asymp.UCL|UCL")
  ) %>%
  filter(y.level != ref_cluster) %>% 
  mutate(river_grp = "All Rivers") %>%
  select(river_grp, y.level, estimate, conf.low, conf.high)

plot_slopes_data <- bind_rows(
  slope_all_df,
  river_slopes_global %>% 
    filter(river_grp != "Multi") %>% 
    select(river_grp, y.level, estimate, conf.low, conf.high)
)

plot_slopes_data$river_grp <- factor(plot_slopes_data$river_grp, levels = custom_order)

# =============================================================================
# Slope Comparison PDF Output (Horizontal)
# =============================================================================
pdf_slope_path <- file.path(SCRIPT_DIR, paste0("cluster_slope_bktk", num_clusters, ".pdf"))
pdf(pdf_slope_path, width = 15, height = 4)

p_slopes_horizontal <- ggplot(plot_slopes_data, aes(x = y.level, y = estimate, color = y.level)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.3, linewidth = 0.8) +
  geom_point(size = 4) +
  facet_wrap(~ river_grp, nrow = 1, labeller = river_labeller) + 
  coord_cartesian(ylim = c(-0.25, 0.25)) + 
  labs(title = paste0("Annual Change Comparison (Ref: Cluster ", ref_cluster, ")"), x = "Cluster", y = "Slope") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p_slopes_horizontal)
dev.off()

cat(paste0("\n✅ Saved successfully: ", pdf_file_path, "\n✅ Saved successfully: ", pdf_slope_path, "\n"))

