# ===================================================================
# Westbrook Portfolio Effect, Synchrony & Temporal Change Analysis
# ===================================================================
rm(list = ls(all.names = TRUE))

# ===== 0. User Configuration =====
CLUSTER_NUM <- 3
WINDOW_SIZE <- 7

TARGET_CLUSTER_COL <- paste0("cluster_k", CLUSTER_NUM)
CLUSTER_IDS <- c("1", "2", "3") # Adjust according to CLUSTER_NUM

# Colors
CLUSTER_COLORS <- c("steelblue", "aquamarine3", "salmon3")
OUTPUT_PDF <- paste0("Westbrook_", WINDOW_SIZE, "yr_window_k", CLUSTER_NUM, ".pdf")
INPUT_FILE <- "../csv/WB_combined_bkt_rev2.csv"

# ===== 1. Package Loading =====
packages <- c("readr", "dplyr", "lubridate", "tidyr", "ggplot2")
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
    if (current_path != "") setwd(dirname(current_path))
  }, silent = TRUE)
}
message("Working directory: ", getwd())

# ===== 3. Data Loading and Preprocessing =====
data <- read_csv(INPUT_FILE, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

if (!TARGET_CLUSTER_COL %in% names(data)) {
  stop(paste("Error: Column", TARGET_CLUSTER_COL, "not found in the data."))
}

# Rename and standardize columns
data <- data %>%
  rename(cluster_id = all_of(TARGET_CLUSTER_COL)) %>%
  mutate(
    year = year(date),
    river = case_when(
      river == "west brook"  ~ "WB",
      river == "wb jimmy"    ~ "OL",
      river == "wb mitchell" ~ "OS",
      river == "wb obear"    ~ "IL",
      TRUE                   ~ river
    )
  ) %>%
  filter(!is.na(cluster_id))

# Expand residence period to seasonal units and filter S3 (2002-2015)
alive_by_season_filtered <- data %>%
  arrange(tag, date) %>%
  group_by(tag) %>%
  mutate(
    next_date = lead(date),
    end_date = if_else(!is.na(next_date), next_date, date)
  ) %>%
  rowwise() %>%
  mutate(season_seq = list({
    dseq <- seq(date, end_date, by = "day")
    unique(paste0(year(dseq), "-S", quarter(dseq)))
  })) %>%
  unnest(season_seq) %>%
  ungroup() %>%
  select(tag, river, cluster_id, season_seq) %>%
  distinct() %>%
  group_by(river, cluster_id, season_seq) %>%
  summarise(n_alive = n_distinct(tag), .groups = "drop") %>%
  mutate(year = as.numeric(substr(season_seq, 1, 4))) %>%
  filter(year >= 2002 & year <= 2015) %>%
  filter(grepl("-S3$", season_seq))

# ===== 4. Define Sliding Windows =====
min_year <- min(alive_by_season_filtered$year)
max_year <- max(alive_by_season_filtered$year)

start_years <- min_year:(max_year - WINDOW_SIZE + 1)
windows <- lapply(start_years, function(y) y:(y + WINDOW_SIZE - 1))
names(windows) <- paste0(start_years, "-", start_years + WINDOW_SIZE - 1)

cat(sprintf("\nDefined %d rolling windows of size %d years.\n", length(windows), WINDOW_SIZE))

# ===== 5. Base calculation for Synchrony Analysis =====
sync_base <- alive_by_season_filtered %>%
  group_by(cluster_id, year) %>%
  summarise(n_alive = sum(n_alive), .groups = "drop")

r_one_vs_others_all <- bind_rows(lapply(CLUSTER_IDS, function(cid) {
  sync_base %>%
    mutate(is_target = if_else(cluster_id == cid, "Target", "Others")) %>%
    group_by(is_target, year) %>%
    summarise(n_alive = sum(n_alive), .groups = "drop") %>%
    arrange(is_target, year) %>%
    group_by(is_target) %>%
    mutate(
      n_prev = lag(n_alive, order_by = year),
      year_diff = year - lag(year, order_by = year),
      r_rate = if_else(year_diff == 1 & n_prev > 0 & n_alive > 0, log(n_alive / n_prev), NA_real_)
    ) %>%
    filter(!is.na(r_rate)) %>%
    mutate(target_cluster = cid) %>%
    ungroup()
}))

# ===== 6. Master Loop for All Analyses =====
contribution_results <- list()
sync_results <- list()
pe_results <- list()

for (w_name in names(windows)) {
  target_years <- windows[[w_name]]
  window_data <- alive_by_season_filtered %>% filter(year %in% target_years)
  
  # Skip if window data has less than 2 seasons (needed for sd)
  if (n_distinct(window_data$season_seq) < 2) next
  
  # ---------------------------------------------------------
  # Analysis A: Portfolio Effect (Overall, River, Cluster)
  # ---------------------------------------------------------
  # CV by River x Cluster
  cv_rc <- window_data %>%
    group_by(river, cluster_id) %>%
    summarise(m = mean(n_alive), cv = sd(n_alive) / mean(n_alive), .groups = 'drop')
  
  # CV by River
  cv_r <- window_data %>%
    group_by(river, season_seq) %>%
    summarise(n = sum(n_alive), .groups = "drop") %>%
    group_by(river) %>%
    summarise(m = mean(n), cv = sd(n) / mean(n), .groups = 'drop')
  
  # CV by Cluster
  cv_c <- window_data %>%
    group_by(cluster_id, season_seq) %>%
    summarise(n = sum(n_alive), .groups = "drop") %>%
    group_by(cluster_id) %>%
    summarise(m = mean(n), cv = sd(n) / mean(n), .groups = 'drop')
  
  # Total CV
  cv_t <- window_data %>%
    group_by(season_seq) %>%
    summarise(tot = sum(n_alive), .groups = 'drop') %>%
    summarise(cv = sd(tot) / mean(tot)) %>% pull(cv)
  
  wm_rc <- weighted.mean(cv_rc$cv, w = cv_rc$m, na.rm = TRUE)
  wm_r  <- weighted.mean(cv_r$cv, w = cv_r$m, na.rm = TRUE)
  wm_c  <- weighted.mean(cv_c$cv, w = cv_c$m, na.rm = TRUE)
  
  pe_results[[w_name]] <- data.frame(
    Window = w_name,
    `Cluster Effect` = wm_rc - wm_r,
    `River Effect`   = wm_rc - wm_c,
    `Overall Effect` = wm_rc - cv_t,
    check.names = FALSE
  )
  
  # ---------------------------------------------------------
  # Analysis B: Cluster Contribution Effect
  # ---------------------------------------------------------
  for (cid in CLUSTER_IDS) {
    temp_cv <- window_data %>%
      mutate(merged_cluster = if_else(cluster_id == cid, as.character(cid), "Others")) %>%
      group_by(river, merged_cluster, season_seq) %>%
      summarise(n = sum(n_alive), .groups = "drop") %>%
      group_by(river, merged_cluster) %>%
      summarise(m = mean(n), cv = sd(n) / mean(n), .groups = "drop")
    
    wm_temp <- weighted.mean(temp_cv$cv, w = temp_cv$m, na.rm = TRUE)
    
    contribution_results[[paste(w_name, cid)]] <- data.frame(
      Window = w_name, cluster_id = cid, contribution_effect = wm_temp - wm_r
    )
  }
  
  # ---------------------------------------------------------
  # Analysis C: Synchrony
  # ---------------------------------------------------------
  window_r_data <- r_one_vs_others_all %>% filter(year %in% target_years)
  for (cid in CLUSTER_IDS) {
    df_wide <- window_r_data %>% filter(target_cluster == cid) %>%
      pivot_wider(id_cols = year, names_from = is_target, values_from = r_rate)
    
    r_val <- NA_real_; p_val <- NA_real_
    if (all(c("Target", "Others") %in% names(df_wide))) {
      valid_data <- df_wide %>% filter(is.finite(Target) & is.finite(Others))
      if (nrow(valid_data) >= 3) {
        test_res <- cor.test(valid_data$Target, valid_data$Others, method = "pearson")
        r_val <- test_res$estimate
        p_val <- test_res$p.value
      }
    }
    
    sig_mark <- case_when(
      is.na(p_val) ~ "", p_val < 0.001 ~ "***", p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*", p_val < 0.1 ~ ".", TRUE ~ ""
    )
    
    sync_results[[paste(w_name, cid)]] <- data.frame(
      Window = w_name, cluster_id = cid, correlation = r_val, p_value = p_val, sig = sig_mark
    )
  }
}

# Bind results
df_contrib <- bind_rows(contribution_results) %>% mutate(Window = factor(Window, levels = names(windows)))
df_sync    <- bind_rows(sync_results) %>% mutate(Window = factor(Window, levels = names(windows)))
df_pe      <- bind_rows(pe_results) %>%
  mutate(Window = factor(Window, levels = names(windows))) %>%
  pivot_longer(cols = c("Cluster Effect", "River Effect", "Overall Effect"), 
               names_to = "effect_type", values_to = "pe_value")

# ===== 7. Plotting (Unified Theme) =====
custom_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(linewidth = 1, color = "black"), 
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Plot 1: Contribution Effect
p_contrib <- ggplot(df_contrib, aes(x = Window, y = contribution_effect, group = cluster_id, color = as.factor(cluster_id))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  scale_color_manual(values = CLUSTER_COLORS) +
  custom_theme +
  labs(title = "Decomposition of Portfolio Effect (Contribution)", x = paste0(WINDOW_SIZE, "-Year Window"), y = "Contribution Effect (CV Reduction)", color = "Cluster")

# Plot 2: Synchrony
p_sync <- ggplot(df_sync, aes(x = Window, y = correlation, group = cluster_id, color = as.factor(cluster_id))) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  geom_text(aes(label = sig, y = correlation + 0.05), show.legend = FALSE, size = 6, color = "black") +
  scale_color_manual(values = CLUSTER_COLORS) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.0, by = 0.2)) +
  custom_theme +
  labs(title = "One vs Others Synchrony (Pearson's r)", x = paste0(WINDOW_SIZE, "-Year Window"), y = "Pearson's Correlation Coefficient (r)", color = "Cluster")

# Plot 3: Temporal Change in Portfolio Effects (Black & White with Shapes)
p_pe <- ggplot(df_pe, aes(x = Window, y = pe_value, group = effect_type, shape = effect_type, linetype = effect_type)) +
  geom_line(linewidth = 1.2, color = "black") + 
  geom_point(size = 4, color = "black", fill = "black", stroke = 1.5) +
  scale_shape_manual(values = c("Cluster Effect" = 22,   
                                "River Effect" = 21,     
                                "Overall Effect" = 24)) + 
  scale_linetype_manual(values = c("Cluster Effect" = "solid", 
                                   "River Effect" = "solid", 
                                   "Overall Effect" = "solid")) +
  custom_theme +
  labs(
    title = "Temporal Change in Portfolio Effects", 
    x = paste0(WINDOW_SIZE, "-Year Window"), 
    y = "Portfolio Effect (CV Reduction)", 
    shape = "Effect Type",
    linetype = "Effect Type"
  )

# ===== 8. Save Plots to a Single PDF =====
cat("\nSaving plots to:", OUTPUT_PDF, "\n")
pdf(file = OUTPUT_PDF, width = 10, height = 6, onefile = TRUE)
print(p_contrib)
print(p_sync)
print(p_pe)
invisible(dev.off())

cat("All Analysis Completed Successfully.\n")

