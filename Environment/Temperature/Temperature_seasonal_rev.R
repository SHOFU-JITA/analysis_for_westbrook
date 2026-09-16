# =============================================================================
# Westbrook Project: Temperature Analysis & Stats
# Output: 
#   1. Graph: seasonal_temperature_change_detailed.pdf
#   2. Graph: temperature_by_season_panel.pdf (B&W compatible)
#   3. Stats: LM/GLS Models, ANOVA, Emmeans, AIC Grouping (CSV outputs)
# =============================================================================
                                
# ===== 0. Reset Environment =====
rm(list = ls(all.names = TRUE))

# ===== 1. Load Packages =====
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(tidyr)
library(broom)   # Tidy regression output
library(car)     # Type III ANOVA
library(purrr)   # Map functions
library(nlme)    # GLS models
library(emmeans) # Post-hoc tests
library(tibble)  # Rownames to columns

# ===== 2. Directory and Data Setup =====
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}
message(paste("Save directory:", SCRIPT_DIR))

csv_file_path <- "../../rawdata/envDataWB.csv" 
if (!file.exists(csv_file_path)) stop("Data file not found. Please check the path.")
df <- read_csv(csv_file_path)

# ===== 3. Data Preprocessing =====
if (!"temperature" %in% names(df)) stop("'temperature' column not found in the dataset.")

df_filtered <- df %>%
  mutate(date = as.Date(date), year = year(date)) %>%
  filter(year >= 2000 & year <= 2015) %>%
  # Rename rivers to abbreviations
  mutate(
    river = case_when(
      grepl("(?i)west brook", river) ~ "WB",
      grepl("(?i)wb mitchel", river) ~ "OS",
      grepl("(?i)wb jimmy", river) ~ "OL",
      grepl("(?i)wb obear", river) ~ "IL",
      TRUE ~ river
    ),
    river = factor(river, levels = c("WB", "OS", "OL", "IL"))
  ) %>%
  mutate(
    month_num = month(date),
    season = case_when(
      month_num %in% 3:5 ~ 1,
      month_num %in% 6:8 ~ 2,
      month_num %in% 9:11 ~ 3,
      TRUE ~ 4
    ),
    season_label = factor(season, levels = c(1, 2, 3, 4),
                          labels = c("Spring", "Summer", "Autumn", "Winter"))
  )

# ===== 4. Aggregate Summary Stats =====
seasonal_summary_stats <- df_filtered %>%
  filter(!is.na(temperature)) %>%
  group_by(river, year, season_label) %>%
  summarise(
    avg_temp = mean(temperature, na.rm = TRUE),
    n_obs = n(),
    sd_temp = sd(temperature, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    se_temp = ifelse(n_obs > 1, sd_temp / sqrt(n_obs), NA_real_),
    ci_margin = ifelse(n_obs > 1, qt(0.975, df = n_obs - 1) * se_temp, NA_real_),
    ci_lower = avg_temp - ci_margin,
    ci_upper = avg_temp + ci_margin
  )

# =============================================================================
# 5. Graph Creation
# =============================================================================

# --- Graph 1: Seasonal Average Temperature Change (Detailed) ---
if(nrow(seasonal_summary_stats) > 0) {
  message("Creating Plot 1: seasonal_temperature_change_detailed.pdf ...")
  p1 <- ggplot(seasonal_summary_stats, aes(x = year, y = avg_temp, color = season_label)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.5, na.rm = TRUE) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, na.rm = TRUE) +
    facet_grid(rows = vars(season_label), cols = vars(river), scales = "free_y") +
    labs(title = "Seasonal Average Temperature Change", x = "Year", y = "Avg Temp (°C)") +
    theme_bw()
  ggsave(file.path(SCRIPT_DIR, "seasonal_temperature_change_detailed.pdf"), plot = p1, width = 14, height = 10)
}

# --- Graph 2: Temperature by Season Panel (B&W Compatible, No Borders) ---
if(nrow(seasonal_summary_stats) > 0) {
  message("Creating Plot 2: temperature_by_season_panel.pdf ...")
  
  p3 <- ggplot(seasonal_summary_stats, aes(x = year, y = avg_temp, group = river, 
                                           color = river, shape = river, linetype = river)) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.3, size = 0.6, alpha = 0.7, linetype = "solid", show.legend = FALSE) +
    geom_point(size = 4, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    facet_wrap(~ season_label, scales = "free_y", ncol = 2) +
    labs(title = "Annual Average Temperature Trends by Season", 
         x = "Year", 
         y = "Average Temperature (°C)",
         color = "River", shape = "River", linetype = "River") +
    scale_color_grey(start = 0, end = 0.6) +
    theme_bw(base_size = 30) + 
    theme(
      panel.border = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = 30),     
      axis.title = element_text(face = "bold", size = 30),     
      axis.text = element_text(size = 30, color = "black"),    
      legend.title = element_text(face = "bold", size = 30),   
      legend.text = element_text(size = 30),                   
      legend.position = "bottom",                              
      legend.key.width = unit(2.5, "cm")
    )
  
  ggsave(file.path(SCRIPT_DIR, "temperature_by_season_panel.pdf"), plot = p3, width = 16, height = 12)
}


# =============================================================================
# 10. Statistical Model Analysis & CSV Exports
# =============================================================================

# -------------------------------------------------------
# 10.1 Linear Model (lm) & Type III ANOVA
# -------------------------------------------------------
message("--- LM Type III ANOVA Running ---")

season_models <- seasonal_summary_stats %>%
  split(.$season_label) %>%
  map(~ lm(avg_temp ~ year * river, data = ., weights = n_obs))

season_results <- map(season_models, ~ Anova(.x, type = "III"))

lm_anova_df <- map2_dfr(season_results, names(season_results), function(dat, name) {
  dat %>%
    as.data.frame() %>%
    rownames_to_column(var = "Term") %>%
    mutate(Season = name) %>%
    select(Season, Term, everything())
})

file_lm <- file.path(SCRIPT_DIR, "stats_LM_ANOVA_results.csv")
write_csv(lm_anova_df, file_lm)
message(paste("Saved LM ANOVA results to:", file_lm))


# -------------------------------------------------------
# 10.2 Generalized Least Squares (GLS)
# -------------------------------------------------------
message("--- GLS Running ---")

gls_season_models <- seasonal_summary_stats %>%
  split(.$season_label) %>%
  map(~ tryCatch({
    gls(
      avg_temp ~ year * river, 
      data = ., 
      weights = varIdent(form = ~1 | river),
      correlation = corAR1(form = ~ year | river),
      na.action = na.omit
    )
  }, error = function(e) {
    message("Error in model fitting: ", e$message)
    return(NULL)
  })) %>%
  compact() # Remove seasons that resulted in errors (NULL)

gls_anova_results <- map(gls_season_models, ~ anova(.))

gls_anova_df <- map2_dfr(gls_anova_results, names(gls_anova_results), function(dat, name) {
  dat %>%
    as.data.frame() %>%
    rownames_to_column(var = "Term") %>%
    mutate(Season = name) %>%
    select(Season, Term, everything())
})

file_gls <- file.path(SCRIPT_DIR, "stats_GLS_ANOVA_results.csv")
write_csv(gls_anova_df, file_gls)
message(paste("Saved GLS ANOVA results to:", file_gls))


# -------------------------------------------------------
# 10.3 Emmeans Analysis (GLS Model)
# -------------------------------------------------------
message("--- Emmeans Analysis Running (GLS Model) ---")

data_list <- seasonal_summary_stats %>%
  split(.$season_label)

valid_seasons <- names(gls_season_models)
data_list_valid <- data_list[valid_seasons]

# Run emtrends (avoid data reference errors)
em_trends_gls <- map2(gls_season_models, data_list_valid, function(model, data) {
  model$call$data <- quote(data)
  emtrends(model, specs = "river", var = "year", data = data)
})

# (A) Slope Estimates (GLS)
emmeans_slopes_gls_df <- map2_dfr(em_trends_gls, names(em_trends_gls), function(emm, name) {
  summary(emm) %>%
    as.data.frame() %>%
    mutate(Season = name) %>%
    select(Season, everything())
})

file_em_slopes_gls <- file.path(SCRIPT_DIR, "stats_emmeans_slopes_GLS.csv")
write_csv(emmeans_slopes_gls_df, file_em_slopes_gls)
message(paste("Saved GLS slope estimates to:", file_em_slopes_gls))

# (B) Pairwise Comparisons (GLS)
emmeans_pairs_gls_df <- map2_dfr(em_trends_gls, names(em_trends_gls), function(emm, name) {
  pairs(emm) %>%
    as.data.frame() %>%
    mutate(Season = name) %>%
    select(Season, everything())
})

file_em_pairs_gls <- file.path(SCRIPT_DIR, "stats_emmeans_pairwise_GLS.csv")
write_csv(emmeans_pairs_gls_df, file_em_pairs_gls)
message(paste("Saved GLS pairwise comparisons to:", file_em_pairs_gls))


# ===== 11. AIC Slope Grouping Analysis =====
message("--- AIC Slope Grouping Analysis Running ---")

analyze_slope_grouping <- function(target_season_name, df_all) {
  
  message(paste("Processing:", target_season_name))
  
  season_data <- df_all %>% filter(season_label == target_season_name)
  
  if(nrow(season_data) == 0) {
    warning(paste("No data for", target_season_name, "- skipping."))
    return(NULL)
  }
  
  river_list <- unique(season_data$river)
  
  target_combinations <- c(
    combn(river_list, 1, simplify = FALSE),
    combn(river_list, 2, simplify = FALSE)
  )
  
  # 1. Calculate Grouped Models (Combinations)
  aic_grouped_results <- map_dfr(target_combinations, function(targets) {
    target_name <- paste(targets, collapse = "+")
    
    # Target group (Group A) vs Others (Group B)
    temp_data <- season_data %>%
      mutate(river_group = if_else(river %in% targets, "Group A", "Group B"))
    
    # Grouped GLS model (ML method for AIC comparison)
    mod <- tryCatch({
      gls(
        avg_temp ~ river + year : river_group + year, # Intercept per river, slope per group
        data = temp_data,
        weights = varIdent(form = ~1 | river),
        correlation = corAR1(form = ~ year | river),
        na.action = na.omit,
        method = "ML" 
      )
    }, error = function(e) NULL)
    
    if (!is.null(mod)) {
      return(data.frame(
        Season = target_season_name,
        Grouping = target_name,
        Type = "Grouped Model",
        AIC = AIC(mod),
        LogLik = as.numeric(logLik(mod))
      ))
    } else {
      return(NULL)
    }
  })
  
  # 2. Baseline Models (Separate & Unified)
  # (A) Separate: Different slopes for each river
  mod_separate <- tryCatch({
    gls(avg_temp ~ year * river, data = season_data,
        weights = varIdent(form = ~1 | river), correlation = corAR1(form = ~ year | river),
        na.action = na.omit, method = "ML")
  }, error = function(e) NULL)
  
  # (B) Unified: Same slope for all rivers
  mod_unified <- tryCatch({
    gls(avg_temp ~ year + river, data = season_data,
        weights = varIdent(form = ~1 | river), correlation = corAR1(form = ~ year | river),
        na.action = na.omit, method = "ML")
  }, error = function(e) NULL)
  
  baseline_results <- data.frame(
    Season = target_season_name,
    Grouping = c("ALL_SEPARATE (4 slopes)", "ALL_UNIFIED (1 slope)"),
    Type = "Baseline",
    AIC = c(if(!is.null(mod_separate)) AIC(mod_separate) else NA, 
            if(!is.null(mod_unified)) AIC(mod_unified) else NA),
    LogLik = c(if(!is.null(mod_separate)) as.numeric(logLik(mod_separate)) else NA, 
               if(!is.null(mod_unified)) as.numeric(logLik(mod_unified)) else NA)
  )
  
  # Combine and sort by AIC
  final_ranking <- bind_rows(aic_grouped_results, baseline_results) %>%
    filter(!is.na(AIC)) %>%
    arrange(AIC) %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    select(Season, Grouping, Type, AIC, delta_AIC, LogLik)
  
  return(final_ranking)
}

# Run analysis for all seasons
seasons_to_analyze <- c("Spring", "Summer", "Autumn", "Winter")

aic_analysis_results <- map_dfr(seasons_to_analyze, ~ analyze_slope_grouping(.x, seasonal_summary_stats))

file_aic <- file.path(SCRIPT_DIR, "stats_AIC_slope_grouping.csv")
write_csv(aic_analysis_results, file_aic)

message("AIC Analysis Results:")
print(aic_analysis_results)
message(paste("Saved AIC grouping analysis results to:", file_aic))

message("===== All analyses completed =====")

