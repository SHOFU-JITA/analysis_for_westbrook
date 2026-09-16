# reset environment
rm(list=ls(all.names=T))

# load packages
library(dplyr)
library(ggplot2)
library(flexmix)
library(RColorBrewer)
library(readr)

# plot setting
theme_set(
  theme_bw() +
    theme(title = element_text(size=20),
          axis.title = element_text(size=20),
          axis.text = element_text(size=15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.position = 'right',
          plot.caption = element_text(size = 12))
)

# read and format data
# ファイルパスは環境に合わせて適宜修正してください
# read and format data
# ファイルパスは環境に合わせて適宜修正してください
# data <- read_csv("F:/藤田バックアップ/運用/csv/WB_all_ats.csv")
data <- read_csv("~/Westbrook_カワマスP/code_SF/csv/WB_electro.csv", show_col_types = FALSE)
#data <- read_csv("/media/sato-lab/HDD_Noda/code_SF/csv/WB_electro.csv", show_col_types = FALSE)
data$date <- as.Date(data$date)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  SCRIPT_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  SCRIPT_DIR <- getwd()
}

data <- data %>%
  filter(!is.na(observedLength)) %>%
  filter(species == "bnt") %>%
  mutate(id = tag,
         fl = observedLength,
         hatching_date = as.Date(paste0(cohort, "-01-01")),
         days_from_hatching = as.numeric(difftime(date, hatching_date, units = "days")))

data <- as.data.frame(data)
data <- data %>%
  filter(!is.na(days_from_hatching))

# エラーデータの排除（成長が著しく減少した個体）
tag_decrease <- data %>%
  arrange(tag, detectionDate) %>%
  group_by(tag) %>%
  mutate(prev_length = lag(observedLength)) %>%
  filter(!is.na(prev_length) & (observedLength - prev_length) < -10) %>%
  distinct(tag) %>%
  pull(tag)

data <- data %>%
  filter(!tag %in% tag_decrease)

# Somers型von Bertalanffy成長関数
sbgm <- function(lmax, lmin, a1, theta, t0 = 0, t) {
  lmax - (lmax - lmin) * exp(-a1 * (t - t0) - 180 * a1 / pi *
                               (cos(pi * (t0 + theta) / 180) - cos(pi * (t + theta) / 180)))
}

# 負の対数尤度関数
llsbgm <- function(para, data) {
  lmax <- exp(para[1]); lmin <- exp(para[2]); a1 <- exp(para[3]); theta <- para[4]
  t <- data[, 1]; lobs <- data[, 2]
  pred_l <- sbgm(lmax, lmin, a1, theta, t = t)
  # 予測値が0以下やNAになる場合のエラー処理
  pred_l[pred_l <= 0 | is.na(pred_l)] <- 1e-6 
  sd_val <- sd(log(lobs), na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) sd_val <- 1e-6 # SDが0になるのを防ぐ
  -sum(dnorm(log(lobs), log(pred_l), sd = sd_val, log = TRUE))
}

llsbgm2 <- function(para, data) {
  lmax <- exp(para[1]); lmin <- exp(para[2]); a1 <- exp(para[3]); theta <- para[4]
  t <- data[, 1]; lobs <- data[, 2]
  pred_l <- sbgm(lmax, lmin, a1, theta, t = t)
  pred_l[pred_l <= 0 | is.na(pred_l)] <- 1e-6
  sd_val <- sd(log(lobs), na.rm = TRUE)
  if (is.na(sd_val) || sd_val == 0) sd_val <- 1e-6
  dnorm(log(lobs), log(pred_l), sd = sd_val, log = TRUE)
}

# 結果保持用
aic_bic <- data.frame(clu = c(), aic = c(), bic = c())
plot_list <- list()

# ★★★ パラメータ保存用のデータフレームを初期化 ★★★
parameter_results <- data.frame()

# 全クラスタで共通の処理
t <- seq(0, max(data$days_from_hatching, na.rm = TRUE), by = 1)
uid <- sort(unique(data$id))

# 全データで初期パラメータを推定
ipara_all <- c(log(max(data$fl)), log(min(data$fl)), log(1.7E-4), 0)
res_all <- optim(ipara_all, llsbgm, data = data[, c("days_from_hatching", "fl")], method = "BFGS")
data$predfl <- sbgm(exp(res_all$par[1]), exp(res_all$par[2]), exp(res_all$par[3]), res_all$par[4], t = data$days_from_hatching)
data$resid <- log(data$fl) - log(data$predfl)

# 全クラスタリング結果を保存するためのデータフレームを初期化
cluster_assignments_df <- data.frame(id = uid)


# ==============================================================================================
#  クラスタリングループ
# ==============================================================================================
for (kk in 1:10) {
  
  message(sprintf("\n\n=== クラスタ数: %d の処理を開始 ===", kk))
  
  # ---------- 1. 初期クラスタ分類 ----------
  m1 <- FLXMRglm(resid ~ 1, family = "gaussian")
  set.seed(123)
  # flexmixの反復回数を増やし、より安定した初期分類を得る
  mp2 <- flexmix(. ~ . | id, data = data, model = list(m1), k = kk,
                 control = list(iter.max = 300, minprior = 0.01, nrep = 10, verbose = 0))
  
  data$cluster <- clusters(mp2)
  k <- sort(unique(data$cluster))
  
  # クラスタ数がkkに満たない場合は次のループへ
  if (length(k) < kk) {
    message(sprintf(">>> 生成されたクラスタ数が%d未満のためスキップします。", kk))
    plot_list[[kk]] <- NULL 
    
    col_name <- paste0("cluster_k", kk)
    cluster_assignments_df[[col_name]] <- NA
    next
  }
  
  # ---------- 2. EMアルゴリズムによるパラメータとクラスタの交互最適化 ----------
  
  # パラメータの初期値を設定 (各クラスタに全データでの推定値を割り当て)
  estpara <- t(replicate(kk, res_all$par))
  
  max_iter <- 100
  loglik_hist <- numeric()
  
  set.seed(123)
  for (iter in 1:max_iter) {
    
    # E-step: パラメータを固定し、クラスタを再割り当て
    if (iter > 1) {
      for (cid in uid) {
        dmy <- data[data$id == cid, ]
        if (nrow(dmy) > 0) {
          cll <- sapply(k, function(cl) {
            loglik_col <- paste0("ll-", cl)
            if (loglik_col %in% colnames(data)) {
              sum(data[data$id == cid, loglik_col])
            } else {
              Inf 
            }
          })
          if (all(is.infinite(cll))) {
            data[data$id == cid, "cluster"] <- prev_cluster_assignment[data$id == cid][1]
          } else {
            data[data$id == cid, "cluster"] <- which.min(cll)
          }
        }
      }
    }
    
    # M-step: クラスタを固定し、パラメータを再推定
    prev_estpara <- estpara
    current_loglik <- numeric(length(k))
    names(current_loglik) <- k
    
    for (cl in k) {
      cluster_data <- data[data$cluster == cl, c("days_from_hatching", "fl")]
      
      if (nrow(cluster_data) < 10) {
        message(sprintf("[iter %d, cl %d] データ不足のためパラメータ更新をスキップ", iter, cl))
        estpara[which(k==cl), ] <- prev_estpara[which(k==cl), ]
        current_loglik[as.character(cl)] <- llsbgm(prev_estpara[which(k==cl), ], cluster_data)
        next
      }
      
      res <- tryCatch({
        optim(prev_estpara[which(k==cl), ], llsbgm, data = cluster_data,
              method = "BFGS", control = list(maxit = 2000))
      }, error = function(e) {
        message(sprintf("[iter %d, cl %d] optimでエラー発生: %s", iter, cl, e$message))
        return(NULL)
      })
      
      if (!is.null(res) && res$convergence == 0) {
        estpara[which(k==cl), ] <- res$par
        current_loglik[as.character(cl)] <- res$value
      } else {
        estpara[which(k==cl), ] <- prev_estpara[which(k==cl), ]
        current_loglik[as.character(cl)] <- llsbgm(prev_estpara[which(k==cl), ], cluster_data)
      }
      
      data[[paste0("ll-", cl)]] <- -llsbgm2(estpara[which(k==cl), ], data[, c("days_from_hatching", "fl")])
    }
    
    # ---------- 3. 収束判定 ----------
    loglik_hist[iter] <- sum(current_loglik)
    
    if (iter > 1) {
      loglik_diff_rate <- abs(loglik_hist[iter] - loglik_hist[iter-1]) / (abs(loglik_hist[iter-1]) + 1e-6)
      cluster_change_rate <- mean(data$cluster != prev_cluster_assignment)
      
      message(sprintf("[iter %d] LogLik: %.2f | Rel_diff: %.6f | Cluster_change: %.4f",
                      iter, loglik_hist[iter], loglik_diff_rate, cluster_change_rate))
      
      if (loglik_diff_rate < 1e-5 && cluster_change_rate < 0.01) {
        message(">>> 収束判定条件を満たしました。ループを終了します。")
        break
      }
    } else {
      message(sprintf("[iter %d] LogLik: %.2f", iter, loglik_hist[iter]))
    }
    prev_cluster_assignment <- data$cluster
  } # EMアルゴリズムのループ終了
  
  # ---------- 4. 結果の整理と保存 ----------
  data$ocluster <- data$cluster
  
  # クラスタをサイズ（最大体長）順にリネーム
  clu_max_len <- data %>%
    filter(days_from_hatching < 300) %>%
    group_by(ocluster) %>%
    summarise(max_len = max(fl, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_len)) %>%
    mutate(ocluster_rename = row_number())
  
  data <- data %>%
    select(-matches("ocluster_rename|max_len")) %>%
    left_join(clu_max_len, by = "ocluster")
  
  # このkkでのクラスタ分け結果を抽出し、全体の結果に結合する
  temp_result <- data %>%
    distinct(id, ocluster_rename) %>%
    rename(!!paste0("cluster_k", kk) := ocluster_rename)
  
  cluster_assignments_df <- cluster_assignments_df %>%
    left_join(temp_result, by = "id")
  
  # パラメータのマトリックスを整理（サイズ順に並び替え）
  estpara2 <- cbind(estpara, ocluster = k) %>%
    as.data.frame() %>%
    left_join(select(clu_max_len, ocluster, ocluster_rename), by = "ocluster") %>%
    arrange(ocluster_rename) %>%
    select(-starts_with("ocluster")) %>%
    as.matrix()
  
  # ★★★ 詳細統計量（N, Sigma）の計算と保存 ★★★
  
  # 各クラスタごとの統計量を計算するためのループ
  stats_list <- list()
  for(i in 1:nrow(estpara2)){
    # このクラスタ（リネーム後）に属するデータを抽出
    sub_dat <- data %>% filter(ocluster_rename == i)
    
    # このクラスタの予測値を計算（estpara2は1-3列が対数、4列が実数であることに注意）
    # パラメータは estpara2[i, ] を使用
    pred_vals <- sbgm(exp(estpara2[i, 1]), exp(estpara2[i, 2]), exp(estpara2[i, 3]), estpara2[i, 4], t = sub_dat$days_from_hatching)
    
    # 対数スケールでの残差標準偏差 (Sigma) を計算
    # ※尤度関数内で sd(log(lobs)) を計算しているのと合わせます
    resid_log <- log(sub_dat$fl) - log(pred_vals)
    sigma_val <- sd(resid_log, na.rm = TRUE)
    
    stats_list[[i]] <- data.frame(
      n = nrow(sub_dat),
      sigma_log = sigma_val
    )
  }
  stats_df <- do.call(rbind, stats_list)
  
  # パラメータ保存用データフレームの作成
  temp_params <- data.frame(
    num_clusters = kk,                 # その時のクラスタ総数 (k)
    cluster_id = 1:kk,                 # クラスタID (サイズ順)
    lmax = exp(estpara2[, 1]),         # 最大体長 (mm)
    lmin = exp(estpara2[, 2]),         # 最小体長 (mm)
    a1 = exp(estpara2[, 3]),           # 成長係数
    theta = estpara2[, 4],             # 位相パラメータ
    n = stats_df$n,                    # クラスタ内の個体数 (N)
    proportion = stats_df$n / nrow(data), # 全体に占める割合
    sigma_log = stats_df$sigma_log     # 対数スケールでのばらつき (Sigma)
  )
  
  parameter_results <- rbind(parameter_results, temp_params)
  # ★★★ 詳細統計量の保存完了 ★★★
  
  
  # AIC/BIC 計算
  n_params <- 5 * kk - 1 
  aic <- 2 * sum(current_loglik) + 2 * n_params
  bic <- 2 * sum(current_loglik) + n_params * log(nrow(data))
  aic_bic <- rbind(aic_bic, data.frame(clu = kk, aic = aic, bic = bic))
  
  # プロット作成
  cb <- brewer.pal(max(3, kk), "Set1")
  g <- ggplot()
  
  for (i in 1:nrow(estpara2)) {
    growth_curve <- tibble(
      t = t,
      fl = sbgm(exp(estpara2[i, 1]), exp(estpara2[i, 2]), exp(estpara2[i, 3]), exp(estpara2[i, 4]), t = t)
    ) %>% filter(t >= min(data$days_from_hatching))
    
    if("ocluster_rename" %in% names(data)){
      plot_data <- data %>% filter(ocluster_rename == i)
    } else {
      plot_data <- data %>% filter(FALSE) 
    }
    
    g <- g +
      geom_point(data = plot_data, aes(x = days_from_hatching, y = fl),
                 size = 2, alpha = 0.5, fill = cb[i], shape = 21, colour = "black") +
      geom_line(data = growth_curve, aes(x = t, y = fl),
                linewidth = 2, col = cb[i])
  }
  
  g <- g + labs(
    title = paste0(kk, " clusters (bkt)"), 
    x = "No. of days since January-1",
    y = "Fork length (mm)"
  )
  plot_list[[kk]] <- g
  
  message(sprintf("==== クラスタ数: %d 完了 ====", kk))
}

# ==============================================================================================
#  結果の出力
# ==============================================================================================

# 1. PDFへのプロット書き出し
pdf("all_growth_bnt_rev2.pdf", width = 8, height = 6)
for (i in seq_along(plot_list)) {
  if (!is.null(plot_list[[i]])) print(plot_list[[i]])
}
dev.off()

# 2. AIC/BICの表示とプロット
print(aic_bic)
write_csv(aic_bic, "AIC_BIC_bnt_rev2.csv")

g_aicbic <- ggplot(aic_bic) +
  geom_line(aes(x = clu, y = aic), linewidth = 2, color = "gray66") +
  geom_point(aes(x = clu, y = aic), size = 5, color = "gray66") +
  geom_line(aes(x = clu, y = bic), linewidth = 2) +
  geom_point(aes(x = clu, y = bic), size = 5) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  labs(x = "Number of clusters", y = "Value", caption = "AIC: gray, BIC: black")

print(g_aicbic)
ggsave("AIC_BIC_bnt_rev2.pdf", plot = g_aicbic, width = 8, height = 6)

# 3. ★★★ 推定パラメータ一覧（N, Proportion, Sigma付き）をCSVに出力 ★★★
write_csv(parameter_results, "estimated_parameters_bnt.csv")

# 4. 全クラスタリング結果を結合してCSVに書き出す
data_clean <- data %>%
  select(-any_of(c("predfl", "resid", "cluster", "ocluster", "ocluster_rename", "max_len"))) %>%
  select(-starts_with("ll-"))

final_data <- left_join(data_clean, cluster_assignments_df, by = "id")

write_csv(final_data, "WB_bnt_clustered_rev2.csv", na = "")

message("\n全ての処理が完了しました。")

