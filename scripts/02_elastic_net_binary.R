## Script 02: Elastic Net - Binary Classification (ALS vs Healthy, iMN day 28)
##
## Unbiased gene discovery: penalized logistic regression sweeps 17,000+ genes
## simultaneously to find which subset best discriminates ALS from healthy
## iPSC-derived motor neurons.
##
## Method: glmnet elastic net (alpha=0.5), 10-fold CV for lambda selection.
## Note: with 3 replicates/line and 13 lines, CV folds may split replicates
##   from the same donor into train/test - treat AUC as an estimate of
##   discriminative power rather than unbiased generalization error.
##
## Outputs (plots/02_elastic_net_binary/):
##   01_cv_curve.pdf          - CV deviance / AUC vs log(lambda)
##   02_coef_path.pdf         - regularization path: coefficients vs log(lambda)
##   03_selected_heatmap.pdf  - top selected genes x samples heatmap
##   04_roc_curve.pdf         - ROC with AUC

REPO_DIR  <- Sys.getenv("REPO_DIR", unset = "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/ALS-iMN-ML")
DATA_DIR  <- file.path(REPO_DIR, "data")
PLOT_DIR  <- file.path(REPO_DIR, "plots", "02_elastic_net_binary")
RES_DIR   <- file.path(REPO_DIR, "results")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR,  showWarnings = FALSE, recursive = TRUE)

local_lib <- "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq/Rlib_tmp"
.libPaths(c(local_lib, .libPaths()))
for (pkg in c("glmnet", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, lib = local_lib, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(glmnet); library(pROC)
  library(ggplot2); library(dplyr); library(tidyr)
})

GENOTYPE_COLORS <- c(Healthy = "#4dac26", SOD1 = "#d01c8b",
                      FUS = "#f1b6da", TARDBP = "#b8e186", C9orf72 = "#7b3294")

# ── Load features ─────────────────────────────────────────────────────────────
feat <- readRDS(file.path(DATA_DIR, "features_iMN_day28.RDS"))
X <- feat$X              # samples x genes, z-scored
y <- feat$y_binary       # Healthy / ALS

cat(sprintf("Input: %d samples x %d genes\n", nrow(X), ncol(X)))
cat("Labels:\n"); print(table(y))

set.seed(42)

# ── 10-fold CV elastic net (alpha = 0.5) ─────────────────────────────────────
cat("\n=== Running 10-fold CV elastic net (alpha = 0.5) ===\n")
cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 0.5,
                    nfolds = 10, type.measure = "auc", standardize = FALSE)

lambda_best <- cv_fit$lambda.1se
cat(sprintf("Best lambda (1se): %.4f\n", lambda_best))

coefs <- coef(cv_fit, s = "lambda.1se")
sel   <- coefs[-1, 1]               # drop intercept
sel   <- sel[sel != 0]
cat(sprintf("Genes selected at lambda.1se: %d\n", length(sel)))

# Save results
res_df <- data.frame(
  gene   = names(sel),
  coefficient = as.numeric(sel),
  direction = ifelse(sel > 0, "up_in_ALS", "down_in_ALS")
) %>% arrange(desc(abs(coefficient)))
write.csv(res_df, file.path(RES_DIR, "02_elastic_net_binary_coefficients.csv"),
          row.names = FALSE)

# ── Plot 1: CV AUC curve ──────────────────────────────────────────────────────
cv_df <- data.frame(
  log_lambda = log(cv_fit$lambda),
  mean_auc   = cv_fit$cvm,
  se         = cv_fit$cvsd,
  nzero      = cv_fit$nzero
)
lambda_min_log <- log(cv_fit$lambda.min)
lambda_1se_log <- log(cv_fit$lambda.1se)

p1 <- ggplot(cv_df, aes(log_lambda, mean_auc)) +
  geom_ribbon(aes(ymin = mean_auc - se, ymax = mean_auc + se), fill = "grey80") +
  geom_line(color = "#2c7bb6", linewidth = 0.8) +
  geom_vline(xintercept = lambda_min_log, linetype = "dashed", color = "#d7191c") +
  geom_vline(xintercept = lambda_1se_log, linetype = "dotted", color = "#d7191c") +
  annotate("text", x = lambda_min_log, y = min(cv_df$mean_auc) + 0.01,
           label = "lambda.min", hjust = -0.1, size = 3, color = "#d7191c") +
  annotate("text", x = lambda_1se_log, y = min(cv_df$mean_auc) + 0.01,
           label = "lambda.1se", hjust = 1.1, size = 3, color = "#d7191c") +
  labs(title = "Elastic Net - 10-fold CV (alpha = 0.5)",
       subtitle = "ALS vs Healthy | iMN day 28",
       x = "log(lambda)", y = "Mean AUC (10-fold CV)") +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "01_cv_auc_curve.pdf"), p1, width = 6, height = 4)

# ── Plot 2: Regularization path (top 20 genes by |final coef|) ───────────────
top_genes <- head(names(sel)[order(abs(sel), decreasing = TRUE)], 20)
path_mat  <- as.matrix(coef(cv_fit$glmnet.fit)[-1, ])    # genes x lambdas
path_df   <- as.data.frame(path_mat[top_genes, , drop = FALSE])
path_df$gene <- rownames(path_df)
path_long <- pivot_longer(path_df, -gene, names_to = "lambda_idx", values_to = "coef") %>%
  mutate(log_lambda = rep(log(cv_fit$glmnet.fit$lambda), each = length(top_genes)))

p2 <- ggplot(path_long, aes(log_lambda, coef, color = gene, group = gene)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_vline(xintercept = lambda_1se_log, linetype = "dotted", color = "grey40") +
  scale_color_viridis_d(option = "turbo") +
  labs(title = "Regularization path - top 20 genes by |coefficient|",
       subtitle = "Dotted = lambda.1se used for feature selection",
       x = "log(lambda)", y = "Coefficient", color = NULL) +
  guides(color = guide_legend(ncol = 1, key.height = unit(0.3, "cm"), label.theme = element_text(size = 7))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")
ggsave(file.path(PLOT_DIR, "02_coef_path.pdf"), p2, width = 9, height = 5)

# ── Plot 3: Selected genes heatmap ───────────────────────────────────────────
top30 <- head(res_df$gene, min(30, nrow(res_df)))
if (length(top30) > 0) {
  heat_mat <- t(X[, top30, drop = FALSE])
  heat_df  <- as.data.frame(heat_mat)
  heat_df$gene <- rownames(heat_df)
  heat_long <- pivot_longer(heat_df, -gene, names_to = "sample", values_to = "zscore")
  meta_ann  <- feat$meta[, c("sample_id", "genotype", "cell_line")]
  heat_long <- left_join(heat_long, meta_ann, by = c("sample" = "sample_id"))
  heat_long$gene <- factor(heat_long$gene,
                            levels = top30[length(top30):1])
  heat_long$sample <- factor(heat_long$sample,
                              levels = rownames(feat$meta)[order(feat$meta$genotype)])

  p3 <- ggplot(heat_long, aes(sample, gene, fill = zscore)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, limits = c(-3, 3), oob = scales::squish,
                         name = "z-score") +
    facet_grid(. ~ genotype, scales = "free_x", space = "free_x") +
    labs(title = sprintf("Top %d elastic-net selected genes (iMN day 28)", length(top30)),
         x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.text = element_text(size = 8, face = "bold"),
          panel.spacing = unit(0.15, "lines"))
  ggsave(file.path(PLOT_DIR, "03_selected_heatmap.pdf"), p3,
         width = 9, height = max(5, length(top30) * 0.22 + 2))
}

# ── Plot 4: ROC curve (in-sample, for illustration) ──────────────────────────
pred_prob <- predict(cv_fit, X, s = "lambda.1se", type = "response")[, 1]
roc_obj   <- roc(y, pred_prob, levels = c("Healthy", "ALS"), direction = "<",
                 quiet = TRUE)
auc_val   <- round(as.numeric(auc(roc_obj)), 3)
roc_df    <- data.frame(
  sensitivity = rev(roc_obj$sensitivities),
  specificity = rev(1 - roc_obj$specificities)
)

p4 <- ggplot(roc_df, aes(specificity, sensitivity)) +
  geom_line(color = "#2c7bb6", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC = %s\n(in-sample, lambda.1se)", auc_val),
           size = 4) +
  labs(title = "ROC - Elastic Net, ALS vs Healthy",
       subtitle = "iMN day 28 | Note: in-sample AUC overestimates generalization",
       x = "1 - Specificity", y = "Sensitivity") +
  coord_equal() + theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "04_roc_curve.pdf"), p4, width = 5, height = 5)

cat(sprintf("\nDone. %d genes selected | In-sample AUC = %s\n", length(sel), auc_val))
cat("Results -> results/02_elastic_net_binary_coefficients.csv\n")
cat("Plots   -> plots/02_elastic_net_binary/\n")
