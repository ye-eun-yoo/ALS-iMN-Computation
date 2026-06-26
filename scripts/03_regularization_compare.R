## Script 03: Regularization Comparison - Ridge vs Elastic Net vs LASSO
##
## Compares three regularization strategies for classifying ALS vs Healthy
## in iMN day 28 RNA-seq:
##   alpha = 0.00  (Ridge)         - shrinks all genes, none dropped
##   alpha = 0.50  (Elastic Net)   - balance of Ridge and LASSO behavior
##   alpha = 1.00  (LASSO)         - aggressively sparse, some genes zeroed out
##
## Outputs (plots/03_regularization/):
##   01_cv_auc_comparison.pdf  - AUC vs log(lambda) for all three
##   02_sparsity.pdf           - # features selected vs alpha at best lambda
##   03_top_gene_stability.pdf - do the same top genes appear across methods?

REPO_DIR  <- Sys.getenv("REPO_DIR", unset = "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/ALS-iMN-ML")
DATA_DIR  <- file.path(REPO_DIR, "data")
PLOT_DIR  <- file.path(REPO_DIR, "plots", "03_regularization")
RES_DIR   <- file.path(REPO_DIR, "results")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

local_lib <- "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq/Rlib_tmp"
.libPaths(c(local_lib, .libPaths()))
for (pkg in c("glmnet")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, lib = local_lib, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2); library(dplyr); library(tidyr)
})

feat <- readRDS(file.path(DATA_DIR, "features_iMN_day28.RDS"))
X <- feat$X
y <- feat$y_binary
set.seed(42)

ALPHAS      <- c(0, 0.5, 1)
ALPHA_NAMES <- c("0" = "Ridge (alpha=0)", "0.5" = "Elastic Net (alpha=0.5)", "1" = "LASSO (alpha=1)")
ALPHA_COLS  <- c("0" = "#2c7bb6", "0.5" = "#1a9641", "1" = "#d7191c")

# ── Run CV for each alpha ─────────────────────────────────────────────────────
cat("Running 10-fold CV for each alpha...\n")
cv_list <- lapply(ALPHAS, function(a) {
  cat(sprintf("  alpha = %s\n", a))
  cv.glmnet(X, y, family = "binomial", alpha = a,
            nfolds = 10, type.measure = "auc", standardize = FALSE)
})
names(cv_list) <- as.character(ALPHAS)

# ── Plot 1: CV AUC curves side by side ───────────────────────────────────────
auc_df <- bind_rows(lapply(names(cv_list), function(a) {
  cv <- cv_list[[a]]
  data.frame(
    method     = a,
    log_lambda = log(cv$lambda),
    mean_auc   = cv$cvm,
    lo         = cv$cvm - cv$cvsd,
    hi         = cv$cvm + cv$cvsd
  )
}))
auc_df$alpha_label <- ALPHA_NAMES[auc_df$method]

vlines <- bind_rows(lapply(names(cv_list), function(a) {
  data.frame(method = a, alpha_label = ALPHA_NAMES[a],
             lambda_1se = log(cv_list[[a]]$lambda.1se),
             best_auc   = max(cv_list[[a]]$cvm))
}))

p1 <- ggplot(auc_df, aes(log_lambda, mean_auc, color = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = method), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = vlines, aes(xintercept = lambda_1se, color = method),
             linetype = "dotted") +
  scale_color_manual(values = ALPHA_COLS, labels = ALPHA_NAMES, name = "Method") +
  scale_fill_manual(values  = ALPHA_COLS, labels = ALPHA_NAMES, name = "Method") +
  facet_wrap(~ alpha_label, nrow = 1) +
  labs(title = "Regularization comparison - 10-fold CV AUC",
       subtitle = "ALS vs Healthy | iMN day 28 | Dotted = lambda.1se",
       x = "log(lambda)", y = "CV AUC") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave(file.path(PLOT_DIR, "01_cv_auc_comparison.pdf"), p1, width = 10, height = 4)

# ── Plot 2: Sparsity - # features selected vs alpha ──────────────────────────
sparsity_df <- bind_rows(lapply(names(cv_list), function(a) {
  cv <- cv_list[[a]]
  bind_rows(lapply(c("lambda.min", "lambda.1se"), function(lam) {
    coefs <- coef(cv, s = lam)[-1, 1]
    data.frame(alpha_val  = as.numeric(a),
               lambda_rule = lam,
               n_selected  = sum(coefs != 0),
               best_auc    = cv$cvm[which.min(abs(cv$lambda - cv[[lam]]))])
  }))
}))

p2 <- ggplot(sparsity_df, aes(alpha_val, n_selected, color = lambda_rule,
                               group = lambda_rule)) +
  geom_line(linewidth = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c(lambda.min = "#2c7bb6", lambda.1se = "#d7191c"),
                     labels = c("lambda.min", "lambda.1se (sparser)"),
                     name = "Lambda rule") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("Ridge\n(alpha=0)", "Elastic Net\n(alpha=0.5)", "LASSO\n(alpha=1)")) +
  labs(title = "Number of genes selected vs regularization strength",
       subtitle = "Ridge retains all genes; LASSO aggressively prunes",
       x = "Alpha", y = "Genes with non-zero coefficient") +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "02_sparsity.pdf"), p2, width = 6, height = 4)

# ── Plot 3: Top-gene overlap across methods ───────────────────────────────────
# Extract top 50 genes by |coefficient| for each method (at lambda.1se)
top_per_method <- lapply(names(cv_list), function(a) {
  cv    <- cv_list[[a]]
  coefs <- coef(cv, s = "lambda.1se")[-1, 1]
  coefs <- coefs[coefs != 0]
  if (length(coefs) == 0) return(character(0))
  head(names(sort(abs(coefs), decreasing = TRUE)), 50)
})
names(top_per_method) <- names(cv_list)

# Build presence/absence matrix
all_genes  <- unique(unlist(top_per_method))
overlap_df <- data.frame(
  gene  = all_genes,
  Ridge = all_genes %in% top_per_method[["0"]],
  EN    = all_genes %in% top_per_method[["0.5"]],
  LASSO = all_genes %in% top_per_method[["1"]]
)
overlap_df$n_methods <- rowSums(overlap_df[, 2:4])
overlap_df <- overlap_df %>% filter(n_methods >= 2) %>%
  arrange(desc(n_methods), gene)

stab_long <- pivot_longer(overlap_df, c(Ridge, EN, LASSO),
                           names_to = "method", values_to = "selected")
stab_long$method <- factor(stab_long$method, levels = c("Ridge", "EN", "LASSO"))
stab_long$gene   <- factor(stab_long$gene,
                            levels = overlap_df$gene[nrow(overlap_df):1])

p3 <- ggplot(stab_long, aes(method, gene, fill = selected)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "#f0f0f0", "TRUE" = "#2c7bb6"),
                    labels = c("Not selected", "Selected"), name = NULL) +
  labs(title = "Genes selected by >=2 regularization methods (top 50 each)",
       subtitle = "Shared genes = robust signal independent of penalty type",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(PLOT_DIR, "03_top_gene_stability.pdf"), p3,
       width = 5, height = max(4, nrow(overlap_df) * 0.2 + 2))

# Summary table
summary_df <- bind_rows(lapply(names(cv_list), function(a) {
  cv    <- cv_list[[a]]
  coefs <- coef(cv, s = "lambda.1se")[-1, 1]
  data.frame(
    alpha        = as.numeric(a),
    method       = ALPHA_NAMES[a],
    best_cv_auc  = round(max(cv$cvm), 4),
    n_genes_1se  = sum(coefs != 0),
    lambda_1se   = round(cv$lambda.1se, 5)
  )
}))
cat("\n=== Summary ===\n"); print(summary_df)
write.csv(summary_df, file.path(RES_DIR, "03_regularization_summary.csv"), row.names = FALSE)
cat("\nPlots -> plots/03_regularization/\n")
