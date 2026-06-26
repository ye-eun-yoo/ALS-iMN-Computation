## Script 04: Random Forest — Multi-class Genotype Classification
##
## Classifies motor neuron samples into 5 genotype classes:
##   Healthy / SOD1 / FUS / TARDBP / C9orf72
##
## Unlike the binary elastic net, RF handles multi-class natively and
## provides gene importance scores without assuming linearity.
## Training: iMN day 28 | Evaluation: OOB error + held-out iMN day 50
##
## Outputs (plots/04_random_forest/):
##   01_oob_error.pdf        — OOB classification error vs n trees
##   02_importance.pdf       — top 30 genes by mean decrease in accuracy
##   03_confusion_matrix.pdf — predicted vs true genotype (iMN day 50)
##   04_class_accuracy.pdf   — per-class precision / recall

REPO_DIR  <- Sys.getenv("REPO_DIR", unset = ".")
DATA_DIR  <- file.path(REPO_DIR, "data")
PLOT_DIR  <- file.path(REPO_DIR, "plots", "04_random_forest")
RES_DIR   <- file.path(REPO_DIR, "results")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

local_lib <- "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq/Rlib_tmp"
.libPaths(c(local_lib, .libPaths()))
for (pkg in c("ranger", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, lib = local_lib, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ranger); library(ggplot2); library(dplyr); library(tidyr)
})

GENOTYPE_COLORS <- c(Healthy = "#4dac26", SOD1 = "#d01c8b",
                      FUS = "#f1b6da", TARDBP = "#b8e186", C9orf72 = "#7b3294")

# ── Load features ─────────────────────────────────────────────────────────────
feat28 <- readRDS(file.path(DATA_DIR, "features_iMN_day28.RDS"))
feat50 <- readRDS(file.path(DATA_DIR, "features_iMN_day50.RDS"))

X28 <- feat28$X;  y28 <- feat28$y_multi
X50 <- feat50$X;  y50 <- feat50$y_multi

cat("Training set (day28):\n"); print(table(y28))
cat("Test set (day50):\n");     print(table(y50))

# Align columns between day28 and day50 (same genes, same order)
shared_genes <- intersect(colnames(X28), colnames(X50))
X28 <- X28[, shared_genes];  X50 <- X50[, shared_genes]

set.seed(42)

# ── Train RF (ranger for speed, ~17k features) ────────────────────────────────
cat("\n=== Training Random Forest (500 trees, iMN day28) ===\n")
train_df      <- as.data.frame(X28)
train_df$y    <- y28

rf_fit <- ranger(
  y ~ ., data = train_df,
  num.trees          = 500,
  importance         = "permutation",
  probability        = TRUE,
  num.threads        = 4,
  seed               = 42
)
cat(sprintf("OOB prediction error: %.3f\n", rf_fit$prediction.error))

# ── Plot 1: OOB error convergence ─────────────────────────────────────────────
# Refit tracking OOB at each tree step (ranger stores this)
rf_oob <- ranger(y ~ ., data = train_df, num.trees = 500,
                 importance = "none", keep.inbag = FALSE,
                 num.threads = 4, seed = 42)

# ranger doesn't expose per-tree OOB natively; refit at intervals
tree_seq <- seq(10, 500, by = 10)
oob_seq  <- vapply(tree_seq, function(nt) {
  rf_tmp <- ranger(y ~ ., data = train_df, num.trees = nt,
                   probability = TRUE, num.threads = 4, seed = 42)
  rf_tmp$prediction.error
}, numeric(1))

oob_df <- data.frame(n_trees = tree_seq, oob_error = oob_seq)
p1 <- ggplot(oob_df, aes(n_trees, oob_error)) +
  geom_line(color = "#2c7bb6", linewidth = 0.8) +
  geom_point(data = oob_df[nrow(oob_df), ], size = 3, color = "#d7191c") +
  labs(title = "Random Forest — OOB error convergence",
       subtitle = "Multi-class: Healthy / SOD1 / FUS / TARDBP / C9orf72 | iMN day 28",
       x = "Number of trees", y = "OOB classification error") +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "01_oob_error.pdf"), p1, width = 6, height = 4)

# ── Plot 2: Variable importance (top 30) ─────────────────────────────────────
imp_df <- data.frame(
  gene       = names(rf_fit$variable.importance),
  importance = rf_fit$variable.importance
) %>% arrange(desc(importance)) %>% head(30) %>%
  mutate(gene = factor(gene, levels = rev(gene)))

p2 <- ggplot(imp_df, aes(importance, gene)) +
  geom_col(fill = "#2c7bb6", alpha = 0.8) +
  labs(title = "Random Forest — Top 30 genes by permutation importance",
       subtitle = "Mean decrease in classification accuracy when gene is permuted",
       x = "Mean decrease in accuracy", y = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(PLOT_DIR, "02_importance.pdf"), p2, width = 7, height = 6)

write.csv(data.frame(gene = names(rf_fit$variable.importance),
                     importance = rf_fit$variable.importance) %>%
            arrange(desc(importance)),
          file.path(RES_DIR, "04_rf_gene_importance.csv"), row.names = FALSE)

# ── Predict on iMN day 50 ────────────────────────────────────────────────────
cat("\n=== Predicting on iMN day 50 (temporal transfer) ===\n")
test_df  <- as.data.frame(X50)
pred50   <- predict(rf_fit, data = test_df)
pred_cls <- factor(colnames(pred50$predictions)[apply(pred50$predictions, 1, which.max)],
                   levels = levels(y28))
true_cls <- factor(y50, levels = levels(y28))

test_acc <- mean(pred_cls == true_cls)
cat(sprintf("Day50 test accuracy: %.3f\n", test_acc))

# ── Plot 3: Confusion matrix (day50 test) ────────────────────────────────────
cm_df <- as.data.frame(table(Predicted = pred_cls, True = true_cls))
cm_df <- cm_df %>% group_by(True) %>% mutate(pct = 100 * Freq / sum(Freq))

p3 <- ggplot(cm_df, aes(True, Predicted, fill = pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", Freq, pct)), size = 3.2) +
  scale_fill_gradient(low = "white", high = "#2c7bb6", name = "% of true class") +
  labs(title = "Random Forest — Confusion matrix (iMN day 50, held-out)",
       subtitle = sprintf("Train: day28 | Test: day50 | Accuracy = %.1f%%", 100 * test_acc),
       x = "True genotype", y = "Predicted genotype") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(PLOT_DIR, "03_confusion_matrix.pdf"), p3, width = 6, height = 5)

# ── Plot 4: Per-class accuracy ────────────────────────────────────────────────
# Compute precision and recall per class
class_stats <- lapply(levels(y28), function(cls) {
  TP <- sum(pred_cls == cls & true_cls == cls)
  FP <- sum(pred_cls == cls & true_cls != cls)
  FN <- sum(pred_cls != cls & true_cls == cls)
  data.frame(
    genotype  = cls,
    precision = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
    recall    = ifelse(TP + FN == 0, NA, TP / (TP + FN)),
    n_test    = sum(true_cls == cls)
  )
})
class_df  <- bind_rows(class_stats) %>%
  pivot_longer(c(precision, recall), names_to = "metric", values_to = "value")

p4 <- ggplot(class_df, aes(genotype, value, fill = metric)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.65) +
  scale_fill_manual(values = c(precision = "#2c7bb6", recall = "#d7191c"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Per-class precision and recall (day50 test set)",
       subtitle = "C9orf72: n=1 line — interpret with caution",
       x = "Genotype", y = "Score") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(PLOT_DIR, "04_class_accuracy.pdf"), p4, width = 6, height = 4)

cat("\nDone. Plots -> plots/04_random_forest/\n")
cat("Gene importance -> results/04_rf_gene_importance.csv\n")
