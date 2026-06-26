## Script 06: Temporal Analysis - Disease Signatures Across Differentiation Stages
##
## iPSC-derived neurons go through: iPSC -> NPC -> iMN day28 -> iMN day50
## Questions:
##   1. Does the ALS transcriptional signature emerge or strengthen over time?
##      (Train on day28, test on day50 - compare AUC)
##   2. Do the same genes drive ALS classification at both timepoints?
##   3. PCA trajectory: how do genotypes diverge across differentiation?
##
## Outputs (plots/06_temporal/):
##   01_pca_trajectory.pdf   - PCA colored by stage + genotype
##   02_auc_temporal.pdf     - binary ALS AUC at day28 vs day50
##   03_feature_overlap.pdf  - Venn: genes selected at day28 vs day50
##   04_logfc_heatmap.pdf    - top genes: logFC at day28 vs day50

REPO_DIR  <- Sys.getenv("REPO_DIR", unset = "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/ALS-iMN-ML")
DATA_DIR  <- file.path(REPO_DIR, "data")
PLOT_DIR  <- file.path(REPO_DIR, "plots", "06_temporal")
RES_DIR   <- file.path(REPO_DIR, "results")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

local_lib <- "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq/Rlib_tmp"
.libPaths(c(local_lib, .libPaths()))
for (pkg in c("glmnet", "pROC", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, lib = local_lib, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(edgeR)
  library(ggplot2); library(dplyr); library(tidyr); library(ggrepel)
})

GENOTYPE_COLORS <- c(Healthy = "#4dac26", SOD1 = "#d01c8b",
                      FUS = "#f1b6da", TARDBP = "#b8e186", C9orf72 = "#7b3294")
STAGE_COLORS    <- c(iPSC = "#ffffcc", NPC = "#a1dab4",
                     iMN_day28 = "#41b6c4", iMN_day50 = "#225ea8")

feat28 <- readRDS(file.path(DATA_DIR, "features_iMN_day28.RDS"))
feat50 <- readRDS(file.path(DATA_DIR, "features_iMN_day50.RDS"))
dge    <- readRDS(file.path(DATA_DIR, "GSE299997_DGEList.RDS"))
meta   <- dge$samples
meta$genotype  <- factor(sub("-[0-9]+$", "", meta$cell_line),
                          levels = c("Healthy", "SOD1", "FUS", "TARDBP", "C9orf72"))
meta$condition <- ifelse(meta$genotype == "Healthy", "Healthy", "ALS")

# ── Plot 1: PCA trajectory across all 4 stages ───────────────────────────────
cat("=== PCA trajectory across all stages ===\n")
logcpm  <- cpm(dge, log = TRUE, prior.count = 1)
pca_all <- prcomp(t(logcpm), scale. = TRUE)
pve     <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)

pca_df <- data.frame(
  PC1      = pca_all$x[, 1],
  PC2      = pca_all$x[, 2],
  Stage    = meta$Stage_label,
  Genotype = meta$genotype,
  Condition= meta$condition
)
stage_order <- c("iPSC", "NPC", "iMN_day28", "iMN_day50")
pca_df$Stage <- factor(pca_df$Stage, levels = stage_order)

# Centroids per Stage x Genotype for trajectory arrows
centroids <- pca_df %>%
  group_by(Stage, Genotype) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop") %>%
  arrange(Genotype, Stage)

p1 <- ggplot(pca_df, aes(PC1, PC2, color = Genotype, shape = Stage)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_path(data = centroids, aes(group = Genotype), linewidth = 0.6,
            arrow = arrow(length = unit(0.15, "cm"), type = "open"),
            linetype = "dashed") +
  geom_point(data = centroids, size = 4, stroke = 1.2, shape = 21,
             aes(fill = Genotype), color = "black") +
  scale_color_manual(values = GENOTYPE_COLORS) +
  scale_fill_manual(values  = GENOTYPE_COLORS) +
  scale_shape_manual(values = c(iPSC = 16, NPC = 17, iMN_day28 = 15, iMN_day50 = 18),
                     name = "Stage") +
  labs(title = "Differentiation trajectory - iPSC -> NPC -> iMN day28 -> day50",
       subtitle = "Arrows connect stage centroids per genotype",
       x = sprintf("PC1 (%s%%)", pve[1]), y = sprintf("PC2 (%s%%)", pve[2])) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_bw(base_size = 11)
ggsave(file.path(PLOT_DIR, "01_pca_trajectory.pdf"), p1, width = 8, height = 6)

# ── Elastic net at each timepoint ─────────────────────────────────────────────
cat("\n=== Elastic net at day28 and day50 ===\n")
set.seed(42)

fit_enet <- function(feat) {
  X  <- feat$X
  y  <- feat$y_binary
  shared <- which(colnames(X) %in% colnames(feat28$X) &
                  colnames(X) %in% colnames(feat50$X))
  X  <- X[, shared]
  cv <- cv.glmnet(X, y, family = "binomial", alpha = 0.5,
                   nfolds = 10, type.measure = "auc", standardize = FALSE)
  coefs <- coef(cv, s = "lambda.1se")[-1, 1]
  sel   <- coefs[coefs != 0]
  list(cv = cv, selected = sel, best_auc = min(max(cv$cvm), 1.0))
}

res28 <- fit_enet(feat28)
res50 <- fit_enet(feat50)

cat(sprintf("Day28: AUC=%.3f | genes=%d\n", res28$best_auc, length(res28$selected)))
cat(sprintf("Day50: AUC=%.3f | genes=%d\n", res50$best_auc, length(res50$selected)))

# ── Plot 2: AUC at day28 vs day50 ─────────────────────────────────────────────
auc_df <- data.frame(
  Stage     = c("iMN day 28", "iMN day 50"),
  CV_AUC    = c(res28$best_auc, res50$best_auc),
  N_genes   = c(length(res28$selected), length(res50$selected))
)

p2 <- ggplot(auc_df, aes(Stage, CV_AUC, fill = Stage)) +
  geom_col(width = 0.45, alpha = 0.85) +
  geom_text(aes(label = sprintf("AUC = %.3f\n%d genes selected",
                                CV_AUC, N_genes)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("iMN day 28" = "#41b6c4", "iMN day 50" = "#225ea8")) +
  scale_y_continuous(limits = c(0, 1.15), oob = scales::squish) +
  labs(title = "Does the ALS signature strengthen over differentiation?",
       subtitle = "Elastic net (alpha=0.5), 10-fold CV AUC | ALS vs Healthy",
       x = NULL, y = "Best CV AUC") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "02_auc_day28_vs_day50.pdf"), p2, width = 5, height = 4)

# ── Plot 3: Feature overlap between day28 and day50 ──────────────────────────
g28  <- names(res28$selected)
g50  <- names(res50$selected)
both <- intersect(g28, g50)
only28 <- setdiff(g28, g50)
only50 <- setdiff(g50, g28)

cat(sprintf("\nGenes selected: day28=%d, day50=%d, shared=%d\n",
            length(g28), length(g50), length(both)))

# Show shared genes: coefficient at day28 vs day50
if (length(both) > 0) {
  shared_df <- data.frame(
    gene   = both,
    coef28 = res28$selected[both],
    coef50 = res50$selected[both]
  ) %>% arrange(desc(abs(coef28)))

  p3 <- ggplot(shared_df, aes(coef28, coef50, label = gene)) +
    geom_point(color = "#2c7bb6", size = 3, alpha = 0.8) +
    geom_text_repel(size = 2.8, max.overlaps = 20) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    labs(title = sprintf("Shared elastic net genes: %d (day28 ∩ day50)", length(both)),
         subtitle = "Coefficient direction should be consistent across timepoints",
         x = "Coefficient at day 28", y = "Coefficient at day 50") +
    theme_bw(base_size = 11)
  ggsave(file.path(PLOT_DIR, "03_shared_gene_coefficients.pdf"), p3, width = 6, height = 5)
}

# ── Plot 4: logFC heatmap at day28 and day50 for shared top genes ─────────────
run_edger <- function(dge_in, idx) {
  dge_s <- dge_in[, idx]
  dge_s <- calcNormFactors(dge_s, method = "TMM")
  m <- dge_s$samples
  m$genotype  <- factor(sub("-[0-9]+$", "", m$cell_line),
                         levels = c("Healthy", "SOD1", "FUS", "TARDBP", "C9orf72"))
  m$condition <- factor(ifelse(m$genotype == "Healthy", "Healthy", "ALS"),
                         levels = c("Healthy", "ALS"))
  dge_s$samples <- m
  des <- model.matrix(~ condition, data = m)
  dge_s <- estimateDisp(dge_s, des, robust = TRUE)
  fit   <- glmQLFit(dge_s, des, robust = TRUE)
  qlf   <- glmQLFTest(fit, coef = "conditionALS")
  topTags(qlf, n = Inf, sort.by = "none")$table
}

de28 <- run_edger(dge, meta$Stage_label == "iMN_day28")
de50 <- run_edger(dge, meta$Stage_label == "iMN_day50")

# Focus on union of selected genes
focus_genes <- union(g28, g50)
focus_genes <- focus_genes[focus_genes %in% rownames(de28) &
                           focus_genes %in% rownames(de50)]
focus_genes <- head(focus_genes, 40)

heat_df <- data.frame(
  gene  = focus_genes,
  logFC_day28 = de28[focus_genes, "logFC"],
  logFC_day50 = de50[focus_genes, "logFC"],
  fdr28 = de28[focus_genes, "FDR"],
  fdr50 = de50[focus_genes, "FDR"],
  selected = case_when(
    focus_genes %in% both   ~ "Both",
    focus_genes %in% only28 ~ "Day28 only",
    focus_genes %in% only50 ~ "Day50 only",
    TRUE                    ~ "Neither"
  )
) %>% pivot_longer(c(logFC_day28, logFC_day50),
                    names_to = "Timepoint", values_to = "logFC") %>%
  mutate(Timepoint = recode(Timepoint,
                             logFC_day28 = "iMN day 28",
                             logFC_day50 = "iMN day 50"),
         gene = factor(gene, levels = rev(focus_genes)))

p4 <- ggplot(heat_df, aes(Timepoint, gene, fill = logFC)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                       midpoint = 0, name = "logFC\n(ALS/Ctrl)") +
  facet_grid(selected ~ ., scales = "free_y", space = "free_y") +
  labs(title = "ALS vs Healthy logFC at day28 and day50",
       subtitle = "Genes selected by elastic net (one or both timepoints)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.text.y = element_text(size = 7, face = "bold"),
        panel.spacing = unit(0.2, "lines"))
ggsave(file.path(PLOT_DIR, "04_logfc_heatmap_temporal.pdf"), p4,
       width = 5, height = max(6, length(focus_genes) * 0.22 + 3))

write.csv(heat_df, file.path(RES_DIR, "06_temporal_logfc.csv"), row.names = FALSE)
cat("\nDone. Plots -> plots/06_temporal/\n")
