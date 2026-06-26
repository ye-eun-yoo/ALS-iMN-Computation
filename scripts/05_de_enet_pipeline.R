## Script 05: DE → Elastic Net — Two-Step Discovery Pipeline (Tier 2)
##
## Combines classical differential expression with penalized regression:
##   Step 1: edgeR ALS vs Healthy (iMN day 28) — identifies DE candidates
##   Step 2: Elastic net on DE-filtered gene set vs full genome
##
## Rationale: DE pre-filtering reduces dimensionality and focuses ML on
## biologically relevant variance, which can improve or hurt generalization
## depending on whether DE genes encode genotype structure beyond ALS/Ctrl.
##
## Outputs (plots/05_de_enet/):
##   01_volcano.pdf          — edgeR ALS vs Healthy volcano plot
##   02_de_summary.pdf       — # DE genes by FDR threshold
##   03_enet_comparison.pdf  — AUC: full genome vs DE-filtered elastic net
##   04_top_features.pdf     — top selected genes from DE-filtered model

REPO_DIR  <- Sys.getenv("REPO_DIR", unset = ".")
DATA_DIR  <- file.path(REPO_DIR, "data")
PLOT_DIR  <- file.path(REPO_DIR, "plots", "05_de_enet")
RES_DIR   <- file.path(REPO_DIR, "results")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

local_lib <- "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq/Rlib_tmp"
.libPaths(c(local_lib, .libPaths()))
for (pkg in c("glmnet", "pROC")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, lib = local_lib, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(edgeR); library(glmnet); library(pROC)
  library(ggplot2); library(dplyr); library(tidyr)
})

feat28 <- readRDS(file.path(DATA_DIR, "features_iMN_day28.RDS"))
dge_all <- readRDS(file.path(DATA_DIR, "GSE299997_DGEList.RDS"))
meta    <- dge_all$samples

# ── Step 1: edgeR DE — ALS vs Healthy, iMN day 28 ────────────────────────────
cat("=== Step 1: edgeR DE (ALS vs Healthy, iMN day28) ===\n")
idx28  <- meta$Stage_label == "iMN_day28"
dge28  <- dge_all[, idx28]
dge28  <- calcNormFactors(dge28, method = "TMM")

meta28 <- dge28$samples
meta28$genotype  <- factor(sub("-[0-9]+$", "", meta28$cell_line),
                            levels = c("Healthy", "SOD1", "FUS", "TARDBP", "C9orf72"))
meta28$condition <- factor(ifelse(meta28$genotype == "Healthy", "Healthy", "ALS"),
                            levels = c("Healthy", "ALS"))
dge28$samples <- meta28

des <- model.matrix(~ condition, data = meta28)
dge28 <- estimateDisp(dge28, des, robust = TRUE)
fit   <- glmQLFit(dge28, des, robust = TRUE)
qlf   <- glmQLFTest(fit, coef = "conditionALS")
de    <- topTags(qlf, n = Inf, sort.by = "none")$table
de$gene <- rownames(de)

cat(sprintf("FDR<0.05: %d genes | FDR<0.20: %d genes\n",
            sum(de$FDR < 0.05), sum(de$FDR < 0.20)))
write.csv(de, file.path(RES_DIR, "05_de_als_vs_healthy_imn28.csv"), row.names = FALSE)

# ── Plot 1: Volcano plot ──────────────────────────────────────────────────────
de$sig <- case_when(
  de$FDR < 0.05 & de$logFC >  1 ~ "Up in ALS (FDR<0.05)",
  de$FDR < 0.05 & de$logFC < -1 ~ "Down in ALS (FDR<0.05)",
  de$FDR < 0.20                  ~ "FDR<0.20",
  TRUE                           ~ "NS"
)
top_labels <- de %>% filter(FDR < 0.10) %>%
  arrange(FDR) %>% head(20)

p1 <- ggplot(de, aes(logFC, -log10(FDR), color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  ggrepel::geom_text_repel(data = top_labels, aes(label = gene),
                            size = 2.5, max.overlaps = 20,
                            box.padding = 0.3, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("Up in ALS (FDR<0.05)"   = "#d7191c",
               "Down in ALS (FDR<0.05)" = "#2c7bb6",
               "FDR<0.20"               = "#fdae61",
               "NS"                     = "grey70"),
    name = NULL) +
  labs(title = "edgeR — ALS vs Healthy (iMN day 28)",
       subtitle = "All ALS genotypes pooled: SOD1 / FUS / TARDBP / C9orf72 vs Healthy",
       x = "log2 fold change (ALS / Healthy)", y = "-log10(FDR)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(PLOT_DIR, "01_volcano.pdf"), p1, width = 7, height = 6)

# ── Plot 2: DE gene counts at multiple FDR thresholds ─────────────────────────
fdr_cuts <- c(0.01, 0.05, 0.10, 0.20, 0.30)
de_counts <- bind_rows(lapply(fdr_cuts, function(fdr) {
  data.frame(fdr_threshold = fdr,
             direction = c("Up in ALS", "Down in ALS"),
             n_genes   = c(sum(de$FDR < fdr & de$logFC > 0),
                           sum(de$FDR < fdr & de$logFC < 0)))
}))

p2 <- ggplot(de_counts, aes(factor(fdr_threshold), n_genes, fill = direction)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Up in ALS" = "#d7191c", "Down in ALS" = "#2c7bb6"),
                    name = NULL) +
  labs(title = "DE gene counts at multiple FDR thresholds",
       subtitle = "edgeR ALS vs Healthy | iMN day 28",
       x = "FDR threshold", y = "Number of genes") +
  theme_bw(base_size = 12)
ggsave(file.path(PLOT_DIR, "02_de_gene_counts.pdf"), p2, width = 6, height = 4)

# ── Step 2: Elastic net on DE-filtered vs full-genome ─────────────────────────
cat("\n=== Step 2: Elastic net — full genome vs DE-filtered ===\n")
X_full <- feat28$X
y      <- feat28$y_binary

# Three gene sets: all genes, DE FDR<0.20, DE FDR<0.05
gene_sets <- list(
  "All genes"     = colnames(X_full),
  "DE FDR<0.20"  = de$gene[de$FDR < 0.20],
  "DE FDR<0.05"  = de$gene[de$FDR < 0.05]
)
cat("Gene set sizes:\n")
for (nm in names(gene_sets)) cat(sprintf("  %s: %d genes\n", nm, length(gene_sets[[nm]])))

set.seed(42)
cv_results <- lapply(names(gene_sets), function(nm) {
  gs <- intersect(gene_sets[[nm]], colnames(X_full))
  if (length(gs) < 10) {
    cat(sprintf("  Skipping '%s': too few genes (%d)\n", nm, length(gs)); return(NULL))
  }
  X_sub <- X_full[, gs, drop = FALSE]
  cat(sprintf("  CV for '%s' (%d genes)...\n", nm, ncol(X_sub)))
  cv <- cv.glmnet(X_sub, y, family = "binomial", alpha = 0.5,
                   nfolds = 10, type.measure = "auc", standardize = FALSE)
  coefs <- coef(cv, s = "lambda.1se")[-1, 1]
  list(cv = cv, gene_set = nm, n_input = length(gs),
       n_selected = sum(coefs != 0), best_auc = max(cv$cvm),
       coefs = coefs[coefs != 0])
})
names(cv_results) <- names(gene_sets)
cv_results <- Filter(Negate(is.null), cv_results)

# ── Plot 3: AUC comparison across gene sets ───────────────────────────────────
auc_df <- bind_rows(lapply(cv_results, function(r) {
  data.frame(gene_set   = r$gene_set,
             n_input    = r$n_input,
             n_selected = r$n_selected,
             best_auc   = r$best_auc)
}))
auc_df$gene_set <- factor(auc_df$gene_set, levels = names(gene_sets))

p3 <- ggplot(auc_df, aes(gene_set, best_auc, fill = gene_set)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_text(aes(label = sprintf("AUC=%.3f\nn=%d→%d", best_auc, n_input, n_selected)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("All genes" = "#4dac26",
                                "DE FDR<0.20" = "#f1b6da",
                                "DE FDR<0.05" = "#d01c8b")) +
  scale_y_continuous(limits = c(0, 1.08), labels = scales::percent) +
  labs(title = "DE pre-filtering — does it improve elastic net?",
       subtitle = "10-fold CV AUC | Elastic Net (α=0.5) | iMN day 28",
       x = "Gene set", y = "Best CV AUC") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "03_enet_auc_comparison.pdf"), p3, width = 6, height = 4)

# ── Plot 4: Top features from DE-filtered model ───────────────────────────────
de_model_key <- if ("DE FDR<0.05" %in% names(cv_results) &&
                    cv_results[["DE FDR<0.05"]]$n_selected > 0) "DE FDR<0.05"
               else "DE FDR<0.20"

if (!is.null(cv_results[[de_model_key]])) {
  top_de_coefs <- sort(abs(cv_results[[de_model_key]]$coefs), decreasing = TRUE)
  top_de_genes <- head(names(top_de_coefs), 30)
  top_df <- data.frame(
    gene      = top_de_genes,
    coef      = cv_results[[de_model_key]]$coefs[top_de_genes],
    de_logfc  = de$logFC[match(top_de_genes, de$gene)],
    de_fdr    = de$FDR[match(top_de_genes, de$gene)]
  ) %>% mutate(gene = factor(gene, levels = rev(gene)))

  p4 <- ggplot(top_df, aes(coef, gene, fill = de_logfc)) +
    geom_col(alpha = 0.85) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, name = "DE logFC\n(ALS/Ctrl)") +
    geom_vline(xintercept = 0, color = "grey40") +
    labs(title = sprintf("Top genes — %s elastic net model", de_model_key),
         subtitle = "Bar = elastic net coefficient | Color = DE logFC",
         x = "Elastic net coefficient (ALS direction = positive)", y = NULL) +
    theme_bw(base_size = 10)
  ggsave(file.path(PLOT_DIR, "04_top_features_de_filtered.pdf"), p4, width = 7, height = 6)
}

write.csv(auc_df, file.path(RES_DIR, "05_de_enet_auc_comparison.csv"), row.names = FALSE)
cat("\nDone. Plots -> plots/05_de_enet/\n")
