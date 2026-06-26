## Script 01: Prepare feature matrices for ML
##
## Reads the DGEList, computes TMM-normalized logCPM, and builds
## ready-to-use feature matrices for each stage (iMN_day28, iMN_day50).
##
## Outputs:
##   data/features_iMN_day28.RDS  — list(X, y_binary, y_multi, meta)
##   data/features_iMN_day50.RDS
##   plots/01_qc/01_libsize.pdf
##   plots/01_qc/02_pca_stages.pdf
##   plots/01_qc/03_pca_genotypes.pdf

REPO_DIR  <- Sys.getenv("REPO_DIR", unset = ".")
DATA_DIR  <- file.path(REPO_DIR, "data")
PLOT_DIR  <- file.path(REPO_DIR, "plots", "01_qc")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

local_lib <- "/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq/Rlib_tmp"
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(edgeR); library(ggplot2); library(dplyr)
})

GENOTYPE_COLORS <- c(Healthy = "#4dac26", SOD1 = "#d01c8b",
                      FUS = "#f1b6da", TARDBP = "#b8e186", C9orf72 = "#7b3294")
STAGE_COLORS    <- c(iPSC = "#ffffcc", NPC = "#a1dab4",
                     iMN_day28 = "#41b6c4", iMN_day50 = "#225ea8")

# ── Load DGEList ──────────────────────────────────────────────────────────────
dge  <- readRDS(file.path(DATA_DIR, "GSE299997_DGEList.RDS"))
meta <- dge$samples
# Ensure genotype is a factor with Healthy first
meta$genotype <- factor(sub("-[0-9]+$", "", meta$cell_line),
                        levels = c("Healthy", "SOD1", "FUS", "TARDBP", "C9orf72"))
meta$condition <- ifelse(meta$genotype == "Healthy", "Healthy", "ALS")
dge$samples    <- meta

logcpm <- cpm(dge, log = TRUE, prior.count = 1)

cat(sprintf("Loaded: %d genes x %d samples across %d stages\n",
            nrow(dge), ncol(dge), length(unique(meta$Stage_label))))
print(table(meta$Stage_label, meta$genotype))

# ── QC Plot 1: Library size distribution ─────────────────────────────────────
p1 <- ggplot(meta, aes(x = reorder(sample_id, lib.size), y = lib.size / 1e6,
                        fill = Stage_label)) +
  geom_col() +
  scale_fill_manual(values = STAGE_COLORS, name = "Stage") +
  labs(title = "GSE299997 — Library sizes",
       x = NULL, y = "Mapped reads (M)") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(PLOT_DIR, "01_libsize.pdf"), p1, width = 10, height = 4)

# ── QC Plot 2: PCA across all stages ─────────────────────────────────────────
pca <- prcomp(t(logcpm), scale. = TRUE)
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  Stage    = meta$Stage_label,
  Genotype = meta$genotype,
  row.names = rownames(meta)
)

p2 <- ggplot(pca_df, aes(PC1, PC2, color = Stage, shape = Genotype)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = STAGE_COLORS) +
  scale_shape_manual(values = c(16, 15, 17, 18, 8)) +
  labs(title = "PCA — all samples, all stages",
       x = sprintf("PC1 (%s%%)", pve[1]), y = sprintf("PC2 (%s%%)", pve[2])) +
  theme_bw(base_size = 11)
ggsave(file.path(PLOT_DIR, "02_pca_stages.pdf"), p2, width = 7, height = 5)

# ── QC Plot 3: PCA within iMN_day28 only ─────────────────────────────────────
meta28  <- meta[meta$Stage_label == "iMN_day28", ]
lc28    <- logcpm[, rownames(meta28)]
pca28   <- prcomp(t(lc28), scale. = TRUE)
pve28   <- round(100 * pca28$sdev^2 / sum(pca28$sdev^2), 1)
df28    <- data.frame(PC1 = pca28$x[, 1], PC2 = pca28$x[, 2],
                      Genotype = meta28$genotype, Line = meta28$cell_line)

p3 <- ggplot(df28, aes(PC1, PC2, color = Genotype, label = Line)) +
  geom_point(size = 4, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_color_manual(values = GENOTYPE_COLORS) +
  labs(title = "PCA — iMN day 28 by genotype",
       x = sprintf("PC1 (%s%%)", pve28[1]), y = sprintf("PC2 (%s%%)", pve28[2])) +
  theme_bw(base_size = 11)
ggsave(file.path(PLOT_DIR, "03_pca_genotype_day28.pdf"), p3, width = 6, height = 5)

# ── Build and save feature objects per stage ──────────────────────────────────
build_features <- function(stage) {
  idx    <- meta$Stage_label == stage
  m      <- meta[idx, ]
  X      <- t(logcpm[, idx])   # samples × genes
  X      <- scale(X)            # z-score across samples per gene
  list(
    X        = X,
    y_binary = factor(m$condition, levels = c("Healthy", "ALS")),
    y_multi  = m$genotype,
    meta     = m
  )
}

feat28 <- build_features("iMN_day28")
feat50 <- build_features("iMN_day50")
saveRDS(feat28, file.path(DATA_DIR, "features_iMN_day28.RDS"))
saveRDS(feat50, file.path(DATA_DIR, "features_iMN_day50.RDS"))

cat("\niMN_day28 feature matrix:", nrow(feat28$X), "samples x", ncol(feat28$X), "genes\n")
cat("Class balance:\n"); print(table(feat28$y_binary))
cat("Genotype balance:\n"); print(table(feat28$y_multi))
cat("\nAll feature files saved to data/\n")
