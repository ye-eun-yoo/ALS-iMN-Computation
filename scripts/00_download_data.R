## GSE299997 - Download and prepare count data from GEO
##
## Dataset: "Integrated temporal profiling of iPSC-derived motor neurons from
##   ALS patients carrying C9orf72, FUS, TARDBP, and SOD1 mutations"
##   GEO: GSE299997 | License: CC BY 4.0
##
## Run this script once to download raw counts from GEO.
## Outputs saved to data/:
##   GSE299997_DGEList.RDS  - filtered/TMM-normalized DGEList
##   GSE299997_metadata.csv - sample metadata with genotype/stage labels

REPO_DIR <- Sys.getenv("REPO_DIR", unset = ".")
DATA_DIR <- file.path(REPO_DIR, "data")
SUPP_DIR <- file.path(DATA_DIR, "supp")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

# Install packages to a local lib if not available
local_lib <- file.path(REPO_DIR, "rlib")
dir.create(local_lib, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("GEOquery", "edgeR")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, lib = local_lib, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({ library(GEOquery); library(edgeR) })

GEO_ID     <- "GSE299997"
COUNT_FILE <- file.path(SUPP_DIR, "GSE299997_CountfileMND.txt.gz")

cat("=== Fetching GEO metadata ===\n")
gse     <- getGEO(GEO_ID, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
meta_raw <- pData(gse)

char_cols <- grep("^characteristics_ch1", colnames(meta_raw), value = TRUE)
parsed    <- lapply(char_cols, function(col) {
  vals <- as.character(meta_raw[[col]])
  key  <- sub(":.*", "", vals[1])
  val  <- sub("^[^:]+:\\s*", "", vals)
  setNames(data.frame(val, stringsAsFactors = FALSE), gsub("\\s+", "_", key))
})
meta_geo <- cbind(
  meta_raw[, c("geo_accession", "title", "source_name_ch1")],
  do.call(cbind, parsed)
)
rownames(meta_geo) <- meta_geo$title

if (!file.exists(COUNT_FILE)) {
  cat("=== Downloading count file ===\n")
  getGEOSuppFiles(GEO_ID, makeDirectory = FALSE, baseDir = SUPP_DIR)
}

cat("=== Reading count matrix ===\n")
counts_raw <- read.table(gzfile(COUNT_FILE), header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE, comment.char = "#")
ANNO_COLS  <- c("Chr", "Start", "End", "Strand", "Length")
counts_mat <- counts_raw[, setdiff(colnames(counts_raw), ANNO_COLS), drop = FALSE]
colnames(counts_mat) <- sub("_hisat2\\.bam$|^\\.bam$", "", colnames(counts_mat))

sample_df <- data.frame(
  sample_id  = colnames(counts_mat),
  patient_id = sub("^(.*)-([^-]+)-([0-9]+)$", "\\1", colnames(counts_mat)),
  Stage      = sub("^(.*)-([^-]+)-([0-9]+)$", "\\2", colnames(counts_mat)),
  Replicate  = as.integer(sub("^(.*)-([^-]+)-([0-9]+)$", "\\3", colnames(counts_mat))),
  stringsAsFactors = FALSE, row.names = colnames(counts_mat)
)
stage_map <- c(IPS = "iPSC", NPC = "NPC", MND28 = "iMN_day28", MND50 = "iMN_day50")
sample_df$Stage_label <- stage_map[sample_df$Stage]

shared <- intersect(colnames(counts_mat), rownames(meta_geo))
counts_mat <- counts_mat[, shared, drop = FALSE]
sample_df  <- cbind(sample_df[shared, ], meta_geo[shared, c("geo_accession", "cell_line")])

# Parse genotype from cell_line ("SOD1-2" -> "SOD1")
sample_df$genotype <- sub("-[0-9]+$", "", sample_df$cell_line)
sample_df$condition <- ifelse(sample_df$genotype == "Healthy", "Healthy", "ALS")
sample_df$genotype  <- factor(sample_df$genotype,
                               levels = c("Healthy", "SOD1", "FUS", "TARDBP", "C9orf72"))

dge  <- DGEList(counts = counts_mat, samples = sample_df)
dge  <- calcNormFactors(dge, method = "TMM")
keep <- filterByExpr(dge, group = sample_df$condition)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")

saveRDS(dge, file.path(DATA_DIR, "GSE299997_DGEList.RDS"))
write.csv(dge$samples, file.path(DATA_DIR, "GSE299997_metadata.csv"))
cat(sprintf("\nSaved: %d genes x %d samples\n", nrow(dge), ncol(dge)))
print(table(dge$samples$Stage_label, dge$samples$genotype))
