# ALS iPSC Motor Neuron — ML Gene Discovery

Supervised machine learning pipeline for unbiased gene discovery in ALS, applied to a public iPSC-derived motor neuron RNA-seq dataset (GSE299997, CC BY 4.0).

## Dataset

**GSE299997** — *Integrated temporal profiling of iPSC-derived motor neurons from ALS patients*  
4 ALS genotypes: **SOD1**, **FUS**, **TARDBP**, **C9orf72**  
4 differentiation stages: iPSC → NPC → iMN day 28 → iMN day 50  
13 iPSC lines × 3 replicates = ~157 bulk RNA-seq samples

## Workflows

| Script | Method | Question |
|--------|--------|----------|
| `01_prepare_features.R` | TMM normalization, PCA | QC and feature matrix construction |
| `02_elastic_net_binary.R` | glmnet, binomial, 10-fold CV | Which genes classify ALS vs Healthy? |
| `03_regularization_compare.R` | Ridge / Elastic Net / LASSO | How does regularization type affect gene selection? |
| `04_random_forest.R` | ranger, permutation importance | Multi-class genotype prediction; train day28 → test day50 |
| `05_de_enet_pipeline.R` | edgeR + glmnet | Does DE pre-filtering improve elastic net AUC? |
| `06_temporal_analysis.R` | glmnet + edgeR per timepoint | Does the ALS signature strengthen during differentiation? |

## Quick start

```r
# 1. Download data from GEO (requires internet + GEOquery)
REPO_DIR="." Rscript scripts/00_download_data.R

# 2. Build feature matrices
REPO_DIR="." Rscript scripts/01_prepare_features.R

# 3. Run analyses
REPO_DIR="." Rscript scripts/02_elastic_net_binary.R
REPO_DIR="." Rscript scripts/03_regularization_compare.R
REPO_DIR="." Rscript scripts/04_random_forest.R
REPO_DIR="." Rscript scripts/05_de_enet_pipeline.R
REPO_DIR="." Rscript scripts/06_temporal_analysis.R
```

## Dependencies

R packages: `edgeR`, `glmnet`, `ranger`, `pROC`, `ggplot2`, `ggrepel`, `dplyr`, `tidyr`  
All available from Bioconductor 3.20 / CRAN.

## Key outputs

- `results/02_elastic_net_binary_coefficients.csv` — genes selected for ALS classification
- `results/04_rf_gene_importance.csv` — permutation importance ranking (~17k genes)
- `results/05_de_enet_auc_comparison.csv` — AUC: full genome vs DE-filtered models
- `plots/*/` — publication-quality PDF figures per analysis

## Data citation

Ma GM, Xia CC, Lyu BY, Liu J et al. Integrated profiling of iPSC-derived motor neurons carrying C9orf72, FUS, TARDBP, or SOD1 mutations. *Stem Cell Reports* 2025 Oct 14;20(10):102649. PMID: 41043426.  
GEO: [GSE299997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299997) | License: CC BY 4.0
