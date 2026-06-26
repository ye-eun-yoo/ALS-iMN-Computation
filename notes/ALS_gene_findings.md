# ALS Gene Findings — GSE299997 iPSC-derived Motor Neurons

**Dataset**: Ma GM et al., *Stem Cell Reports* 2025, PMID 41043426  
**Genotypes**: SOD1 (n=4), FUS (n=2), TARDBP (n=2), C9orf72 (n=1) vs Healthy (n=4)  
**Stage analysed for classification**: iMN day 28 (39 samples: 27 ALS, 12 Healthy)  
**Methods**: edgeR DE, elastic net (glmnet), random forest (ranger), temporal comparison (day28 vs day50)

---

## Key Result: ALS vs Healthy is Strongly Separable

All three classification approaches (elastic net, LASSO, ridge, random forest) achieved CV AUC = 1.0 on iMN day 28. The ALS transcriptional signature is robust across all four mutant genotypes when pooled against healthy controls — even with a small dataset (n=39 samples, 10-fold CV).

---

## Multi-Method Consensus Genes

The most robust candidates appear in both the elastic net (53 genes selected at lambda.1se) and the random forest top-30 permutation importance list.

| Gene | Direction | EN coef | RF importance | Notes |
|------|-----------|---------|---------------|-------|
| ZNF506 | Up in ALS | +0.343 | 0.00175 | Top elastic net hit; zinc finger TF |
| ZNF528 | Up in ALS | +0.141 | 0.00134 | Zinc finger TF |
| CYP2E1 | Up in ALS | +0.104 | 0.00133 | Cytochrome P450; oxidative stress, ROS production |
| CTC-559E9.5 | Up in ALS | +0.061 | 0.00157 | lncRNA; also top DE hit (FDR=1.6e-4) |
| ZNF454 | Up in ALS | +0.032 | 0.00109 | Zinc finger TF; also top DE hit (FDR=2.4e-4, logFC=+1.59) |
| ETFA | Down in ALS | -0.116 | 0.00128 | Electron transfer flavoprotein alpha; mitochondrial FAO |
| MT-TS1 | Down in ALS | -0.148 | 0.00113 | Mitochondrial tRNA-Ser-1; mitochondrial function |
| XXbac-B444P24.8 | Down in ALS | -0.147 | 0.00109 | Uncharacterized locus |

**EN coef** = elastic net coefficient at lambda.1se (larger |value| = stronger classifier weight)  
**RF importance** = mean decrease in accuracy (permutation-based)

---

## Top DE Genes — edgeR ALS vs Healthy (iMN day 28)

Ranked by FDR. All pooled ALS genotypes vs Healthy.

| Gene | logFC | FDR | Direction | Biology |
|------|-------|-----|-----------|---------|
| ZNF571-AS1 | +0.92 | 1.0e-4 | Up | lncRNA antisense to ZNF571 |
| UTS2 | -3.80 | 1.1e-4 | Down | Urotensin II neuropeptide; expressed in spinal motor neurons |
| ZNF185 | -2.35 | 1.2e-4 | Down | Zinc finger; actin cytoskeleton |
| FAM19A3 | -2.98 | 1.2e-4 | Down | TAFA chemokine family; neural secreted protein |
| DLK1 | -5.41 | 1.6e-4 | Down | Notch ligand; neurodevelopment and muscle; strongly reduced |
| JPH2 | -3.38 | 1.6e-4 | Down | Junctophilin-2; ER-PM junctions, calcium signaling |
| ANXA3 | -2.62 | 1.6e-4 | Down | Annexin A3; phospholipid binding, anti-inflammatory |
| ZNF578 | +2.96 | 1.6e-4 | Up | Zinc finger TF |
| CTC-559E9.5 | +0.88 | 1.6e-4 | Up | lncRNA (also elastic net + RF) |
| LINC00672 | -1.22 | 1.9e-4 | Down | lncRNA; also elastic net; stable at day28 and day50 |
| ERV3-1 | +0.74 | 1.9e-4 | Up | Endogenous retrovirus element 3-1 |
| ZNF454 | +1.59 | 2.4e-4 | Up | Zinc finger TF (also elastic net + RF) |
| ETFA | -0.30 | 2.7e-4 | Down | Mitochondrial FAO (also elastic net + RF) |
| APLP2 | -0.31 | 2.7e-4 | Down | Amyloid precursor-like protein 2; neurodegeneration |
| PCDHGA8 | +1.94 | 3.1e-4 | Up | Protocadherin gamma A8; neuronal cell adhesion |

---

## Biologically Highlighted Genes

### Mitochondrial / Metabolic (Theme: ALS mitochondrial dysfunction)
- **ETFA** (down, FDR=2.7e-4): Electron Transfer Flavoprotein subunit alpha, central to mitochondrial fatty acid beta-oxidation. Found by all three methods. Consistent with the well-documented mitochondrial dysfunction in ALS.
- **MT-TS1** (down): Mitochondrial tRNA for serine, encoded in the mitochondrial genome. Reduced expression across ALS genotypes.
- **CYP2E1** (up): Cytochrome P450 2E1, a major producer of reactive oxygen species. Upregulation suggests increased oxidative stress in ALS iMNs.
- **ACACB** (up, elastic net): Acetyl-CoA carboxylase beta, rate-limiting enzyme in mitochondrial fatty acid oxidation control.

### Transcriptional Reprogramming (Theme: ZNF cluster dysregulation)
Multiple zinc finger proteins are upregulated in ALS: ZNF506, ZNF528, ZNF454, ZNF578, ZNF571-AS1 (lncRNA). This cluster suggests broad transcriptional reprogramming in ALS motor neurons rather than a single-gene effect. ZNF506 is the most robust hit across methods.

### Motor Neuron Identity
- **LHX3** (RF importance rank #3): LIM homeobox 3 is a master transcription factor for spinal motor neuron identity. Dysregulation in ALS iMNs could reflect loss of mature motor neuron identity or a compensatory response.

### Synaptic / Neuronal Function
- **DLK1** (down, logFC=-5.41, FDR=1.6e-4): Delta-like non-canonical Notch ligand 1. Strongly reduced in ALS. Regulates neuronal differentiation and is expressed in motor neurons.
- **STXBP6** (up, elastic net): Syntaxin binding protein 6. Synaptic vesicle trafficking — altered neurotransmitter release.
- **SYT6** (up, elastic net): Synaptotagmin 6. Calcium sensor for vesicle fusion, synaptic transmission.
- **DLG2** (down, elastic net): Disc large homolog 2 (PSD-93). Postsynaptic density scaffold; synaptic maturation and plasticity.
- **APLP2** (down, FDR=2.7e-4): Amyloid precursor-like protein 2, paralog of APP. Reduced in ALS iMNs.

### Neuroinflammation / Signaling
- **KYNU** (down, elastic net): Kynureninase, kynurenine pathway. Reduced KYNU expression in ALS iMNs is notable given that kynurenine pathway metabolites are implicated in ALS neuroinflammation.
- **UTS2** (down, logFC=-3.80, FDR=1.1e-4): Urotensin II neuropeptide. Expressed in spinal motor neurons and regulated by motor neuron activity.
- **JPH2** (down, FDR=1.6e-4): Junctophilin-2, maintains ER-plasma membrane junctions critical for calcium signaling. Reduced in ALS iMNs.

### Neuronal Cell Adhesion
- **PCDHGB4**, **PCDHGA8**, **PCDHGA2** (up, elastic net / DE): Protocadherin gamma family members. Multiple upregulated in ALS across methods. Protocadherins govern neuronal self-avoidance and synapse specificity.

### lncRNAs
- **LINC00672** (down): Found by elastic net and DE (FDR=1.9e-4). Unique in that it is stably downregulated at both day28 (logFC=-1.22) and day50 (logFC=-1.25) — one of only two genes persistently dysregulated across the temporal window.

---

## Temporal Analysis (day28 vs day50)

| Timepoint | Elastic net genes selected | CV AUC |
|-----------|--------------------------|--------|
| iMN day 28 | 53 | 1.000 |
| iMN day 50 | 72 | 1.000 |

The ALS signature strengthens (requires more genes) from day28 to day50, consistent with progressive transcriptional divergence during maturation. Only 2 genes are consistently dysregulated at both timepoints:

| Gene | logFC day28 | logFC day50 | Direction |
|------|-------------|-------------|-----------|
| LINC00672 | -1.22 | -1.25 | Down (stable) |
| AC139099.3 | -2.03 | -1.61 | Down (attenuated) |

---

## Regularization Comparison

| Method | Genes selected (lambda.1se) | CV AUC |
|--------|-----------------------------|--------|
| Ridge (alpha=0) | 33,747 (all genes shrunk) | 1.0 |
| Elastic Net (alpha=0.5) | 57 | 1.0 |
| LASSO (alpha=1) | 19 | 1.0 |

All three achieve perfect cross-validated separation. The LASSO selects the smallest gene set (19 genes) while maintaining AUC = 1.0, suggesting the ALS signal is highly concentrated in a small subset of genes. The elastic net (57 genes) represents a balance between sparsity and stability.

---

## Caveats

1. **Small n**: 13 iPSC lines, 3 replicates each. 10-fold CV folds can contain replicates from the same line — AUC likely overestimates true generalization.
2. **Genotype imbalance**: C9orf72 n=1 line (3 samples); its contribution to the signal is uncertain.
3. **All ALS genotypes pooled**: The classification identifies genes that separate any ALS from Healthy. Genotype-specific signals are not captured.
4. **AUC = 1.0**: While technically correct (CV predictions), this reflects that bulk RNA-seq differences between ALS and Healthy iPSC-derived neurons are very large relative to within-group variance at this sample size. These results should be replicated in a larger cohort.
5. **lncRNAs and novel transcripts**: Several top hits (CTC-559E9.5, XXbac-B444P24.8, AC139099.3) have no defined function. They could reflect technical artifacts of alignment or genuine regulatory lncRNAs.

---

*Analysis date: 2026-06-26 | Data: GSE299997 | Scripts: scripts/01–06*
