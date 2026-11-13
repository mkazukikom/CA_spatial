# CA_spatial

R code for figure generation and statistical analysis in:

> Matsuda KM, Chen Y-Y, Morita A, Kotani H, Toyama S, Omatsu J, Lan C-CE, Yu H-S, Shibata S.
> **Spatial transcriptomics reveals shared artery-specific interferon program and anti-CD320 antibody–defined molecular subsets in cutaneous arteritis.**

## Overview

The study integrates three transcriptomic datasets:

| Dataset | Description | Samples | DDBJ accession |
|---|---|---|---|
| Human GeoMx DSP | Whole Transcriptome Atlas spatial profiling of FFPE skin biopsies from patients with cutaneous arteritis (CA) | 6 patients (3 anti-CD320 seropositive, 3 seronegative) × 4 ROI types = 24 AOIs | DRR1271807–DRR1271830 |
| Mouse passive transfer | Bulk RNA-seq of back skin after subcutaneous injection of patient-derived IgG (anti-CD320+ vs anti-CD320−) in WT, *Cd320*-KO and WT + methylcobalamin mice | n = 3 per condition per group | DRR1085606–DRR1085623 |
| HDBEC in vitro | Bulk RNA-seq of human dermal blood endothelial cells treated with anti-CD320 mAb or isotype control for 6 h | n = 2 per condition | DRR1085624–DRR1085627 |

## Repository structure

Each folder contains one standalone R script that reproduces the corresponding figure panel.

| Script | Panel | Analysis |
|---|---|---|
| `Figure1/Figure1B` | Fig. 1B | PCA of GeoMx AOIs (`prcomp` on log2 Q3-normalized expression) |
| `Figure1/Figure1C` | Fig. 1C | Canonical B-cell marker expression across ROI types |
| `Figure1/Figure1D` | Fig. 1D | Volcano plot, artery vs vessel (limma-voom results) |
| `Figure1/Figure1E` | Fig. 1E | Hallmark GSEA, artery vs vessel (`fgsea`) |
| `Figure2/Figure2C` | Fig. 2C | Volcano plots, anti-CD320+ vs anti-CD320− IgG in WT, *Cd320*-KO and WT + MetCbl mice |
| `Figure2/Figure2D` | Fig. 2D | GO Biological Process GSEA, WT mice (`clusterProfiler::gseGO`) |
| `Figure3/Figure3B` | Fig. 3B | Cross-species log2FC comparison of the top 10 conserved IFN core genes |
| `Figure3/Figure3C` | Fig. 3C | Hierarchical clustering of the 56-gene conserved IFN core across all AOIs |
| `Figure4/Figure4B` | Fig. 4B | Volcano plot, HDBEC anti-CD320 mAb vs isotype control |
| `Figure4/Figure4C` | Fig. 4C | Hallmark GSEA, HDBEC (`fgsea`) |
| `Figure4/Figure4D` | Fig. 4D | cGAS–STING and RIG-I/MDA5 sensor gene expression in HDBECs |
| `Figure5/Figure5A` | Fig. 5A | Volcano plot, seropositive vs seronegative within artery ROIs, coloured by gene family |
| `Figure5/Figure5B` | Fig. 5B | Representative histone and IFN core genes by anti-CD320 status |
| `Figure5/Figure5C` | Fig. 5C | Direction-of-effect binomial sign test |
| `FigureS2/FigureS2` | Supplementary Fig. 2 | Human–mouse CD320 global alignment and dot plot |
| `FigureS3/FigureS3` | Supplementary Fig. 3 | Cross-species log2FC comparison of the remaining conserved IFN core genes |
| `FigureS4S5/FigureS4S5` | Supplementary Figs. 4–5 | Additional histone and IFN core genes by anti-CD320 status |

Supplementary Figure 1 (histology) has no associated code.

## Key algorithms

**Conserved IFN core gene set** (used in Figs. 3, 5 and Supplementary Figs. 3–5)

```
human_LE <- leading-edge genes of HALLMARK_INTERFERON_GAMMA_RESPONSE ∪ HALLMARK_INTERFERON_ALPHA_RESPONSE
            from the GeoMx artery-vs-vessel GSEA                                 # n = 125
mouse_LE <- the same, from the WT anti-CD320+ vs anti-CD320− IgG GSEA, toupper()  # n = 84
IFN_core <- intersect(human_LE, mouse_LE)                                         # n = 56
```

**Direction-of-effect sign test** (`Figure5/Figure5C`)

```
for gene_set in {IFN_core (n = 56), histone genes (n = 73; HGNC histone nomenclature)}:
  for region in {artery, vein/vasculature}:
    d[g] <- mean(log2 expr | seropositive AOIs) − mean(log2 expr | seronegative AOIs)
    k    <- number of genes with d[g] > 0
    p    <- binom.test(k, length(gene_set), p = 0.5)   # two-sided exact test
```

**CD320 sequence comparison** (`FigureS2/FigureS2`)

```
alignment <- Needleman–Wunsch global alignment of human Q9NPF0 vs mouse Q9Z1P5
             (match +1, mismatch −1, gap −1)
identity  <- identical positions / alignment length
dot plot  <- every (i, j) at which ≥ 5 of 9 residues in a sliding 9-residue window are identical
epitope   <- human aa 189–198 mapped onto the alignment
```

Thresholds used throughout: Benjamini–Hochberg adjusted *P* < 0.05 and |log2FC| > 1 for DEGs. In volcano plots, the y-axis shows nominal −log10 *P* and the dashed line marks nominal *P* = 0.05. For HDBECs (n = 2), points are coloured by nominal *P* < 0.05. GSEA ranks genes by log2FC × −log10(*P*) and uses 10,000 permutations.

## Input files

The scripts read the following files from the working directory:

| File | Content | Used by |
|---|---|---|
| `02_se_explored.rds` | `SummarizedExperiment` of Q3-normalized GeoMx expression (11,148 genes × 24 AOIs) with AOI metadata (patient, ROI type, anti-CD320 status) | 1B, 1C, 3C, 5B, 5C, S4S5 |
| `DGE_artery_combined_vs_vessel.csv` | limma-voom results, artery vs vessel | 1D, 1E, 3B, S3 |
| `DGE_antibody_pos_vs_neg_artery.csv` | limma-voom results, seropositive vs seronegative within artery ROIs | 5A |
| `DEG_BF_CA_vs_HC.csv` | Mouse DEG, WT, anti-CD320+ vs anti-CD320− IgG | 2C, 2D, 3B, S3 |
| `DEG_ko_cd320_vs_ko_ctrl.csv` | Mouse DEG, *Cd320*-KO | 2C |
| `DEG_wt_cd320metcbl_vs_ctrl_metcbl.csv` | Mouse DEG, WT + MetCbl | 2C |
| `DEG_HDBEC_ab_vs_ig.csv` | HDBEC DEG, anti-CD320 mAb vs isotype control | 4B, 4C, 4D |
| `IFN_core_genes.rds` | Character vector of the 56 conserved IFN core genes | 3C, 5A, 5C, S3, S4S5 |
| `h.all.v2023.2.Hs.entrez.gmt` | MSigDB Hallmark gene sets (download from https://www.gsea-msigdb.org) | 1E, 4C |
| `Q9NPF0.fasta`, `Q9Z1P5.fasta` | Human and mouse CD320 sequences (download from https://www.uniprot.org) | S2 |

Upstream processing: GeoMx counts were Q3-normalized with the GeoMx NGS pipeline (Bruker Spatial Biology). Mouse reads were aligned with STAR (GRCm38, GENCODE M18). HDBEC reads were aligned with HISAT2 v2.1.0 (GRCh38) and quantified with StringTie v2.2.1 on the Subio Platform v1.24.5853. Differential expression for GeoMx used limma-voom with patient as a random effect.

## Requirements

R 4.4.3 with the following packages:

```r
install.packages(c("BiocManager", "dplyr", "tidyr", "tibble", "stringr",
                   "ggplot2", "ggrepel", "patchwork", "pheatmap", "RColorBrewer"))
BiocManager::install(c("SummarizedExperiment", "AnnotationDbi", "fgsea",
                       "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db"))
```

## Usage

Place the required input files in the same folder as the script, then run it:

```bash
git clone https://github.com/mkazukikom/CA_spatial.git
cd CA_spatial/Figure1
Rscript Figure1B
```

Output panels (PDF/PNG) and result tables (CSV) are written to the working directory.

## License

Released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the article above. Citation details will be updated upon publication.

## Contact

Kazuki M. Matsuda, Department of Dermatology, The University of Tokyo Hospital.
For questions about the code, please open a GitHub issue.
