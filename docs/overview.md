# Overview

The Nutkin software implements a linear-based variability metric and a corresponding
differential testing framework for quantifying and testing differences in overall transcriptional
cell-to-cell variability between groups of cells using single-cell RNA-seq data.

Nutkin is designed to help researchers gain a deeper understanding of cellular heterogeneity
in various biological and biomedical contexts by providing interpretable statistical insights
into overall transcriptional variations.

## Pipeline Summary

Nutkin consists of two main steps:

**Step 1 — Data Preprocessing**

The preprocessing module (`preprocessing.py`) loads single-cell RNA-seq data in `.h5ad` format
and automatically detects whether the data has already been preprocessed. If not, it runs a
standard Scanpy pipeline including filtering, normalization, log-transformation, highly variable
gene selection, and scaling. Optionally, it generates a PCA scree plot and a UMAP embedding.

**Step 2 — Variability Quantification and Differential Testing**

The main Nutkin class (`nutkin.py`) applies a PCA-based variability metric to each cell group,
producing a summary table of variability scores. It then performs pairwise bootstrap-based
differential tests between groups, reporting test statistics and p-values.

## Example Data

As example data, this manual uses scRNA-seq data from GEO (GSM7872694), which corresponds
to 10X Genomics 3' scRNA v3 of SPLINTR-barcoded MLL-AF9 cells cultured in vitro. The
processed dataset (`test_data.h5ad`) is provided in the `data/` directory of the repository.

To use the provided test data directly, proceed to [Preprocessing](preprocessing.md).
To generate `test_data.h5ad` from raw GEO files, see [Preparing the Test Data](prepare_data.md).
