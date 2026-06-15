# Step 1: Data Preprocessing

Nutkin takes single-cell RNA-seq data in the form of an AnnData object (`.h5ad` format) as
input. The preprocessing module `nutkin/preprocessing.py` provides a complete pipeline for
loading and preparing scRNA-seq data, implemented through the `load_data` function.

## What `load_data` Does

`load_data` starts by loading the `.h5ad` file and automatically detects whether the data
has already been preprocessed by checking for `log1p` in `adata.uns` or `highly_variable`
in `adata.var.columns`. Depending on the result:

- **Raw data**: executes a core preprocessing routine (see below)
- **Preprocessed data**: skips preprocessing and proceeds directly to PCA

After preprocessing, `load_data` optionally runs PCA (storing results in the AnnData object
and saving a scree plot) and generates low-dimensional embeddings (UMAP).

## Core Preprocessing Steps (Raw Data Only)

1. Filter low-quality cells (min 200 genes) and genes (min 3 cells)
2. Normalize total counts per cell to 10,000
3. Log-transform (`log1p`)
4. Select highly variable genes
5. Subset to highly variable genes
6. Scale gene expression values (max value = 10)

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_path` | str | — | Path to the input `.h5ad` file |
| `n_comps` | int | 50 | Number of PCA components to compute |
| `output_path` | str / None | None | Directory to save output files |
| `run_pca` | bool | True | Whether to run PCA and generate the scree plot |
| `run_embedding` | bool | False | Whether to generate a UMAP embedding |
| `group_col` | str | `"group"` | Column in `adata.obs` specifying group labels |
| `embedding_method` | str | `"umap"` | Embedding method (`"umap"` currently supported) |

## Example

Open a Python session and import the preprocessing module:

```python
import sys
import os
import scanpy as sc

# Append the parent directory of Nutkin to the Python path
project_root = os.path.abspath("..")
sys.path.append(project_root)

import nutkin.preprocessing
```

Read the test data:

```python
data_path = "../Nutkin/data/test_data.h5ad"
adata = sc.read_h5ad(data_path)
print(adata)
# AnnData object with n_obs × n_vars = 1932 × 31053
#   obs: 'barcode', 'bc.umi.count', 'num_barcodes', 'max_umi',
#        'lineage_barcodes_filtered', 'umi_counts_filtered',
#        'clone_id', 'louvain', 'clone_plot_group'
#   var: 'gene_ids'
#   uns: 'log1p', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'umap'
#   obsm: 'X_pca', 'X_umap'
#   varm: 'PCs'
#   obsp: 'connectivities', 'distances'
```

Run `load_data` with PCA and UMAP embedding enabled, grouping by `clone_id`:

```python
output_path = "../Nutkin/test/result/"

adata = nutkin.preprocessing.load_data(
    data_path=data_path,
    output_path=output_path,
    n_comps=50,
    run_pca=True,
    run_embedding=True,
    group_col='clone_id'
)
# Successfully loaded data from: ../Nutkin/data/test_data.h5ad
# Data detected as PREPROCESSED. Skipping core preprocessing steps.
# Generating UMAP embedding...
```

## Output Files

After running `load_data`, the following files are saved to `output_path`:

| File | Description |
|---|---|
| `preprocessed_data.h5ad` | Processed AnnData object |
| `pca_scree_plot.png` | PCA scree plot (if `run_pca=True`) |
| `umapumap_plot.png` | UMAP embedding colored by `group_col` (if `run_embedding=True`) |

```{figure} _static/pca_scree_plot.png
:alt: PCA Scree Plot
:width: 500px

Example PCA scree plot showing variance explained by each principal component.
```
