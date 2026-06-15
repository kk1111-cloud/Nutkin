# Tutorial: Full Analysis Walkthrough
test
This tutorial walks through a complete Nutkin analysis using the provided test dataset
(`test_data.h5ad`). The same steps are available as a Jupyter notebook in `test/test.ipynb`.

## Test Dataset

The test dataset contains scRNA-seq data for 1,932 cells from four MLL-AF9 cell clones:

- `mCHERRY_Barcode_1614`
- `mCHERRY_Barcode_5774`
- `mCHERRY_Barcode_1755`
- `mCHERRY_Barcode_65372`

Group membership is stored in the `clone_id` column of `adata.obs`.

---

## Step 1: Load and Preprocess Data

```python
import sys
import os
import scanpy as sc

project_root = os.path.abspath("..")
sys.path.append(project_root)

import nutkin.preprocessing

data_path = "../Nutkin/data/test_data.h5ad"
output_path = "../Nutkin/test/result/"

adata = nutkin.preprocessing.load_data(
    data_path=data_path,
    output_path=output_path,
    n_comps=50,
    run_pca=True,
    run_embedding=True,
    group_col='clone_id'
)
```

Since `test_data.h5ad` is already preprocessed, Nutkin will skip the core preprocessing steps
and proceed directly to PCA and UMAP embedding.

**Output files:**

- `preprocessed_data.h5ad`
- `pca_scree_plot.png`
- `umapumap_plot.png`

---

## Step 2: Initialize Nutkin

```python
from nutkin.nutkin import Nutkin

adata = sc.read("../Nutkin/test/result/preprocessed_data.h5ad")

nut = Nutkin(
    adata=adata,
    group_col='clone_id',
    output_path=output_path,
    detail=True
)
```

---

## Step 3: Quantify Variability

```python
var_df = nut.measure()
print(var_df)
```

This saves `var_measurement.csv` to the output directory, listing the variability score for
each clone.

---

## Step 4: Differential Testing

### Compare specific group pairs

```python
nut.differential_test(
    group1=['mCHERRY_Barcode_1614', 'mCHERRY_Barcode_1614'],
    group2=['mCHERRY_Barcode_5774', 'mCHERRY_Barcode_1755']
)
```

### Automatic pairwise comparisons

```python
nut.differential_test()
```

Both calls save results to `differential_variability.csv`.

---

## Step 5: Visualize PCA Results

```python
nut.visualize_pca_results(detail=True, plot_type='both')
```

This generates:

- `pairwise_pca_both.png` — scatter + contour plots for each group pair and PC combination
- `summary_specific.csv` — per-PC variability table for each group

---

## Output Directory Summary

After the full tutorial, your `output_path` will contain:

| File | Description |
|---|---|
| `preprocessed_data.h5ad` | Preprocessed AnnData object |
| `pca_scree_plot.png` | PCA scree plot |
| `umapumap_plot.png` | UMAP embedding |
| `var_measurement.csv` | Variability scores per group |
| `differential_variability.csv` | Pairwise test results |
| `summary_specific.csv` | Per-PC variability summary |
| `pairwise_pca_scatter.png` | PCA scatter plots |
| `pairwise_pca_contour.png` | PCA contour plots |
| `pairwise_pca_both.png` | PCA scatter + contour overlaid |
