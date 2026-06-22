# Step 2: Nutkin Analysis

After preprocessing, the main Nutkin analysis is performed using the `Nutkin` class in
`nutkin/nutkin.py`. This class handles variability quantification, differential testing,
and optional PCA visualization.

## Initializing Nutkin

Load the preprocessed data and create a `Nutkin` instance:

```python
import scanpy as sc
from nutkin.nutkin import Nutkin

adata = sc.read("../test/result/preprocessed_data.h5ad")

output_path = "../Nutkin/test/result/"

nut = Nutkin(
    adata=adata,
    group_col='clone_id',
    output_path=output_path,
    detail=True
)
```

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adata` | AnnData | — | Preprocessed single-cell RNA-seq data |
| `num_pc` | int | 15 | Number of PCs to use in variability computation |
| `group_col` | str | `"group"` | Column in `adata.obs` specifying group labels |
| `group1` | str or list | — | Name(s) of the first group(s) to compare |
| `group2` | str or list | — | Name(s) of the second group(s) to compare |
| `metric_type` | str | `"sum"` | How to summarize variability: `"sum"`, `"product"`, or `"both"` |
| `variability_method` | str | `"mad"` | Method to compute variability: `"sd"`, `"mad"`, or `"both"` |
| `verbal` | bool | True | Whether to print progress messages |
| `detail` | bool | False | Whether to generate detailed PCA outputs and visualizations |
| `plot_pc_num` | int | 3 | Number of PCs for visualization |
| `plot_type` | str | `"scatter"` | PCA plot type: `"scatter"`, `"contour"`, or `"both"` |

---

## 2.1 Quantify Variability per Group

Use `measure()` to compute the variability metric for each cell group:

```python
var_df = nut.measure()
# Quantifying variability using MAD method
# Running dimensional reduction (PCA)...
# Saving results
```

The result is saved as `var_measurement.csv` in the output directory. Each row contains:

- Group name
- Number of cells in the group
- Computed variability statistic(s)

```{figure} output/measure.png
:alt: var_measurement.csv
:width: 500px

Example the variability metric for each cell group.
```

---

## 2.2 Differential Testing Between Groups

Use `differential_test()` to perform pairwise bootstrap-based statistical comparisons.
The null hypothesis assumes no difference in variability between the two groups.
A two-sided p-value below 0.05 is considered statistically significant.

### Comparing specific groups

Pass explicit lists to `group1` and `group2`. Each element in `group1` is compared to
the element at the same index in `group2`:

```python
nut.differential_test(
    group1=['mCHERRY_Barcode_1614', 'mCHERRY_Barcode_1614'],
    group2=['mCHERRY_Barcode_5774', 'mCHERRY_Barcode_1755']
)
# Conducting differential analysis using MAD method
# Comparing mCHERRY_Barcode_1614 and mCHERRY_Barcode_5774
# Running Bootstrapping
# Comparing mCHERRY_Barcode_1614 and mCHERRY_Barcode_1755
# Running Bootstrapping
# Saving results
```
```{figure} output/compare1.png
:alt: differential_variability.csv
:width: 1000px

Example differential analysis results.
```

### Automatic pairwise comparisons

If `group1` and `group2` are both `None`, Nutkin automatically compares adjacent groups
in the order of variability from `var_measurement.csv`:

```python
nut.differential_test()
# Conducting differential analysis using MAD method
# Comparing mCHERRY_Barcode_1614 and mCHERRY_Barcode_65372
# Running Bootstrapping
# Comparing mCHERRY_Barcode_65372 and mCHERRY_Barcode_1755
# Running Bootstrapping
# Comparing mCHERRY_Barcode_1755 and mCHERRY_Barcode_5774
# Running Bootstrapping
# Saving results
```

The output is saved as `differential_variability.csv`. Each row contains:

| Column | Description |
|---|---|
| `groupA` | Name of the first group |
| `number of cells in groupA` | Cell count for groupA |
| `sum/product statistics for groupA` | Variability metric for groupA |
| `groupB` | Name of the second group |
| `number of cells in groupB` | Cell count for groupB |
| `sum/product statistics for groupB` | Variability metric for groupB |
| `null hypothesis` | The null hypothesis tested (`A = B`, `A < B`, or `A > B`) |
| `test statistics` | Observed difference in variability (groupA − groupB) |
| `p value` | Two-sided bootstrap p-value |

```{figure} output/compare2.png
:alt: differential_variability.csv
:width: 1000px

Example differential analysis results.
```
---

## 2.3 Visualize PCA Results (Optional)

When `detail=True`, call `visualize_pca_results()` to generate PCA scatter and/or contour
plots for each pair of compared groups:

```python
nut.visualize_pca_results(detail=True, plot_type='both')
# PCA visualization (both) saved to ../Nutkin/test/result
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `detail` | False | Must be `True` to generate plots |
| `plot_pc_num` | 3 | Number of PCs to show (plots all pairwise combinations) |
| `plot_type` | `"scatter"` | `"scatter"`, `"contour"`, or `"both"` |

### Output Files (when `detail=True`)

| File | Description |
|---|---|
| `summary_specific.csv` | Per-group variability for each PC (PC1–PC15) and overall sum |
| `pairwise_pca_scatter.png` | Scatter plots in pairwise PC space |
| `pairwise_pca_contour.png` | KDE contour plots in pairwise PC space |
| `pairwise_pca_both.png` | Scatter + contour overlaid |


```{figure} output/summary.png
:alt: summary_specific.csv
:width: 1000px

Example Per-group variability for each PC (PC1–PC15) and overall statistics.
```
```{figure} output/pairwise_pca_both.png
:alt: pairwise_pca_both.png
:width: 800px

Example Scatter and KDE contour plots in pairwise PC space.
```
