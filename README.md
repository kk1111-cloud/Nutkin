# Nutkin

Nutkin is a Python package for quantifying and testing differences in overall
transcriptional cell-to-cell variability between groups of cells using
single-cell RNA-seq data, based on a linear embedding–based statistical framework.

## Features

- Supports sum, product, standard deviation, median absolute deviation combination statistics
- Works for both balanced and unbalanced group designs
- Bootstrap-based inference with parallel computing support
- Designed for single-cell RNA-seq data analysis

## Installation

### Requirements

- Python >= 3.9
- NumPy
- ScanPy
- pandas

See `pyproject.toml` for the full list of dependencies.

### Install from PyPI

```bash
pip install nutkin
```

### Install from source

```bash
git clone https://github.com/kk1111-cloud/Nutkin.git
cd nutkin
pip install -e .
```

## Quick Start

```python
from nutkin import run_test

pval = run_test(
    data=expr_matrix,
    group=group_labels,
    method="sum",
    n_bootstrap=1000
)

print(pval)
```

## Usage

### Basic workflow

1. Prepare single-cell RNA-seq data
2. Specify group labels and analysis parameters
3. Run Nutkin
4. Interpret the output statistics

### Required input

- **scRNA-seq data**: Data can be provided either as a preprocessed or raw dataset. If raw data are supplied, Nutkin will perform preprocessing automatically.

### Running Nutkin from the command line

```bash
python main.py \
    --input_dataset <path_to_input> \
    --output_path <path_to_output> \
    --group_col <Metadata column(s) in adata.obs> \
    --group1 <group_name_1> \
    --group2 <group_name_2> \
    --detail
```

## Nutkin pipeline and outputs

Nutkin consists of two main steps.

### Step 1: Data preprocessing

The preprocessing module (`preprocessing.py`) performs the following:

- Reads the input data and checks its format
- If the input is an AnnData object:
  - Automatically detects whether preprocessing has been performed
  - If not, preprocessing and normalization are conducted using Scanpy
- If the input is not an AnnData object, an error is raised

**Input**
- Single-cell RNA-seq data (`.h5ad`)

**Output**
- A cleaned and processed AnnData object ready for Nutkin analysis

**Optional outputs**
- Scree plot (`--run_pca`): Shows the variance explained by each principal component
- Embedding plot (`--run_embedding`): UMAP, PHATE, or den-SNE visualizations for exploratory analysis

**Parameters**
| Parameter       | Type     | Default   | Description                                                                 |
|-----------------|----------|-----------|-----------------------------------------------------------------------------|
| `data_path`      | str      | —         | Path to the input single-cell RNA-seq `.h5ad` file.                        |
| `n_comps`        | int      | 50        | Number of principal components to compute during PCA.                      |
| `output_path`    | str/None | None      | Directory path to save optional outputs (scree plot, embedding plot).      |
| `run_pca`        | bool     | True      | Whether to run PCA and generate the scree plot.                            |
| `run_embedding`  | bool     | False     | Whether to generate a low-dimensional embedding (UMAP, PHATE, or den-SNE). |
| `method`         | str      | "umap"    | Embedding method to use if `run_embedding` is True. Only "umap" implemented.|

**Notes**
- The function automatically detects whether preprocessing has already been done by checking for `log1p` in `adata.uns` or `highly_variable` in `adata.var`.  
- If preprocessing is required, it performs filtering, normalization, log transformation, HVG selection, and scaling.  
- Scree plot and embedding plot are optional outputs controlled by `run_pca` and `run_embedding`.

### Step 2: Nutkin main analysis

**Input**
- Preprocessed AnnData object generated in Step 1

**Outputs**
- `var_measurement.csv`: Quantified variability for each defined group
- `differential_variability.csv`: Pairwise statistical test results between specified groups or all groups

**Optional outputs (`--detail`)**
- `summary_specific.csv`: Variability measures for individual principal components
- `all_pairwise_scatter.png`, `all_pairwise_kde.png`: PCA scatter plots and contour plots (default: top 3 PCs)

**Parameters**
| Parameter            | Type     | Default           | Description                                                                                      |
|----------------------|----------|-----------------|--------------------------------------------------------------------------------------------------|
| `adata`              | AnnData  | —               | Preprocessed single-cell RNA-seq data from Step 1                                                |
| `num_pc`             | int      | DEFAULT_NUM_PC  | Number of principal components to use in PCA                                                     |
| `group_col`          | str      | "group"         | Column in `adata.obs` specifying group labels                                                   |
| `metric_type`        | str      | "sum"           | How to summarize variability: `"sum"`, `"product"`, or `"both"`                                 |
| `variability_method` | str      | "mad"           | Method to compute variability: `"sd"`, `"mad"`, or `"both"`                                     |
| `verbal`             | bool     | True            | Whether to print progress messages                                                               |
| `detail`             | bool     | False           | Whether to generate detailed PCA outputs and visualizations                                      |

**Notes**
- `measure()` computes variability metrics for each group and saves `var_measurement.csv`.  
- `differential_test()` performs pairwise statistical comparisons between groups and saves `differential_variability.csv`.  
- Bootstrapping is used internally for hypothesis testing.  
- PCA-based visualization can be generated using `visualize_pca_results(detail=True)`.  
- Optional outputs (`summary_specific.csv` and pairwise PCA plots) are only generated if `detail=True`.  
- The class automatically handles combinations of `metric_type` and `variability_method`, e.g., `"sum"` with `"sd"`, `"both"` with `"both"`, etc.

## Documentation

Full documentation will be available at:

https://nutkin.readthedocs.io

The documentation includes:
- Installation guide
- Detailed usage examples
- API reference
- Methodological background

## License

MIT License

