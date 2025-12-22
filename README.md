# Nutkin

Nutkin is a Python package for quantifying and testing differences in overall
transcriptional cell-to-cell variability between groups of cells using
single-cell RNA-seq data, based on a linear embedding–based statistical framework.

## Features

- Supports sum, product, variance, and median combination statistics
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
git clone https://github.com/yourname/nutkin.git
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
    --group_col clone_id \
    --group1 <group_name_1> \
    --group2 <group_name_2> \
    --n_comps 30 \
    --skip_preprocessing \
    --save_scree_plot \
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
- Scree plot (`--save_scree_plot`): Shows the variance explained by each principal component
- Embedding plots (`--save_embedding_plot`): UMAP, PHATE, or den-SNE visualizations for exploratory analysis

### Step 2: Nutkin main analysis

**Input**
- Preprocessed AnnData object generated in Step 1

**Outputs**
- `var_measurement.csv`: Quantified variability for each defined group
- Pairwise statistical test results between specified groups

**Optional outputs (`--detail`)**
- Variability measures for individual principal components
- PCA scatter plots and contour plots (default: top 3 PCs)
- Summary tables and visualizations:
  - `summary_specific.csv`
  - `all_pairwise_scatter.png`
  - `all_pairwise_kde.png`

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

