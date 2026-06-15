# Command-Line Interface

Nutkin can be run directly from the command line. The CLI runs the same preprocessing,
variability measurement, differential testing, and visualization steps as the Python API.

## Usage

Navigate to the project root directory first:

```bash
cd ../Nutkin
```

Then run the main script:

```bash
python -m nutkin.main \
    --input_path PATH_TO_YOUR_DATA \
    --output_path PATH_TO_YOUR_OUTPUT_FOLDER \
    --group_col NAME_OF_GROUP_COLUMN \
    --group1 NAME_OF_FIRST_GROUP \
    --group2 NAME_OF_SECOND_GROUP
```

## Example

```bash
python -m nutkin.main \
    --input_path ../Nutkin/data/test_data.h5ad \
    --output_path ../Nutkin/test/result/ \
    --group_col clone_id \
    --group1 mCHERRY_Barcode_1614 mCHERRY_Barcode_1614 \
    --group2 mCHERRY_Barcode_5774 mCHERRY_Barcode_1755 \
    --detail \
    --plot_type both
```

The outputs produced by this command are identical to those obtained using the Python API.

## All Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `--input_path` | str | *required* | Path to input `.h5ad` file |
| `--output_path` | str | *required* | Path to save output files |
| `--n_comps` | int | 50 | Number of PCA components for the scree plot |
| `--run_embedding` | flag | False | Compute and plot a low-dimensional embedding |
| `--embedding_method` | str | `umap` | Embedding method: `umap`, `phate`, or `densmap` |
| `--group_col` | str | `group` | Column in `adata.obs` for group labels |
| `--group1` | str(s) | None | Name(s) of the first group(s); space-separated for multiple |
| `--group2` | str(s) | None | Name(s) of the second group(s); space-separated for multiple |
| `--metric_type` | str | `sum` | Variability summary metric: `sum`, `product`, or `both` |
| `--variability_method` | str | `mad` | Variability measure: `sd`, `mad`, or `both` |
| `--detail` | flag | False | Enable detailed PCA outputs and visualizations |
| `--plot_pc_num` | int | 3 | Number of PCs to visualize |
| `--plot_type` | str | `scatter` | PCA plot type: `scatter`, `contour`, or `both` |
