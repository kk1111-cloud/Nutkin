\# Nutkin



Nutkin is a Python package for quantifying and testing differences in overall

transcriptional cell-to-cell variability between groups of cells using

single-cell RNA-seq data, based on a linear embedding–based statistical framework.



\## Features



\- Supports sum, product, variance, and median combination statistics

\- Works for both balanced and unbalanced group designs

\- Bootstrap-based inference with parallel computing support

\- Designed for single-cell RNA-seq data analysis



\## Installation



\### Requirements

\- Python >= 3.9

\- NumPy

\- ScanPy

\- pandas



See `pyproject.toml` for the full list of dependencies.



\### Install from PyPI



```bash

pip install nutkin

```



\### Install from source



```bash

git clone https://github.com/yourname/nutkin.git

cd nutkin

pip install -e .

```



\## Quick Start



```python

from nutkin import run\_test



pval = run\_test(

&nbsp;   data=expr\_matrix,

&nbsp;   group=group\_labels,

&nbsp;   method="sum",

&nbsp;   n\_bootstrap=1000

)



print(pval)

```



\## Usage



\### Basic workflow



1\. Prepare single-cell RNA-seq data

2\. Specify group labels and analysis parameters

3\. Run Nutkin

4\. Interpret the output statistics



\### Required input



\- \*\*scRNA-seq data\*\*

&nbsp; Data can be provided either as a preprocessed or raw dataset.

&nbsp; If raw data are supplied, Nutkin will perform preprocessing automatically.



\### Running Nutkin from the command line



```bash

python main.py \\

&nbsp; --input\_dataset <path\_to\_input> \\

&nbsp; --output\_path <path\_to\_output> \\

&nbsp; --group\_col clone\_id \\

&nbsp; --group1 <group\_name\_1> \\

&nbsp; --group2 <group\_name\_2> \\

&nbsp; --n\_comps 30 \\

&nbsp; --skip\_preprocessing \\

&nbsp; --save\_scree\_plot \\

&nbsp; --detail

```



\## Nutkin pipeline and outputs



Nutkin consists of two main steps.



\### Step 1: Data preprocessing



The preprocessing module (`preprocessing.py`) performs the following:



\- Reads the input data and checks its format

\- If the input is an AnnData object:

&nbsp; - Automatically detects whether preprocessing has been performed

&nbsp; - If not, preprocessing and normalization are conducted using Scanpy

\- If the input is not an AnnData object, an error is raised



\*\*Input\*\*

\- Single-cell RNA-seq data (`.h5ad`)



\*\*Output\*\*

\- A cleaned and processed AnnData object ready for Nutkin analysis



\*\*Optional outputs\*\*

\- Scree plot (`--save\_scree\_plot`):

&nbsp; Shows the variance explained by each principal component

\- Embedding plots (`--save\_embedding\_plot`):

&nbsp; UMAP, PHATE, or den-SNE visualizations for exploratory analysis



\### Step 2: Nutkin main analysis



\*\*Input\*\*

\- Preprocessed AnnData object generated in Step 1



\*\*Outputs\*\*

\- `var\_measurement.csv`:

&nbsp; Quantified variability for each defined group

\- Pairwise statistical test results between specified groups



\*\*Optional outputs (`--detail`)\*\*

\- Variability measures for individual principal components

\- PCA scatter plots and contour plots (default: top 3 PCs)

\- Summary tables and visualizations:

&nbsp; - `summary\_specific.csv`

&nbsp; - `all\_pairwise\_scatter.png`

&nbsp; - `all\_pairwise\_kde.png`



\## Documentation



Full documentation will be available at:



https://nutkin.readthedocs.io



The documentation includes:

\- Installation guide

\- Detailed usage examples

\- API reference

\- Methodological background



\## License



MIT License

\# Nutkin



Nutkin is a Python package for quantifying and testing differences in overall

transcriptional cell-to-cell variability between groups of cells using

single-cell RNA-seq data, based on a linear embedding–based statistical framework.



\## Features



\- Supports sum, product, variance, and median combination statistics

\- Works for both balanced and unbalanced group designs

\- Bootstrap-based inference with parallel computing support

\- Designed for single-cell RNA-seq data analysis



\## Installation



\### Requirements

\- Python >= 3.9

\- NumPy

\- SciPy

\- pandas



See `pyproject.toml` for the full list of dependencies.



\### Install from PyPI



```bash

pip install nutkin

```



\### Install from source



```bash

git clone https://github.com/yourname/nutkin.git

cd nutkin

pip install -e .

```



\## Quick Start



```python

from nutkin import run\_test



pval = run\_test(

&nbsp;   data=expr\_matrix,

&nbsp;   group=group\_labels,

&nbsp;   method="sum",

&nbsp;   n\_bootstrap=1000

)



print(pval)

```



\## Usage



\### Basic workflow



1\. Prepare single-cell RNA-seq data

2\. Specify group labels and analysis parameters

3\. Run Nutkin

4\. Interpret the output statistics



\### Required input



\- \*\*scRNA-seq data\*\*

&nbsp; Data can be provided either as a preprocessed or raw dataset.

&nbsp; If raw data are supplied, Nutkin will perform preprocessing automatically.



\### Running Nutkin from the command line



```bash

python main.py \\

&nbsp; --input\_dataset <path\_to\_input> \\

&nbsp; --output\_path <path\_to\_output> \\

&nbsp; --group\_col clone\_id \\

&nbsp; --group1 <group\_name\_1> \\

&nbsp; --group2 <group\_name\_2> \\

&nbsp; --n\_comps 30 \\

&nbsp; --skip\_preprocessing \\

&nbsp; --save\_scree\_plot \\

&nbsp; --detail

```



\## Nutkin pipeline and outputs



Nutkin consists of two main steps.



\### Step 1: Data preprocessing



The preprocessing module (`preprocessing.py`) performs the following:



\- Reads the input data and checks its format

\- If the input is an AnnData object:

&nbsp; - Automatically detects whether preprocessing has been performed

&nbsp; - If not, preprocessing and normalization are conducted using Scanpy

\- If the input is not an AnnData object, an error is raised



\*\*Input\*\*

\- Single-cell RNA-seq data (`.h5ad`)



\*\*Output\*\*

\- A cleaned and processed AnnData object ready for Nutkin analysis



\*\*Optional outputs\*\*

\- Scree plot (`--save\_scree\_plot`):

&nbsp; Shows the variance explained by each principal component

\- Embedding plots (`--save\_embedding\_plot`):

&nbsp; UMAP, PHATE, or den-SNE visualizations for exploratory analysis



\### Step 2: Nutkin main analysis



\*\*Input\*\*

\- Preprocessed AnnData object generated in Step 1



\*\*Outputs\*\*

\- `var\_measurement.csv`:

&nbsp; Quantified variability for each defined group

\- Pairwise statistical test results between specified groups



\*\*Optional outputs (`--detail`)\*\*

\- Variability measures for individual principal components

\- PCA scatter plots and contour plots (default: top 3 PCs)

\- Summary tables and visualizations:

&nbsp; - `summary\_specific.csv`

&nbsp; - `all\_pairwise\_scatter.png`

&nbsp; - `all\_pairwise\_kde.png`



\## Documentation



Full documentation will be available at:



https://nutkin.readthedocs.io



The documentation includes:

\- Installation guide

\- Detailed usage examples

\- API reference

\- Methodological background



\## License



MIT License



