# Installation

## Requirements

- Python >= 3.9
- NumPy >= 1.21
- pandas >= 1.4
- SciPy >= 1.7
- Scanpy >= 1.9
- anndata >= 0.8
- scikit-learn >= 1.1

Optional (for plotting):

- umap-learn >= 0.5
- matplotlib >= 3.5
- seaborn >= 0.12

## Install from PyPI

```bash
pip install nutkin
```

## Install from Source

```bash
git clone https://github.com/shimlab/Nutkin.git
cd Nutkin
pip install -e .
```

To include optional plotting dependencies:

```bash
pip install -e ".[plot]"
```
