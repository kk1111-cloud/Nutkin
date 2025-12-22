import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

import umap


def load_data(data_path):
    # --- Step 1: Load input data based on format ---
    if data_path.endswith(".h5ad"):
        adata = sc.read_h5ad(data_path)

    elif os.path.isdir(data_path) and os.path.exists(os.path.join(data_path, "matrix.mtx")):
        adata = sc.read_10x_mtx(data_path, var_names="gene_symbols", cache=True)

    elif data_path.endswith(".csv"):
        df = pd.read_csv(data_path)
        group_col = "group" if "group" in df.columns else ("Group" if "Group" in df.columns else None)
        if group_col:
            adata = ad.AnnData(df.drop(columns=[group_col]))
            adata.obs['group'] = df[group_col].values
        else:
            adata = ad.AnnData(df)

    else:
        raise ValueError("Unsupported input format. Use .h5ad, 10X folder, or .csv")
    return adata
 
def preprocess_adata(adata,n_comps,save_scree_plot=True, skip_preprocessing=False, output_path=None):
    # --- Step 1: Standard preprocessing steps ---
    if not skip_preprocessing:
        sc.pp.filter_cells(adata, min_genes=200)  # Filter out low-quality cells
        sc.pp.filter_genes(adata, min_cells=3)    # Filter out rarely expressed genes
        sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize total counts per cell
        sc.pp.log1p(adata)                        # Log-transform the data
        sc.pp.highly_variable_genes(
            adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )  # Identify highly variable genes
        adata = adata[:, adata.var.highly_variable]  # Subset to HVGs
        sc.pp.scale(adata, max_value=10)         # Scale each gene to unit variance and clip values
    else:
        print("Skipping preprocessing... assuming data is already processed.")

    # --- Step 2: PCA ---
    sc.tl.pca(adata, svd_solver='arpack',n_comps=n_comps)    # Perform PCA
    
    # --- Step 3: Variance explained by PCs ---
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    cumulative = np.cumsum(variance_ratio)

    print("\nExplained variance by top PCs:")
    for i in range(min(10, len(variance_ratio))):
        print(f"  PC{i+1}: {variance_ratio[i]:.4f} (Cumulative: {cumulative[i]:.4f})")

    if save_scree_plot == True:
        plt.figure(figsize=(8, 5))
        plt.plot(range(1, len(variance_ratio)+1), variance_ratio, marker='o', label='Individual PC')
        plt.plot(range(1, len(cumulative)+1), cumulative, marker='s', label='Cumulative')
        plt.axhline(0.01, color='gray', linestyle='--', label='1% Threshold')
        plt.xlabel("Principal Component")
        plt.ylabel("Explained Variance Ratio")
        plt.title("PCA Variance Explained")
        plt.legend()
        plt.tight_layout()
        if output_path is not None:
            scree_plot_file = os.path.join(output_path, "pca_scree_plot.png")
            plt.savefig(scree_plot_file, dpi=300)
            print(f"PCA scree plot saved to {scree_plot_file}")
    
        plt.close()

    return adata


def plot_embedding(
    adata,
    method="umap",
    color="group",
    use_rep="X_pca",
    save_embedding_plot=False,
    output_path=None):
    """
    Generate and plot low-dimensional embeddings for AnnData.
    
    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object with .X or .obsm["X_pca"].
    method : str
        One of "umap", "phate", "densmap".
    color : str or list of str
        Column(s) in adata.obs to color the plot by.
    use_rep : str
        Representation to use for embedding, e.g., "X_pca" or "X".
    save : bool or str
        If True, saves the plot to a file. If str, saves with that filename.
    show : bool
        Whether to display the plot.
    **kwargs : additional arguments passed to the embedding function
    """

    method = method.lower()
    assert method in ["umap", "phate", "densmap"], f"Unsupported method: {method}"
    rep = adata.obsm[use_rep] if use_rep in adata.obsm else adata.X

    # Compute embedding
    if method == "umap":
        sc.pp.neighbors(adata, use_rep=use_rep)
        sc.tl.umap(adata)
        basis = "umap"

    if save_embedding_plot and output_path is not None:
        save_dir = os.path.dirname(output_path)
        os.makedirs(save_dir, exist_ok=True)

        sc.settings.figdir = save_dir
        save_filename = os.path.basename(output_path)  
    else:
        save_filename = None



    # Plot
    sc.pl.embedding(
        adata,
        basis=basis,
        color=color,
        show=False,
        save=save_filename 
    )
