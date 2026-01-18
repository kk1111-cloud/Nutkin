import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import umap



def load_data(data_path, n_comps=50, output_path=None, run_pca=True, run_embedding=False, method="umap"):
    """
    Complete pipeline: Load .h5ad, auto-preprocess if needed, run PCA, and generate embeddings.
    """
    # --- Step 1: Restrict input to .h5ad only ---
    if not data_path.endswith(".h5ad"):
        raise ValueError("Unsupported input format. Only .h5ad files are allowed.")
    
    adata = sc.read_h5ad(data_path)
    print(f"Successfully loaded data from: {data_path}")

    # --- Step 2: Automatic detection of preprocessing status ---
    # Check for 'log1p' or 'highly_variable' as indicators of prior processing
    is_processed = 'log1p' in adata.uns or 'highly_variable' in adata.var.columns
    
    if not is_processed:
        print("Data detected as RAW. Executing preprocessing pipeline...")
        adata = run_core_preprocessing(adata)
    else:
        print("Data detected as PREPROCESSED. Skipping core preprocessing steps.")

    # --- Step 3: Run PCA and generate scree visualization ---
    if run_pca:
        # Assign back to adata to store PCA results (X_pca, uns['pca'], etc.)
        run_pca_and_save_plot(adata, n_comps=n_comps, output_path=output_path)

    # --- Step 4: Generate Low-Dimensional Embedding ---
    if run_embedding:
        print(f"Step 4: Generating {method.upper()} embedding...")
        emb_filename = f"{method}_plot.png" if output_path else None
        emb_full_path = os.path.join(output_path, emb_filename) if output_path else None
        
        # plot_embedding internally modifies adata (adds X_umap etc.)
        plot_embedding(
            adata, 
            method=method, 
            output_path=emb_full_path
        )
        
    return adata

def run_core_preprocessing(adata):
    """Standard QC, normalization, and scaling."""
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Subset to HVGs and copy to avoid view/fragmentation issues
    adata = adata[:, adata.var.highly_variable].copy() 
    sc.pp.scale(adata, max_value=10)
    return adata

def run_pca_and_save_plot(adata, n_comps=50, output_path=None):
    """Computes PCA and saves Scree Plot."""
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
    
    if output_path is not None:
        os.makedirs(output_path, exist_ok=True)
        variance_ratio = adata.uns["pca"]["variance_ratio"]
        plt.figure(figsize=(8, 5))
        plt.plot(range(1, len(variance_ratio)+1), variance_ratio, marker='o', label='Individual')
        plt.title("PCA Scree Plot")
        plt.savefig(os.path.join(output_path, "pca_scree_plot.png"), dpi=300)
        plt.close()


def plot_embedding(adata, method="umap", group_col="group", use_rep="X_pca", output_path=None):
    """Computes and plots embeddings like UMAP."""
    method = method.lower()
    
    if method == "umap":
        sc.pp.neighbors(adata, use_rep=use_rep)
        sc.tl.umap(adata)
        basis = "umap"
    else:
        raise ValueError(f"Method {method} not currently implemented in this block.")

    # Handle file saving logic
    save_filename = None
    if output_path:
        save_dir = os.path.dirname(output_path)
        if save_dir: # Ensure directory string isn't empty
            os.makedirs(save_dir, exist_ok=True)
        sc.settings.figdir = save_dir
        save_filename = os.path.basename(output_path)  

    sc.pl.embedding(
        adata,
        basis=basis,
        color=group_col,
        show=False,
        save=save_filename 
    )
