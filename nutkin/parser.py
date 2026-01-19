import textwrap
import logging
import argparse
import numpy as np
from config import *
import sys

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def print_logo(args):
    if not args.minimal_stdout:
        print(textwrap.dedent(
            f'''\n\nWelcome to \n
                {NUTKIN_LOGO}
        '''))
    return

def parse_arg():
    parser = argparse.ArgumentParser(
    description=textwrap.dedent(
    '''NUTKIN is a tool to quantify and testing diffences of cell-to-cell variability using normalized single cell RNA-seq (scRNA-Seq) data. 
    It takes a normalized scRNA-Seq file as an input and output two pandas dataframes, 
    one is the quantification result, the other is dimensional reduction scRNA-seq data which can be later used as a required input for differential analysis. 
    \n
    '''),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser = argparse.ArgumentParser(description="Run Nutkin variability analysis.")
    parser.add_argument("--input_path", type=str, required=True, help="Path to input file (.h5ad, .csv, or 10X folder)")
    parser.add_argument("--output_path", type=str, required=True,help="Path to save processed AnnData file (.h5ad)")
    parser.add_argument("--n_comps", type=int, default=50, help="Number of PCA components in scree plot")
    parser.add_argument("--run_pca", action="store_true", help="Whether to compute and plot PCA scree plot")
    parser.add_argument("--run_embedding", action="store_true", help="Whether to compute and plot low-dimensional embedding (e.g. UMAP)")
    parser.add_argument("--embedding_method", choices=["umap", "phate", "densmap"], default="umap",
                        help="Embedding method to use")
    parser.add_argument("--group_col", type=str, default="group",
                        help="Metadata column(s) in adata.obs used to color the plot and group")
    parser.add_argument("--use_rep", default="X_pca", help="Representation to use for embedding (e.g., X_pca or X)")


    
    parser.add_argument("--group1",nargs='+', default=None, help="Group 1 name (optional)")
    parser.add_argument("--group2",nargs='+', default=None, help="Group 2 name (optional)")
    parser.add_argument("--metric_type", choices=["sum", "product", "both"], default='sum', help="Type of metric to test: 'sum', 'product', or 'both")
    parser.add_argument("--variability_method", choices=["sd", "mad", "both"], default='mad', help="Type of metric to calculate variability: 'standard deviation', 'median absolute deviation', or 'both")
    parser.add_argument("--detail", action="store_true", help="Detailed PCA visualization")
    parser.add_argument("--plot_pc_num", default=3, help="Number of PCs to plot")
    parser.add_argument("--plot_type", choices=["scatter", "contour","both"], default='scatter', help="Type of PCA visualization: 'scatter', 'contour','both'")
    return parser.parse_args()
    
    
