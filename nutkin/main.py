from matplotlib import pyplot as plt
import logging
import textwrap
from preprocessing import *
from config import *
from parser import parse_arg
from nutkin import Nutkin
import pandas as pd

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def print_logo(args):
  
    print(textwrap.dedent(
            f'''\n\nWelcome to \n
                {NUTKIN_LOGO}
        '''))
    return



def main():
    args = parse_arg()
    print_logo(args)

    # === Step 1: Load and preprocess the data ===
    logger.info("Loading and preprocessing data...")
    adata = load_data(
        data_path=args.input_path,
        n_comps=args.n_comps,
        output_path=args.output_path,
        run_pca=True,
        run_embedding=args.run_embedding,
        method=args.embedding_method
    )

    from matplotlib import pyplot as plt
import logging
import textwrap
from preprocessing import *
from config import *
from parser import parse_arg
from nutkin import Nutkin
import pandas as pd

# setup logging
LOG_FORMAT = \
'(%(asctime)s) %(message)s'
DATE_FORMATE = '%d/%m/%Y %H:%M:%S' #'%a, %d %b %Y %H:%M:%S'
logging.basicConfig(format=LOG_FORMAT, datefmt=DATE_FORMATE)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def print_logo(args):
  
    print(textwrap.dedent(
            f'''\n\nWelcome to \n
                {NUTKIN_LOGO}
        '''))
    return



def main():
    args = parse_arg()
    print_logo(args)

    # === Step 1: Load and preprocess the data ===
    logger.info("Loading and preprocessing data...")
    adata = load_data(
        data_path=args.input_path,
        n_comps=args.n_comps,
        output_path=args.output_path,
        run_pca=True,
        run_embedding=args.run_embedding,
        method=args.embedding_method
    )

    
   
   
    # === Step 2: Initialize Nutkin ===
    logger.info("Initializing Nutkin...")
    processor = Nutkin(
        adata=adata,
        group_col=args.group_col,
        output_path=args.output_path,
        num_pc=DEFAULT_NUM_PC,
        metric_type=args.metric_type,
        variability_method=args.variability_method,
        verbal=True,
        detail=args.detail
    )

    # === Step 3: Measure variability within each cell group ===
    logger.info("Measuring variability...")
    var_df = processor.measure()
    logger.info("Variability summary:\n%s", var_df)

    # === Step 4: Perform differential testing ===
    logger.info("Performing pairwise differential tests and plotting results...")
    processor.differential_test(group1=args.group1,group2=args.group2)
    processor.visualize_pca_results(detail=args.detail,plot_pc_num=args.plot_pc_num,plot_type=args.plot_type)
   

if __name__ == "__main__":
    main()
