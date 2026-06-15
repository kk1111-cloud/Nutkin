import os
import numpy as np


# Default number of PC components
DEFAULT_NUM_PC = 15
# Optimal testing parameters 
BOOTSTRAP_SAMPLE_SIZE_UPPER_BOUND = 80
NUMBER_OF_BOOTSTRAP_SAMPLES = 1000
MAXIMUM_BOOTSTRAP_PROPORTION = 0.8
DEFAULT_SIGNIFICANCE = 0.05
RANDOM_SEED = 42
DECIMAL = 6


ARG_NOT_GIVEN = "arg_was_not_given"
# text logo 
NUTKIN_LOGO = \
"""
 N     N  U     U  TTTTTTTT  K     K  II  N     N
 N N   N  U     U     TT     K    K   II  N N   N
 N  N  N  U     U     TT     K  K     II  N  N  N
 N   N N  U     U     TT     KKK      II  N   N N
 N    NN  U     U     TT     K  K     II  N    NN
 N     N  U     U     TT     K    K   II  N     N
 N     N  UUUUUUU     TT     K     K  II  N     N

"""