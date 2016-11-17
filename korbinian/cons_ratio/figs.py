from scipy.stats import ttest_ind
import ast
import csv
import itertools
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def save_figures_describing_proteins_in_list(pathdict, s, logging):
     logging.info("~~~~~~~~~~~~         starting run_save_figures_describing_proteins_in_list          ~~~~~~~~~~~~")

     logging.info("~~~~~~~~~~~~        run_save_figures_describing_proteins_in_list is finished        ~~~~~~~~~~~~")