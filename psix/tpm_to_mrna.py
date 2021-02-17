import numpy as np
import pandas as pd
import subprocess as sp
import os
from scipy.stats import gaussian_kde
from numpy import random as r
from tqdm import tqdm

from numba import jit

@jit
def normalize_equation(cell, moda):
    n = np.sum(np.log10(cell) <= moda)
    interval = np.sum(cell.loc[np.log10(cell) <= moda])/np.sum(cell)
    return n/interval

@jit
def get_transcripts_per_cell(cell, remove_outliers, bw_method, adjust_high):
    z = gaussian_kde(np.log10(cell), bw_method)
    
    moda = np.arange(-1, 10, 0.1)[z.pdf(np.arange(-1, 10, 0.1)).argmax()]
    
    molecules_in_cell = normalize_equation(cell, moda)
    
    
    if (molecules_in_cell > 10**5.5) and adjust_high:
        moda = np.arange(1, 10, 0.1)[z.pdf(np.arange(1, 10, 0.1)).argmax()]
    
        molecules_in_cell = normalize_equation(cell, moda)
            
    
    if remove_outliers and (molecules_in_cell > 10**5.5):
        
        return 0
        
    return molecules_in_cell
        
@jit    
def transform_cell(cell, remove_outliers, bw_method, adjust_high):
    cell_filtered = cell.loc[cell > 0.1]
    molecules_in_cell = get_transcripts_per_cell(cell_filtered, remove_outliers, bw_method, adjust_high)
    cell_remove_zeros = cell * (cell > 0.1)
    normalize_lysate = molecules_in_cell / 10**5.5
    cell_transcript_counts = cell_remove_zeros * normalize_lysate
    
    return cell_transcript_counts


def tpm2mrna(tpm_file, bw_method='scott', adjust_high = True, remove_outliers=True):
    tpm_dataset = pd.read_csv(tpm_file, sep='\t', index_col=0)
    mrna_counts = pd.DataFrame()
    mrna_counts_per_cell = []
    cells = tpm_dataset.columns
    tpm_dataset_filtered = tpm_dataset.loc[tpm_dataset.max(axis=1) > 0.1]
    
    for cell in tqdm(cells, position=0, leave=True):
        cell_mrna = transform_cell(tpm_dataset_filtered[cell], remove_outliers, bw_method, adjust_high)
        if all([x == 0 for x in cell_mrna]):
            continue
        mrna_counts[cell] = cell_mrna
        
    return mrna_counts


