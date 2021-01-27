import numpy as np
import pandas as pd
import subprocess as sp
import os
from tqdm import tqdm

def get_psi_table(SJ_table_name, minJR=5, minCell=20, drop_duplicates = False):
    '''
    
    This functions splits this table into one individual for
    each junction type. It additionally creates a PSI table
    based on PSI = (I1 + I2) / ((I1 + I2) + 2*SE)
    
    Input:
    - SJ_counts_table is a pandas dataframe, with rows corresponding
      to individual splicing junctions, and columns corresponding to
      single cells.
    
      The format of the index is:
      Gene_X_SJ
      Where Gene correspond to the gene where the splicing event
      is present. X is a number assigned to one particulat ASE. 
      NMD events have additionally _nmdSE_ between Gene and X.
      SJ corresponds to I1 (included 1), I2 (included 2) or SE
      (skipped exon).
    
    - minCell (int) is the minimum number of cells that are required
      to have at least minJR reads.
    
    - minJR (int) determines the minimum number of reads required on
      at least minCell cells for a junction to be acceptable.
    
    Output:
    - I1_counts (DF) Counts for one Included SJ
    - I2_counts (DF) Counts for the other Included SJ 
    - SE_counts (DF) Counts for Skipped Exon SJ 
    - PSI_table (DF) As in name.
    - total_counts (DF) SE + (I1 + I2)
    
    '''
    
    if drop_duplicates:
        SJ_counts_table = pd.read_csv(SJ_table_name, sep='\t', index_col=0).drop_duplicates('last')
    else:
        SJ_counts_table = pd.read_csv(SJ_table_name, sep='\t', index_col=0)
    
    
    events = sorted(set([x[:-3] for x in SJ_counts_table.index]))
    
    i1_events = [x + '_I1' for x in events]
    I1_table = SJ_counts_table.loc[i1_events]
    I1_table.index = events
    
    i2_events = [x + '_I2' for x in events]
    I2_table = SJ_counts_table.loc[i2_events]
    I2_table.index = events
    
    se_events = [x + '_SE' for x in events]
    SE_table = SJ_counts_table.loc[se_events]
    SE_table.index = events
    
    I1_filt = I1_table.index[(I1_table > minJR).sum(axis=1) > minCell]
    I2_filt = I2_table.index[(I2_table > minJR).sum(axis=1) > minCell]
    SE_filt = SE_table.index[(SE_table > minJR).sum(axis=1) > minCell]
    
    filtered_events = I1_filt & I2_filt & SE_filt
    
    I1_table = I1_table.loc[filtered_events]
    I2_table = I2_table.loc[filtered_events]
    SE_table = SE_table.loc[filtered_events]
    
    PSI_table = (I1_table + I2_table) /(2*SE_table + I1_table + I2_table)
    total_counts = SE_table + I1_table + I2_table

    return PSI_table, total_counts