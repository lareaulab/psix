import numpy as np
import pandas as pd
import subprocess as sp   
import os
from tqdm import tqdm
import gzip



def process_SJ_dir(rnaseq_dir,
                       intron_file,
                       save_files_in = '',
                       cell_list = []
                      ):
    
    splice_junction_table = get_SJ_table(rnaseq_dir, intron_file, cell_list = cell_list)
    
    if os.path.isdir(save_files_in):
        splice_junction_table.to_csv(save_files_in + '/splice_junctions.tab.gz', sep='\t', header=True, index=True)
    
    return splice_junction_table


def read_intron_file(intron_file):
    intron_table = pd.read_csv(intron_file, sep='\t', index_col=0)
    return intron_table[['intron']]


def process_SJ_table(cell_sj_file, cell):
    
    idx = []
    sj_counts = []
    with gzip.open(cell_sj_file, 'rb') as fh:#open(sj_file, 'r') as fh:
        for line in fh:
            line = line.decode().rstrip().split('\t')
            
            if line[3] == '1':
                strand = '+'
            elif line[3] == '2':
                strand = '-'
            else:
                continue
            
            intron_idx = line[0] + ':' + line[1] + '-' + line[2] + ':' + strand
            idx.append(intron_idx)
            sj_counts.append(int(line[6]))
            
    cell_series = pd.Series(sj_counts, index=idx, name=cell)
    
    return cell_series


def get_SJ_table(rnaseq_dir, intron_file, cell_list):
    
    if len(cell_list) > 0:
        cells = cell_list
        
    else:
        cells = sorted([x.split('.')[0] for x in os.listdir(rnaseq_dir)])
    
    series_list = []
    
    for cell in tqdm(cells, position = 0, leave=True):
        
        cell_sj_file = sj_file = rnaseq_dir + '/' + cell + '.SJ.out.tab.gz'
        cell_series = process_SJ_table(cell_sj_file, cell)
        series_list.append(cell_series)
        
    sj_table = pd.concat(series_list, axis=1)
    
    intron_table = read_intron_file(intron_file)
    
    sj_table = pd.merge(intron_table, sj_table, left_on='intron', how='left', right_index=True).fillna(0).drop('intron', axis=1)
    
    return sj_table


def get_psi_table(SJ_counts_table, minJR=1, minCell=1, tenX = False):
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
    
#     if drop_duplicates:
#         SJ_counts_table = pd.read_csv(SJ_table_name, sep='\t', index_col=0).drop_duplicates('last')
#     else:
#         SJ_counts_table = pd.read_csv(SJ_table_name, sep='\t', index_col=0)
        
    events_i1 = pd.Index([x[:-3] for x in SJ_counts_table.index if '_I1' in x])
    events_i2 = pd.Index([x[:-3] for x in SJ_counts_table.index if '_I2' in x])
    events_se = pd.Index([x[:-3] for x in SJ_counts_table.index if '_SE' in x])
    
    constitutive_intron_idx = pd.Index([x for x in SJ_counts_table.index if '_CI' in x])
    constitutive_sj = SJ_counts_table.loc[constitutive_intron_idx]
    constitutive_sj.index = pd.Index([x[:-3] for x in SJ_counts_table.index if '_CI' in x])
    
    events = events_i1 & events_i2 & events_se
    
    
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
    
    if tenX:
        I_table = pd.concat([I1_table, I2_table]).max(level=0)
        PSI_table = I_table /(SE_table + I_table)
        total_counts = SE_table + I2_table
        
    else:
        PSI_table = (I1_table + I2_table) /(2*SE_table + I1_table + I2_table)
        total_counts = SE_table + I1_table + I2_table

    return PSI_table, total_counts, constitutive_sj
