import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import csv
import gzip
import scipy.io

def read_solo_features(solo_features_path):
    
    print('Reading solo splice junctions')
    
    features_solo = pd.read_csv(solo_features_path, sep='\t',
                           names = ['chrom', 'start', 'end', 'strand', 'motif', 'annot', 'unique', 'multimap', 'overlap'])
    intron_list = []
    for idx, row in features_solo.iterrows():
        intron = row.chrom + ':' + str(row.start) + '-' + str(row.end) + ':'
        if row.strand == 1:
            strand = '+'
        else:
            strand = '-'
        intron += strand

        intron_list.append(intron)
        
    return intron_list

def read_solo_barcodes(solo_barcodes_path):
    print('getting barcodes')
    barcodes = [row[0] for row in csv.reader(gzip.open(solo_barcodes_path, 'rt'), delimiter="\t")]
    return barcodes

def read_solo_matrix(solo_matrix_path, intron_list, barcodes, cell_list):
    print('getting matrix')
    matrix = pd.DataFrame.sparse.from_spmatrix(scipy.io.mmread(solo_matrix_path),
                                               index=intron_list, columns=barcodes)
    
    matrix = matrix[cell_list]
    matrix = matrix.loc[[x.split(':')[0][:3]=='chr' for x in matrix.index]]
    matrix = matrix.loc[(matrix.sum(axis=1) > 0)]
    
    return matrix

def process_solo(solo_dir, intron_tab, cell_list):
    
    solo_features_path = solo_dir + '/features.tsv.gz'
    solo_barcodes_path = solo_dir + '/barcodes.tsv.gz'
    solo_matrix_path = solo_dir + '/matrix.mtx.gz'
    
    intron_list = read_solo_features(solo_features_path)
    barcodes = read_solo_barcodes(solo_barcodes_path)
    
    if len(cell_list) == 0:
        cell_list = barcodes
        
    matrix = read_solo_matrix(solo_matrix_path, intron_list, barcodes, cell_list)
    
    intron_events = pd.read_csv(intron_tab, sep='\t', index_col=0)
    
    intron_mtx = intron_events.drop_duplicates().merge(matrix.drop_duplicates(), left_on='intron', right_index=True)[matrix.columns]
    
    intron_mtx_CI = intron_mtx.loc[[x[-3:] == '_CI' for x in intron_mtx.index]]
    intron_mtx_exons = intron_mtx.loc[[x[-3:] != '_CI' for x in intron_mtx.index]]
    
    return intron_mtx_CI, intron_mtx_exons


def get_psi_table_solo(intron_mtx_exons, minJR=5, minCell=20, tenX = False):
    print('getting psi table')
    '''
    
    This functions splits this table into one individual for
    each junction type. It additionally creates a PSI table
    based on PSI = (I1 + I2) / ((I1 + I2) + 2*SE)
    
    Input:
    - intron_mtx_exons is a pandas dataframe, with rows corresponding
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
        
    events_i1 = pd.Index([x[:-3] for x in intron_mtx_exons.index if '_I1' in x])
    events_i2 = pd.Index([x[:-3] for x in intron_mtx_exons.index if '_I2' in x])
    events_se = pd.Index([x[:-3] for x in intron_mtx_exons.index if '_SE' in x])
    
    events = events_i1 & events_i2 & events_se
    
    
    i1_events = [x + '_I1' for x in events]
    I1_table = intron_mtx_exons.loc[i1_events]
    I1_table.index = events
    
    i2_events = [x + '_I2' for x in events]
    I2_table = intron_mtx_exons.loc[i2_events]
    I2_table.index = events
    
    se_events = [x + '_SE' for x in events]
    SE_table = intron_mtx_exons.loc[se_events]
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

    return PSI_table, total_counts

    
    