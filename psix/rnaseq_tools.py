import numpy as np
import pandas as pd
import subprocess as sp
import os
import anndata
from tqdm import tqdm
from scipy.stats import gaussian_kde
from numpy import random as r
from tpm_to_mrna import *
from solo_tools import *


def get_mrna_per_event(mrna, psi, reads, constitutive_sj, solo):#constitutive_sj_file):

    #constitutive_sj = pd.read_csv(constitutive_sj_file, sep='\t', index_col=0)
    obs_mrna = mrna.index[mrna.median(axis=1) >= 1]
    obs_junctions = [x for x in constitutive_sj.index if x.split('_')[0] in obs_mrna]

    mrna_per_junction = mrna.loc[[x.split('_')[0] for x in obs_junctions]]
    mrna_per_junction.index = obs_junctions

    
    if solo:
        reads_per_junction = (constitutive_sj.loc[obs_junctions] / mrna_per_junction).sparse.to_dense().replace([np.inf, -np.inf], np.nan)
    else:
        reads_per_junction = (constitutive_sj.loc[obs_junctions] / mrna_per_junction).replace([np.inf, -np.inf], np.nan)
    SJ_mean = reads_per_junction.mean()
    
    mrna_events = (reads/(SJ_mean * (1+psi)))
    
    return mrna_events
    

def get_psi_table(SJ_table_name, minJR=5, minCell=20, drop_duplicates = False, tenX = False):
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
        
    events_i1 = pd.Index([x[:-3] for x in SJ_counts_table.index if '_I1' in x])
    events_i2 = pd.Index([x[:-3] for x in SJ_counts_table.index if '_I2' in x])
    events_se = pd.Index([x[:-3] for x in SJ_counts_table.index if '_SE' in x])
    
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

    return PSI_table, total_counts


def normalize_equation(cell, moda):
    n = np.sum(np.log10(cell) <= moda)
    interval = np.sum(cell.loc[np.log10(cell) <= moda])/np.sum(cell)
    return n/interval


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
        
    
def transform_cell(cell, remove_outliers, bw_method, adjust_high):
    cell_filtered = cell.loc[cell > 0.1]
    molecules_in_cell = get_transcripts_per_cell(cell_filtered, remove_outliers, bw_method, adjust_high)
    cell_remove_zeros = cell * (cell > 0.1)
    normalize_lysate = molecules_in_cell / 10**6
    cell_transcript_counts = cell_remove_zeros * normalize_lysate
    
    return cell_transcript_counts


def tpm2mrna(tpm_file, bw_method='scott', adjust_high = True, remove_outliers=True):
    tpm_dataset = pd.read_csv(tpm_file, sep='\t', index_col=0)
    mrna_counts = pd.DataFrame()
    mrna_counts_per_cell = []
    cells = tpm_dataset.columns
    tpm_dataset_filtered = tpm_dataset.loc[tpm_dataset.max(axis=1) > 0.1]
    
    for cell in tqdm(cells):
        cell_mrna = transform_cell(tpm_dataset_filtered[cell], remove_outliers, bw_method, adjust_high)
        if all([x == 0 for x in cell_mrna]):
            continue
        mrna_counts[cell] = cell_mrna
        
    return mrna_counts


def process_rnaseq_files(
    self,
    exon_sj_file,
    constitutive_sj_file = '',
    tpm_file = '',
    minJR = 1,
    minCell = 1,
    drop_duplicates = False,
    min_psi = 0.05,
    min_observed = 0.25,
    tenX = False
):

    print('Obtaining PSI tables...')

    psi, reads = get_psi_table(exon_sj_file, minJR, minCell, drop_duplicates, tenX=tenX)

    alt_exons = psi.index[np.abs(0.5 - psi.mean(axis=1)) <= (0.5-min_psi)]
    obs_exons = psi.index[psi.isna().mean(axis=1) <= 1-min_observed]
    selected_exons = alt_exons & obs_exons

    psi = psi.loc[selected_exons]
    reads = reads.loc[selected_exons]

    if tenX:
        mrna_per_event = reads
    else:

        print('Reading TPM and transforming to mRNA counts...')

        tpm_exists = os.path.isfile(tpm_file)
        constitutive_sj_exists = os.path.isfile(constitutive_sj_file)

        if not (tpm_exists and constitutive_sj_exists):
            raise Exception('TPM file and constitutive junctions are required when processing smart-seq data')

        mrna = tpm_to_mrna2(tpm_file)
        ##### New thing
        cells = psi.columns & mrna.columns
        mrna = mrna[cells]
        psi = psi[cells]
        
        constitutive_sj = pd.read_csv(constitutive_sj_file, sep='\t', index_col=0)
        mrna_per_event = get_mrna_per_event(mrna, psi, reads, constitutive_sj, solo=False) #constitutive_sj_file)

    if len(self.adata.obs) > 0:
        idx = self.adata.obs.index & mrna_per_event.index
    else:
        idx = mrna_per_event.index

    self.adata.uns['psi'] = psi.loc[idx].T
    self.adata.uns['mrna_per_event'] = mrna_per_event.loc[idx].T

    print('Successfully processed RNA-seq data')

    
    
