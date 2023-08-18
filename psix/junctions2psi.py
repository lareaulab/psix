import numpy as np
import pandas as pd
import subprocess as sp   
import os
from tqdm import tqdm
import gzip
from .mrna_census import *



def process_SJ_dir(rnaseq_dir, intron_file, save_files_in = '', cell_list = [], dtype=np.float64):
    
    splice_junction_table = get_SJ_table(rnaseq_dir, intron_file, cell_list = cell_list, dtype=dtype)
    
    if os.path.isdir(save_files_in):
        splice_junction_table.to_csv(save_files_in + '/splice_junctions.tab.gz', sep='\t', header=True, index=True)
    
    return splice_junction_table


def read_intron_file(intron_file):
    intron_table = pd.read_csv(intron_file, sep='\t', index_col=0)
    return intron_table[['intron']]


def process_SJ_line(line, idx, sj_counts):
    
    try:
        line = line.decode().rstrip().split('\t')
    except:
        line = line.rstrip().split('\t')

    if line[3] == '1':
        strand = '+'
    elif line[3] == '2':
        strand = '-'
    else:
        return idx, sj_counts

    intron_idx = line[0] + ':' + line[1] + '-' + line[2] + ':' + strand
    idx.append(intron_idx)
    sj_counts.append(int(line[6]))
    
    return idx, sj_counts


def process_SJ_table(cell_sj_file, cell, dtype):
    
    idx = []
    sj_counts = []
    
    if cell_sj_file[-3:] == '.gz':
        with gzip.open(cell_sj_file, 'rb') as fh:
            for line in fh:
                idx, sj_counts = process_SJ_line(line, idx, sj_counts)
    else:
        with open(cell_sj_file, 'r') as fh:
            for line in fh:
                idx, sj_counts = process_SJ_line(line, idx, sj_counts)
            
    cell_series = pd.Series(sj_counts, index=idx, name=cell, dtype=dtype)
    
    return cell_series


def get_SJ_table(rnaseq_dir, intron_file, cell_list, dtype):
    
    if len(cell_list) > 0:
        cells = cell_list
        
    else:
        cells = sorted([x.split('.')[0] for x in os.listdir(rnaseq_dir)])
    
    series_list = []
    
    for cell in tqdm(cells, position = 0, leave=True):
        
        
        if os.path.isfile(rnaseq_dir + '/' + cell + '.SJ.out.tab.gz'):
        
            cell_sj_file = sj_file = rnaseq_dir + '/' + cell + '.SJ.out.tab.gz'
        elif os.path.isfile(rnaseq_dir + '/' + cell + '.SJ.out.tab'):
            cell_sj_file = sj_file = rnaseq_dir + '/' + cell + '.SJ.out.tab'
        else:
            raise Exception('Missing ' + cell + ' SJ.out.tab file')
            
        cell_series = process_SJ_table(cell_sj_file, cell, dtype)
        series_list.append(cell_series)
        
    sj_table = pd.concat(series_list, axis=1)
    
    intron_table = read_intron_file(intron_file)
    
    sj_table = pd.merge(intron_table, sj_table, left_on='intron', 
                        how='left', right_index=True).fillna(0).drop('intron', axis=1)
    
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


def junctions_dir_to_psi(
        self,
        sj_dir,
        intron_file,
        tpm_file,
        cell_list = [],
        minJR = 1,
        minCell = 1,
        minPsi = 0.05,
        min_observed = 0.25,
        tenX = False,
        save_files_in = '',
        dtype=np.float64
    ):
        
    if dtype not in [np.float64, np.float32, np.float16]:
        raise Exception('dtype has to be numpy float type')
    

    print('Collecting splice junctions....', flush=True)


    if os.path.isfile(sj_dir):
        sj_file = pd.read_csv(sj_dir, sep='\t', index_col=0)
        
    else:
        sj_file = process_SJ_dir(sj_dir,
                             intron_file,
                             save_files_in = save_files_in,
                             cell_list = cell_list,
                             dtype=dtype
                            )

    print('Obtaining PSI tables...')

    psi, reads, constitutive_sj = get_psi_table(sj_file, minJR, minCell, tenX=tenX)

    alt_exons = psi.index[np.abs(0.5 - psi.mean(axis=1)) <= (0.5-minPsi)]
    obs_exons = psi.index[psi.isna().mean(axis=1) <= 1-min_observed]
    selected_exons = alt_exons & obs_exons

    psi = psi.loc[selected_exons]
    reads = reads.loc[selected_exons]
    

    if tenX:
        mrna_per_event = reads
    else:

        print('Reading TPM and transforming to mRNA counts...')

        tpm_exists = os.path.isfile(tpm_file)

        if len(cell_list) == 0:
            cell_list = psi.columns
        mrna = tpm2mrna(tpm_file, cell_list, dtype=dtype)
        ##### New thing
        cells = psi.columns & mrna.columns
        mrna = mrna[cells]
        psi = psi[cells]
        constitutive_sj = constitutive_sj[cells]

        mrna_per_event = get_mrna_per_event(mrna, psi, reads, constitutive_sj, solo=False) #constitutive_sj_file)

    if len(self.adata.obs) > 0:
        idx = self.adata.obs.index & mrna_per_event.index
    else:
        idx = mrna_per_event.index
        
    if dtype == np.float64:
        psi = psi.loc[idx]
        mrna_per_event = mrna_per_event.loc[idx]
    else:
        # setting to np.float32 for compatibility with numba
        psi = psi.loc[idx].astype(np.float32)
        mrna_per_event = mrna_per_event.loc[idx].astype(np.float32)

    if os.path.isdir(save_files_in):
        psi.to_csv(save_files_in + '/psi.tab.gz', sep='\t', 
                   index=True, header=True)
        mrna_per_event.to_csv(save_files_in + '/mrna.tab.gz', sep='\t', 
                   index=True, header=True)

    psi = psi.T
    mrna_per_event = mrna_per_event.T
    ncells_former = mrna_per_event.shape[0]

    mrna_per_event = mrna_per_event[(mrna_per_event.sum(axis=1) > 0.1) & (mrna_per_event.sum(axis=1) < np.inf)]
    ncells_current = mrna_per_event.shape[0]
    if ncells_former > ncells_current:
        n_diff = str(ncells_former - ncells_current)
        print('removed ' + n_diff + 'cells with missing of "inf" mRNA values.')
        print('This can be the consequence of very shallow coverage in the cell.')
        psi = psi.loc[mrna_per_event.index]
    
    self.adata.uns['psi'] = psi.T
    self.adata.uns['mrna_per_event'] = mrna_per_event.T

    print('Successfully processed RNA-seq data')
