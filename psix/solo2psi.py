import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import csv
import gzip
from scipy.io import mmread
from scipy.sparse import csr_matrix
from .mrna_census import *

def read_solo_features(solo_features_path):
        
    features_solo = pd.read_csv(solo_features_path, sep='\t',
                           names = ['chrom', 'start', 'end', 'strand', 'motif', 'annot', 'unique', 'multimap', 'overlap'])
    intron_list = []
    for idx, row in features_solo.iterrows():
        intron = str(row.chrom) + ':' + str(row.start) + '-' + str(row.end) + ':'
        if row.strand == 1:
            strand = '+'
        else:
            strand = '-'
        intron += strand

        intron_list.append(intron)
        
    return intron_list

def read_solo_barcodes(solo_barcodes_path):
    if solo_barcodes_path[-3:] == '.gz':
        barcodes = [row[0] for row in csv.reader(gzip.open(solo_barcodes_path, 'rt'), delimiter="\t")]
    else:
        barcodes = []
        with open(solo_barcodes_path, 'r') as fh:
            for row in fh:
                barcodes.append(row.rstrip())
    return barcodes

def read_solo_matrix(solo_matrix_path, intron_list, barcodes, cell_list, intron_set):
    print('loading sparse matrix')
    m = mmread(solo_matrix_path)
    m = csr_matrix(m)
    print('sub-setting sparse matrix')
    Z = [i for i, x in enumerate(intron_list) if ((x in intron_set) and (x[:3]=='chr'))]
    Y = [i for i, x in enumerate(barcodes) if x in cell_list]

    Z_names = [x for i, x in enumerate(intron_list) if  ((x in intron_set) and (x[:3]=='chr'))]
    Y_names = [x for i, x in enumerate(barcodes) if x in cell_list]

    m = m.transpose()[Y].transpose()[Z]
    m = pd.DataFrame.sparse.from_spmatrix(m, index=Z_names, columns=Y_names)
    #matrix = pd.DataFrame.sparse.from_spmatrix(scipy.io.mmread(solo_matrix_path),
    #                                           index=intron_list, columns=barcodes)
    
    m = m[cell_list]
    #m = m.loc[[x.split(':')[0][:3]=='chr' for x in m.index]]
    m = m.loc[(m.sum(axis=1) > 0)]
    print('intron matrix shape: ' + str(m.shape))
    return m

def process_solo(solo_dir, intron_file, cell_list):
    
    if os.path.isfile(solo_dir + '/features.tsv.gz'):
        solo_features_path = solo_dir + '/features.tsv.gz'
    elif os.path.isfile(solo_dir + '/features.tsv'):
        solo_features_path = solo_dir + '/features.tsv'
    else:
        raise Exception('features file not found in solo directory')
        
    if os.path.isfile(solo_dir + '/barcodes.tsv.gz'):
        solo_barcodes_path = solo_dir + '/barcodes.tsv.gz'
    elif os.path.isfile(solo_dir + '/barcodes.tsv'):
        solo_barcodes_path = solo_dir + '/barcodes.tsv'
    else:
        raise Exception('barcodes file not found in solo directory')
        
    if os.path.isfile(solo_dir + '/matrix.mtx.gz'):
        solo_matrix_path = solo_dir + '/matrix.mtx.gz'
    elif os.path.isfile(solo_dir + '/matrix.mtx'):
        solo_matrix_path = solo_dir + '/matrix.mtx'
    else:
        raise Exception('matrix file not found in solo directory')
    
    print('reading barcodes, features')
    intron_list = read_solo_features(solo_features_path)
    barcodes = read_solo_barcodes(solo_barcodes_path)
    
    if len(cell_list) == 0:
        cell_list = barcodes
        
    
    intron_events = pd.read_csv(intron_file, sep='\t', index_col=0).drop_duplicates()
    intron_set = set(intron_events.intron)
    print('reading matrix')
    matrix = read_solo_matrix(solo_matrix_path, intron_list, barcodes, cell_list, intron_set)
    
    print('Matrix shape:')
    print(matrix.shape)
    # intron_events = pd.read_csv(intron_file, sep='\t', index_col=0)
    
    intron_mtx = intron_events.drop_duplicates().merge(matrix.drop_duplicates(), left_on='intron', right_index=True)[matrix.columns]
    
    intron_mtx_CI = intron_mtx.loc[[x[-3:] == '_CI' for x in intron_mtx.index]]
    intron_mtx_exons = intron_mtx.loc[[x[-3:] != '_CI' for x in intron_mtx.index]]
    
    return intron_mtx_CI, intron_mtx_exons


def get_psi_table_solo(intron_mtx_exons, minJR=5, minCell=20, tenX = False):
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
    
    events = events_i1.intersection(events_i2).intersection(events_se)
    
    
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
    
    filtered_events = I1_filt.intersection(I2_filt).intersection(SE_filt)
    
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

    
def solo_to_psi(
    self,
    solo_dir,
    intron_file,
    tpm_file,
    cell_list = [],
    minJR = 1,
    minCell = 1,
    minPsi = 0.05,
    min_observed = 0.25,
    tenX = False,
    save_files_in = ''
):

    print('Processing STARsolo output. This might take a few minutes...')
    intron_mtx_CI, intron_mtx_exons = process_solo(solo_dir, intron_file, cell_list)

    print('Obtaining PSI tables...')

    psi, reads = get_psi_table_solo(intron_mtx_exons, minJR, minCell, tenX=tenX)

    alt_exons = psi.index[np.abs(0.5 - psi.mean(axis=1)) <= (0.5-minPsi)]
    obs_exons = psi.index[psi.isna().mean(axis=1) <= 1-min_observed]
    selected_exons = alt_exons.intersection(obs_exons)

    psi = psi.loc[selected_exons]
    reads = reads.loc[selected_exons]

    if tenX:
        mrna_per_event = reads
    else:

        print('Reading TPM and transforming to mRNA counts...')

        tpm_exists = os.path.isfile(tpm_file)

        if not tpm_exists:
            raise Exception('TPM file is required when processing smart-seq data')

        if len(cell_list) == 0:
            cell_list = psi.columns
        mrna = tpm2mrna(tpm_file, cell_list)
        ##### New thing
        cells = psi.columns.intersection(mrna.columns)
        mrna = mrna[cells]
        psi = psi[cells]

        mrna_per_event = get_mrna_per_event(mrna, psi, reads, intron_mtx_CI, solo=True) #constitutive_sj_file)

    if len(self.adata.obs) > 0:
        idx = self.adata.obs.index.intersection(mrna_per_event.index)
    else:
        idx = mrna_per_event.index

    psi = psi.loc[idx].T
    mrna_per_event = mrna_per_event.loc[idx].T
    ncells_former = mrna_per_event.shape[0]

    mrna_per_event = mrna_per_event[(mrna_per_event.sum(axis=1) > 0.1) & (mrna_per_event.sum(axis=1) < np.inf)]
    ncells_current = mrna_per_event.shape[0]
    if ncells_former > ncells_current:
        n_diff = str(ncells_former - ncells_current)
        print('removed ' + n_diff + ' cells with all missing or "inf" mRNA values.')
        print('This can be the consequence of very shallow coverage in the cell.')
        psi = psi.loc[mrna_per_event.index]
        
    if os.path.isdir(save_files_in):
        psi.T.to_csv(save_files_in + '/psi.tab.gz', sep='\t', 
                   index=True, header=True)
        mrna_per_event.T.to_csv(save_files_in + '/mrna.tab.gz', sep='\t', 
                   index=True, header=True)

    self.adata.uns['psi'] = psi
    self.adata.uns['mrna_per_event'] = mrna_per_event

    print('Successfully processed RNA-seq data')
    
