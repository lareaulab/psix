import numpy as np
import pandas as pd
import subprocess as sp   
import os
from tqdm import tqdm
import gzip



def process_rnaseq_dir(rnaseq_dir,
                       intron_file,
                       get_tpm = False,
                       ensembl2gene = '',
                       save_files_in = ''):
    
    splice_junction_table = get_SJ_table(rnaseq_dir, intron_file)
    
    if os.path.isdir(save_files_in):
        splice_junction_table.to_csv(save_files_in + '/splice_junctions.tab.gz', sep='\t', header=True, index=True)
    
    if get_tpm:
        tpm_table = get_rsem_table(rnaseq_dir, ensembl2gene)
        
        if os.path.isdir(save_files_in):
            tpm_table.to_csv(save_files_in + '/gene_tpm.tab.gz', sep='\t', header=True, index=True)
            
        return splice_junction_table, tpm_table
    
    return splice_junction_table



def read_intron_file(intron_file):
    intron_table = pd.read_csv(intron_file, sep='\t', index_col=0)
    return intron_table[['intron']]


def process_SJ_table(rnaseq_dir, cell):
    
    sj_file = rnaseq_dir + '/' + cell + '/star_output/SJ.out.tab.gz'
    
    idx = []
    sj_counts = []
    with gzip.open(sj_file, 'rb') as fh:#open(sj_file, 'r') as fh:
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


def get_SJ_table(rnaseq_dir, intron_file):
    cells = os.listdir(rnaseq_dir)
    series_list = []
    for cell in tqdm(cells, position = 0, leave=True):
        cell_series = process_SJ_table(rnaseq_dir, cell)
        series_list.append(cell_series)
        
    splice_junction_table = pd.concat(series_list, axis=1) 
    intron_table = read_intron_file(intron_file)
    splice_junction_table = pd.merge(intron_table, splice_junction_table, 
                                     left_on='intron', how='left', right_index=True).fillna(0).drop('intron', axis=1)
    
    return splice_junction_table


def get_rsem_table(rnaseq_dir, ensembl2gene):
    cells = os.listdir(rnaseq_dir)
    
    series_list = []
    for cell in tqdm(cells, position = 0, leave=True):
        cell_series = process_rsem_table(rnaseq_dir, cell)
        series_list.append(cell_series)
    
    tpm_table = pd.concat(series_list, axis=1)
    
    if os.path.isfile(ensembl2gene):
        
        ensembl2gene_dir = pd.read_csv(ensembl2gene, 
                                   sep='\t').dropna().drop_duplicates()

        ensembl2gene_dir = ensembl2gene_dir.set_index('ensembl_gene_id').drop_duplicates()
        ensembl2gene_dir = ensembl2gene_dir[~ensembl2gene_dir.index.duplicated(keep='first')]

        shared_idx = ensembl2gene_dir.index & tpm_table.index

        tpm_table = tpm_table.loc[shared_idx]
        tpm_table.index = ensembl2gene_dir.loc[shared_idx].	mgi_symbol
    
    tpm_table.sort_index(inplace=True)
    
    return tpm_table
    
        
def process_rsem_table(rnaseq_dir, cell):
    tpm_file = rnaseq_dir + '/' + cell + '/rsem_output/rsem_output.genes.results.gz'
    tpm_tab = pd.read_csv(tpm_file, sep='\t', index_col=0)
    cell_series = pd.Series(np.array(tpm_tab.TPM), index=[x.split('.')[0] for x in tpm_tab.index], name=cell)
    return cell_series