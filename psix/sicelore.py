import numpy as np
import pandas as pd
from tqdm import tqdm

def read_sicelore(isomatrix, min_cell_counts = 1e3, max_nan_psi = 0.75):
    isoforms = pd.read_csv(isomatrix, sep='\t')

    isoforms.index = isoforms.geneId + ':' + isoforms.transcriptId
    X = isoforms[isoforms.columns[3:]]
    cell_names = isoforms.columns[3:]
    suma = isoforms.groupby('geneId')[cell_names].sum()
    
    gene_sums = isoforms.groupby('geneId')[isoforms.columns[3:]].transform('sum')
    psi_tab = isoforms[isoforms.columns[3:]]/gene_sums
    
    _isoforms = isoforms.groupby('geneId').size().index[isoforms.groupby('geneId').size() > 1]
    
    isoforms = isoforms.loc[isoforms.geneId.isin(_isoforms)]
    _isoforms_table = isoforms.loc[isoforms.geneId.isin(_isoforms)].sort_values('geneId')
    
    gene_list = []
    transcript_id = []
    
    for idx, df in tqdm(_isoforms_table.groupby('geneId')):
        if (df.transcriptId == 'undef').any():
            transcript = list(df.loc[(df.transcriptId != 'undef')].index)
        else:
            transcript = list(df.index[0:-1])
            
        gene_list.extend([idx]*len(transcript))
        transcript_id.extend(transcript)
    
    gene_list = pd.Index(gene_list)
    transcript_id = pd.Index(transcript_id)
    
    cells = suma.columns[suma.loc[list(set(gene_list))].sum(axis=0) > min_cell_counts]
    
    psi_selected = psi_tab.loc[transcript_id, cells]
    
    psi_selected = psi_selected.loc[psi_selected.isna().mean(axis=1) < max_nan_psi]

    mrna_counts = gene_sums.loc[psi_selected.index, psi_selected.columns]

    return psi_selected, mrna_counts