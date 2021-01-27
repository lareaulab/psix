import numpy as np
import pandas as pd
import os
from cell_metric import compute_cell_metric
from model_functions import psix_score
from tpm_to_mrna import tpm2mrna
import anndata
from smartseq_tools import *


class Psix:
    
    def __init__(
        self,
        adata = anndata.AnnData(),
        reads_file = ''
    ):
        self.adata = adata
        
        if os.path.isfile(reads_file):
            self.adata = anndata.read_csv(reads_file, delimiter='\t', first_column_names=True)
        
    def process_smartseq(
        self,
        exon_sj_file,
        constitutive_sj_file,
        tpm_file,
        minJR = 5,
        minCell=20,
        drop_duplicates = False,
        min_psi = 0.05,
        min_observed = 0.1
    ):
        
        print('Obtaining psi tables...')
            
        psi, reads = get_psi_table(exon_sj_file, minJR, minCell, drop_duplicates)
        
        alt_exons = psi.index[np.abs(0.5 - psi.mean(axis=1)) <= (0.5-min_psi)]
        obs_exons = psi.index[psi.isna().mean(axis=1) < 1-min_observed]
        selected_exons = alt_exons & obs_exons
        
        psi = psi.loc[selected_exons]
        reads = reads.loc[selected_exons]
        
        print('Reading TPM and transforming to mRNA counts...')
        
        mrna = tpm2mrna(tpm_file)
        mrna_per_event = get_mrna_per_event(mrna, psi, reads, constitutive_sj_file)
        
        self.adata.uns['psi'] = psi.T
        self.adata.uns['mrna_per_event'] = mrna_per_event.T
        
        print('Successfully processed smart-seq data')

