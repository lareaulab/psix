import numpy as np
from functools import partial
import pandas as pd
import time as t
from multiprocessing import Pool
import scipy.special as special
from sklearn.utils import shuffle
import multiprocessing as mp
from itertools import combinations
from tqdm import tqdm
import argparse
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
import os
        
import psix_functions as pr


import warnings
warnings.filterwarnings('ignore')


if __name__ == '__main__':
    
    data='/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'
    
    psi_table = data + 'tiklova_neurogenesis/skipped_exons_psi.tab'
    mrna_table = data + 'tiklova_neurogenesis/mrna_per_event.tab'
    rd = data + 'tiklova_neurogenesis/rd_pc2.tab'
    
    rd = pd.read_csv(rd, sep='\t', index_col=0)
    psi_table = pd.read_csv(psi_table, sep='\t', index_col=0)[rd.index]
    mrna_table = pd.read_csv(mrna_table, sep='\t', index_col=0)[rd.index]
    
    times = 1000
    c = 0.1
    m = 0.05
    n = 0.25
    min_probability = 0.01
    k = 100
    pv = 1
    capture_var =  0
    
    if k == 0:
        k = int(round(np.sqrt(len(rd.index))))
    
    # exons with PSI between m and 1-m
    exons = psi_table.index[np.abs(0.5 - psi_table.mean(axis=1)) <= (0.5-m)]
    # exons with PSI observations in at least n % of the cells
    exons = exons & psi_table.index[psi_table.isna().mean(axis=1) <= (1-n)]
    
    
    ### This should be just to speed things up
    psi_table = psi_table.loc[exons]
    mrna_table = mrna_table.loc[exons]
    ###
    
    
    W = pr.get_distance_matrix(rd, k=k)
    
    L_df = pd.DataFrame()
    
    
    weight_observations = False
    
    print('psix_v1')
    
    exon_score_array_v1 = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v1') for exon in tqdm(exons)]
    
    L_df['psix_v1'] = exon_score_array_v1
    
    L_df.index = exons
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    print('psix_v1_discrete')
    
    exon_score_array_v1_discrete = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v1_discrete') for exon in tqdm(exons)]
    
    L_df['psix_v1_discrete'] = exon_score_array_v1_discrete
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    print('psix_v2')
    
    exon_score_array_v2 = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v2') for exon in tqdm(exons)]
    
    L_df['psix_v2'] = exon_score_array_v2
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    print('psix_v2_discrete')
    
    exon_score_array_v2_discrete = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v2_discrete') for exon in tqdm(exons)]
    
    L_df['psix_v2_discrete'] = exon_score_array_v2_discrete
    L_df.index = exons
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    
    ####
    
    weight_observations = True
    
    print('psix_v1_wo')
    
    exon_score_array_v1_wo = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v1') for exon in tqdm(exons)]
    
    L_df['psix_v1_wo'] = exon_score_array_v1_wo
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    print('psix_v1_discrete_wo')
    
    exon_score_array_v1_discrete_wo = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v1_discrete') for exon in tqdm(exons)]
    
    L_df['psix_v1_discrete_wo'] = exon_score_array_v1_discrete_wo
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    print('psix_v2_wo')
    
    exon_score_array_v2_wo = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v2') for exon in tqdm(exons)]
    
    L_df['psix_v2_wo'] = exon_score_array_v2_wo
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)
    
    print('psix_v2_discrete_wo')
    
    exon_score_array_v2_discrete_wo = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                            exon, 0, c, True, False, 0, min_probability, pv, capture_var, times, 
                                                weight_observations, 'psix_v2_discrete') for exon in tqdm(exons)]
    
    L_df['psix_v2_discrete_wo'] = exon_score_array_v2_discrete_wo
    
    ####
    
    

    
    L_df.to_csv('test/psix.scores.txt', sep='\t', index=True, header=True)


    print('Successfully finished running Psix.')

