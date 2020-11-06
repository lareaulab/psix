import numpy as np
from functools import partial
import pandas as pd
import time as t
from multiprocessing import Pool
import scipy.special as special
from sklearn.utils import shuffle
# import sys
import probability_functions as pr
import multiprocessing as mp
from itertools import combinations
from tqdm import tqdm
import argparse
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
import os
        

import warnings
warnings.filterwarnings('ignore')



parser = argparse.ArgumentParser(description='Psix: autocorrelated exon discovery from scRNA-seq data.')

parser.add_argument('-psi', '--psi_table', type=str, required=True,
                   help='Table with exon PSI values per cell. Right now Psix supports only skipped exon events.')

parser.add_argument('-mrna', '--mrna_table', type=str, required=True,
                   help='Table with mRNA counts per splicing event.')

parser.add_argument('-rd', '--reduced_dimensions', type=str, required=True,
                   help='Table with a low-dimensionality projection of the cells (e.g., first 5 principal components).')

parser.add_argument('-p', '--psix_output', type=str, required=True,
                   help='Psix output.')

parser.add_argument('-o', '--output_name', type=str, required=False, default='psix_cross_output',
                   help='Name for output files.')

parser.add_argument('-c', '--capture_efficiency', type=float, required=False, default=0.1,
                   help='Average capture efficiency for the probabilistic model. Default: c = 0.1.')

parser.add_argument('-m', '--min_psi', type=float, required=False, default=0.05,
                   help='Consider only exons with avergae global PSI between m and 1-m. Default: m = 0.05.')

parser.add_argument('-n', '--max_missing', type=float, required=False, default = 0,
                   help='Maximum percent of missing values per exon.')

parser.add_argument('-t', '--threads', type=int, required=False, default = 1,
                   help='Threads (to run in parallel). Default: 1 thread.')

parser.add_argument('-k', '--k_nearest_neighbors', type=int, required=False, default = 0,
                   help='K nearest neighbors. Set k==0 to use the squared root of the number of cells.')

parser.add_argument('-mp', '--min_probability', type=float, required=False, default = 0.01,
                   help='Min probability for very unlikely observations.')

parser.add_argument('-s', '--sum_times', type=float, required=False, default = 10,
                   help='How many times r/c molecules will be counter for the probability.')

parser.add_argument('-l', '--l_min', type=float, required=False, default = 0,
                   help='How many times r/c molecules will be counter for the probability.')

parser.add_argument('-a', '--approximate', type=float, required=False, default = 0,
                   help='If > 0, it will approximate the probabilities of observations by sampling a times.')


if __name__ == '__main__':
    
    args = parser.parse_args()
    
    # Required arguments
    rd = pd.read_csv(args.reduced_dimensions, sep='\t', index_col=0)
    psi_table = pd.read_csv(args.psi_table, sep='\t', index_col=0)[rd.index]
    mrna_table = pd.read_csv(args.mrna_table, sep='\t', index_col=0)[rd.index]
    
    psix_table = pd.read_csv(args.psix_output + '.scores.txt', sep='\t', index_col=0)
    
    
    # Optional arguments
    output_name = args.output_name
    c = args.capture_efficiency
    m = args.min_psi
    t = args.threads
    n = args.max_missing
    min_probability = args.min_probability
    s = args.sum_times
    k = args.k_nearest_neighbors
    l_min = args.l_min
    a = args.approximate
    
    if k == 0:
        k = int(round(np.sqrt(len(rd.index))))
    
    exons = psix_table.loc[(psix_table.L_score > l_min) & (psix_table.qvals <= 0.05)].index
    
    ### This should be just to speed things up
    psi_table = psi_table.loc[exons]
    mrna_table = mrna_table.loc[exons]
    ########################################################
    
    print('Running Psix on:')
    print(str(len(psi_table.columns)) + ' cells.')
    print(str(len(exons)) + ' exons.')
    
    print('...')
    
    print('Calculating distance matrix.')
    W = pr.get_distance_matrix(rd, k=k)
    
    print('...')
    
    print('Running on {t} threads.'.format(t=str(t)))
    
    
    interactions_df = pd.DataFrame(np.identity(len(exons))*np.array(round(psix_table.loc[exons].L_score, 3)))
    interactions_df.columns = exons
    interactions_df.index = exons
    
    exon_comb = list(combinations(exons, r=2))
    
    if t > 1:
    
        pool = mp.Pool(t)

        exon_score_array_mp = [pool.apply_async(pr.calculate_cross_L, args=(psi_table, W, mrna_table,          
                            exon[0], exon[1], 0, c, True, False, 0, min_probability, s, n, a)) for exon in exon_comb]

        exon_score_array = [p.get() for p in tqdm(exon_score_array_mp)]

        pool.terminate()
        
    else:
        
        exon_score_array = [pr.calculate_cross_L(psi_table, W, mrna_table,          
                            exon[0], exon[1], 0, c, True, False, 0, min_probability, s, n, a) for exon in tqdm(exon_comb)]

    print('finished running; storing data...')
    
    for i in range(len(exon_comb)):
        exon = exon_comb[i]
        interactions_df.loc[exon[0], exon[1]] = round(exon_score_array[i][0], 3)
        interactions_df.loc[exon[1], exon[0]] = round(exon_score_array[i][1], 3)
    
    print('saving file...')
    
    interactions_df.to_csv(output_name + '.cross_scores.tab', sep='\t', index=True, header=True)
    
    print('Successfully finished running Psix cross.')

