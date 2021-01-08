import numpy as np
from functools import partial
import pandas as pd
import time as t
from multiprocessing import Pool
import scipy.special as special
from sklearn.utils import shuffle
# import sys
import psix_functions as pr
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

parser.add_argument('-sp', '--skip_pvalues', action='store_true', required=False,
                   help='Skip permutation steps to calculate p-values. Default: do permutations.')

parser.add_argument('-o', '--output_name', type=str, required=False, default='psix_output',
                   help='Name for output files. Files will be <output_name>.scores.txt; and a directory <output_name>_pvals/ (only if permutations are performed). Careful: Psix will overwrite any files with the same name.')

parser.add_argument('-c', '--capture_efficiency', type=float, required=False, default=0.1,
                   help='Average capture efficiency for the probabilistic model. Default: c = 0.1.')

parser.add_argument('-m', '--min_psi', type=float, required=False, default=0.05,
                   help='Consider only exons with avergae global PSI between m and 1-m. Default: m = 0.05.')

parser.add_argument('-p', '--permutations', type=int, required=False, default = 1000,
                   help='How many permutations per bin. Default: 1000 permutations.')

parser.add_argument('-b', '--bins', type=int, required=False, default = 5,
                   help='How many bins to split exons for permutations. Default: 5 bins.')

parser.add_argument('-t', '--threads', type=int, required=False, default = 1,
                   help='Threads (to run in parallel). Default: 1 thread.')

parser.add_argument('-k', '--k_nearest_neighbors', type=int, required=False, default = 0,
                   help='K nearest neighbors. Set k==0 to use the squared root of the number of cells.')

parser.add_argument('-n', '--max_missing', type=float, required=False, default = 0.25,
                   help='Maximum percent of missing values per exon.')

parser.add_argument('-mp', '--min_probability', type=float, required=False, default = 0.01,
                   help='Min probability for very unlikely observations.')

parser.add_argument('-s', '--sum_times', type=float, required=False, default = 10,
                   help='How many times r/c molecules will be counter for the probability.')

parser.add_argument('--return_cell_scores', action='store_true', required=False,
                   help='For each exon, return score for each cell, instead of a single score per exon. Default: single score per exon.')

# parser.add_argument('-a', '--approximate', type=float, required=False, default = 0,
#                    help='If > 0, it will approximate the probabilities of observations by sampling a times.')


if __name__ == '__main__':
    
    args = parser.parse_args()
    
    # Required arguments
    rd = pd.read_csv(args.reduced_dimensions, sep='\t', index_col=0)
    psi_table = pd.read_csv(args.psi_table, sep='\t', index_col=0)[rd.index]
    mrna_table = pd.read_csv(args.mrna_table, sep='\t', index_col=0)[rd.index]
    
    # Optional arguments
    skip_pvalues = args.skip_pvalues
    output_name = args.output_name
    c = args.capture_efficiency
    m = args.min_psi
    p = args.permutations
    b = args.bins
    t = args.threads
    n = args.max_missing
    min_probability = args.min_probability
    s = args.sum_times
    k = args.k_nearest_neighbors
    return_cell_scores = args.return_cell_scores
#     a = args.approximate
    
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
    
    print('Running Psix on:')
    print(str(len(psi_table.columns)) + ' cells.')
    print(str(len(exons)) + ' exons.')
    
    print('...')
    
    print('Calculating distance matrix.')
    W = pr.get_distance_matrix(rd, k=k)
    
    print('...')
    
    print('Running on {t} threads.'.format(t=str(t)))
    
    if return_cell_scores:
        ## This is a temporal fix. If this works for the correlation analysis I'll make it return everything in the future 
        ## when return_cell_scores == True
        
        psix_df = pd.DataFrame()
        
        exon_score_array = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                                exon, 0, c, True, False, 0, min_probability, s, True) for exon in tqdm(exons)]
        
        print('retrieving exon scores.')
        for i in tqdm(range(len(exons))):
            exon_df = pd.DataFrame()
            exon_df[exons[i]] = exon_score_array[i][0]
            exon_df.index = exon_score_array[i][1]
            
            psix_df = pd.concat([psix_df, exon_df], axis=1)
            
        #psix_df.index = W.index
        
        psix_df.T.to_csv(output_name + '.scores_tab.txt', sep='\t', index=True, header=True)
        
    else:
    
        if t > 1:
            pool = mp.Pool(t)

            exon_score_array_mp = [pool.apply_async(pr.calculate_exon_L, args=(psi_table, W, mrna_table,          
                                exon, 0, c, True, False, 0, min_probability, s, False)) for exon in exons]

            exon_score_array = [p.get() for p in tqdm(exon_score_array_mp)]


            pool.terminate()

        else:
            exon_score_array = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                                exon, 0, c, True, False, 0, min_probability, s, False) for exon in tqdm(exons)]

        L_df = pd.DataFrame()
        L_df['L_score'] = exon_score_array
        L_df.index = exons

        print('Successfully ran Psix on exons.')
        print('...')
        print('')
        print('Estimating p-values')

        express_mean = np.log10(mrna_table.loc[exons]+1).mean(axis=1)
        x_step = (express_mean.max() - express_mean.min())/b
        psi_var = psi_table.loc[exons].var(axis=1)
        p_step = (psi_var.max() - psi_var.min())/b

        bins = {}
        for i in range(b):
            mean_bins = {}

            exons_x_small = express_mean.index[express_mean >= express_mean.min() + i*x_step]
            exons_x_large = express_mean.index[express_mean <= express_mean.min() + (i+1)*x_step]
            exons_x = exons_x_small.intersection(exons_x_large)

            if i == (b-1):
                exons_x = exons_x_small

            for j in range(b):
                exons_var_small = psi_var.index[psi_var >= psi_var.min() + j*p_step]
                exons_var_large = psi_var.index[psi_var <= psi_var.min() + (j+1)*p_step]
                exons_var = exons_var_small.intersection(exons_var_large).intersection(exons_x)

                if j == (b-1):
                    exons_var = exons_var_small.intersection(exons_x)

                mean_bins.update({'var_' + str(j+1):exons_var})

            bins.update({'mean_' + str(i+1):mean_bins})

        bin_dir = pd.DataFrame()
        for mean in bins.keys():
            for var in bins[mean].keys():
                for exon in bins[mean][var]:
                    bin_dir[exon] = [mean + '_' + var]

        assert len([x for x in exons if x not in bin_dir.T.index]) == 0

        buckets_dir = output_name + '_pvals/'

        if not skip_pvalues:
            print('Running {p} permutations to obtain p-values. This will take a while.'.format(p=str(p)))

            if not os.path.isdir(buckets_dir):
                print('Making dir')
                os.mkdir(buckets_dir)

            np.random.seed(0)
            print(str(p) + ' randomized')
            for i in range(b):
                for j in range(b):
                    exon_list = bins['mean_'+str(i+1)]['var_'+str(j+1)]
                    if len(exon_list) > 0:


                        print('mean_'+str(i+1) + ', ' + 'var_'+str(j+1))

                        r_choice = np.random.choice(exon_list, p, replace=True)

                        if t > 1:

                            pool = mp.Pool(t)



                            random_score_array_mp = [pool.apply_async(pr.calculate_exon_L, args=(psi_table, W, mrna_table, 
                                             r_choice[exon], 0, c, True, True, exon, min_probability, s)) for exon in range(p)]

                            random_score_array = [p.get() for p in tqdm(random_score_array_mp)]


                            pool.terminate()

                        else:

                            random_score_array = [pr.calculate_exon_L(psi_table, W, mrna_table,          
                                r_choice[exon], 0, c, True, True, 0, min_probability, s) for exon in tqdm(range(p))]


                        fh = open(buckets_dir + 'mean_'+str(i+1)+'_var_'+str(j+1)+'.txt', 'w')
                        fh.writelines([str(x) + '\n' for x in random_score_array])
                        fh.close()



        else:
            print('Skipping permutations. Retrieving permutated scores from ' + buckets_dir)

        print('...')
        print('')

        print('Estimating p-values.')
        print('...')
        print('')

        L_dir = pd.DataFrame()
        for mean in bins.keys():
            for var in bins[mean].keys():
                for exon in bins[mean][var]:
                    L_dir[exon] = [mean + '_' + var]


        L_dir = L_dir.T
        L_dir.columns = ['bin']
        L_dir['L_score'] = np.array(L_df.loc[L_dir.index].L_score)


        bin_randoms = pd.DataFrame()

        for mv_bin in L_dir.bin.unique():
            fh = open(buckets_dir + mv_bin + '.txt', 'r')
            mv = [float(x.rstrip()) for x in fh.readlines()]
            bin_randoms[mv_bin] = mv


        pvals = []
        for idx, row in tqdm(L_dir.iterrows(), position=0, leave=True):
            mv = row.bin
            L = row.L_score
            p = (np.sum(bin_randoms[mv] > L) + 1)/(len(bin_randoms[mv])+1)
            pvals.append(p)

        L_dir['pvals'] = pvals

        L_dir['qvals'] = multipletests(L_dir.pvals, method='fdr_bh')[1]
        pvals = []

        for idx, row in tqdm(L_dir.iterrows(), position=0, leave=True):
            mv = row.bin
            L = row.L_score
            norm_mean, norm_var = norm.fit(bin_randoms[mv])
            p = norm.sf(L, norm_mean, norm_var)
            pvals.append(p)

        L_dir['norm_pvals'] = pvals

        L_dir['norm_qvals'] = multipletests(L_dir.norm_pvals, method='fdr_bh')[1]

        L_dir.to_csv(output_name + '.scores.txt', sep='\t', index=True, header=True)


    print('Successfully finished running Psix.')

