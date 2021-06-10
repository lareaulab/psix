import numpy as np
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests
from model_functions import psix_score
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
import itertools  

def compute_psix_scores(self,
                        n_jobs=1,
                        capture_efficiency = 0.1, 
                        min_probability = 0.01,
                        seed=0,
                        pvals_bins=5,
                        n_random_exons = 2000,
                        latent='latent', 
                        n_neighbors = 100, 
                        weight_metric=True,
                        turbo = False
                       ):

    self.capture_efficiency = capture_efficiency
    self.min_probability = min_probability
    self.seed = seed
    self.pvals_bins = pvals_bins
    self.n_random_exons = n_random_exons
    self.n_jobs = n_jobs

    if turbo:
        self.turbo = load_turbo(turbo_dir = turbo)

    else:
        self.turbo = False

    if not 'cell_metric' in self.adata.uns:
        print('Computing cell-cell metric...', flush=True)
        self.get_cell_metric(latent=latent, n_neighbors = n_neighbors, weight_metric=weight_metric)

    print('Computing Psix score in ' + str(len(self.adata.uns['psi'].columns)) + ' exons')


    exon_list = self.adata.uns['psi'].columns

    if self.n_jobs == 1:

        exon_score_array = [run_psix_exon(exon, self) for exon in tqdm(self.adata.uns['psi'].columns, position=0, leave=True)]


    else:
        
        psix_score_parallel=partial(run_psix_exon, self=self)



        with Pool(
            processes=self.n_jobs
        ) as pool:

            chunksize = np.int(np.ceil(len(exon_list)/self.n_jobs))

            exon_score_array = list(
                tqdm(
                    pool.imap(psix_score_parallel, exon_list, chunksize=chunksize),
#                     pool.imap(self.psix_score_parallel, exon_list, chunksize=chunksize),
                    total=len(exon_list), position=0, leave=True
                )
            )
#              = results

    self.psix_results = pd.DataFrame()
    self.psix_results['psix_score'] = exon_score_array
    self.psix_results.index = exon_list
    self.psix_results = self.psix_results.sort_values('psix_score', ascending=False)

    print('Successfully computed Psix score of exons.')
    print('Estimating p-values. This might take a while...')

    compute_pvalues(self)

    print('Successfully estimated p-values')


    
    
def run_psix_exon(exon, self):
    ##############

    return psix_score(np.array(self.adata.uns['psi'][exon]), 
                                       np.array(self.adata.uns['mrna_per_event'][exon]), 
                                       self.metric, 
                                       capture_efficiency = self.capture_efficiency, 
                                       randomize = False,  
                                       min_probability = self.min_probability,
                                       seed=self.seed,
                                       turbo = self.turbo
                           )




def compute_pvalues(self):
    get_bins(self)

    np.random.seed(self.seed)

    mean_list = ['mean_'+str(i+1) for i in range(self.pvals_bins)]
    var_list = ['var_'+str(i+1) for i in range(self.pvals_bins)]
    self.random_scores = {}
    for m in mean_list:
        m_dict = {}
        for v in var_list:
            m_dict.update({v:[]})
        self.random_scores.update({m:m_dict})

    all_buckets = list(itertools.product(mean_list, var_list))

    if self.n_jobs == 1:
        for bucket in all_buckets:
            buckets_scores = compute_random_exons(bucket, self)
            self.random_scores[bucket[0]][bucket[1]] = buckets_scores

    else:
        
        random_exons_parallel=partial(compute_random_exons, self=self)

        with Pool(
            processes=self.n_jobs
        ) as pool:

            chunksize = np.int(np.ceil((self.pvals_bins**2)/self.n_jobs))

            results = list(
                tqdm(
                    pool.imap(random_exons_parallel, all_buckets, chunksize=chunksize),
                    total=len(all_buckets), position=0, leave=True
                )
            )
        self.exon_score_array = results

        for i in range(self.pvals_bins**2):
            bucket = all_buckets[i]
            self.random_scores[bucket[0]][bucket[1]] = list(results[i])

    self.psix_results['pvals'] = [1]*len(self.psix_results)

    for mean in self.bins.keys():
        for var in self.bins[mean].keys():
            for exon in self.bins[mean][var]:
                pval = calculate_exon_pvalue(self, exon, mean, var)
                self.psix_results.loc[exon, 'pvals'] = pval

    self.psix_results['qvals'] = multipletests(self.psix_results.pvals, method='fdr_bh')[1]




def calculate_exon_pvalue(self, exon, mean, var):
    exon_score = self.psix_results.loc[exon, 'psix_score']
    random_scores_array = np.array(self.random_scores[mean][var])
    total_random_larger = np.sum(random_scores_array >= exon_score)
    total_random = len(random_scores_array)

    return (total_random_larger+1)/(total_random+1)



def compute_random_exons(bucket, self):

    exon_list = self.bins[bucket[0]][bucket[1]]

    if len(exon_list)==0:
        return np.array([])

    else:
        r_choice = np.random.choice(exon_list, self.n_random_exons, replace=True)

        if self.n_jobs == 1:


            ####################
            output = [psix_score(np.array(self.adata.uns['psi'][exon]), 
                       np.array(self.adata.uns['mrna_per_event'][exon]), 
                       self.metric, 
                       capture_efficiency = self.capture_efficiency, 
                       randomize = True,  
                       min_probability = self.min_probability,
                       seed=self.seed,
                       turbo = self.turbo) for exon in tqdm(r_choice, position=0, leave=True)]
            ######################


        else:
            ####################
            output = [psix_score(np.array(self.adata.uns['psi'][exon]), 
                       np.array(self.adata.uns['mrna_per_event'][exon]), 
                       self.metric, 
                       capture_efficiency = self.capture_efficiency, 
                       randomize = True,  
                       min_probability = self.min_probability,
                       seed=self.seed,
                       turbo = self.turbo) for exon in r_choice]
            ####################

        return np.array(output)

def get_bins(self):
    mrna_table = self.adata.uns['mrna_per_event']
    psi_table = self.adata.uns['psi']
    b = self.pvals_bins

    express_mean = np.log10(mrna_table+1).mean(axis=0)
    x_step = (express_mean.max() - express_mean.min())/b
    psi_var = psi_table.var(axis=0)
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

    assert len([x for x in self.adata.uns['psi'].columns if x not in bin_dir.T.index]) == 0
    self.exon_bins = bin_dir
    self.bins = bins

