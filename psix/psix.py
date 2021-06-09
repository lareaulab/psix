import numpy as np
import pandas as pd
import os

from cell_metric import *
from model_functions import psix_score
from mrna_census import *
# from tpm_to_mrna import *
import anndata
# from rnaseq_tools import *

################# from mrna_census import *
from junctions2psi import *

from turbo_tools import *
from solo2psi import *
from modules import local_correlation_plot, compute_modules_function

from statsmodels.stats.multitest import multipletests


from multiprocessing import Pool
from tqdm import tqdm
import itertools  

from matplotlib import pyplot as plt
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["axes.facecolor"] = 'white'
import matplotlib as mpl
import numpy as np
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams['pdf.fonttype'] = 42


import warnings
warnings.filterwarnings("ignore")


class Psix:
    
    def __init__(self, 
                 psi_table = '', 
                 mrna_table = ''
                ):
        self.adata = anndata.AnnData()
        
        if os.path.isfile(psi_table) and os.path.isfile(mrna_table):
            psi = pd.read_csv(psi_table, sep='\t', index_col=0).T
            mrna_per_event = pd.read_csv(mrna_table, sep='\t', index_col=0).T
            
            self.adata.uns['psi'] = psi
            self.adata.uns['mrna_per_event'] = mrna_per_event
            
            
    def junctions2psi(
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
        solo = False,
        save_files_in = ''
       ):
        
        if (len(save_files_in) > 0) and not os.path.isdir(save_files_in):
            os.mkdir(save_files_in)
            
        if solo:
            solo_to_psi(
                self,
                solo_dir = sj_dir,
                intron_file = intron_file,
                tpm_file = tpm_file,
                cell_list = cell_list,
                minJR = minJR,
                minCell = minCell,
                minPsi = minPsi,
                min_observed = min_observed,
                tenX = tenX
               )
            
        else:
            junctions_dir_to_psi(
                self,
                sj_dir = sj_dir,
                intron_file = intron_file,
                tpm_file = tpm_file,
                cell_list = cell_list,
                minJR = minJR,
                minCell = minCell,
                minPsi = minPsi,
                min_observed = min_observed,
                tenX = tenX,
                save_files_in = save_files_in
               )
        
#         if (len(save_files_in) > 0) and not os.path.isdir(save_files_in):
#             os.mkdir(save_files_in)

#         print('Collecting splice junctions...')

#     #     if os.path.isfile(sj_file):
#     #         sj_file = pd.read_csv(sj_file, sep='\t', index_col=0)

#     #     else:
#         sj_file = process_SJ_dir(sj_dir,
#                                  intron_file,
#                                  save_files_in = save_files_in,
#                                  cell_list = cell_list
#                                 )

#         print('Obtaining PSI tables...')

#         psi, reads, constitutive_sj = get_psi_table(sj_file, minJR, minCell, tenX=tenX)

#         alt_exons = psi.index[np.abs(0.5 - psi.mean(axis=1)) <= (0.5-minPsi)]
#         obs_exons = psi.index[psi.isna().mean(axis=1) <= 1-min_observed]
#         selected_exons = alt_exons & obs_exons

#         psi = psi.loc[selected_exons]
#         reads = reads.loc[selected_exons]

#         if tenX:
#             mrna_per_event = reads
#         else:

#             print('Reading TPM and transforming to mRNA counts...')

#             tpm_exists = os.path.isfile(tpm_file)
#     #         constitutive_sj_exists = os.path.isfile(constitutive_sj_file)

#     #         if not (tpm_exists and constitutive_sj_exists):
#     #             raise Exception('TPM file and constitutive junctions are required when processing smart-seq data')

#             if len(cell_list) == 0:
#                 cell_list = psi.columns
#             mrna = tpm2mrna(tpm_file, cell_list)
#             ##### New thing
#             cells = psi.columns & mrna.columns
#             mrna = mrna[cells]
#             psi = psi[cells]
#             constitutive_sj = constitutive_sj[cells]

#     #         constitutive_sj = pd.read_csv(constitutive_sj_file, sep='\t', index_col=0)
#             mrna_per_event = get_mrna_per_event(mrna, psi, reads, constitutive_sj, solo=False) #constitutive_sj_file)

#         if len(self.adata.obs) > 0:
#             idx = self.adata.obs.index & mrna_per_event.index
#         else:
#             idx = mrna_per_event.index

#         if os.path.isdir(save_files_in):
#             psi.loc[idx].to_csv(save_files_in + '/psi.tab.gz', sep='\t', 
#                        index=True, header=True)
#             mrna_per_event.loc[idx].to_csv(save_files_in + '/mrna.tab.gz', sep='\t', 
#                        index=True, header=True)

#         self.adata.uns['psi'] = psi.loc[idx].T
#         self.adata.uns['mrna_per_event'] = mrna_per_event.loc[idx].T

#         print('Successfully processed RNA-seq data')

        
#     def process_rnaseq(
#         self,
#         solo_dir = '',
#         intron_tab = '',
#         cell_list = [],
#         sj_file = '',
#         tpm_file = '',
#         minJR = 1,
#         minCell = 1,
#         drop_duplicates = False,
#         minPsi = 0.05,
#         min_observed = 0.25,
#         tenX = False,
#     ):
#         process_rnaseq_files(
#             self,
#             sj_file = sj_file,
#             tpm_file = tpm_file,
#             minJR = minJR,
#             minCell = minCell,
#             drop_duplicates = drop_duplicates,
#             minPsi = minPsi,
#             min_observed = min_observed,
#             tenX = tenX
#         )
        
#     def process_rnaseq_solo(
#         self,
#         solo_dir,
#         intron_tab,
#         cell_list = [],
#         tpm_file = '',
#         minJR = 1,
#         minCell = 1,
#         minPsi = 0.05,
#         min_observed = 0.25,
#         tenX = False
#     ):
#         process_rnaseq_solo(
#                 self,
#                 solo_dir,
#                 intron_tab,
#                 cell_list,
#                 tpm_file = tpm_file,
#                 minJR = minJR,
#                 minCell = minCell,
#                 minPsi = minPsi,
#                 min_observed = min_observed,
#                 tenX = tenX
#             )
        
#     def process_rnaseq_solo(
#         self,
#         solo_dir,
#         intron_tab,
#         cell_list = [],
#         tpm_file = '',
#         minJR = 1,
#         minCell = 1,
#         minPsi = 0.05,
#         min_observed = 0.25,
#         tenX = False
#     ):

#         print('Processing STARsolo output. This might take a few minutes...')
#         intron_mtx_CI, intron_mtx_exons = process_solo(solo_dir, intron_tab, cell_list)

#         print('Obtaining PSI tables...')

#         psi, reads = get_psi_table_solo(intron_mtx_exons, minJR, minCell, tenX=tenX)

#         alt_exons = psi.index[np.abs(0.5 - psi.mean(axis=1)) <= (0.5-minPsi)]
#         obs_exons = psi.index[psi.isna().mean(axis=1) <= 1-min_observed]
#         selected_exons = alt_exons & obs_exons

#         psi = psi.loc[selected_exons]
#         reads = reads.loc[selected_exons]

#         if tenX:
#             mrna_per_event = reads
#         else:

#             print('Reading TPM and transforming to mRNA counts...')

#             tpm_exists = os.path.isfile(tpm_file)

#             if not tpm_exists:
#                 raise Exception('TPM file is required when processing smart-seq data')

#             mrna = tpm_to_mrna2(tpm_file)
#             ##### New thing
#             cells = psi.columns & mrna.columns
#             mrna = mrna[cells]
#             psi = psi[cells]

#             mrna_per_event = get_mrna_per_event(mrna, psi, reads, intron_mtx_CI, solo=True) #constitutive_sj_file)

#         if len(self.adata.obs) > 0:
#             idx = self.adata.obs.index & mrna_per_event.index
#         else:
#             idx = mrna_per_event.index

#         self.adata.uns['psi'] = psi.loc[idx].T
#         self.adata.uns['mrna_per_event'] = mrna_per_event.loc[idx].T

#         print('Successfully processed RNA-seq data')
        
    
        
    def get_cell_metric(self, latent, n_neighbors = 100, weight_metric=True):
        
        if type(latent) == str:
            self.latent = pd.read_csv(latent, sep='\t', index_col=0).loc[self.adata.uns['psi'].index]
        else:
            self.latent = latent
        
        self.metric = compute_cell_metric(self.latent, 
                                          n_neighbors = n_neighbors, 
                                          weight_metric = weight_metric
                                         )        
        print('Successfully computed cell-cell metric')
            
                
    def compute_neighbors_psi(self, latent, n_neighbors=100, weight_metric=True):
        if not hasattr(self, 'metric'):
            try:
                self.get_cell_metric(latent=latent, n_neighbors = n_neighbors, weight_metric=weight_metric)
            except:
                raise Exception('Latent space "' + latent +'" does not exist.')
        get_background(self)
                
        print('Successfully computed neighbors')
        
    
            
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
            
            exon_score_array = [psix_score(np.array(self.adata.uns['psi'][exon]), 
                                           np.array(self.adata.uns['mrna_per_event'][exon]), 
                                           self.metric, 
                                           capture_efficiency = self.capture_efficiency, 
                                           randomize = False,  
                                           min_probability = self.min_probability,
                                           seed=0,
                                           turbo = self.turbo
                                        ) for exon in tqdm(self.adata.uns['psi'].columns, position=0, leave=True)]

            
        else:
            
            with Pool(
                processes=self.n_jobs
            ) as pool:
                
                chunksize = np.int(np.ceil(len(exon_list)/self.n_jobs))

                exon_score_array = list(
                    tqdm(
                        pool.imap(self.psix_score_parallel, exon_list, chunksize=chunksize),
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
        
        self.compute_pvalues()
        
        print('Successfully estimated p-values')
        
    
              
    def psix_score_parallel(self, exon):
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
        self.get_bins()

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
                buckets_scores = self.compute_random_exons(bucket)
                self.random_scores[bucket[0]][bucket[1]] = buckets_scores

        else:

            with Pool(
                processes=self.n_jobs
            ) as pool:

                chunksize = np.int(np.ceil((self.pvals_bins**2)/self.n_jobs))

                results = list(
                    tqdm(
                        pool.imap(self.compute_random_exons, all_buckets, chunksize=chunksize),
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
                    pval = self.calculate_exon_pvalue(exon, mean, var)
                    self.psix_results.loc[exon, 'pvals'] = pval
                    
        self.psix_results['qvals'] = multipletests(self.psix_results.pvals, method='fdr_bh')[1]
                   

                
        
    def calculate_exon_pvalue(self, exon, mean, var):
        exon_score = self.psix_results.loc[exon, 'psix_score']
        random_scores_array = np.array(self.random_scores[mean][var])
        total_random_larger = np.sum(random_scores_array >= exon_score)
        total_random = len(random_scores_array)
        
        return (total_random_larger+1)/(total_random+1)
        
        

    def compute_random_exons(self, bucket):

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


    def read_psi(self, psi_file):
        self.adata.uns['psi'] = pd.read_csv(psi_file, sep='\t', index_col=0).T
        
    def read_mrna(self, mrna_file):
        self.adata.uns['mrna_per_event'] = pd.read_csv(mrna_file, sep='\t', index_col=0).T
        
    def save_psix_object(self, psix_dir = 'psix_object/', overwrite = False):
        if os.path.isdir(psix_dir):
            if overwrite:
                sp.run('rm -r ' + psix_dir, shell=True)
            else:
                raise Exception('Directory ' + psix_dir +' already exists')
        os.mkdir(psix_dir)
        self.adata.write(psix_dir+'adata_object', compression='gzip')
        try:
            self.psix_results.to_csv(psix_dir+'psix_results.tab.gz', sep='\t', index=True, header=True)
        except:
            print('No scores to save.')
            
    def read_psix_object(self, psix_dir = 'psix_object/'):
        self.adata = anndata.read_h5ad(psix_dir+'adata_object.gz')
        self.psix_results = pd.read_csv(psix_dir+'psix_results.tab.gz', sep='\t', index_col=0)
        
    
    def compute_modules(self,
                        min_gene_threshold=20, 
                        fdr_threshold=None, 
                        z_threshold=0.3, 
                        core_only=False
                        ):
        
        background_psi = self.adata.uns['neighbors_psi'].mask(self.adata.uns['psi'].isna()).T
        self.sig_exons = self.psix_results.loc[(self.psix_results.psix_score > 0) & (self.psix_results.qvals <= 0.05)].index
        
        self.exon_correlation = self.adata.uns['neighbors_psi'].mask(self.adata.uns['psi'].isna())[self.sig_exons].corr().fillna(0)
        
        self.modules, self.linkage = compute_modules_function(self.exon_correlation,
                                     min_gene_threshold=min_gene_threshold, 
                                     fdr_threshold=fdr_threshold, 
                                     z_threshold=z_threshold, 
                                     core_only=core_only
                                     )
        
    def plot_correlation_modules(self,
                                 z_cmap='RdBu_r', yticklabels=False,
                                 plot_name = '',
                                ):
        
        return local_correlation_plot(self.exon_correlation, 
                               self.modules, self.linkage,
                               z_cmap=z_cmap, yticklabels=yticklabels,
                               plot_name = plot_name
                              )
        
        
