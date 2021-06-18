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
from score_functions import *
from turbo_tools import *
from solo2psi import *
from modules import local_correlation_plot, compute_modules_function



from deprecated_rnaseq_tools import *




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
                 psix_object = '',
                 psi_table = '', 
                 mrna_table = ''
                ):
        
        
        if os.path.isdir(psix_object):
            self.adata = anndata.read_h5ad(psix_object+'/adata')
            self.psix_results = pd.read_csv(psix_object+'/psix_results.tab.gz', sep='\t', index_col=0)
            
        else:
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
            solo_to_psi(self, solo_dir = sj_dir, intron_file = intron_file, tpm_file = tpm_file,
                        cell_list = cell_list, minJR = minJR, minCell = minCell, minPsi = minPsi,
                        min_observed = min_observed, tenX = tenX)
            
        else:
            junctions_dir_to_psi(self, sj_dir = sj_dir, intron_file = intron_file, tpm_file = tpm_file,
                cell_list = cell_list, minJR = minJR, minCell = minCell, minPsi = minPsi,
                min_observed = min_observed, tenX = tenX, save_files_in = save_files_in)
        
        
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
        
    
    def run_psix(self,
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
        
        self.n_neighbors = n_neighbors
        
        compute_psix_scores(self,
                        n_jobs=n_jobs,
                        capture_efficiency = capture_efficiency, 
                        min_probability = min_probability,
                        seed=seed,
                        pvals_bins = pvals_bins,
                        n_random_exons = n_random_exons,
                        latent = latent, 
                        n_neighbors = n_neighbors, 
                        weight_metric = weight_metric,
                        turbo = turbo
                       )
        

    def save_psix_object(self, psix_dir = 'psix_object', overwrite = False):
        if os.path.isdir(psix_dir):
            if overwrite:
                sp.run('rm -r ' + psix_dir, shell=True)
            else:
                raise Exception('Directory ' + psix_dir +' already exists. If you want to overwrite it, set overwrite = True')
        os.mkdir(psix_dir)
        self.adata.write(psix_dir+'/adata', compression='gzip')
        try:
            self.psix_results.to_csv(psix_dir+'/psix_results.tab.gz', sep='\t', index=True, header=True)
        except:
            print('No scores to save.')
            
    
    def compute_modules(self,
                        min_gene_threshold=20, 
                        fdr_threshold=None, 
                        z_threshold=0.3, 
                        core_only=False,
                        n_neighbors = 100,
                        weight_metric=True,
                        plot = False,
                        z_cmap='RdBu_r', 
                        yticklabels=False,
                        plot_name = ''
                        ):
        
        if not 'neighbors_psi' in self.adata.uns:
            self.compute_neighbors_psi(self.latent, n_neighbors=self.n_neighbors, weight_metric=True)
        
        background_psi = self.adata.uns['neighbors_psi'].mask(self.adata.uns['psi'].isna()).T
        self.sig_exons = self.psix_results.loc[(self.psix_results.psix_score > 0) & (self.psix_results.qvals <= 0.05)].index
        
        self.exon_correlation = self.adata.uns['neighbors_psi'].mask(self.adata.uns['psi'].isna())[self.sig_exons].corr().fillna(0)
        
        self.modules, self.linkage = compute_modules_function(self.exon_correlation,min_gene_threshold=min_gene_threshold, 
                                     fdr_threshold=fdr_threshold, z_threshold=z_threshold, core_only=core_only)
        
        if plot:
             local_correlation_plot(self.exon_correlation, self.modules, self.linkage,
                                    z_cmap=z_cmap, yticklabels=yticklabels, plot_name = plot_name)
        
        
    #### DEPRECATED
        
    def process_rnaseq(
        self,
        solo_dir = '',
        intron_tab = '',
        cell_list = [],
        exon_sj_file = '',
        constitutive_sj_file = '',
        tpm_file = '',
        minJR = 1,
        minCell = 1,
        drop_duplicates = False,
        min_psi = 0.05,
        min_observed = 0.25,
        tenX = False,
    ):
        process_rnaseq_files(
            self,
            exon_sj_file,
            constitutive_sj_file = constitutive_sj_file,
            tpm_file = tpm_file,
            minJR = minJR,
            minCell = minCell,
            drop_duplicates = drop_duplicates,
            min_psi = min_psi,
            min_observed = min_observed,
            tenX = tenX
        )
#     def plot_correlation_modules(self,
#                                  z_cmap='RdBu_r', 
#                                  yticklabels=False,
#                                  plot_name = ''
#                                 ):
        
#         return local_correlation_plot(self.exon_correlation, 
#                                self.modules, self.linkage,
#                                z_cmap=z_cmap, yticklabels=yticklabels,
#                                plot_name = plot_name
#                               )
        
        
