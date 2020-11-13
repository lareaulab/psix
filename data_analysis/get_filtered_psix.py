import numpy as np
import os
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from scipy.stats import zscore
from scipy.special import logit
import seaborn as sns
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from IPython.core.pylabtools import figsize
from sklearn.neighbors import NearestNeighbors

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from scipy.stats import zscore
from scipy.special import logit
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

import seaborn as sns
from tqdm import tqdm

import sys
sys.path.append('/mnt/lareaulab/cfbuenabadn/psix/utils/')
import psix_functions as pr

from statsmodels.stats.multitest import multipletests

def get_distance_matrix(pca, k=100):
    '''
    This is the same function used for Psix.
    '''
    
    nbrs = NearestNeighbors(n_neighbors=k).fit(pca)
    distances, indices = nbrs.kneighbors(pca)
    
    cells = list(pca.index)
    
    W = pd.DataFrame(np.zeros((len(cells), len(cells))))
    W.columns = cells
    W.index = cells
    
    for i in tqdm(range(len(cells))):
        cell_i = cells[i]
        sigma = np.max(distances[i])
        for j in range(1, len(distances[i])):
            cell_j = cells[indices[i][j]]
            d = distances[i][j]
            w = np.exp(-(d**2)/(sigma**2))        
            W.loc[cell_i, cell_j] = w
    
    return W

### Rank comparison functions

def get_ranks(scores_df, score_name, bins=50):
    score_series = scores_df.sort_values(score_name).sort_values(score_name)
    ranks_dict = dict()
    total_exons = len(scores_df.index)
    exons_per_rank = int(total_exons/bins)
    for i in range(bins):
        bin_str = i+1
        bin_lims = -exons_per_rank*i
        if i == bins - 1:
            rank_exons = scores_df.sort_values(score_name).index[:bin_lims]
        elif i == 0:
            rank_exons = scores_df.sort_values(score_name).index[bin_lims-exons_per_rank:]
        else:
            
            rank_exons = scores_df.sort_values(score_name).index[bin_lims-exons_per_rank:bin_lims]
        ranks_dict.update({bin_str:rank_exons})
    return ranks_dict
    

def make_comparison(df1, df2, score1, score2, bins=50):
    
    rank_df = pd.DataFrame()
    rank_df_exons = pd.DataFrame()
    
    ranks1 = get_ranks(df1, score1, bins)
    ranks2 = get_ranks(df2, score2, bins)
    
    for i in range(1, bins+1):
        shared_list = []
        shared_exons_list = []
        r1_i = ranks1[i]
        for j in range(1, bins+1):
            r2_j = ranks2[j]
            
            shared = len(r1_i & r2_j)
            
            shared_list.append(len(r1_i & r2_j))
            shared_exons_list.append(r1_i & r2_j)
            
        rank_df['bin_' + str(i)] = shared_list
        rank_df_exons['bin_' + str(i)] = shared_exons_list
        
    rank_df.index = ['bin_' + str(i) for i in range(1, bins+1)]
    rank_df_exons.index = ['bin_' + str(i) for i in range(1, bins+1)]
    
    return rank_df, rank_df_exons

def plot_ranks(df1, df2, score1, score2, bins=20, name1='', name2='', title=''):
    
    shared_idx = df1.index & df2.index
    
    if (name1=='') or (name1==''):
        name1 = score1
        name2 = score2
    
    comparison, comparison_exons = make_comparison(df1.loc[shared_idx], df2.loc[shared_idx],
                             score1, score2, bins=bins)
    mask = comparison == 0

    g = sns.heatmap(comparison, mask=mask, cmap='viridis', cbar_kws={'label': 'shared exons'}, 
                    yticklabels=False, xticklabels=False)
    g.set_facecolor('lightgray')

    # g.tick_params(labelsize=0)
    g.figure.axes[-1].tick_params(labelsize=18)
    g.figure.axes[-1].yaxis.label.set_size(18)

    plt.title(title, fontsize=18)
    plt.xlabel(name1, fontsize=18)
    plt.ylabel(name2, fontsize=18)
    plt.show()
    
    return comparison, comparison_exons


def plot_event(rd, psi, exon, size = (5.5, 3)):
    figsize(size[0],size[1])
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    sc = ax.scatter(rd.PC_1, rd.PC_2, c=psi.loc[exon])
    cb = plt.colorbar(sc)
    
    cb.set_label(label='observed $\hat{\Psi}$',size=18)
    cb.ax.tick_params(labelsize=18, length=5)
    cb.outline.set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


    ax.tick_params(labelsize=18, length=5)

    plt.ylabel('PC2', fontsize=18)
    plt.xlabel('PC1', fontsize=18)
    plt.title(exon, fontsize=22)


    plt.show()
    
    
def plot_event_average(rd, psi, exon, W, size = (5.5, 3)):
    
    
    x = psi.loc[exon].dropna().index
    avg = pd.DataFrame(np.array(W.loc[x, x]) *np.array(psi.loc[exon, x])).sum(axis=1)/np.array(W.loc[x, x].sum(axis=1))
    
    figsize(size[0],size[1])
    fig = plt.figure()
    ax  = plt.subplot(1,1,1)
    sc = ax.scatter(rd.loc[x].PC_1, rd.loc[x].PC_2, c=avg, vmin=0, vmax=1)
    cb = plt.colorbar(sc)
    
    cb.set_label(label='observed $\hat{\Psi}$',size=18)
    cb.ax.tick_params(labelsize=18, length=5)
    cb.outline.set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


    ax.tick_params(labelsize=18, length=5)

    plt.ylabel('PC2', fontsize=18)
    plt.xlabel('PC1', fontsize=18)
    plt.title(exon, fontsize=22)


    plt.show()


def get_averages_bulk(bulk_tab):
    samples = bulk_tab.columns
    conditions = []
    for s in samples:
        if not s.split('_')[0] in conditions:
            conditions.append(s.split('_')[0])
        
    diff_df = pd.DataFrame()
    condition_1 = [x for x in samples if x.split('_')[0] == conditions[0]]
    for c in conditions[1:]:
        condition_c = [x for x in samples if x.split('_')[0] == c]
        
        diff_df[c] = np.abs(np.array(bulk_tab[condition_1].mean(axis=1)) - np.array(bulk_tab[condition_c].mean(axis=1)))
        
    diff_df.index = bulk_tab.index
        
    return diff_df


#### UPGMA functions


def get_discrepant_exons(comparison_exons, top = 10):
    x = len(comparison_exons.index)
    rows_good = pd.Index([])
    cols_good = pd.Index([])
    for bin_good in ['bin_'+str(x) for x in range(1, top+1)]:
        for bin_bad in ['bin_'+str(x) for x in range(x-top+1, x+1)]:

            rows_good = rows_good | comparison_exons.loc[bin_good, bin_bad]
            cols_good = cols_good | comparison_exons.loc[bin_bad, bin_good]

    
    return rows_good, cols_good
            
    
def make_bins(psi_table, mrna_table, exons, b=5):
    express_mean = np.log10(mrna_table.loc[exons]+1).mean(axis=1)
    x_step = (express_mean.max() - express_mean.min())/b
    psi_var = psi_table.loc[exons].var(axis=1)
    p_step = (psi_var.max() - psi_var.min())/b

    return express_mean.min(), x_step, psi_var.min(), p_step

def get_mean_var_bin(express, psi_var, limits, b=5):
    mean_bin = b
    var_bin = b
    for i in range(1, b)[::-1]:
        if express <= (limits[0] + (limits[1]*i)):
            mean_bin = i
            
        if psi_var <= (limits[2] + (limits[3]*i)):
            var_bin = i
            
    return 'mean_' + str(mean_bin) + '_var_' + str(var_bin)


def process_pval_permutations(pval_dir):
    permutation_df = pd.DataFrame()
    list_permutations = os.listdir(pval_dir)
    for lista in list_permutations:
        fh = [float(x.rstrip()) for x in open(pval_dir + '/' + lista).readlines()]
        permutation_df[lista.split('.')[0]] = fh
        
    return permutation_df

def calculate_pval(score, score_list):
    return (np.sum([float(x) > float(score) for x in score_list])+1)/np.sum(len(score_list)+1)


def compare_psix(psi_table, mrna_event, exon_list, all_exons, mrna_min_list, W, pval_dir, b=5):
    
    limits = make_bins(psi_table, mrna_event, all_exons, b=b)
    pvals_df = process_pval_permutations(pval_dir)
    
    total_exons = len(exon_list)
    total_mrnas = len(mrna_min_list)
    
    psix_comparison_df = pd.DataFrame(np.zeros((total_exons, total_mrnas)))
    psix_comparison_df.index = exon_list
    psix_comparison_df.columns = mrna_min_list
    
    psix_pval_df = pd.DataFrame(np.zeros((total_exons, total_mrnas)))
    psix_pval_df.index = exon_list
    psix_pval_df.columns = mrna_min_list
    
    for exon in tqdm(exon_list):
        for mrna in mrna_min_list:
            
            if np.sum(mrna_event.loc[exon] >= mrna) < 100:
                psix_score = np.nan
                psix_pval = np.nan
                
            else:
            
                bin_name = get_mean_var_bin(np.log10(mrna_event.mask(mrna_event < mrna).loc[exon]+1).mean(), 
                                 psi_table.mask(mrna_event < mrna).loc[exon].var(), limits, b=b)
                
                if bin_name not in pvals_df.columns:
                    print('this should rarely happen')
                    bin_name = get_mean_var_bin(np.log10(mrna_event.loc[exon]+1).mean(), 
                                 psi_table.loc[exon].var(), limits, b=b)
                
                psix_score = pr.calculate_exon_L(psi_table.mask(mrna_event < mrna), W, 
                    mrna_event.mask(mrna_event < mrna), exon, k = 0, c = 0.1, 
                    weight_distance=True, randomize = False, seed=0, 
                    min_probability = 0.01, approximate=10)
                
                psix_pval = calculate_pval(psix_score, pvals_df[bin_name])
            
            psix_comparison_df.loc[exon, mrna] = psix_score
            psix_pval_df.loc[exon, mrna] = psix_pval
            
    return psix_comparison_df, psix_pval_df


def compare_kw(psi_table, mrna_event, exon_list, mrna_min_list, clusters):
    
    n_clusters = len(clusters.unique())
    
    total_exons = len(exon_list)
    total_mrnas = len(mrna_min_list)
    
    kw_comparison_df = pd.DataFrame(np.zeros((total_exons, total_mrnas)))
    kw_comparison_df.index = exon_list
    kw_comparison_df.columns = mrna_min_list
    
    kw_comparison_df_p = pd.DataFrame(np.zeros((total_exons, total_mrnas)))
    kw_comparison_df_p.index = exon_list
    kw_comparison_df_p.columns = mrna_min_list
    
    for exon in tqdm(exon_list):
        for mrna in mrna_min_list:
            
            if np.sum(mrna_event.loc[exon] >= mrna) < 100:
                kw_score = (np.nan, np.nan)
            else:
                kw_score = test_exon_bimodal_anova(psi_table.mask(mrna_event < mrna), 
                                                   exon, clusters, obs_min=0, linearize=False, n=n_clusters)

                kw_comparison_df.loc[exon, mrna] = kw_score[0]
                kw_comparison_df_p.loc[exon, mrna] = kw_score[1]
            
    return kw_comparison_df, kw_comparison_df_p



def cluster_rd(rd, n_clusters):
    ac = AgglomerativeClustering(n_clusters=n_clusters)
    ac_clusters = ac.fit_predict(rd[['PC_1', 'PC_2']])
    rd['clusters'] = ac_clusters
    return rd





if __name__ == '__main__':
    
    
    data_dir = '/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'
    tiklova_mrna_event = pd.read_csv(data_dir + 'tiklova_neurogenesis/mrna_per_event.tab', sep='\t', index_col=0)
    tiklova_rd = pd.read_csv(data_dir + 'tiklova_neurogenesis/rd_pc2.tab', sep='\t', index_col=0)
    tiklova_PSI = pd.read_csv(data_dir + 'tiklova_neurogenesis/skipped_exons_psi.tab', sep='\t', index_col=0)
    
    tiklova_psix = pd.read_csv('/mnt/lareaulab/cfbuenabadn/psix/psix_runs/tiklova_neurogenesis.scores.txt', 
                               sep='\t', index_col=0)
    
    tiklova_psix_reduced_var = pd.read_csv(
        '/mnt/lareaulab/cfbuenabadn/psix/development/psix_runs/tiklova_neurogenesis_reduced_var.scores.txt', 
                               sep='\t', index_col=0)
    tiklova_cross_psix = pd.read_csv('../psix_runs/tiklova_neurogenesis.cross_scores.tab', sep='\t', index_col=0)
    tiklova_kw = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/tiklova_neurogenesis_kw.tab', sep='\t', index_col=0)
    tiklova_kw = tiklova_kw.loc[tiklova_psix.index]
    tiklova_kw['qvals'] = multipletests(tiklova_kw.pvals, method='fdr_bh')[1]
    
    W_tiklova = get_distance_matrix(tiklova_rd[['PC_1', 'PC_2']], k=100)


    
    
    psix_comparison, psix_comparison_pvals = compare_psix(tiklova_PSI.loc[tiklova_psix.index], 
                                                      tiklova_mrna_event.loc[tiklova_psix.index], 
                                                      tiklova_psix.index,  tiklova_psix.index, range(21), W_tiklova,
                          '/mnt/lareaulab/cfbuenabadn/psix/psix_runs/tiklova_neurogenesis_pvals/',
                                                        b=5)

    psix_comparison.to_csv('tables/tiklova_psix_comparison.tab', sep='\t', index=True, header=True)
    psix_comparison_pvals.to_csv('tables/tiklova_psix_comparison_pvals.tab', sep='\t', index=True, header=True)
    
    
#     tiklova_rd = cluster_rd(tiklova_rd, 5)
    
#     kw_comparison, kw_comparison_pvals = compare_kw(tiklova_PSI, tiklova_mrna_event, tiklova_psix.index, 
#                            range(21), tiklova_rd.clusters)
    
    
#     kw_comparison.to_csv('tables/tiklova_kw_comparison.tab', sep='\t', index=True, header=True)
#     kw_comparison_pvals.to_csv('tables/tiklova_kw_comparison_pvals.tab', sep='\t', index=True, header=True)
    