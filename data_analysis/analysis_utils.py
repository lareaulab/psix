import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from scipy.stats import zscore
from scipy.special import logit
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

import seaborn as sns
from tqdm import tqdm

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

def L_upgma_step(cross_L_matrix, score_min = 0.1):
#     score_min = 0.1
    counter = 1
    for exon in cross_L_matrix.index:
        exon_max = cross_L_matrix.loc[exon].dropna().sort_values().index[-1]
        exon_pair = (exon, exon_max)
        pair_score = cross_L_matrix.loc[exon, exon_max]
        pair_score_r = cross_L_matrix.loc[exon_max, exon]
        if counter == 1:
            max_pair = exon_pair
            max_pair_score = pair_score
            max_pair_score_r = pair_score_r
        else:
            token_1 = (pair_score > max_pair_score) and (pair_score_r >= score_min)
            token_2 = (pair_score >= score_min) and (pair_score_r >= score_min) and (max_pair_score_r < score_min)
            if token_1 or token_2:
                max_pair = (exon, exon_max)
                max_pair_score = cross_L_matrix.loc[exon, exon_max]
                max_pair_score_r = cross_L_matrix.loc[exon_max, exon]
        counter += 1

    if (max_pair_score >= score_min) and (max_pair_score_r >= score_min):
        new_idx = [x for x in cross_L_matrix.index if x not in max_pair]
        new_df = cross_L_matrix.loc[new_idx, new_idx]

        combined_cross_1 = cross_L_matrix.loc[new_idx, list(max_pair)].mean(axis=1)
        combined_cross_2 = cross_L_matrix.loc[list(max_pair), new_idx].mean(axis=0)

        new_df[','.join(max_pair)] = combined_cross_2
        new_df = new_df.T
        new_df[','.join(max_pair)] = list(combined_cross_1) + [0]
        new_df = new_df.T
        
        return new_df
    else:
        print(max_pair_score)
        return [0]
            
            
def get_modules(cross_L_matrix, score_min = 0.0):
    counter = 1
    for i in tqdm(range(len(cross_L_matrix.index))):
        combined_cross = L_upgma_step(cross_L_matrix, score_min)
        if len(combined_cross) == 1:
            return cross_L_matrix
        else:
            cross_L_matrix = combined_cross
            
    return cross_L_matrix

def get_discrepant_exons(comparison_exons, top = 10):
    x = len(comparison_exons.index)
    rows_good = pd.Index([])
    cols_good = pd.Index([])
    for bin_good in ['bin_'+str(x) for x in range(1, top+1)]:
        for bin_bad in ['bin_'+str(x) for x in range(x-top+1, x+1)]:

            rows_good = rows_good | comparison_exons.loc[bin_good, bin_bad]
            cols_good = cols_good | comparison_exons.loc[bin_bad, bin_good]

    
    return rows_good, cols_good
            
        
        

if __name__ == '__main__':
    
    
    data_dir = '/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'
    tiklova_mrna_event = pd.read_csv(data_dir + 'tiklova_neurogenesis/mrna_per_event.tab', sep='\t', index_col=0)
    tiklova_rd = pd.read_csv(data_dir + 'tiklova_neurogenesis/rd_pc2.tab', sep='\t', index_col=0)
    tiklova_PSI = pd.read_csv(data_dir + 'tiklova_neurogenesis/skipped_exons_psi.tab', sep='\t', index_col=0)
    tiklova_psix = pd.read_csv('../psix_runs/tiklova_neurogenesis.scores.txt', sep='\t', index_col=0)
    tiklova_cross_psix = pd.read_csv('../psix_runs/tiklova_neurogenesis.cross_scores.tab', sep='\t', index_col=0)
    tiklova_kw = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/tiklova_neurogenesis_kw.tab', sep='\t', index_col=0)
    tiklova_kw = tiklova_kw.loc[tiklova_psix.index]
    tiklova_kw['qvals'] = multipletests(tiklova_kw.pvals, method='fdr_bh')[1]
    tiklova_geary_C = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/tiklova_autocor_one_matrix/tiklova_GearyC_k100.tab',
                         sep='\t', index_col=0)
    tiklova_geary_C = tiklova_geary_C.loc[tiklova_psix.index & tiklova_geary_C.index]
    tiklova_geary_C.columns = ['C_score', 'pvals']
    
    
    data_dir = '/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'
    chen_mrna_event = pd.read_csv(data_dir + 'chen/mrna_per_event.tab', sep='\t', index_col=0)
    chen_rd = pd.read_csv(data_dir + 'chen/rd_pc2.tab', sep='\t', index_col=0)
    chen_PSI = pd.read_csv(data_dir + 'chen/skipped_exons_psi.tab', sep='\t', index_col=0)
    chen_psix = pd.read_csv('../psix_runs/chen.scores.txt', sep='\t', index_col=0)
    chen_cross_psix = pd.read_csv('../psix_runs/chen.cross_scores.tab', sep='\t', index_col=0)
    chen_kw = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/chen_kw.tab', sep='\t', index_col=0)
    chen_kw = chen_kw.loc[chen_psix.index]
    chen_kw['qvals'] = multipletests(chen_kw.pvals, method='fdr_bh')[1]
    chen_geary_C = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/chen_autocor/chen_GearyC_k30.tab',
                         sep='\t', index_col=0)
    chen_geary_C = chen_geary_C.loc[chen_psix.index & chen_geary_C.index]
    chen_geary_C.columns = ['C_score', 'pvals']
    
    
    data_dir = '/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'
    tiklova_all_mrna_event = pd.read_csv(data_dir + 'tiklova/mrna_per_event.tab', sep='\t', index_col=0)
    tiklova_all_rd = pd.read_csv(data_dir + 'tiklova/rd_pc2.tab', sep='\t', index_col=0)
    tiklova_all_PSI = pd.read_csv(data_dir + 'tiklova/skipped_exons_psi.tab', sep='\t', index_col=0)
    tiklova_all_psix = pd.read_csv('../psix_runs/tiklova.scores.txt', sep='\t', index_col=0)
    tiklova_all_cross_psix = pd.read_csv('../psix_runs/tiklova.cross_scores.tab', sep='\t', index_col=0)
    tiklova_all_kw = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/tiklova_kw.tab', sep='\t', index_col=0)
    tiklova_all_kw = tiklova_all_kw.loc[tiklova_all_psix.index]
    tiklova_all_kw['qvals'] = multipletests(tiklova_all_kw.pvals, method='fdr_bh')[1]
    tiklova_all_geary_C = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/tiklova_autocor_one_matrix/tiklova_GearyC_k100.tab',
                         sep='\t', index_col=0)
    tiklova_all_geary_C = tiklova_all_geary_C.loc[tiklova_all_psix.index & tiklova_all_geary_C.index]
    tiklova_all_geary_C.columns = ['C_score', 'pvals']
    
    
    
    data_dir = '/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'
    song_mrna_event = pd.read_csv(data_dir + 'song/mrna_per_event.tab', sep='\t', index_col=0)
    song_rd = pd.read_csv(data_dir + 'song/rd_pc2.tab', sep='\t', index_col=0)
    song_PSI = pd.read_csv(data_dir + 'song/skipped_exons_psi.tab', sep='\t', index_col=0)
    song_psix = pd.read_csv('../psix_runs/song.scores.txt', sep='\t', index_col=0)
    song_cross_psix = pd.read_csv('../psix_runs/song.cross_scores.tab', sep='\t', index_col=0)
    song_kw = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/song_kw.tab', sep='\t', index_col=0)
    song_kw = song_kw.loc[song_psix.index]
    song_kw['qvals'] = multipletests(song_kw.pvals, method='fdr_bh')[1]
    song_geary_C = pd.read_csv('~/sc_splicing_regulation/sc_neurogenesis/song_autocor/song_GearyC_k30.tab',
                         sep='\t', index_col=0)
    song_geary_C = song_geary_C.loc[song_psix.index & song_geary_C.index]
    song_geary_C.columns = ['C_score', 'pvals']
    
    
    hubbard_pvals = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/Hubbard_pvals.tab', sep='\t', index_col=0)
    hubbard_fdr = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/Hubbard_fdr.tab', sep='\t', index_col=0)
    hubbard_psi = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/Hubbard_psi.tab', sep='\t', index_col=0)

    hubbard_ds = get_averages_bulk(hubbard_psi)
    
    song_bulk_pvals = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/song_bulk_pvals.tab', sep='\t', index_col=0)
    song_bulk_fdr = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/song_bulk_fdr.tab', sep='\t', index_col=0)
    song_bulk_psi = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/song_bulk_psi.tab', sep='\t', index_col=0)

    song_bulk_ds = get_averages_bulk(song_bulk_psi)
    
    weyn_pvals = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/Weyn_pvals.tab', sep='\t', index_col=0)
    weyn_fdr = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/Weyn_fdr.tab', sep='\t', index_col=0)
    weyn_psi = pd.read_csv('~/sc_splicing_regulation/timeseries_neurogenesis/Weyn_psi.tab', sep='\t', index_col=0)

    weyn_ds = get_averages_bulk(weyn_psi)
    
    
    
    x = tiklova_PSI.index[tiklova_PSI[tiklova_rd.index].isna().mean(axis=1) <= 0.75] & tiklova_psix.index
    y = x & weyn_fdr.index
    weyn_rmats_fdr = -pd.DataFrame(weyn_fdr.loc[y].min(axis=1).sort_values())
    weyn_rmats_fdr.columns = ['fdr']
    
    x = chen_PSI.index[chen_PSI[chen_rd.index].isna().mean(axis=1) <= 0.75] & chen_psix.index
    y = x & hubbard_fdr.index
    hubbard_rmats_fdr = -pd.DataFrame(hubbard_fdr.loc[y].min(axis=1).sort_values())
    hubbard_rmats_fdr.columns = ['fdr']
    
    x = song_PSI.index[song_PSI[song_rd.index].isna().mean(axis=1) <= 0.75] & song_psix.index
    y = x & song_bulk_fdr.index
    song_bulk_rmats_fdr = -pd.DataFrame(song_bulk_fdr.loc[y].min(axis=1).sort_values())
    song_bulk_rmats_fdr.columns = ['fdr']
    
    W_tiklova = get_distance_matrix(tiklova_rd[['PC_1', 'PC_2']], k=100)
    W_chen = get_distance_matrix(chen_rd[['PC_1', 'PC_2']], k=30)
    W_song = get_distance_matrix(song_rd[['PC_1', 'PC_2']], k=30)
    
    pseudotime = pd.read_csv('~/data_sc_regulation/tiklova/pseudotime.tab', sep='\t', index_col=0)
    ordered_cells = pseudotime.loc[tiklova_rd.index].lineage_1_pseudotime.dropna().sort_values().index