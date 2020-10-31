from tqdm import tqdm 
import numpy as np
from scipy.stats import friedmanchisquare
from scipy.stats import zscore
from scipy.special import logit
from scipy.stats import kruskal
import pandas as pd
from statsmodels.stats.multitest import multipletests
from sklearn.cluster import AgglomerativeClustering



def get_int_events(psi_table, m = 0.05, n=0.25):
    exons = psi_table.index[np.abs(0.5 - psi_table.mean(axis=1)) <= (0.5-m)]
    exons = exons & psi_table.index[psi_table.isna().mean(axis=1) <= (1-n)]
    return exons
    
    
def apply_kruskal_wallis(psi_table, mrna_per_event, pca_clusters, clusters = 5,
                     psi_min = 0.05, obs_min = 0.25, mrna_min = 10, filter_obs = False, 
                    linearize=False):
    '''
    Wrapper function that manages the run of the Kruskal-Wallis test in the dataset, in addition to
    running basic filtering and exon selection functions. At the moment of writing this note, many parts 
    of the code are vestigial and will be removed.
    
    Input:
      psi_table: Matrix of PSI
      mrna_counts: matrix of mRNA molecules per gene
      mrna_per_event: mrna_counts with psi_table index; extended for multiple exons per gene
      read_counts: SJ counts used to estimate observations in psi_table
      coverage_tab: splice junction coverage rate
      pca_clust: metadata matrix with cluster information for cells
      clusters: column in pca_clust with clusters (default AC, but cell_type can also be used)
      psi_min: consider only exons with PSI average between [psi_min, 1-psi_min]
      obs_min: minimum % of cells in cluster that have an observation to run the test 
               (default: 50%; for three clusters min)
      mrna_min: minimum number of mRNAs in a quality observation (default 10)
      mrna_read_min: set an additional baseline minimum of reads for the mRNA filter (default 0)
      read_min: flat minimum of informative reads per event for the read filter (default 10)
      filter_obs: vestigial; should not be run as True
      dset_name: vestigial; this function used to make plots as well
      correct_multitest: vestigial
      linearize: whether if linearize the PSI before running KW test. We've found that it has
                 very little effect in the results.
    Output
      change_tab: table with KW p-values for exons with minimum observations
      mrna_selected: set of exons that pass the mRNA filter
      mrna_only_selected: vestigial; set of exons that pass the 10 mRNA in gene minimum, but not
                          necessarily the minimum reads expected for 10 mRNAs
      read_selected: set of exons that pass the flat-read minimum filter
    '''

#     int_exons = get_int_events(psi_table, psi_min, obs_min)
    
#     good_observations_table = (mrna_per_event.loc[int_exons].fillna(0) >= mrna_min)
    
#     selected_exons = int_exons[good_observations_table.mean(axis=1) >= obs_min]
    
#     psi_table = psi_table.loc[selected_exons]
    
    if filter_obs:
        psi_table = psi_table.mask(~good_observations_table)
        
    kw_tab = cluster_anova_test(psi_table, pca_clusters, obs_min=obs_min, linearize=linearize)
    
    
    
    int_exons = get_int_events(psi_table, psi_min, obs_min)
    
    good_observations_table = (mrna_per_event.loc[int_exons].fillna(0) >= mrna_min)
    
    selected_exons = int_exons[good_observations_table.mean(axis=1) >= obs_min]
    
    for exon in kw_tab.index:
        if exon not in selected_exons:
            kw_tab.loc[exon, 'KW_score'] = 0
            kw_tab.loc[exon, 'pvals'] = 1
    
    return kw_tab


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


def elife_test(psi_table, mrna_table, clusters, obs_min, mrna_min = 10, mask_diff = False):
    kw = apply_kruskal_wallis(psi_table, mrna_table, clusters, 
                          obs_min = obs_min, mrna_min=mrna_min, linearize=False)
    
    kw['KW_score'] = kw.KW_score.fillna(0)
    kw['pvals'] = kw.pvals.fillna(1)
    
    return kw


def filtered_kruskal_wallis(psi_table, mrna_per_event, rd, clusters = 5,
                     psi_min = 0.05, obs_min = 0.25, mrna_min = 10,
                    linearize=False):
    '''
    Wrapper function that manages the run of the Kruskal-Wallis test in the dataset, in addition to
    running basic filtering and exon selection functions. At the moment of writing this note, many parts 
    of the code are vestigial and will be removed.
    
    Input:
      psi_table: Matrix of PSI
      mrna_counts: matrix of mRNA molecules per gene
      mrna_per_event: mrna_counts with psi_table index; extended for multiple exons per gene
      read_counts: SJ counts used to estimate observations in psi_table
      coverage_tab: splice junction coverage rate
      pca_clust: metadata matrix with cluster information for cells
      clusters: column in pca_clust with clusters (default AC, but cell_type can also be used)
      psi_min: consider only exons with PSI average between [psi_min, 1-psi_min]
      obs_min: minimum % of cells in cluster that have an observation to run the test 
               (default: 50%; for three clusters min)
      mrna_min: minimum number of mRNAs in a quality observation (default 10)
      mrna_read_min: set an additional baseline minimum of reads for the mRNA filter (default 0)
      read_min: flat minimum of informative reads per event for the read filter (default 10)
      filter_obs: vestigial; should not be run as True
      dset_name: vestigial; this function used to make plots as well
      correct_multitest: vestigial
      linearize: whether if linearize the PSI before running KW test. We've found that it has
                 very little effect in the results.
    Output
      change_tab: table with KW p-values for exons with minimum observations
      mrna_selected: set of exons that pass the mRNA filter
      mrna_only_selected: vestigial; set of exons that pass the 10 mRNA in gene minimum, but not
                          necessarily the minimum reads expected for 10 mRNAs
      read_selected: set of exons that pass the flat-read minimum filter
    '''
    
    int_exons = get_int_events(psi_table, psi_min, obs_min)
    
    good_observations_table = (mrna_per_event.loc[int_exons].fillna(0) >= mrna_min)
    
    selected_exons = int_exons[good_observations_table.mean(axis=1) >= obs_min]
    
    
    psi_table = psi_table.loc[selected_exons]
    
    psi_table = psi_table.mask(~good_observations_table)
    
    kwstats = []
    pvals = []
    exon_pass = []
    not_pass = 0
    
    for exon in selected_exons:
        cells = psi_table.loc[exon].dropna().index
        
        rd_c = rd.loc[cells]
        
        ac = AgglomerativeClustering(n_clusters=clusters)
        ac_clusters = ac.fit_predict(rd_c)
        
        rd_c['clusters'] = ac_clusters
        
        
        n = len(rd_c.clusters.unique())
        
        anova_stat, anova_p = test_exon_bimodal_anova(psi_table[cells], exon, rd_c.clusters, obs_min=obs_min, linearize=linearize, n=n)
        
        if not np.isnan(anova_p):
            kwstats.append(anova_stat)
            pvals.append(anova_p)
            exon_pass.append(exon)
            
        else:
            not_pass += 1
    
    cluster_df = pd.DataFrame()
    cluster_df['KW_score'] = kwstats
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df
    
    
def cluster_anova_test(psi_table, pca_clusters, correction = 'fdr_bh', 
                          obs_min = 0.5, linearize=False):
    '''
    Runs the Kruskal-Wallis test for a PSI matrix, and a given set of clusters.
    It wraps the test_exon_bimodal_anova function for all exons.
    
    Input
      psi_table: Matrix of PSI
      pca_clust: metadata dataframe with cluster information
      clusters: column of pca_clust with cluster information
      correction: vestigial; used to include an option to correct p-values
      correct_multitest: vestigial. p-values are not corrected anymore, because the significance of
                         individual exons is not the focus of this test.
      obs_min: minimum % of cells in cluster that have an observation to run the test 
               (default: 50%; for three clusters min)
      linearize: if calculate logit of the PSI (default: False)
    Output
      cluster_df: dataframe with p-values for each exon that meets the observed minimum
    '''
    
    n = len(pca_clusters.unique())
    
    kwstats = []
    pvals = []
    exon_pass = []
    not_pass = 0
        
    for i in range(len(psi_table.index)):
        
        exon = psi_table.index[i]
        anova_stat, anova_p = test_exon_bimodal_anova(psi_table, exon, pca_clusters, obs_min=obs_min, linearize=linearize, n=n)

        kwstats.append(anova_stat)
        pvals.append(anova_p)
        exon_pass.append(exon)
    
    cluster_df = pd.DataFrame()
    cluster_df['KW_score'] = kwstats
    cluster_df['pvals'] = pvals
    cluster_df.index = exon_pass
    
    return cluster_df


def test_exon_bimodal_anova(PSI_tab, exon, pca_clusters, obs_min=0.5, linearize=False, n=3):
    '''
    Run Kruskal-Wallis test for one exon, to get significance in the differences 
    in PSI between multiple clusters. 
    
    Input
      PSI_tab: Matrix of PSI
      exon: name of the exon to test
      pca_clust: metadata dataframe with cluster information
      clusters: column of pca_clust with cluster information
      obs_min: minimum % of cells in cluster that have an observation to run the test (default: 50%)
      linearize: if calculate logit of the PSI (default: False)
    Output:
      anova_p: p-value of the KW test on the input exon
      10: vestigial. Downstream code expects two outputs
    '''
    
    obs_cells = PSI_tab.loc[exon].dropna().index # select cells that have a PSI observation
    cluster_psi = [] # list of lists of PSI observations per cluster
    
    for i in pca_clusters.unique(): # iterate through each cluster
        clust_cells = pca_clusters.index[pca_clusters == i] # select the cells in the cluster
        c_cells = [x for x in obs_cells if x in clust_cells] # cells in cluster with PSI observation of the target exon
        if len(c_cells)/len(clust_cells) >= obs_min: # Make sure that minimum % of cells in the cluster is met
            psi = list(PSI_tab.loc[exon, c_cells]) # list of PSI observations of the exon in cluster
            cluster_psi.append(psi)
    if len(cluster_psi) >= n: # run the test only if the exon is observed in at least n different clusters
        try:
            
            anova_stat, anova_p = run_anova(cluster_psi, linearize) # stat and p-value of the KW test
        except:
            anova_stat = np.nan
            anova_p = np.nan # if the test crashes for some reason
    else:
        anova_stat = np.nan
        anova_p = np.nan # not enough observations to run test
    return anova_stat, anova_p

def run_anova(samples, linearize = False):
    '''
    Runs the Kruskal-Wallis analysis of variance for data between 2 up to 6 clusters.
    This function is necessary because for some reason, the Python implementation of the
    Kruskal-Wallis test can take mutiple observations, but not as an array, thus I have 
    to resort to this function to test different number of clusters.
    
    Input:
      samples: list of list of PSI observations
    Output:
      returns Kruskal-Wallis test results over the list of lists
    '''
    
    if linearize:
        samples = [linearize_psi(x) for x in samples]
        samples = [x - np.mean(x) for x in samples]
        
    if len(samples) == 2:
        return kruskal(samples[0], samples[1])
    elif len(samples) == 3:
        return kruskal(samples[0], samples[1], samples[2])
    elif len(samples) == 4:
        return kruskal(samples[0], samples[1], samples[2], samples[3])
    elif len(samples) == 5:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4])
    elif len(samples) == 6:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5])
    elif len(samples) == 7:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6])
    elif len(samples) == 8:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7])
    elif len(samples) == 9:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8])
    elif len(samples) == 10:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8], samples[9])
    elif len(samples) == 11:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8], samples[9], samples[10])
    elif len(samples) == 12:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8], samples[9], samples[10], samples[11])
    elif len(samples) == 13:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8], samples[9], samples[10], samples[11], samples[12])
    elif len(samples) == 14:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8], samples[9], samples[10], samples[11], samples[12], samples[13])
    elif len(samples) == 15:
        return kruskal(samples[0], samples[1], samples[2], samples[3], samples[4], samples[5], samples[6], samples[7], samples[8], samples[9], samples[10], samples[11], samples[12], samples[13], samples[14])
    else:
        raise Exception('Too many clusters. You can just edit this function and add compatibility for more clusters.')


def get_subsampled_exons(exon_list, index, tipo='subsampled1'):

    exon_list_subsampled = [x + '_{tipo}'.format(tipo=tipo) for x in exon_list]
    exon_list_subsampled = [x for x in exon_list_subsampled if x in index]
    return exon_list_subsampled


def get_good_bad_exons(L_dir, L_test_dir, kw_dir, kw_test_dir, bulk_fdr, bulk_ds, tipo):
    good_dset_L = L_dir.index[(L_dir.pvals <= 0.05) & (L_dir.L_score > 0)]
    good_dset_L_subsampled = get_subsampled_exons(good_dset_L, L_test_dir.index, tipo='subsampled'+tipo)
    #
    good_dset_kw = kw_dir.index[(kw_dir.pvals <= 0.05)]
    good_dset_kw_subsampled = get_subsampled_exons(good_dset_kw, kw_test_dir.index, tipo='subsampled'+tipo)
    #
    good_bulk = bulk_fdr.index[(bulk_fdr.min(axis=1) <= 0.05)] & bulk_ds.index[(bulk_ds.max(axis=1) >= 0.2)]
    #
    good_bulk_subsampled = get_subsampled_exons(good_bulk, L_test_dir.index, tipo='subsampled'+tipo)
    good_bulk_subsampled = [x for x in good_bulk_subsampled if x in L_test_dir.index]
    #
    good_for_all_dset = good_dset_L & good_dset_kw & good_bulk
    good_for_all_dset_subsample = get_subsampled_exons(good_for_all_dset, 
                                                          L_test_dir.index & kw_test_dir.index, 
                                                          tipo='subsampled'+tipo)
    
    good_dict = {'L':good_dset_L,
                 'KW':good_dset_kw,
                 'bulk':good_bulk,
                 'L_sub':good_dset_L_subsampled,
                 'KW_sub':good_dset_kw_subsampled,
                 'bulk_sub':good_bulk_subsampled}
    
    ###
    bad_dset_L = [x for x in L_dir.index if x not in good_dset_L]
    bad_dset_L_subsampled = get_subsampled_exons(bad_dset_L, L_test_dir.index, tipo='subsampled'+tipo)
    #
    bad_dset_kw = [x for x in kw_dir.index if x not in good_dset_kw]
    bad_dset_kw_subsampled = get_subsampled_exons(bad_dset_kw, kw_test_dir.index, tipo='subsampled'+tipo)
    #
    bad_bulk =bulk_fdr.index[(bulk_fdr.min(axis=1) > 0.05)] & bulk_ds.index[(bulk_ds.max(axis=1) < 0.2)]
    bad_bulk_subsampled = get_subsampled_exons(bad_bulk, L_test_dir.index, tipo='subsampled'+tipo)
    bad_bulk_subsampled = [x for x in bad_bulk_subsampled if x in L_test_dir.index]
    #
    bad_for_all_dset = sorted(set(bad_dset_L) & set(bad_dset_kw) & set(bad_bulk))
    bad_for_all_dset_subsample = get_subsampled_exons(bad_for_all_dset, 
                                                         L_test_dir.index & kw_test_dir.index, tipo='subsampled'+tipo)


    bad_dict = {'L':bad_dset_L,
                 'KW':bad_dset_kw,
                 'bulk':bad_bulk,
                 'L_sub':bad_dset_L_subsampled,
                 'KW_sub':bad_dset_kw_subsampled,
                 'bulk_sub':bad_bulk_subsampled}
    
    random_dset_L = [x for x in L_test_dir.index if 'random'+tipo in x]
    random_dset_kw = [x for x in kw_test_dir.index if 'random'+tipo in x]


    dset_random = [x for x in random_dset_L if x in random_dset_kw]
    dset_negatives = dset_random + bad_for_all_dset_subsample
    dset_positives = good_for_all_dset_subsample
    
    return dset_positives, dset_negatives, good_dict, bad_dict


def get_elife_tests(rd, PSI, mrna_event, L_score, n_clusters = 5):
    
    ac = AgglomerativeClustering(n_clusters=n_clusters)
    ac_clusters = ac.fit_predict(rd[['PC_1', 'PC_2']])
    rd['clusters'] = ac_clusters

    kw = elife_test(PSI.loc[L_score.index], 
                            mrna_event.loc[L_score.index], rd.clusters, obs_min=0, mrna_min=0)
    kw_1 = elife_test(PSI.loc[L_score.index], 
                                  mrna_event.loc[L_score.index], rd.clusters, obs_min=0.1, mrna_min=0)
    kw_25 = elife_test(PSI.loc[L_score.index], 
                                  mrna_event.loc[L_score.index], rd.clusters, obs_min=0.25, mrna_min=0)
    elife = elife_test(PSI.loc[L_score.index], 
                               mrna_event.loc[L_score.index], rd.clusters, obs_min=0.5, mrna_min=10)

    return kw, kw_1, kw_25, elife





# from metrics.precision_recall_fscore_support(…)
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import balanced_accuracy_score

def psix_predictions(psix_dir):
    predictions = ((psix_dir.pvals <= 0.05) & (psix_dir.L_score > 0)).astype(int)
    return predictions


def kw_predictions(kw_dir):
    predictions = (kw_dir.pvals <= 0.05).astype(int)
    return predictions
    
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve

def plot_sensitivity_v_specificity(
    pred_list, pos_exons, neg_exons, ax, label, color='navy', marker='o', reps=100, pred_score=''
):
    
    true_pos = sorted(set(pred_list.index) & set(pos_exons))
    true_neg = sorted(set(pred_list.index) & set(neg_exons))
    
    
    
    x_pos = int(len(true_pos)*0.1)
    x_neg = int(len(true_neg)*0.1)

    spe_list = []
    sen_list = []
    acc_list = []
    
    for i in tqdm(range(reps), leave=True, position=0):
        
        pos_select = np.random.choice(true_pos, x_pos, replace=False)
        neg_select = np.random.choice(true_neg, x_neg, replace=False)
        exon_select = list(pos_select) + list(neg_select)
        exon_y = [1]*len(pos_select) + [0]*len(neg_select)
        
        predictions_select = np.array(pred_list.loc[exon_select])
        
        sen_list.append(recall_score(exon_y, predictions_select))
        exon_y_inv = np.abs(1-np.array(exon_y))
        predictions_select_inv = np.abs(1-predictions_select)
        spe_list.append(recall_score(exon_y_inv, predictions_select_inv))
        
        acc_list.append(balanced_accuracy_score(exon_y, predictions_select))
        
        
    true_exons = true_pos + true_neg
    y = [1]*len(true_pos) + [0]*len(true_neg)
    predictions = np.array(pred_list.loc[true_exons])
    
    y_inv = np.abs(1-np.array(y))
    predictions_inv = np.abs(1-predictions)    
        
    specificity = recall_score(y_inv, predictions_inv)
    sensitivity = recall_score(y, predictions)
    balanced_accuracy = balanced_accuracy_score(y, predictions)
    
    if len(pred_score) != 0:
        y1 = list(pred_score.loc[true_pos]) + list(pred_score.loc[true_neg])
        
        fpr, tpr, thres = roc_curve(y, y1)

        ax.plot(tpr, 1-fpr, c=color)
    
    
    ax.plot([sensitivity-np.nanstd(sen_list), sensitivity+np.nanstd(sen_list)], [specificity, specificity], c=color)
    ax.plot([sensitivity, sensitivity], [specificity-np.nanstd(spe_list), specificity+np.nanstd(spe_list)], c=color)
    ax.scatter([sensitivity], [specificity], c=color, marker=marker, label=label)
  
    return balanced_accuracy, acc_list, specificity, spe_list, sensitivity, sen_list





def plot_sensitivity_v_precision(pred_list, pos_exons, neg_exons, ax, label, color='navy', marker='o', reps=100, pred_score=''):
    
    true_pos = sorted(set(pred_list.index) & set(pos_exons))
    true_neg = sorted(set(pred_list.index) & set(neg_exons))
    
    x_pos = int(len(true_pos)*0.1)
    x_neg = int(len(true_neg)*0.1)

    pre_list = []
    sen_list = []
    f1_list = []
    
    for i in tqdm(range(reps), leave=True, position=0):
        
        pos_select = np.random.choice(true_pos, x_pos, replace=False)
        neg_select = np.random.choice(true_neg, x_neg, replace=False)
        exon_select = list(pos_select) + list(neg_select)
        exon_y = [1]*len(pos_select) + [0]*len(neg_select)
        
        predictions_select = np.array(pred_list.loc[exon_select])
        
        sen_list.append(recall_score(exon_y, predictions_select))
        
        pre_list.append(precision_score(exon_y, predictions_select))
        
        f1_list.append(f1_score(exon_y, predictions_select))
        
        
    true_exons = true_pos + true_neg
    y = [1]*len(true_pos) + [0]*len(true_neg)
    predictions = np.array(pred_list.loc[true_exons])
    
        
    precision = precision_score(y, predictions)
    sensitivity = recall_score(y, predictions)
    f1 = f1_score(y, predictions)
    
    if len(pred_score) != 0:
        y1 = list(pred_score.loc[true_pos]) + list(pred_score.loc[true_neg])
        pre, rec, thres = precision_recall_curve(y, y1)
        ax.plot(rec, pre, c=color)
    
    ax.plot([sensitivity-np.nanstd(sen_list), sensitivity+np.nanstd(sen_list)], [precision, precision], c=color)
    ax.plot([sensitivity, sensitivity], [precision-np.nanstd(pre_list), precision+np.nanstd(pre_list)], c=color)
    ax.scatter([sensitivity], [precision], c=color, marker=marker, label=label)
  
    return f1, f1_list, precision, pre_list, sensitivity, sen_list



# def plot_sensitivity_v_specificity(L_dir, kw_dir, elife_dir, pos_exons, neg_exons, ax, label,
#                                    n_clusters=5, color='navy', obs_min=0.5):
    
#     true_pos = sorted(set(L_dir.index) & set(pos_exons))
# #     print(len(true_pos))
#     true_neg = sorted(set(L_dir.index) & set(neg_exons))

#     x_pos = int(len(true_pos)*0.2)
#     x_neg = int(len(true_neg)*0.2)

#     spe_list_L = []
#     sen_list_L = []
    
#     spe_list_kw = []
#     sen_list_kw = []
    
#     spe_list_elife = []
#     sen_list_elife = []
    
#     spe_list_elife_kw = []
#     sen_list_elife_kw = []
    
#     kw_true_pos = [x for x in true_pos if x in kw_dir.index]
#     kw_true_neg = [x for x in true_neg if x in kw_dir.index]
    
#     elife_true_pos = [x for x in true_pos if x in elife_dir.index]
#     elife_true_neg = [x for x in true_neg if x in elife_dir.index]
    
    
    
#     for i in tqdm(range(100), leave=True, position=0):
#         pos_select = np.random.choice(true_pos, x_pos, replace=False)
#         neg_select = np.random.choice(true_neg, x_neg, replace=False)

#         spe_list_L.append(np.mean((L_dir.loc[neg_select].pvals >= 0.05) | (L_dir.loc[neg_select].L_score <= 0)))
#         sen_list_L.append(np.mean((L_dir.loc[pos_select].pvals <= 0.05) & (L_dir.loc[pos_select].L_score > 0)))
        
#         kw_pos_select = [x for x in pos_select if x in kw_true_pos]
#         kw_neg_select = [x for x in neg_select if x in kw_true_neg]
        
#         spe_list_kw.append(np.mean((kw_dir.loc[kw_neg_select].pvals >= 0.05)))
#         sen_list_kw.append(np.mean((kw_dir.loc[kw_pos_select].pvals <= 0.05)))
        
#         elife_pos_select = [x for x in pos_select if x in elife_true_pos]
#         elife_neg_select = [x for x in neg_select if x in elife_true_neg]
        
#         if len(elife_pos_select) == 0:
#             sen_list_elife.append(0)
#             sen_list_elife_kw.append(0)
#         else:
#             sen_list_elife.append(len(elife_pos_select)/len(pos_select))
#             sen_list_elife_kw.append(np.sum(elife_dir.loc[elife_pos_select].pvals <= 0.05)/len(pos_select))
            
#         if len(elife_neg_select) == 0:
#             spe_list_elife.append(1)
#             spe_list_elife_kw.append(1)
            
#         else:
#             spe_list_elife.append(1-len(elife_neg_select)/len(neg_select))
#             spe_list_elife_kw.append(1-np.sum(elife_dir.loc[elife_neg_select].pvals <= 0.05)/len(neg_select))
            
        
#     accuracy_list_L = [np.mean([spe_list_L[i], sen_list_L[i]]) for i in range(len(spe_list_L))]
#     accuracy_list_kw = [np.mean([spe_list_kw[i], sen_list_kw[i]]) for i in range(len(spe_list_kw))]
#     accuracy_list_elife = [np.mean([spe_list_elife[i], sen_list_elife[i]]) for i in range(len(spe_list_elife))]
#     accuracy_list_elife_kw = [np.mean([spe_list_elife_kw[i], sen_list_elife_kw[i]]) for i in range(len(spe_list_elife_kw))]
        
        
#     specificity_L = np.mean((L_dir.loc[true_neg].pvals >= 0.05) | (L_dir.loc[true_neg].L_score <= 0))
#     sensitivity_L = np.mean((L_dir.loc[true_pos].pvals <= 0.05) & (L_dir.loc[true_pos].L_score > 0))
    
#     ax.plot([sensitivity_L-np.nanstd(sen_list_L), sensitivity_L+np.nanstd(sen_list_L)], [specificity_L, specificity_L], c=color)
#     ax.plot([sensitivity_L, sensitivity_L], [specificity_L-np.nanstd(spe_list_L), specificity_L+np.nanstd(spe_list_L)], c=color)
#     ax.scatter([sensitivity_L], [specificity_L], c=color)
    
#     specificity_kw = np.sum((kw_dir.loc[kw_true_neg].pvals >= 0.05))/len(true_neg)
#     sensitivity_kw = np.sum((kw_dir.loc[kw_true_pos].pvals <= 0.05))/len(true_pos)
#     specificity_elife = 1-len(elife_true_neg)/len(true_neg)
#     sensitivity_elife = len(elife_true_pos)/len(true_pos)
    
#     specificity_elife_kw = 1-np.sum((elife_dir.loc[elife_true_neg].pvals >= 0.05))/len(true_neg)
#     sensitivity_elife_kw = np.sum((elife_dir.loc[elife_true_pos].pvals <= 0.05))/len(true_pos)
    
#     ax.plot([sensitivity_kw-np.nanstd(sen_list_kw), sensitivity_kw+np.nanstd(sen_list_kw)], [specificity_kw, specificity_kw], c=color, linestyle='--')
#     ax.plot([sensitivity_kw, sensitivity_kw], [specificity_kw-np.nanstd(spe_list_kw), specificity_kw+np.nanstd(spe_list_kw)], c=color, linestyle='--')
#     ax.scatter(sensitivity_kw, specificity_kw, c=color, marker='^')
    
#     ax.plot([sensitivity_elife-np.nanstd(sen_list_elife), sensitivity_elife+np.nanstd(sen_list_elife)], [specificity_elife, specificity_elife], c=color, linestyle=':')
#     ax.plot([sensitivity_elife, sensitivity_elife], [specificity_elife-np.nanstd(spe_list_elife), specificity_elife+np.nanstd(spe_list_elife)], c=color, linestyle=':')
#     ax.scatter(sensitivity_elife, specificity_elife, c=color, marker='*')
    
#     ax.plot([sensitivity_elife_kw-np.nanstd(sen_list_elife_kw), sensitivity_elife_kw+np.nanstd(sen_list_elife_kw)], [specificity_elife_kw, specificity_elife_kw], c=color, linestyle=':')
#     ax.plot([sensitivity_elife_kw, sensitivity_elife_kw], [specificity_elife_kw-np.nanstd(spe_list_elife_kw), specificity_elife_kw+np.nanstd(spe_list_elife_kw)], c=color, linestyle=':')
#     ax.scatter(sensitivity_elife_kw, specificity_elife_kw, c=color, marker='v')
    
# #     print(sensitivity_L)
# #     print(specificity_L)
  
#     accuracy_L = np.mean([sensitivity_L, specificity_L])
#     accuracy_kw = np.mean([sensitivity_kw, specificity_kw])
# #     print(accuracy_kw)
#     accuracy_elife = np.mean([sensitivity_elife, specificity_elife])
#     accuracy_elife_kw = np.mean([sensitivity_elife_kw, specificity_elife_kw])
    
#     return accuracy_L, accuracy_list_L, accuracy_kw, accuracy_list_kw, accuracy_elife, accuracy_list_elife, accuracy_elife_kw, accuracy_list_elife_kw 



# def plot_sensitivity_v_precision(L_dir, kw_dir, elife_dir, pos_exons, neg_exons, ax, label,
#                                    n_clusters=5, color='navy', obs_min=0.5):
    
#     true_pos = sorted(set(L_dir.index) & set(pos_exons))
#     true_neg = sorted(set(L_dir.index) & set(neg_exons))
    
#     L_dir = L_dir.loc[true_pos + true_neg]

# #     print(len(true_neg))
    
# #     true_pos = L_dir.index & reference_fdr.loc[(reference_fdr.min(axis=1) <= 0.05) & (reference_ds.max(axis=1) >= 0.2)].index
# #     true_neg = L_dir.index & reference_fdr.loc[(reference_fdr.min(axis=1) > 0.05) | (reference_ds.max(axis=1) < 0.2)].index

#     x_pos = int(len(true_pos)*0.1)
#     x_neg = int(len(true_neg)*0.1)

#     pre_list_L = []
#     sen_list_L = []
    
#     pre_list_kw = []
#     sen_list_kw = []
    
#     pre_list_elife = []
#     sen_list_elife = []
    
#     pre_list_elife_kw = []
#     sen_list_elife_kw = []
    
#     f1_list_L = []
#     f1_list_kw = []
#     f1_list_elife = []
#     f1_list_elife_kw = []
    
#     kw_true_pos = [x for x in true_pos if x in kw_dir.index]
#     kw_true_neg = [x for x in true_neg if x in kw_dir.index]
    
#     kw_dir = kw_dir.loc[kw_true_pos + kw_true_neg]
    
#     elife_true_pos = [x for x in true_pos if x in elife_dir.index]
#     elife_true_neg = [x for x in true_neg if x in elife_dir.index]
    
#     elife_dir = elife_dir.loc[elife_true_pos + elife_true_neg]
    
#     for i in tqdm(range(200), leave=True, position=0):
#         pos_select = np.random.choice(true_pos, x_pos, replace=False)
#         neg_select = np.random.choice(true_neg, x_neg, replace=False)

        
#         L_dir_sample = L_dir.loc[list(pos_select) + list(neg_select)]
#         pred_pos = L_dir_sample.loc[(L_dir_sample.L_score > 0) & (L_dir_sample.pvals < 0.05)].index
#         pre_L = np.mean([x in pos_select for x in pred_pos])

#         sen_L = np.mean((L_dir.loc[pos_select].pvals <= 0.05) & (L_dir.loc[pos_select].L_score > 0))
        
#         kw_pos_select = [x for x in pos_select if x in kw_true_pos]
#         kw_neg_select = [x for x in neg_select if x in kw_true_neg]
        
#         kw_dir_sample = kw_dir.loc[list(kw_pos_select) + list(kw_neg_select)]
#         pred_pos = kw_dir_sample.loc[(kw_dir_sample.pvals < 0.05)].index
#         pre_kw = np.mean([x in pos_select for x in pred_pos])

#         sen_kw = np.mean((kw_dir.loc[kw_pos_select].pvals <= 0.05))
        
#         elife_pos_select = [x for x in pos_select if x in elife_true_pos]
#         elife_neg_select = [x for x in neg_select if x in elife_true_neg]
        
#         if len(elife_pos_select) == 0:
#             sen_elife = 0
#             sen_elife_kw = 0 
#         else:
#             sen_elife = len(elife_pos_select)/len(pos_select)
#             sen_elife_kw = np.sum(elife_dir.loc[elife_pos_select].pvals <= 0.05)/len(pos_select)
            
        
#         elife_dir_sample = elife_dir.loc[list(elife_pos_select) + list(elife_neg_select)]
        
#         if len(elife_dir_sample.index) > 0:

#             pre_elife = np.mean([x in pos_select for x in elife_dir_sample.index])
            
#             pred_pos_elife_kw = elife_dir_sample.loc[(elife_dir_sample.pvals < 0.05)].index
#             pre_elife_kw = np.mean([x in pos_select for x in pred_pos_elife_kw])
            
#             if pre_elife+sen_elife > 0:
#                 f1_list_elife.append(2*(pre_elife*sen_elife)/(pre_elife+sen_elife))
#             if pre_elife_kw+sen_elife_kw > 0:
#                 f1_list_elife_kw.append(2*(pre_elife_kw*sen_elife_kw)/(pre_elife_kw+sen_elife_kw))
                
#             pre_list_elife.append(pre_elife)
#             pre_list_elife.append(pre_elife_kw)
            
            
#         sen_list_L.append(sen_L)
#         pre_list_L.append(pre_L)
#         sen_list_kw.append(sen_kw)
#         pre_list_kw.append(pre_kw)
#         sen_list_elife.append(sen_elife)
#         sen_list_elife.append(sen_elife_kw)
        
            
#         f1_list_L.append(2*(pre_L*sen_L)/(pre_L+sen_L))
#         f1_list_kw.append(2*(pre_kw*sen_kw)/(pre_kw+sen_kw))

        
#     pred_pos = L_dir.loc[(L_dir.L_score > 0) & (L_dir.pvals < 0.05)].index
#     precision_L = np.mean([x in true_pos for x in pred_pos])
#     sensitivity_L = np.mean((L_dir.loc[true_pos].pvals <= 0.05) & (L_dir.loc[true_pos].L_score > 0))
    
#     ax.plot([sensitivity_L-np.nanstd(sen_list_L), sensitivity_L+np.nanstd(sen_list_L)], [precision_L, precision_L], c=color)
#     ax.plot([sensitivity_L, sensitivity_L], [precision_L-np.nanstd(pre_list_L), precision_L+np.nanstd(pre_list_L)], c=color)
#     ax.scatter([sensitivity_L], [precision_L], c=color)
    
#     pred_true_pos = np.sum(kw_dir.loc[kw_true_pos].pvals <= 0.05)
#     pred_pos = np.sum(kw_dir.pvals <= 0.05)
#     precision_kw = pred_true_pos/pred_pos
    
#     sensitivity_kw = np.sum((kw_dir.loc[kw_true_pos].pvals <= 0.05))/len(true_pos)
    
#     pred_pos = elife_dir.index
#     precision_elife = np.mean([x in true_pos for x in pred_pos])
#     sensitivity_elife = len(elife_true_pos)/len(true_pos)
    
#     pred_pos = elife_dir.loc[elife_dir.pvals < 0.05].index
#     precision_elife_kw = np.mean([x in true_pos for x in pred_pos])
#     sensitivity_elife_kw = np.sum((elife_dir.loc[elife_true_pos].pvals <= 0.05))/len(true_pos)
    
#     ax.plot([sensitivity_kw-np.nanstd(sen_list_kw), sensitivity_kw+np.nanstd(sen_list_kw)], 
#             [precision_kw, precision_kw], c=color, linestyle='--')
#     ax.plot([sensitivity_kw, sensitivity_kw], [precision_kw-np.nanstd(pre_list_kw), 
#                                                precision_kw+np.nanstd(pre_list_kw)], c=color, linestyle='--')
#     ax.scatter(sensitivity_kw, precision_kw, c=color, marker='^')
    
#     ax.plot([sensitivity_elife-np.nanstd(sen_list_elife), sensitivity_elife+np.nanstd(sen_list_elife)], 
#             [precision_elife, precision_elife], c=color, linestyle=':')
#     ax.plot([sensitivity_elife, sensitivity_elife], [precision_elife-np.nanstd(pre_list_elife), 
#                                                      precision_elife+np.nanstd(pre_list_elife)], c=color, linestyle=':')
#     ax.scatter(sensitivity_elife, precision_elife, c=color, marker='*')
    
#     ax.plot([sensitivity_elife_kw-np.nanstd(sen_list_elife_kw), sensitivity_elife_kw+np.nanstd(sen_list_elife_kw)], 
#             [precision_elife_kw, precision_elife_kw], c=color, linestyle=':')
#     ax.plot([sensitivity_elife_kw, sensitivity_elife_kw], [precision_elife_kw-np.nanstd(pre_list_elife_kw),
#                                                            precision_elife_kw+np.nanstd(pre_list_elife_kw)], c=color, linestyle=':')
#     ax.scatter(sensitivity_elife_kw, precision_elife_kw, c=color, marker='v')
    


#     f1_L = 2*(precision_L*sensitivity_L)/(precision_L+sensitivity_L)
#     f1_kw = 2*(precision_kw*sensitivity_kw)/(precision_kw+sensitivity_kw)
#     f1_elife = 2*(precision_elife*sensitivity_elife)/(precision_elife+sensitivity_elife)
#     f1_elife_kw = 2*(precision_elife_kw*sensitivity_elife_kw)/(precision_elife_kw+sensitivity_elife_kw)
    
#     return f1_L, f1_list_L, f1_kw, f1_list_kw, f1_elife, f1_list_elife, f1_elife_kw, f1_list_elife_kw 

def plot_summary_statistic(ax, stat, stat_list, color='navy', marker='o', x=1, label=''):
    ax.plot([x, x], [stat-np.nanstd(stat_list), stat+np.nanstd(stat_list)], c=color, linestyle=':')
    ax.plot([x-0.05, x+0.05], [stat-np.nanstd(stat_list), stat-np.nanstd(stat_list)], c=color)
    ax.plot([x-0.05, x+0.05], [stat+np.nanstd(stat_list), stat+np.nanstd(stat_list)], c=color)
    if label != '':
        ax.scatter([x], [stat], c=color, s=100, facecolors='white', label=label)
    else:
        ax.scatter([x], [stat], c=color, s=100, facecolors='white')
    
    
    

def plot_summary(ax, factor_list, color, x=1):
    
    L, list_L, kw, list_kw, elife, list_elife, elife_kw, list_elife_kw = factor_list
    
    ax.scatter([x], [L], c=color, s=100)
    ax.plot([x, x], [L-np.nanstd(list_L), L+np.nanstd(list_L)], c=color)
    ax.plot([x-0.05, x+0.05], [L-np.nanstd(list_L), L-np.nanstd(list_L)], c=color)
    ax.plot([x-0.05, x+0.05], [L+np.nanstd(list_L), L+np.nanstd(list_L)], c=color)

    ax.scatter([x + 0.2], [kw], c=color, s=100, marker='^')
    ax.plot([x + 0.2, x + 0.2], [kw-np.nanstd(list_kw), kw+np.nanstd(list_kw)], c=color, linestyle='--')
    ax.plot([x + 0.2-0.05, x + 0.2+0.05], [kw-np.nanstd(list_kw), kw-np.nanstd(list_kw)], c=color, linestyle='--')
    ax.plot([x + 0.2-0.05, x + 0.2+0.05], [kw+np.nanstd(list_kw), kw+np.nanstd(list_kw)], c=color, linestyle='--')

    ax.scatter([x + 0.4], [elife], c=color, s=100, marker='v')
    ax.plot([x + 0.4, x + 0.4], [elife-np.nanstd(list_elife), elife+np.nanstd(list_elife)], c=color, linestyle=':')
    ax.plot([x + 0.4-0.05, x + 0.4+0.05], [elife-np.nanstd(list_elife), elife-np.nanstd(list_elife)], c=color, linestyle=':')
    ax.plot([x + 0.4-0.05, x + 0.4+0.05], [elife+np.nanstd(list_elife), elife+np.nanstd(list_elife)], c=color, linestyle=':')

    ax.scatter([x + 0.6], [elife_kw], c=color, s=100, marker='*')
    ax.plot([x + 0.6, x + 0.6], [elife_kw-np.nanstd(list_elife_kw), elife_kw+np.nanstd(list_elife_kw)], c=color, linestyle=':')
    ax.plot([x + 0.6-0.05, x + 0.6+0.05], [elife_kw-np.nanstd(list_elife_kw), elife_kw-np.nanstd(list_elife_kw)], c=color, linestyle=':')
    ax.plot([x + 0.6-0.05, x + 0.6+0.05], [elife_kw+np.nanstd(list_elife_kw), elife_kw+np.nanstd(list_elife_kw)], c=color, linestyle=':')

    
    