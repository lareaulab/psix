import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from scipy.stats import spearmanr, pearsonr

import numpy as np
from scipy.stats import friedmanchisquare
from scipy.stats import zscore
from scipy.special import logit
from scipy.stats import kruskal
from IPython.core.pylabtools import figsize

import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, '/mnt/lareaulab/cfbuenabadn/sc_splicing_regulation/utils/')
from utils_functions import *

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams['pdf.fonttype'] = 42


data_dir = '/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/'

tiklova_mrna_event = pd.read_csv(data_dir + 'tiklova/mrna_per_event.tab', sep='\t', index_col=0)
tiklova_rd = pd.read_csv(data_dir + 'tiklova/rd_pc2.tab', sep='\t', index_col=0).loc[tiklova_mrna_event.columns]
tiklova_PSI = pd.read_csv(data_dir + 'tiklova/skipped_exons_psi.tab', sep='\t', index_col=0)[tiklova_mrna_event.columns]
tiklova_L_score = pd.read_csv('../tiklova_mp/tiklova.scores.txt', sep='\t', index_col=0)
tiklova_mrna = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/mrna_counts.tab', sep='\t', index_col=0)[tiklova_rd.index]
tiklova_reads = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/skipped_exons_SJreads.tab', 
                            sep='\t', index_col=0)[tiklova_rd.index]
tiklova_SJ_reads = pd.read_csv('~/data_sc_regulation/tiklova_extended/SE_counts.tab.gz', sep='\t', index_col=0)



from scipy.special import expit
from scipy.special import logit
import numpy as np
from functools import partial
import pandas as pd
import time as t
from multiprocessing import Pool
from tqdm import tqdm
import seaborn as sns
from scipy.stats import zscore
from scipy.special import logit
from scipy.special import comb
import scipy.integrate as integrate
import scipy.special as special
from sklearn.neighbors import NearestNeighbors
from sklearn.utils import shuffle


def probability_psi_observation(psi_o, psi, c, r):
    '''
    Calculate P(observed_PSI | true_PSI, capture efficiency, captured molecules)
    
    Input:
      psi_o (float): observed PSI
      psi (float): underlying PSI
      c (float): capture efficiency
      r (int): caputured gene mRNA molecules
      
    Output:
      P(psi_o | psi, c, r)
    
    '''
    m_array = np.arange(r, np.int(10*r/c))
    
    comb_1 = comb(m_array*psi, r*psi_o)
    comb_2 = comb(m_array*(1-psi), r*(1-psi_o))
    proba_1 = c**(r+1)
    proba_2 = (1-c)**(m_array-r)
    
    prob_array = comb_1*comb_2*proba_1*proba_2
    
    return np.sum(prob_array)


def probability_psi_m_known (psi_o, psi, c, m):
    
    '''
    Calculate P(observed_PSI | true_PSI, capture efficiency, molecules in cell)
    
    Input:
      psi_o (float): observed PSI
      psi (float): underlying PSI
      c (float): capture efficiency
      r (int): gene mRNA molecules in cell (before capture)
      
    Output:
      P(psi_o | psi, c, m)
      
    '''

    r = np.arange(1, m+1)
    
    comb_1 = comb(m*psi, r*psi_o)
    comb_2 = comb(m*(1-psi), r*(1-psi_o))
    proba_1 = c**r
    proba_2 = (1-c)**(m-r)

    prob_array = comb_1*comb_2*proba_1*proba_2/(1-(1-c)**m)
        
    return np.sum(prob_array)


def probability_binary_m_known(psi, c, m):
    '''
    P(psi_o in {0, 1} | psi, c, m)
    
    Input:
      psi_o (float): observed PSI
      psi (float): underlying PSI
      c (float): capture efficiency
      r (int): gene mRNA molecules in cell (before capture)
      
    Output:
      P(psi_o | psi, c, m)
      
    '''
    return probability_psi_m_known(1, psi, c, m) + probability_psi_m_known(0, psi, c, m)


def integrate_probability(psi, min_p, max_p, c, r):
    '''
    Integrate P(psi_o|psi, c, r) with psi_o as a continuous variable.
    
    '''
    return integrate.quad(probability_psi_observation, min_p, max_p, args=(psi, c, r))[0]


def psi_observation_range_probability(psi, min_p, max_p, c, r):
    '''
    r is a discrete variable. For this reason, psi_o has a discrete and limited number of
    possible values. However, realistically speaking r is a rough estimate, unless we model
    directly a correspondence between mRNA molecules and observed reads.
    
    Here we estimate the density of the probability of psi observations with psi_o as a 
    continuos value. P(psi_o|psi, c, r) is NOT a continuous probability distribution.
    
    '''
    total_proba = integrate_probability(psi, 0, 1, c, r)
    range_proba = integrate_probability(psi, min_p, max_p, c, r)
    return range_proba/total_proba


def observation_delta_probability(psi, delta, c, r):
    '''
    Probability that psi_o falls within a certain delta from the underlying psi.
    '''
    return psi_observation_range_probability(psi, psi-delta, psi+delta, c, r)


def L_observation(psi_o, psi_a, psi_null, r, c):
    
    L_a = np.max([0.01, probability_psi_observation(psi_o, psi_a, c, r)])
    L_null = np.max([0.01, probability_psi_observation(psi_o, psi_null, c, r)])
    
    L_a = np.log10(L_a)
    L_null = np.log10(L_null)
        
    if (np.isnan(L_a) or np.isnan(L_null)):
        L = 0
    else:
        L = L_a - L_null
    return L


def L_statistic_vec(psi_o_array, psi_a_array, psi_null, mrna_array, c):
    L = []
    psi_a_array = np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in psi_a_array])

    func = lambda i: L_observation(psi_o_array[i], psi_a_array[i], psi_null, mrna_array[i], c)
    x = range(len(psi_o_array))
    x_func = np.vectorize(func)
    L = x_func(x)

    return L


def calculate_cross_L(PSI_tab, distance_metric, mrna_counts, exon1, exon2, k = 0, c = 0.1, adjust_psi = False, randomize=False, seed=0):
    
    try:
    
        cell_list1 = PSI_tab.loc[exon1].dropna().index
        cell_list2 = PSI_tab.loc[exon2].dropna().index

        cell_list = cell_list1 & cell_list2

        data = pd.DataFrame()

        data['psi_o1'] = PSI_tab.loc[exon1,  cell_list]
        data['mrna1'] = mrna_counts.loc[exon1,  cell_list]

        data['psi_o2'] = PSI_tab.loc[exon2,  cell_list]
        data['mrna2'] = mrna_counts.loc[exon2,  cell_list]

        data.index =  cell_list

        distance_metric = distance_metric.loc[cell_list]
        
        if randomize:
            np.random.seed(seed)
            shuffled_cells = shuffle(cell_list)
            data = data.loc[shuffled_cells]
            data.index = cell_list

        psi_o_array1 = np.array(data.psi_o1)
        psi_o_array2 = np.array(data.psi_o2)
        
        mrna_array1 = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in data.mrna1])
        mrna_array1 = np.round(mrna_array1).astype(int)

        mrna_array2 = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in data.mrna2])
        mrna_array2 = np.round(mrna_array2).astype(int)

        total_cells = [1 if x > 0 else 0 for x in ((mrna_array1)*(mrna_array2))]

        psi_null1 = np.mean(psi_o_array1)
        psi_null2 = np.mean(psi_o_array2)

        if k == 0:
            k = int(round(np.sqrt(len(cell_list))))

        if np.sum(total_cells) <= np.sqrt(len(PSI_tab.index)):
                return np.nan

        nbrs = NearestNeighbors(n_neighbors=k+1).fit(distance_metric.loc[cell_list])
        distances, idxs = nbrs.kneighbors(distance_metric.loc[cell_list])

        func = lambda idx: np.average(data.psi_o1.loc[cell_list[idxs[idx][range(1, len(idxs[idx]))]]], 
                                      weights=1/(distances[idx][range(1, len(idxs[idx]))]))
        x = range(len(idxs))
        x_func = np.vectorize(func)
        psi_a_array1 = x_func(x)

        func = lambda idx: np.average(data.psi_o2.loc[cell_list[idxs[idx][range(1, len(idxs[idx]))]]], 
                                      weights=1/(distances[idx][range(1, len(idxs[idx]))]))
        x = range(len(idxs))
        x_func = np.vectorize(func)
        psi_a_array2 = x_func(x)

        if adjust_psi:
            psi_a_array1 = logit(np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in psi_a_array1]))
            psi_a_array2 = logit(np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in psi_a_array2]))
            
            psi_a_array1_mean = np.mean(psi_a_array1)
            psi_a_array2_mean = np.mean(psi_a_array2)
            
            psi_a_array1 = expit(psi_a_array1 - psi_a_array2_mean)
            
            psi_a_array2 = expit(psi_a_array2 - psi_a_array1_mean)
            

        L_vec1 = L_statistic_vec(psi_o_array1, psi_a_array2, psi_null2, mrna_array1, c)
        L_vec2 = L_statistic_vec(psi_o_array2, psi_a_array1, psi_null1, mrna_array2, c)

        return np.sum(L_vec1)/np.sum(total_cells), np.sum(L_vec2)/np.sum(total_cells)
    
    except:
        print('error')
        return np.nan, np.nan


def calculate_exon_L(PSI_tab, distance_metric, mrna_counts, exon, k = 0, c = 0.1, weight_distance=True, randomize = False, seed=0):
    
    try:

        cell_list = PSI_tab.loc[exon].dropna().index

        data = pd.DataFrame()

        data['psi_o'] = PSI_tab.loc[exon,  cell_list]
        data['mrna'] = mrna_counts.loc[exon,  cell_list]

        data.index =  cell_list

        distance_metric = distance_metric.loc[cell_list]

        if randomize:
#             print('random!')
            np.random.seed(seed)
            shuffled_cells = shuffle(cell_list)
            data = data.loc[shuffled_cells]
            data.index = cell_list

        psi_o_array = np.array(data.psi_o)

        mrna_array = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in data.mrna])
        mrna_array = np.round(mrna_array).astype(int)

        total_cells = round((len(cell_list) - np.sum(mrna_array == 0)))

        if total_cells <= 0:
            return np.nan

        psi_null = np.mean(psi_o_array)

        if k == 0:
            k = int(round(np.sqrt(len(cell_list))))

        nbrs = NearestNeighbors(n_neighbors=k+1).fit(distance_metric.loc[cell_list])
        distances, idxs = nbrs.kneighbors(distance_metric.loc[cell_list])

        if weight_distance:
            func = lambda idx: np.average(data.psi_o.loc[cell_list[idxs[idx][range(1, len(idxs[idx]))]]], 
                                          weights=1/(distances[idx][range(1, len(idxs[idx]))]))
        else:
            func = lambda idx: np.mean(data.psi_o.loc[cell_list[idxs[idx][range(1, len(idxs[idx]))]]])

        x = range(len(idxs))
        x_func = np.vectorize(func)
        psi_a_array = x_func(x)
        
        L_vec = L_statistic_vec(psi_o_array, psi_a_array, psi_null, mrna_array, c)
        
#         return L_vec, cell_list
        ######
        return np.sum(L_vec)/total_cells

    except:
        print('error')
        return np.nan

    
    
def downsample_exon(reads_table, mrna_table, exon, cells, downsample_rate = 0.5):
    exon_mrnas = round(mrna_table.loc[exon, cells].fillna(0)).astype(int)
    sub_mrnas = np.random.binomial(exon_mrnas, downsample_rate)
    ratios = [sub_mrnas[i]/exon_mrnas[i] if exon_mrnas[i]>0 else 0 for i in range(len(exon_mrnas))]
#     print(ratios)
    exon_SJ = [x for x in reads_table.index if exon + '_' in x]
    reads_df = pd.DataFrame(np.random.binomial(reads_table.loc[exon_SJ, cells], ratios))
    reads_df.columns = cells
    reads_df.index = exon_SJ
    
    i_sj = reads_df.loc[[exon + '_I1', exon + '_I2']].sum(axis=0)
    e_sj = reads_df.loc[[exon + '_SE']]
    
    psi = i_sj/(2*e_sj + i_sj)
    
    psi.index = [exon]
    
    sub_mrnas_df = pd.DataFrame()
    sub_mrnas_df[exon] = sub_mrnas
    sub_mrnas_df.index = cells
    
    return psi, sub_mrnas_df.T


def downsample_dataset(psi_table, reads_table, mrna_table, exon_list, cells, downsample_rate = 0.5):
    for exon in exon_list:
        psi, sub_mrnas_df = downsample_exon(reads_table, mrna_table, exon, cells, downsample_rate=downsample_rate)
        psi_table.loc[exon, cells] = psi.loc[exon, cells]
        mrna_table.loc[exon, cells] = sub_mrnas_df.loc[exon, cells]
        
    return psi_table, mrna_table, psi, sub_mrnas_df


def calculate_exon_pval(PSI_tab, distance_metric, mrna_counts, exon, bin_cuts, bin_randoms,
                        k = 0, c = 0.1, weight_distance=True, randomize = False, seed=0, randoms = 100):
    L = calculate_exon_L(PSI_tab, distance_metric, mrna_counts, exon, k = k, c = c, weight_distance=weight_distance, 
                         randomize = randomize, seed=seed)
        
    mean_cut = np.log10(mrna_counts.loc[exon]+1).mean()
    var_cut = PSI_tab.loc[exon].var()
    
    mean_bin = np.sum(bin_cuts.mean_cuts <= mean_cut) + 1
    var_bin = np.sum(bin_cuts.var_cuts <= var_cut) + 1
    
    bin_name = 'mean_{m}_var_{v}---'.format(m=str(mean_bin), v=str(var_bin))

    if bin_name in bin_randoms.columns:
        pval = (np.sum([L < bin_randoms[bin_name]])+1)/(len(bin_randoms[bin_name])+1)
        return L, pval
        
    else:
        L_r = []
        for i in tqdm(range(randoms), position=0, leave=True):
#             print(i)
            L_r.append(calculate_exon_L(PSI_tab, distance_metric, mrna_counts, exon, k = k, c = c, 
                                        weight_distance=weight_distance, 
                             randomize = True, seed=i))
            
        pval = (np.sum([L < np.array(L_r)])+1)/(len(L_r)+1)
        return L, pval#, L_r
        
    
    
    
def make_random(psi_table, mrna_table, read_table, exon):
#     print(exon)
    data = pd.DataFrame()
    shuffled_cells = shuffle(psi_table.columns)
    data['psi'] = list(psi_table.loc[exon, shuffled_cells])
    data['mrna'] = list(mrna_table.loc[exon, shuffled_cells])
    data['reads'] = list(read_table.loc[exon, shuffled_cells])
    data.index = psi_table.columns
    return data


def make_test(psi_table, mrna_table, read_table, exon_list, rd, seed = 0, randomize_exons = True, all_cells=False, 
              prob = 0.5):
    
    np.random.seed(seed)
    
    print(prob)
    
    if randomize_exons:
        shuffled_cells = shuffle(psi_table.columns)
    else:
        shuffled_cells = psi_table.columns
    
    random_psi = psi_table[shuffled_cells].copy()
    random_mrna = mrna_table[shuffled_cells].copy()
    random_read = read_table[shuffled_cells].copy()
    
    random_psi.columns = psi_table.columns
    random_mrna.columns = psi_table.columns
    random_read.columns = psi_table.columns
    
    all_cell_counts = 0
    
    if all_cells:
        print('all cells')
        random_psi, random_mrna, p, s = downsample_dataset(random_psi, 
                                                   random_read, 
                                             random_mrna, 
                                                   exon_list, psi_table.columns, 
                                             downsample_rate = prob)
        
    else:
        
        print('selecting cells')
        
        
    
    
        for exon in tqdm(exon_list, position=0, leave=True):

    #         prob = np.random.uniform(0.05, 0.5)




            lim_pc_1 = 0
            lim_pc_2 = 0
            coin_2_pc1 = 2
            coin_2_pc2 = 2
            cells_1 = []
            cells_2 = []

            coin_toss_1 = np.random.choice(range(4))
            if (coin_toss_1 == 0) or (coin_toss_1 == 2):
                lim_pc_1 = np.random.uniform(rd.PC_1.quantile(0.1), rd.PC_1.quantile(0.9))
                coin_2_pc1 = np.random.binomial(1, 0.5)

                if coin_2_pc1 == 0:
                    cells_1 = rd.index[rd.PC_1 < lim_pc_1]
                elif coin_2_pc1 == 1:
                    cells_1 = rd.index[rd.PC_1 > lim_pc_1]

            if (coin_toss_1 == 1) or (coin_toss_1 == 2):
                lim_pc_2 = np.random.uniform(rd.PC_2.quantile(0.1), rd.PC_2.quantile(0.9))
                coin_2_pc2 = np.random.binomial(1, 0.5)

                if coin_2_pc2 == 0:
                    cells_2 = rd.index[rd.PC_2 < lim_pc_2]
                elif coin_2_pc2 == 1:
                    cells_2 = rd.index[rd.PC_2 > lim_pc_2]

            if coin_toss_1 == 2:
                cells = [x for x in cells_1 if x in cells_2]

            else:
                cells = list(cells_1) + list(cells_2)

            if len(cells) > 100:    

                random_psi, random_mrna, p, s = downsample_dataset(random_psi, 
                                                       random_read, 
                                                 random_mrna, 
                                                       [exon], cells, 
                                                 downsample_rate = prob)
            else:
                random_psi, random_mrna, p, s = downsample_dataset(random_psi, 
                                                       random_read, 
                                                 random_mrna, 
                                                       [exon], random_psi.columns, 
                                                 downsample_rate = prob)
                
                all_cell_counts += 1
                

    print(all_cell_counts)    
    return random_psi, random_mrna
        
    
    
    
random_psi_1, random_mrna_1 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=0)

random_psi_2, random_mrna_2 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=1, prob=0.25)

random_psi_3, random_mrna_3 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=2, prob=0.1)

random_psi_4, random_mrna_4 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=3, prob=0.01)


subsampled_psi_1, subsampled_mrna_1 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], 
                                        seed=4, randomize_exons = False)

subsampled_psi_2, subsampled_mrna_2 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=5, 
                                        prob=0.25, randomize_exons = False)

subsampled_psi_3, subsampled_mrna_3 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=6, 
                                        prob=0.1, randomize_exons = False)

subsampled_psi_4, subsampled_mrna_4 = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
                                        tiklova_mrna_event.loc[tiklova_L_score.index], 
                                    tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=7, 
                                        prob=0.01, randomize_exons = False)



random_psi_1.index = [x + '_random1' for x in random_psi_1.index]
random_mrna_1.index = [x + '_random1' for x in random_mrna_1.index]

random_psi_2.index = [x + '_random2' for x in random_psi_2.index]
random_mrna_2.index = [x + '_random2' for x in random_mrna_2.index]

random_psi_3.index = [x + '_random3' for x in random_psi_3.index]
random_mrna_3.index = [x + '_random3' for x in random_mrna_3.index]

random_psi_4.index = [x + '_random4' for x in random_psi_4.index]
random_mrna_4.index = [x + '_random4' for x in random_mrna_4.index]


subsampled_psi_1.index = [x + '_subsampled1' for x in subsampled_psi_1.index]
subsampled_mrna_1.index = [x + '_subsampled1' for x in subsampled_mrna_1.index]

subsampled_psi_2.index = [x + '_subsampled2' for x in subsampled_psi_2.index]
subsampled_mrna_2.index = [x + '_subsampled2' for x in subsampled_mrna_2.index]

subsampled_psi_3.index = [x + '_subsampled3' for x in subsampled_psi_3.index]
subsampled_mrna_3.index = [x + '_subsampled3' for x in subsampled_mrna_3.index]

subsampled_psi_4.index = [x + '_subsampled4' for x in subsampled_psi_4.index]
subsampled_mrna_4.index = [x + '_subsampled4' for x in subsampled_mrna_4.index]



pd.concat([random_psi_1, random_psi_2, random_psi_3, random_psi_4,
          subsampled_psi_1, subsampled_psi_2, subsampled_psi_3, 
           subsampled_psi_4]).to_csv('tiklova_subsampled_exon_psi.tab', sep='\t', index=True, header=True)

pd.concat([random_mrna_1, random_mrna_2, random_mrna_3, random_mrna_4,
          subsampled_mrna_1, subsampled_mrna_2, subsampled_mrna_3, 
           subsampled_mrna_4]).to_csv('tiklova_subsampled_exon_mrna.tab', sep='\t', index=True, header=True)

