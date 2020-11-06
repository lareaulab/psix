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
from scipy.stats import norm


# from scipy.stats import binom
# from scipy.stats import beta
from scipy.stats import poisson

### For psix; the difference from verum is that is uses a gaussian kernel for the averages


def psi_o_prior(psi_o, psi, r, c, psi_var = 1, capture_var = 0.02, times=100):
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
    
#     c_array = np.random.normal(c, capture_var, times)
    c_array = norm.rvs(c, capture_var, size=times)
    c_array = np.array([x if x > 0.01 else 0.01 for x in c_array])

    m_array = np.random.poisson(r/c_array, times)
#     m_array = poisson.rvs(r/c_array, loc=0, size=1, random_state=None)
    m_array = m_array[~(m_array <= r)]
    times = len(m_array)
    
    try:
        psi_array = expit(norm.rvs(logit(psi), psi_var, size=times))
    except:
        print(psi)
        print(psi_var)
        print(times)
    
    comb_1 = comb(m_array*psi_array, r*psi_o)
    comb_2 = comb(m_array*(1-psi_array), r*(1-psi_o))
    comb_3 = comb(m_array, r)**(-1)

    
    prob_array = comb_1*comb_2*comb_3
    
    return np.mean(prob_array)


def L_observation(psi_o, psi_a, psi_null, r, c, min_probability, psi_var, capture_var, times):
    
    L_a = np.max([min_probability, psi_o_prior(psi_o, psi_a, r, c, psi_var, capture_var, times)])
    L_null = np.max([min_probability, psi_o_prior(psi_o, psi_null, r, c, psi_var, capture_var, times)])
    
    L_a = np.log10(L_a)
    L_null = np.log10(L_null)
        
    if (np.isnan(L_a) or np.isnan(L_null) or np.isnan(L_a - L_null)):
        L = 0
    else:
        L = L_a - L_null
    return L


def L_statistic_vec(psi_o_array, psi_a_array, psi_null, mrna_array, c, min_probability, psi_var, capture_var, times):
    L = []
    psi_a_array = np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in psi_a_array])

    func = lambda i: L_observation(psi_o_array[i], psi_a_array[i], psi_null, mrna_array[i], c, min_probability, psi_var, capture_var, times)
    x = range(len(psi_o_array))
    x_func = np.vectorize(func)
    L = x_func(x)

    return L
    
    
def calculate_exon_L(PSI_tab, W, mrna_counts, exon, k = 0, c = 0.1, weight_distance=True, randomize = False, seed=0, min_probability = 0.01, psi_var = 1, capture_var = 0.02, times=100):
    

#     try:
    cell_list = PSI_tab.loc[exon].dropna().index

    if randomize:
        np.random.seed(seed)
        shuffled_cells = shuffle(cell_list)

    else:
        shuffled_cells = cell_list

    psi_o_array = np.array(PSI_tab.loc[exon, shuffled_cells])

    mrna_array = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in mrna_counts.loc[exon,  shuffled_cells]])
    mrna_array = np.round(mrna_array).astype(int)

    total_cells = round((len(cell_list) - np.sum(mrna_array == 0)))

    if total_cells <= 0:
        return np.nan

    psi_a_array = np.array(
    np.array(pd.DataFrame(
        np.array(W.loc[cell_list, cell_list])*psi_o_array).sum(axis=1)
    )/np.array(W.loc[cell_list, cell_list].sum(axis=1)))

    psi_null = np.mean(psi_o_array)

    L_vec = L_statistic_vec(psi_o_array, psi_a_array, psi_null, mrna_array, c, min_probability, psi_var, capture_var, times)

    ######
    return np.sum(L_vec)/total_cells

#     except:
#         print('error')
#         return np.nan
    
    
def calculate_cross_L(PSI_tab, W, mrna_counts, exon1, exon2, k = 0, c = 0.1, weight_distance=True, randomize = False, seed=0, min_probability = 0.01, sum_times=10, min_obs=0.25, aproximate = 0):
    
    try:

        cell_list1 = PSI_tab.loc[exon1].dropna().index
        cell_list2 = PSI_tab.loc[exon2].dropna().index

        cell_list = cell_list1 & cell_list2

        if len(cell_list)/len(PSI_tab.columns) < min_obs:
            return np.nan, np.nan


        psi_o_array_1 = np.array(PSI_tab.loc[exon1, cell_list])
        psi_o_array_2 = np.array(PSI_tab.loc[exon2, cell_list])

        mrna_array_1 = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in mrna_counts.loc[exon1,  cell_list]])
        mrna_array_1 = np.round(mrna_array_1).astype(int)

        mrna_array_2 = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in mrna_counts.loc[exon2,  cell_list]])
        mrna_array_2 = np.round(mrna_array_2).astype(int)

        total_cells_1 = round((len(cell_list) - np.sum(mrna_array_1 == 0)))
        total_cells_2 = round((len(cell_list) - np.sum(mrna_array_2 == 0)))

        if (total_cells_1 <= 0) or (total_cells_2 <= 0):
            return np.nan, np.nan

        psi_a_array_1 = np.array(
        np.array(pd.DataFrame(
            np.array(W.loc[cell_list, cell_list])*psi_o_array_1).sum(axis=1)
        )/np.array(W.loc[cell_list, cell_list].sum(axis=1)))

        psi_null_1 = np.mean(psi_o_array_1)


        psi_a_array_2 = np.array(
        np.array(pd.DataFrame(
            np.array(W.loc[cell_list, cell_list])*psi_o_array_2).sum(axis=1)
        )/np.array(W.loc[cell_list, cell_list].sum(axis=1)))

        psi_null_2 = np.mean(psi_o_array_2)

        L_vec_1 = L_statistic_vec(psi_o_array_1, psi_a_array_2, psi_null_2, mrna_array_1, c, min_probability, sum_times, aproximate)
        L_vec_2 = L_statistic_vec(psi_o_array_2, psi_a_array_1, psi_null_1, mrna_array_2, c, min_probability, sum_times, aproximate)

        ######
        return np.sum(L_vec_1)/total_cells_1, np.sum(L_vec_2)/total_cells_2

    except:
#         print('error')
        return np.nan, np.nan
    
    
def get_distance_matrix(pca, k=100):
    print('winds of change')
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

