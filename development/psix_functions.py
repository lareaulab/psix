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


def L_observation(psi_o, psi_a, psi_null, r, c, min_probability, psi_var, capture_var, times, version='current'):
    
    if version == 'current':
        L_a = np.max([min_probability, psi_o_prior(psi_o, psi_a, r, c, psi_var, capture_var, times)])
        L_null = np.max([min_probability, psi_o_prior(psi_o, psi_null, r, c, psi_var, capture_var, times)])
        
    elif version == 'psix_v1':
        L_a = np.max([min_probability, psix_v1_continuous(psi_o, psi_a, c, r)])
        L_null = np.max([min_probability, psix_v1_continuous(psi_o, psi_null, c, r)])
        
    elif version == 'psix_v1_discrete':
        L_a = np.max([min_probability, psix_v1_discrete(psi_o, psi_a, c, r)])
        L_null = np.max([min_probability, psix_v1_discrete(psi_o, psi_null, c, r)])
        
    elif version == 'psix_v2':
        
        L_a = np.max([min_probability, psix_v2_continuous(psi_o, psi_a, c, r, zvar = 1, capture_var = 0, times=times)])
        L_null = np.max([min_probability, psix_v2_continuous(psi_o, psi_null, c, r, zvar = 1, capture_var = 0, times=times)])
        
    elif version == 'psix_v2_discrete':
        
        L_a = np.max([min_probability, psix_v2_discrete(psi_o, psi_a, c, r, psi_var=1, times=times)])
        L_null = np.max([min_probability, psix_v2_discrete(psi_o, psi_null, c, r, psi_var=1, times=times)])
        
        
    
    L_a = np.log10(L_a)
    L_null = np.log10(L_null)
        
    if (np.isnan(L_a) or np.isnan(L_null) or np.isnan(L_a - L_null)):
        L = 0
    else:
        L = L_a - L_null
    return L


def L_statistic_vec(psi_o_array, psi_a_array, psi_null, mrna_array, c, min_probability, psi_var, capture_var, times, version='current'):
    L = []
    psi_a_array = np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in psi_a_array])

    func = lambda i: L_observation(psi_o_array[i], psi_a_array[i], psi_null, mrna_array[i], c, min_probability, psi_var, capture_var, times, version=version)
    x = range(len(psi_o_array))
    x_func = np.vectorize(func)
    L = x_func(x)

    return L
    
    
def calculate_exon_L(PSI_tab, W, mrna_counts, exon, k = 0, c = 0.1, weight_distance=True, randomize = False, seed=0, min_probability = 0.01, psi_var = 1, capture_var = 0.02, times=100, weight_observations=False, version='current'):
    

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

    L_vec = L_statistic_vec(psi_o_array, psi_a_array, psi_null, mrna_array, c, min_probability, psi_var, capture_var, times, version=version)

    if weight_observations:
#         if exon == 'Clta_4':
#             print('weight confirmation')
        weights = (W.loc[cell_list, cell_list] > 0).sum(axis=1)
        L = np.average(L_vec, weights=weights)
        
    else:
        L = np.sum(L_vec)/total_cells
    
    ######
    return L

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




##################### Different versions of Psix ########################

##### Psix 1, continuous. This should be the same as the results from Psix regular runs.
##### This is deterministic, meaning that it should be exactly the same as regular Psix runs.

def psix_v1_continuous(psi_o, psi, c, r):
    m = r/c
    if (m*psi < r*psi_o) or (m*(1-psi) < r*(1-psi_o)):
        proba = probability_psi_observation(psi_o, psi, c, r)
    else:
        proba = prob_psi_m_known(psi_o, psi, r, m)
    return proba


def probability_psi_observation(psi_o, psi, c, r, sum_times=10):
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
    m_array = np.arange(r, np.int(sum_times*r/c))
    
    comb_1 = comb(m_array*psi, r*psi_o)
    comb_2 = comb(m_array*(1-psi), r*(1-psi_o))
    proba_1 = c**(r+1)
    proba_2 = (1-c)**(m_array-r)
    
    prob_array = comb_1*comb_2*proba_1*proba_2
    
    return np.sum(prob_array)

    
def prob_psi_m_known(psi_o, psi, r, m):
    comb_1 = comb(m*psi, r*psi_o)
    comb_2 = comb(m*(1-psi), r*(1-psi_o))
    comb_3 = comb(m, r)**(-1)
    return comb_1*comb_2*comb_3
    
    
##### Psix 1, discrete. This is deterministic, meaning that it should be always exactly the same results.
##### In this case, we round the observed mrna molecules and mrna molecules with the exon, in order to
##### make discrete estimates of the probabilities.

def psix_v1_discrete(psi_o, psi, c, r):
    try:
        m = np.int(r/c)
        m_exon = np.int(round(psi*m))

        if m == 0:
            return np.nan

        psi_adj = m_exon/m

        r_cap = np.int(round(r))
        r_exon = np.int(round(psi_o*r))
        psi_o_adj = r_exon/r_cap
        if (m_exon < r_exon) or ((m-m_exon) < (r_cap-r_exon)):
            proba = probability_psi_observation(psi_o_adj, psi_adj, c, r_cap)
        else:
            proba = prob_psi_m_known_discrete(m, m_exon, r_cap, r_exon)
        return proba
    except:
        return np.nan
    
def prob_psi_m_known_discrete(m, m_exon, r_cap, r_exon):
    comb_1 = comb(m_exon, r_exon)
    comb_2 = comb((m-m_exon), (r_cap-r_exon))
    comb_3 = comb(m, r_cap)**(-1)
    return comb_1*comb_2*comb_3


##### Psix 2, continuous
##### Here we are adding a probabilistic prior to Psi

def psix_v2_continuous(psi_o, psi, c, r, zvar = 1, capture_var = 0, times=1000):
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
    
    c_array = np.random.normal(c, capture_var, times)
    c_array = np.array([x if x > 0.01 else 0.01 for x in c_array])

    m_array = np.int(r/c)

    psi_array = expit(np.random.normal(logit(psi), zvar, times))
    
    comb_1 = comb(m_array*psi_array, r*psi_o)
    comb_2 = comb(m_array*(1-psi_array), r*(1-psi_o))
    comb_3 = comb(m_array, r)**(-1)

    
    prob_array = comb_1*comb_2*comb_3
    
    return np.mean(prob_array)


##### Psix 2, discrete. 
##### In this case, we round the observed mrna molecules and mrna molecules with the exon, in order to
##### make discrete estimates of the probabilities.

def psix_v2_discrete(psi_o, psi, c, r, psi_var=1, times=1000):
    try:
        m = np.int(r/c)
        if m == 0:
            return np.nan

        psi = expit(norm.rvs(logit(psi), psi_var, size=times))
        proba_list = []
        m_exon = np.array([np.int(round(x*m)) for x in psi])
        psi_adj = m_exon/m
        r_cap = np.int(round(r))
        r_exon = np.int(round(psi_o*r))
        psi_o_adj = r_exon/r_cap
        proba = prob_psi_psix2(m, m_exon, r_cap, r_exon)
        return proba.mean()
    except:
        return np.nan


def prob_psi_psix2(m, m_exon, r_cap, r_exon):
    comb_1 = comb(m_exon, r_exon)
    comb_2 = comb((m-m_exon), (r_cap-r_exon))
    comb_3 = comb(m, r_cap)**(-1)
    return comb_1*comb_2*comb_3
