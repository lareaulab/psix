import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
from scipy.special import comb
from sklearn.utils import shuffle
import itertools  

from numba import jit, njit

import logging



@jit(nopython=True)
def round_numba(x):
    int_part = int(x)
    frac_part = x - int_part
    if frac_part >= 0.5:
        return int_part + 1
    else:
        return int_part

@jit(nopython=True)
def n_choose_k(n, k):
    if k == 0: return 1
    return (n * n_choose_k(n - 1, k - 1)) / k

@jit(nopython=True)
def probability_psi_observation(observed_psi, true_psi, capture_efficiency, captured_mrna):
    """
    Calculate Pr(observed_PSI | true_PSI, capture_efficiency, captured_mrna)
    
    Input:
      observed_psi (float): observed PSI
      psi (float): underlying PSI
      c (float): capture efficiency
      r (int): caputured gene mRNA molecules
      
    Output:
      Pr(observed_psi | psi, c, captured_mrna)
    """
    cell_molecules = captured_mrna/capture_efficiency
    if (cell_molecules*true_psi < captured_mrna*observed_psi) or (cell_molecules*(1-true_psi) < captured_mrna*(1-observed_psi)):
        proba = probability_psi_m_unknown(observed_psi, true_psi, capture_efficiency, captured_mrna)
    else:
        proba = probability_psi_m_known(observed_psi, true_psi, captured_mrna, cell_molecules)
        
    return proba

@jit(nopython=True)
def probability_psi_m_unknown(observed_psi, true_psi, capture_efficiency, captured_mrna):
    '''
    
    '''
    
    proba_total = 0

    max_iter = 10*captured_mrna/capture_efficiency
    
    for mrna in range(captured_mrna, max_iter):
        mrna_truePsi = round_numba(mrna*true_psi)
        capturedMRNA_obsPsi = round_numba(captured_mrna*observed_psi)

        comb_1 = n_choose_k(mrna_truePsi, capturedMRNA_obsPsi)
        
        mrna_complementPsi = round_numba(mrna*(1-true_psi))
        capturedMRNA_complementPsi = round_numba(captured_mrna*(1-observed_psi))
        
        comb_2 = n_choose_k(mrna_complementPsi, capturedMRNA_complementPsi)

        proba_1 = capture_efficiency**(captured_mrna+1)
        proba_2 = (1-capture_efficiency)**(mrna-captured_mrna)

        proba_total += comb_1*comb_2*proba_1*proba_2
    
    return proba_total

@jit(nopython=True)
def probability_psi_m_known(observed_psi, true_psi, captured_mrna, cell_molecules):
    """
    Pr(observed_PSI | true_PSI, captured_mrna, cell_molecules)
    
    """

    cellMol_truePsi = round_numba(cell_molecules*true_psi)
    capturedmRNA_obsPsi = round_numba(captured_mrna*observed_psi)
    comb_1 = n_choose_k(cellMol_truePsi, capturedmRNA_obsPsi)

    cellMol_inversePsi = round_numba(cell_molecules*(1-true_psi))
    capturedmRNA_inversePsi = round_numba(captured_mrna*(1-observed_psi))
    comb_2 = n_choose_k(cellMol_inversePsi, capturedmRNA_inversePsi)

    comb_3 = 1/n_choose_k(cell_molecules, captured_mrna)
    
    proba_out = comb_1*comb_2*comb_3
    return proba_out



# def probability_psi_observation(observed_psi, true_psi, capture_efficiency, captured_mrna):
#     """
#     Calculate Pr(observed_PSI | true_PSI, capture_efficiency, captured_mrna)
    
#     Input:
#       observed_psi (float): observed PSI
#       psi (float): underlying PSI
#       c (float): capture efficiency
#       r (int): caputured gene mRNA molecules
      
#     Output:
#       Pr(observed_psi | psi, c, captured_mrna)
#     """
#     cell_molecules = captured_mrna/capture_efficiency
#     if (cell_molecules*true_psi < captured_mrna*observed_psi) or (cell_molecules*(1-true_psi) < captured_mrna*(1-observed_psi)):
#         proba = probability_psi_m_unknown(observed_psi, true_psi, capture_efficiency, captured_mrna)
#     else:
#         proba = probability_psi_m_known(observed_psi, true_psi, captured_mrna, cell_molecules)
        
#     return proba

# def probability_psi_m_unknown(observed_psi, true_psi, capture_efficiency, captured_mrna):
#     '''
    
#     '''
#     m_array = np.arange(captured_mrna, int(10*captured_mrna/capture_efficiency))
    
#     comb_1 = comb(m_array*true_psi, captured_mrna*observed_psi)
#     comb_2 = comb(m_array*(1-true_psi), captured_mrna*(1-observed_psi))
#     proba_1 = capture_efficiency**(captured_mrna+1)
#     proba_2 = (1-capture_efficiency)**(m_array-captured_mrna)
    
#     prob_array = comb_1*comb_2*proba_1*proba_2
    
#     return np.sum(prob_array)
    
# def probability_psi_m_known(observed_psi, true_psi, captured_mrna, cell_molecules):
#     """
#     Pr(observed_PSI | true_PSI, captured_mrna, cell_molecules)
    
#     """
#     comb_1 = comb(cell_molecules*true_psi, captured_mrna*observed_psi)
#     comb_2 = comb(cell_molecules*(1-true_psi), captured_mrna*(1-observed_psi))
#     comb_3 = comb(cell_molecules, captured_mrna)**(-1)
#     return comb_1*comb_2*comb_3

# def psi_observation_score(
#     observed_psi, 
#     neighborhood_psi, 
#     global_psi, 
#     captured_mrna, 
#     capture_efficiency, 
#     min_probability
# ):
#     """
#     log Pr(observed_psi | neighborhood_psi) - log Pr(observed_psi | global_psi)
#     """
    
#     if np.isnan(observed_psi):
#         return 0
    
#     probability_given_neighborhood = np.max([
#         min_probability, probability_psi_observation(observed_psi, neighborhood_psi, capture_efficiency, captured_mrna)
#     ])
#     probability_given_global = np.max([
#         min_probability, probability_psi_observation(observed_psi, global_psi, capture_efficiency, captured_mrna)
#     ])
    
#     log_probability_neighborhood = np.log(probability_given_neighborhood)
#     log_probability_global = np.log(probability_given_global)
    
#     any_nan = False
#     if (np.isnan(log_probability_neighborhood)):
#         any_nan=True
#     if np.isnan(log_probability_global):
#         any_nan=True
#     if np.isnan(log_probability_neighborhood - log_probability_global):
#         any_nan=True
    
#     if any_nan:
#         return 0
    
#     else: 
#         return log_probability_neighborhood - log_probability_global

@jit(nopython=True)
def max_pair_numba(a, b):
    if a > b:
        return a
    else:
        return b

@jit(nopython=True)
def min_pair_numba(a, b):
    if a < b:
        return a
    else:
        return b
        
@jit(nopython=True)
def psi_observation_score(
    observed_psi, 
    neighborhood_psi, 
    global_psi, 
    captured_mrna, 
    capture_efficiency, 
    min_probability
):
    """
    log Pr(observed_psi | neighborhood_psi) - log Pr(observed_psi | global_psi)
    """
    
    if np.isnan(observed_psi):
        return 0
    
    probability_given_neighborhood = max_pair_numba(
        min_probability, probability_psi_observation(observed_psi, neighborhood_psi, capture_efficiency, captured_mrna)
    )
    probability_given_global = max_pair_numba(
        min_probability, probability_psi_observation(observed_psi, global_psi, capture_efficiency, captured_mrna)
    )
    
    log_probability_neighborhood = np.log(probability_given_neighborhood)
    log_probability_global = np.log(probability_given_global)
    
    any_nan = False
    if (np.isnan(log_probability_neighborhood)):
        any_nan=True
    if np.isnan(log_probability_global):
        any_nan=True
    if np.isnan(log_probability_neighborhood - log_probability_global):
        any_nan=True
    
    if any_nan:
        return 0
    
    else: 
        return log_probability_neighborhood - log_probability_global
       

# def psi_observations_scores_vec(
#     observed_psi_array, 
#     neighborhood_psi_array, 
#     global_psi, 
#     mrna_array, 
#     capture_efficiency, 
#     min_probability
# ):

#     # neighborhood_psi_array = np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in neighborhood_psi_array])

#     func = lambda i: psi_observation_score(
#         observed_psi_array[i], 
#         neighborhood_psi_array[i], 
#         global_psi, 
#         mrna_array[i], 
#         capture_efficiency, 
#         min_probability
#     )
#     x = range(len(observed_psi_array))
#     x_func = np.vectorize(func)
#     psi_scores_vec = x_func(x)

#     return psi_scores_vec
    

    
##################### Psix light ####################

@jit(nopython=True)
def get_arrays(observed_psi_array, mrna_obs_array, cell_metric):
    psi_o_array = []
    psi_a_array = []
    mrna_array = []
    for i in range(len(observed_psi_array)):
        psi_o = observed_psi_array[i]
        if not np.isnan(psi_o):
            
            psi_o_array.append(psi_o)
            mrna_array.append(mrna_obs_array[i])
            
            neighbors = cell_metric[0][i]
            weights = cell_metric[1][i]

            psi_sum = 0
            weight_sum = 0
            for j in range(1, len(neighbors)):
                psi_n = observed_psi_array[neighbors[j]]
                if not np.isnan(psi_n):
                    psi_sum += (psi_n * weights[j])
                    weight_sum += weights[j]
            if weight_sum > 0:
                psi_a_array.append(psi_sum/weight_sum)
            else:
                psi_a_array.append(np.nan)
                
    mrna_array = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in mrna_array])
    
    return psi_o_array, psi_a_array, mrna_array


# @njit
# def permute_array_fix_nan(idx_array):

#     idx_non_nan = np.array([i for i in range(len(idx_array)) if not np.isnan(idx_array[i])])
#     idx_permute = np.random.permutation(idx_non_nan)

#     new_idx = []
#     j = 0
#     for i in range(len(idx_array)):
#         if np.isnan(idx_array[i]):
#             new_idx.append(i)
#         else:
#             new_idx.append(int(idx_permute[j]))
#             j+=1

#     return np.array(new_idx)
        
    
def psix_score(
    exon_psi_array, 
    exon_mrna_array, 
    cell_metric, 
    capture_efficiency = 0.1, 
    randomize = False,  
    min_probability = 0.01,
    seed=0,
    turbo = 'lookup/',
    cap_mrna = True,
    no_lookup=False
):
    

    ncells = len(exon_psi_array)

    if randomize:
        np.random.seed(seed)
        shuffled_cells = shuffle(range(ncells))
        exon_psi_array = exon_psi_array[shuffled_cells]
        exon_mrna_array = exon_mrna_array[shuffled_cells]


    observed_psi_array, neighborhood_psi_array, mrna_array = get_arrays(exon_psi_array, exon_mrna_array, cell_metric)

    observed_psi_array = list(observed_psi_array)
    neighborhood_psi_array = list(neighborhood_psi_array)
    mrna_array = list(np.round(mrna_array).astype(int))
    

    total_cells = round((len(observed_psi_array) - np.sum(mrna_array == 0)))

    if total_cells <= 0:
        return np.nan

    global_psi = np.nanmean(observed_psi_array)

    max_mrna = len(turbo)

    # L_vec = psi_observations_scores_vec(
    #     observed_psi_array, 
    #     neighborhood_psi_array, 
    #     global_psi, 
    #     mrna_array, 
    #     capture_efficiency, 
    #     turbo, 
    #     max_mrna = max_mrna, 
    #     min_probability=0.01
    # )

    neighborhood_psi_array = [0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in neighborhood_psi_array]

    if cap_mrna:
        mrna_max = len(turbo)
        mrna_array = [mrna_max if x >= mrna_max else x for x in mrna_array]
    
    if no_lookup:

        logging.basicConfig(filename='/project2/yangili1/cfbuenabadn/analysis_psix/debug.log', level=logging.INFO)

        observed_psi_str = ','.join([str(x) for x in observed_psi_array])

        logging.info(type(observed_psi_array[0]))

        logging.info(observed_psi_array)

        logging.info(neighborhood_psi_array)

        logging.info(mrna_array)



        # neighborhood_psi_array = np.array(neighborhood_psi_array)
    
        L_vec = psi_observations_scores_vec(
            observed_psi_array, 
            neighborhood_psi_array, 
            global_psi, 
            mrna_array, 
            capture_efficiency, 
            min_probability
        )

        logging.info(L_vec)
        
    else:
        
        
        L_vec = psi_observations_scores_vec_turbo(
            observed_psi_array, 
            neighborhood_psi_array, 
            global_psi, 
            mrna_array,
            turbo
        )
        

        

    return np.sum(L_vec)/total_cells


##################### lookup functions ########################

@jit(nopython=True)
def psi_observations_scores_vec(
    observed_psi_array, 
    neighborhood_psi_array, 
    global_psi, 
    mrna_array, 
    capture_efficiency, 
    min_probability = 0.01
):

    neighborhood_psi_array = [0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in neighborhood_psi_array]
    
    L_vec = []
    
    for i in range(len(observed_psi_array)):
        psi_o = observed_psi_array[i]
        psi_a = neighborhood_psi_array[i]
        mrna = mrna_array[i]

        L = psi_observation_score(
                psi_o, 
                psi_a, 
                global_psi, 
                mrna, 
                capture_efficiency, 
                min_probability
            )
        L_vec.append(L)
    return L_vec


##################### Turbo functions #######################

@jit(nopython=True)
def psix_turbo(psi_o, psi_a, mrna, turbo_dict):
    
    psi_o_idx = int(round_numba(psi_o*100))
    psi_a_idx = int(round_numba(psi_a*100))-1
    mrna_idx = mrna-1
    
    if psi_a_idx < 0:
        psi_a_idx = 0
    elif psi_a_idx > 98:
        psi_a_idx = 98
    psix_score_out = turbo_dict[mrna_idx][psi_o_idx, psi_a_idx]
    return psix_score_out

@jit(nopython=True)
def L_score_lookup(psi_o, psi_a, psi_n, mrna, turbo_dict):
    L_a = np.log(psix_turbo(psi_o, psi_a, mrna, turbo_dict))
    L_null = np.log(psix_turbo(psi_o, psi_n, mrna, turbo_dict))
    
    return L_a - L_null

@jit(nopython=True)
def psi_observations_scores_vec_turbo(observed_psi_array, neighborhood_psi_array, global_psi, mrna_array, turbo_dict):
    L_vec = []
    
    for i in range(len(observed_psi_array)):
        psi_o = observed_psi_array[i]
        psi_a = neighborhood_psi_array[i]
        mrna = mrna_array[i]
        L = L_score_lookup(psi_o, psi_a, global_psi, mrna, turbo_dict)
        L_vec.append(L)
    return L_vec


