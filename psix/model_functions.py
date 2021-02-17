import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
from scipy.special import comb
from sklearn.utils import shuffle
import itertools  

from numba import jit

@jit
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

@jit
def probability_psi_m_unknown(observed_psi, true_psi, capture_efficiency, captured_mrna):
    '''
    
    
    '''
    m_array = np.arange(captured_mrna, np.int(10*captured_mrna/capture_efficiency))
    
    comb_1 = comb(m_array*true_psi, captured_mrna*observed_psi)
    comb_2 = comb(m_array*(1-true_psi), captured_mrna*(1-observed_psi))
    proba_1 = capture_efficiency**(captured_mrna+1)
    proba_2 = (1-capture_efficiency)**(m_array-captured_mrna)
    
    prob_array = comb_1*comb_2*proba_1*proba_2
    
    return np.sum(prob_array)
    
@jit
def probability_psi_m_known(observed_psi, true_psi, captured_mrna, cell_molecules):
    """
    Pr(observed_PSI | true_PSI, captured_mrna, cell_molecules)
    
    """
    comb_1 = comb(cell_molecules*true_psi, captured_mrna*observed_psi)
    comb_2 = comb(cell_molecules*(1-true_psi), captured_mrna*(1-observed_psi))
    comb_3 = comb(cell_molecules, captured_mrna)**(-1)
    return comb_1*comb_2*comb_3

@jit
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
    
    probability_given_neighborhood = np.max([
        min_probability, probability_psi_observation(observed_psi, neighborhood_psi, capture_efficiency, captured_mrna)
    ])
    probability_given_global = np.max([
        min_probability, probability_psi_observation(observed_psi, global_psi, capture_efficiency, captured_mrna)
    ])
    
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
        

def psi_observations_scores_vec(
    observed_psi_array, 
    neighborhood_psi_array, 
    global_psi, 
    mrna_array, 
    capture_efficiency, 
    min_probability
):

    neighborhood_psi_array = np.array([0.99 if x > 0.99 else 0.01 if x < 0.01 else x for x in neighborhood_psi_array])

    func = lambda i: psi_observation_score(
        observed_psi_array[i], 
        neighborhood_psi_array[i], 
        global_psi, 
        mrna_array[i], 
        capture_efficiency, 
        min_probability
    )
    x = range(len(observed_psi_array))
    x_func = np.vectorize(func)
    psi_scores_vec = x_func(x)

    return psi_scores_vec
    

def psix_score(
    exon_psi_array, 
    exon_mrna_array, 
    cell_metric, 
    capture_efficiency = 0.1, 
    randomize = False,  
    min_probability = 0.01,
    seed=0
):
    
#     try:
    cell_list = exon_psi_array.dropna().index

    if randomize:
        np.random.seed(seed)
        shuffled_cells = shuffle(cell_list)

    else:
        shuffled_cells = cell_list

    observed_psi_array = np.array(exon_psi_array.loc[shuffled_cells])

    mrna_array = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in exon_mrna_array.loc[shuffled_cells]])
    mrna_array = np.round(mrna_array).astype(int)

    total_cells = round((len(cell_list) - np.sum(mrna_array == 0)))

    if total_cells <= 0:
        return np.nan

    neighborhood_psi_array = np.array(
    np.array(pd.DataFrame(
        np.array(cell_metric.loc[cell_list, cell_list])*observed_psi_array).sum(axis=1)
    )/np.array(cell_metric.loc[cell_list, cell_list].sum(axis=1)))

    global_psi = np.mean(observed_psi_array)

    L_vec = psi_observations_scores_vec(
        observed_psi_array, 
        neighborhood_psi_array, 
        global_psi, 
        mrna_array, 
        capture_efficiency, 
        min_probability
    )

    return np.sum(L_vec)/total_cells


    

    
### Deprecated, do not use. Gives the same result, but slower.
def psix_score_precomputed_smooth(
    observed_psi_array, 
    exon_mrna_array, 
    neighborhood_psi_array, 
    capture_efficiency = 0.1, 
    min_probability = 0.01,
):
    
#     try:
    
#     observed_psi_array = exon_psi_array

    mrna_array = np.array([1 if ((x > 0.1) and (x <= 1)) else x for x in exon_mrna_array])
    mrna_array = np.round(mrna_array).astype(int)

    total_cells = np.sum(~np.isnan(observed_psi_array) & (mrna_array > 0))#round((len(cell_list) - np.sum(mrna_array == 0)))

    if total_cells <= 0:
        return np.nan

    

    global_psi = np.nanmean(observed_psi_array)

    L_vec = psi_observations_scores_vec(
        observed_psi_array, 
        neighborhood_psi_array, 
        global_psi, 
        mrna_array, 
        capture_efficiency, 
        min_probability
    )

    return np.sum(L_vec)/total_cells

#     except:
#         return np.nan
