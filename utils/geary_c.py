import importlib
import numpy as np
import pandas as pd
import os
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy import stats as st
import seaborn as sns
import argparse

import numpy.random as r

import numpy as np
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

from scipy.stats import combine_pvalues
# %run -i '../../utils/Kruskal_Wallis_test_functions.py'

from tqdm import tqdm

############

from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import euclidean_distances


def get_distance_matrix(pca, k=10):
    nbrs = NearestNeighbors(n_neighbors=k).fit(pca[['PC_1', 'PC_2']])
    distances, indices = nbrs.kneighbors(pca[['PC_1', 'PC_2']])
    
    cells = list(pca.index)
    
    W = pd.DataFrame(np.zeros((len(cells), len(cells))))
    W.columns = cells
    W.index = cells
    
    for i in tqdm(range(len(cells))):
        cell_i = cells[i]
        sigma = np.max(distances[i])
        for j in range(len(distances[i])):
            cell_j = cells[indices[i][j]]
            d = distances[i][j]
            w = np.exp(-(d**2)/(sigma**2))
            if cell_i == cell_j:
                W.loc[cell_i, cell_j] = 0
            else:
                W.loc[cell_i, cell_j] = w
    
    return W


def get_signature_matrix(PSI_tab):
    return (PSI_tab - PSI_tab.mean())/PSI_tab.std()


def make_mock_C_scores(norm_PSI, Ws, exon_list, total_cells, mock=100000):
    exon_out_list = []
    C_scores = []
    for i in tqdm(range(mock)):
        mock_run = True
        while mock_run:
            exon = r.choice(exon_list, 1)[0]
#             print(exon)
            scramble_cells = r.choice(norm_PSI.columns, total_cells, replace=False)
            mock_PSI = pd.DataFrame(norm_PSI.loc[exon, scramble_cells]).T

#             print(mock_PSI.shape)
#             print(norm_PSI.shape)

            mock_PSI.columns = norm_PSI.columns
            mock_df = mock_PSI.loc[exon]
#             print(type(mock_df))
#             print(mock_df)
            mock_score = get_C(mock_df, Ws)
            
#             print(mock_score)

            if mock_score >= 0:
                C_scores.append(mock_score)
                exon_out_list.append('mock_'+exon+'_'+str(i))
                mock_run = False
    return exon_out_list, C_scores
                    
def get_C(exon_score, W):
    exon_score = exon_score.dropna()
    obs_cells = exon_score.index
    x = (exon_score.values.reshape(-1, 1) - exon_score.values.reshape(1, -1))
    w = W.loc[obs_cells, obs_cells]
    num = (len(obs_cells)-1)*((w*(x**2)).sum().sum())
    den = (2*w.sum().sum())*np.sum((exon_score - exon_score.mean())**2)
    C = num/den
    score = 1 - C
    return score
    
    
##############################


# Ahora si el bueno

def get_mock_dict(PSI_tab, norm_PSI, Ws, mock=200):

    total_cells = len(PSI_tab.columns)

    exons_05_10 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.4) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.45)] 
    exons_10_20 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.3) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.40)] 
    exons_20_30 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.2) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.30)] 
    exons_30_40 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.1) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.20)] 
    exons_40_50 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.1)] 

    exons_obs_10_20 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.9) & (PSI_tab.isna().mean(axis=1) > 0.8)]
    exons_obs_20_30 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.8) & (PSI_tab.isna().mean(axis=1) > 0.7)]
    exons_obs_30_40 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.7) & (PSI_tab.isna().mean(axis=1) > 0.6)]
    exons_obs_40_50 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.6) & (PSI_tab.isna().mean(axis=1) > 0.5)]
    exons_obs_50_60 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.5) & (PSI_tab.isna().mean(axis=1) > 0.4)]
    exons_obs_60_70 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.4) & (PSI_tab.isna().mean(axis=1) > 0.3)]
    exons_obs_70_80 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.3) & (PSI_tab.isna().mean(axis=1) > 0.2)]
    exons_obs_80_90 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.2) & (PSI_tab.isna().mean(axis=1) > 0.1)]
    exons_obs_90_100 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.1)]

    list1 = [exons_05_10, exons_10_20, exons_20_30, exons_30_40, exons_40_50]
    list2 = [exons_obs_10_20, exons_obs_20_30, exons_obs_30_40, exons_obs_40_50,
             exons_obs_50_60, exons_obs_60_70, exons_obs_70_80, exons_obs_80_90, exons_obs_90_100]


    exon_out_list = []
    C_score_list = []

    for lista_1 in list1:
        for lista_2 in list2:
            combination = lista_1 & lista_2
            if len(combination) > 0:
                exon_out, C_scores = make_mock_C_scores(norm_PSI, Ws, lista_1&lista_2, total_cells, mock=mock)
                exon_out_list.append(exon_out)
                C_score_list.append(C_scores)
            else:
                exon_out_list.append([])
                C_score_list.append([])


    psi_key = ['psi_05_10', 'psi_10_20', 'psi_20_30', 'psi_30_40', 'psi_40_50']
    obs_key = ['obs_10_20', 'obs_20_30', 'obs_30_40', 'obs_40_50', 
               'obs_50_60', 'obs_60_70', 'obs_70_80', 'obs_80_90', 'obs_90_100']

    counter = 0
    mock_dict = {}
    for pk in psi_key:
        obs_dict = {}
        for ok in obs_key:
            
            #fit_alpha, fit_loc, fit_beta = st.gamma.fit(C_score_list[counter])
            #random_data = st.gamma.rvs(fit_alpha, loc=fit_loc, scale=fit_beta, size=1000000)
            random_data = C_score_list[counter]
            obs_dict.update({ok:random_data})
            counter += 1
        mock_dict.update({pk:obs_dict})

    return mock_dict



#######################


def get_C_score_pval_gamma(PSI_tab, norm_PSI, Ws, exon_list, total_cells, mock_dict):
    
    
    exon_out_list = []
    C_list = []
    p_list = []
    
    for exon in tqdm(exon_list):
        psi_mean = PSI_tab.loc[exon].mean()
        obs_mean = PSI_tab.loc[exon].isna().mean()
        
        
        
        exon_df = norm_PSI.loc[exon] # to make things faster
        exon_score = get_C(exon_df, Ws)
#         if exon_score >= 0:
        C_list.append(exon_score)
        exon_out_list.append(exon)

        if (np.abs(0.5 - psi_mean) > 0.4) and (np.abs(0.5 - psi_mean) <= 0.45):
            pk = 'psi_05_10'
        elif (np.abs(0.5 - psi_mean) > 0.3) and (np.abs(0.5 - psi_mean) <= 0.4):
            pk = 'psi_10_20'
        elif (np.abs(0.5 - psi_mean) > 0.2) and (np.abs(0.5 - psi_mean) <= 0.3):
            pk = 'psi_20_30'
        elif (np.abs(0.5 - psi_mean) > 0.1) and (np.abs(0.5 - psi_mean) <= 0.2):
            pk = 'psi_30_40'
        elif (np.abs(0.5 - psi_mean) <= 0.1):
            pk = 'psi_40_50'

        if (obs_mean <= 0.9) and (obs_mean > 0.8):
            ok = 'obs_10_20'
        elif (obs_mean <= 0.8) and (obs_mean > 0.7):
            ok = 'obs_20_30'
        elif (obs_mean <= 0.7) and (obs_mean > 0.6):
            ok = 'obs_30_40'
        elif (obs_mean <= 0.6) and (obs_mean > 0.5):
            ok = 'obs_40_50'
        elif (obs_mean <= 0.5) and (obs_mean > 0.4):
            ok = 'obs_50_60'
        elif (obs_mean <= 0.4) and (obs_mean > 0.3):
            ok = 'obs_60_70'
        elif (obs_mean <= 0.3) and (obs_mean > 0.2):
            ok = 'obs_70_80'
        elif (obs_mean <= 0.2) and (obs_mean > 0.1):
            ok = 'obs_80_90'
        elif (obs_mean <= 0.1):
            ok = 'obs_90_100'
        else:
            ok = 'obs_10_20'


        random_data = mock_dict[pk][ok]

        x = np.sum(random_data > exon_score)
        n = len(random_data)
        pv = (x+1)/(n+1)

        p_list.append(pv)
            
    pval_df = pd.DataFrame()
    pval_df['C_score'] = C_list
    pval_df['pvals'] = p_list
    pval_df.index = exon_out_list
    return pval_df

     
    
parser = argparse.ArgumentParser(description='Psix: autocorrelated exon discovery from scRNA-seq data.')

parser.add_argument('-psi', '--psi_table', type=str, required=True,
                   help='Table with exon PSI values per cell. Right now Psix supports only skipped exon events.')

parser.add_argument('-rd', '--reduced_dimensions', type=str, required=True,
                   help='Table with a low-dimensionality projection of the cells (e.g., first 5 principal components).')

parser.add_argument('-o', '--output_name', type=str, required=False, default='gearyc_output',
                   help='Name for output files.')

parser.add_argument('-p', '--permutations', type=int, required=False, default = 1000,
                   help='How many permutations per bin. Default: 1000 permutations.')

parser.add_argument('-k', '--k_nearest_neighbors', type=int, required=False, default = 100,
                   help='K nearest neighbors. Set k==0 to use the squared root of the number of cells.')

parser.add_argument('-n', '--max_missing', type=float, required=False, default = 0.25,
                   help='Maximum percent of missing values per exon.')

if __name__ == '__main__':
    args = parser.parse_args()
    
    psi_table = pd.read_csv(args.psi_table, sep='\t', index_col=0)
    
    rd = pd.read_csv(args.reduced_dimensions, sep='\t', index_col=0)
    
    
    # tiklova_PSI = pd.read_csv('../data/preprocess/tables/psi.tab.gz', sep='\t', index_col=0)

    # tiklova_pca = pd.read_csv('../data/preprocess/tables/pc2_rd.tab.gz', sep='\t', index_col=0)

    exons = psi_table.index#[np.abs(0.5 - psi_table.mean(axis=1)) <= 0.45] & psi_table.index[psi_table.isna().mean(axis=1) <= (1- args.max_missing)]


    #tiklova_exons = tiklova_PSI.index[np.abs(0.5 - tiklova_PSI.mean(axis=1)) <= 0.45] & tiklova_PSI.index[tiklova_PSI.isna().mean(axis=1) < 0.90]

    psi_table = psi_table.loc[exons]
    print(psi_table.shape)

    k = args.k_nearest_neighbors #int(round(np.sqrt(len(tiklova_PSI.columns))))

    W = get_distance_matrix(rd, k=k)
    norm_PSI = get_signature_matrix(psi_table)

    mock_dict = get_mock_dict(psi_table, norm_PSI, W, mock=args.permutations)

    c_table = get_C_score_pval_gamma(psi_table, norm_PSI, 
                                    W, psi_table.index, len(psi_table.columns), mock_dict)

    c_table.to_csv(args.output_name + '.tab.gz', sep='\t', header=True, index=True)

    
    
