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

tiklova_mrna_event = pd.read_csv(data_dir + 'tiklova_neurogenesis/mrna_per_event.tab', sep='\t', index_col=0)
tiklova_rd = pd.read_csv(data_dir + 'tiklova_neurogenesis/rd_pc2.tab', sep='\t', index_col=0).loc[tiklova_mrna_event.columns]
tiklova_PSI = pd.read_csv(data_dir + 'tiklova_neurogenesis/skipped_exons_psi.tab', sep='\t', index_col=0)[tiklova_mrna_event.columns]
tiklova_psix = pd.read_csv('~/psix/data/tiklova_neurogenesis/tiklova_neurogenesis.scores.txt', sep='\t', index_col=0)
tiklova_mrna = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/mrna_counts.tab', sep='\t', index_col=0)[tiklova_rd.index]
tiklova_reads = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/tiklova/skipped_exons_SJreads.tab', 
                            sep='\t', index_col=0)[tiklova_rd.index]
tiklova_SJ_reads = pd.read_csv('~/data_sc_regulation/tiklova_extended/SE_counts.tab.gz', sep='\t', index_col=0)[tiklova_rd.index]



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

def downsample_exon(reads_table, mrna_table, exon, cells, downsample_rate = 0.5):
    exon_mrnas = round(mrna_table.loc[exon, cells].fillna(0)).astype(int)
    sub_mrnas = np.random.binomial(exon_mrnas, downsample_rate)
    ratios = [sub_mrnas[i]/exon_mrnas[i] if exon_mrnas[i]>0 else 0 for i in range(len(exon_mrnas))]
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

def make_random(psi_table, mrna_table, read_table, exon):
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

            lim_pc_1 = 0
            lim_pc_2 = 0
            coin_2_pc1 = 2
            coin_2_pc2 = 2
            cells_1 = []
            cells_2 = []

            coin_toss_1 = np.random.choice(range(4))
            if (coin_toss_1 == 0) or (coin_toss_1 == 2):
                lim_pc_1 = np.random.uniform(rd.PC_1.quantile(0.3), rd.PC_1.quantile(0.7))
                coin_2_pc1 = np.random.binomial(1, 0.5)

                if coin_2_pc1 == 0:
                    cells_1 = rd.index[rd.PC_1 < lim_pc_1]
                elif coin_2_pc1 == 1:
                    cells_1 = rd.index[rd.PC_1 > lim_pc_1]

            if (coin_toss_1 == 1) or (coin_toss_1 == 2):
                lim_pc_2 = np.random.uniform(rd.PC_2.quantile(0.3), rd.PC_2.quantile(0.7))
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
    return random_psi, random_mrna, shuffled_cells
        
    
idx = tiklova_psix.sort_values('L_score').index
random_psi_1, random_mrna_1, shuffled_cells = make_test(tiklova_PSI.loc[idx], 
                                        tiklova_mrna_event.loc[idx], 
                                    tiklova_SJ_reads, idx, tiklova_rd[['PC_1', 'PC_2']], seed=0, prob=0.1)

r_psi_1 = tiklova_PSI[shuffled_cells].copy()
r_psi_1.columns = tiklova_PSI.columns
r_mrna_1 = tiklova_mrna_event[shuffled_cells].copy()
r_mrna_1.columns = tiklova_mrna_event.columns

# random_psi_2, random_mrna_2, shuffled_cells = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
#                                         tiklova_mrna_event.loc[tiklova_L_score.index], 
#                                     tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=1, prob=0.1)

# r_psi_2 = tiklova_PSI[shuffled_cells].copy()
# r_psi_2.columns = tiklova_PSI.columns
# r_mrna_2 = tiklova_mrna_event[shuffled_cells].copy()
# r_mrna_2.columns = tiklova_mrna_event.columns

# random_psi_3, random_mrna_3, shuffled_cells = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
#                                         tiklova_mrna_event.loc[tiklova_L_score.index], 
#                                     tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=2, prob=0.1)

# r_psi_3 = tiklova_PSI[shuffled_cells].copy()
# r_psi_3.columns = tiklova_PSI.columns
# r_mrna_3 = tiklova_mrna_event[shuffled_cells].copy()
# r_mrna_3.columns = tiklova_mrna_event.columns

# random_psi_4, random_mrna_4, shuffled_cells = make_test(tiklova_PSI.loc[tiklova_L_score.index], 
#                                         tiklova_mrna_event.loc[tiklova_L_score.index], 
#                                     tiklova_SJ_reads, tiklova_L_score.index, tiklova_rd[['PC_1', 'PC_2']], seed=3, prob=0.1)

# r_psi_4 = tiklova_PSI[shuffled_cells].copy()
# r_psi_4.columns = tiklova_PSI.columns
# r_mrna_4 = tiklova_mrna_event[shuffled_cells].copy()
# r_mrna_4.columns = tiklova_mrna_event.columns


# random_psi_1.index = [x + '_random1' for x in random_psi_1.index]
# random_mrna_1.index = [x + '_random1' for x in random_mrna_1.index]

# random_psi_2.index = [x + '_random2' for x in random_psi_2.index]
# random_mrna_2.index = [x + '_random2' for x in random_mrna_2.index]

# random_psi_3.index = [x + '_random3' for x in random_psi_3.index]
# random_mrna_3.index = [x + '_random3' for x in random_mrna_3.index]

# random_psi_4.index = [x + '_random4' for x in random_psi_4.index]
# random_mrna_4.index = [x + '_random4' for x in random_mrna_4.index]


# r_psi_1.index = [x + '_random1' for x in r_psi_1.index]
# r_mrna_1.index = [x + '_random1' for x in r_1.index]

# r_psi_2.index = [x + '_random2' for x in r_psi_2.index]
# r_mrna_2.index = [x + '_random2' for x in r_mrna_2.index]

# r_psi_3.index = [x + '_random3' for x in r_psi_3.index]
# r_mrna_3.index = [x + '_random3' for x in r_mrna_3.index]

# r_psi_4.index = [x + '_random4' for x in r_psi_4.index]
# r_mrna_4.index = [x + '_random4' for x in r_mrna_4.index]


random_psi_1.to_csv('tiklova_random_subsampled_psi.tab', sep='\t', index=True, header=True)

random_mrna_1.to_csv('tiklova_random_subsampled_mrna.tab', sep='\t', index=True, header=True)

r_psi_1.to_csv('tiklova_random_psi.tab', sep='\t', index=True, header=True)

r_mrna_1.to_csv('tiklova_random_mrna.tab', sep='\t', index=True, header=True)

