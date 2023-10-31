import os 
import numpy as np
import pandas as pd
from .model_functions import probability_psi_observation
from tqdm import tqdm


def load_turbo(turbo_dir = 'lookup/', max_mrna = 30):
    turbo_files = ['psix_'+str(i)+'.tab.gz' for i in range(1, max_mrna+1)]
    turbo_dict = []
    for x in tqdm(turbo_files, position=0, leave=True):
        mrna_counts = int(x.split('.')[0].split('_')[-1])
        turbo_table = np.array(pd.read_csv(turbo_dir + x, sep='\t', index_col=0))
        turbo_dict.append(turbo_table)
    return turbo_dict


def make_turbo_function(out_dir = 'lookup/', granularity = 0.01, max_mrna = 30, capture_efficiency=0.1, min_probability=0.01):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    observed_range = np.arange(0, 1+granularity, granularity)
    true_range = np.arange(granularity, 1, granularity)
    mrna_range = np.arange(1, max_mrna+1)
    
    for mrna in tqdm(mrna_range, leave=True, position=0):
        mrna_df = pd.DataFrame(np.ones((len(observed_range), len(true_range)))*min_probability, 
                               index=observed_range, columns=true_range)
        for observed_psi in observed_range:
            for true_psi in true_range:
                probability = probability_psi_observation(observed_psi, true_psi, capture_efficiency, mrna)
                mrna_df.loc[observed_psi, true_psi] = np.max((mrna_df.loc[observed_psi, true_psi], probability))

        mrna_df.to_csv(out_dir + 'psix_' + str(mrna)+'.tab.gz', sep='\t', index=True, header=True)
