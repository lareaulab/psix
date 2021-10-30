import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm 
from numba import jit

def compute_cell_metric(
    manifold, 
    n_neighbors=100, 
    weight_metric = True,
):
    """
    Computes cell-cell metric from low-dimensional manifold.
    
    input:
    
        manifold: a pandas dataframe with a N-dimensional manifold. E.g., a PCA projection, or an scVI space.
                  Shape = (cells, dimensions).
                  
        n_neighbors: number of neighbors per neighborhood. It's recommended to be larger than those
                     normally used for gene expression analysis.
                     
        weight_metric: Wether if the influence of neighbours should be weighted by their proximity.
                       Similar as the weighted metric from VISION.If false, all neighbors have the 
                       same influence.
                       Non-neighbors always have 0 weight.
                   
    output: 
    
        cell_metric: a dataframe structure of (cells, cells) dimensions with the weights of
                     each neighbor.
    """
    
    cells = manifold.index
    n_cells = len(cells)
    
    knn_neighbors = NearestNeighbors(n_neighbors=n_neighbors+1).fit(manifold)
    distances, indices = knn_neighbors.kneighbors(manifold)
    
    neighbor_indices = pd.DataFrame(indices, index=cells)
    
    weights = np.ones((len(cells), (n_neighbors+1)))
    
    for i in tqdm(range(len(manifold.index)), position=0, leave=True):
        sigma = np.max(distances[i])
        for j in range(1, len(distances[i])):
            d = distances[i][j]
            w = compute_weight(d, sigma)
            weights[i, j] = w
        
    cell_metric = (indices, weights)
    return cell_metric


@jit(nopython=True)
def compute_weight(d, sigma):
    return np.exp(-(d**2)/(sigma**2)) 

@jit(nopython=True)
def get_exon_neighbors_psi(observed_psi_array, cell_metric):
    psi_a_array = []
    for i in range(len(observed_psi_array)):

        neighbors = cell_metric[0][i]
        weights = cell_metric[1][i]

        psi_sum = 0
        weight_sum = 0
        for j in range(len(neighbors)):
            psi_n = observed_psi_array[neighbors[j]]
            if not np.isnan(psi_n):
                psi_sum += (psi_n * weights[j])
                weight_sum += weights[j]
        if weight_sum > 0:
            psi_a_array.append(psi_sum/weight_sum)
        else:
            psi_a_array.append(np.nan)
                
    return psi_a_array


def get_all_exons_neighbors(psi, cell_metric):
    
    neighbors_psi = []
    
    for exon_psi in psi.T:
        neighbors_psi.append(get_exon_neighbors_psi(exon_psi, cell_metric))
        
    return np.array(neighbors_psi)


def get_background(self):
    psi = np.array(self.adata.uns['psi'])
    
    neighbors_psi = pd.DataFrame(get_all_exons_neighbors(psi, self.metric),
                                 columns = self.adata.uns['psi'].index,
                                 index = self.adata.uns['psi'].columns)
    self.adata.uns['neighbors_psi'] = neighbors_psi.T
    
    

