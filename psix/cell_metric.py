import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm 


def compute_cell_metric(
    manifold, 
    n_neighbors=100, 
    weight_metric = True
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
    
    knn_neighbors = NearestNeighbors(n_neighbors=n_neighbors).fit(manifold)
    distances, indices = knn_neighbors.kneighbors(manifold)
    
    cell_metric = pd.DataFrame(np.zeros((n_cells, n_cells)))
    cell_metric.columns = cells
    cell_metric.index = cells
    
    if weight_metric:
        
        for i in tqdm(range(n_cells)):
            cell_i = cells[i]
            sigma = np.max(distances[i])
            for j in range(1, len(distances[i])):
                cell_j = cells[indices[i][j]]
                d = distances[i][j]
                w = np.exp(-(d**2)/(sigma**2))        
                cell_metric.loc[cell_i, cell_j] = w
                
    else:
        for i in tqdm(range(n_cells)):
            cell_i = cells[i]
            for j in range(1, len(distances[i])):
                cell_j = cells[indices[i][j]]     
                cell_metric.loc[cell_i, cell_j] = 1
    
    return cell_metric