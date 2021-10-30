import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["axes.facecolor"] = 'white'
import matplotlib as mpl
import numpy as np
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams['pdf.fonttype'] = 42

"""
DEPRECATED
"""

def plot_exon(self, exon, axis_list=[], cmap='viridis'):
    
    fig = plt.figure(figsize=(8, 5))
    

    ax.scatter(np.linspace(i-0.25, i+0.25, 26),
#                np.random.uniform(i-0.2, i+0.2, 25), 
               li, s=50, c=color_code, alpha=0.75, linewidth=0)
    

ax.set_xlabel('Modules', size=12)




ax.set_xlim([0, 11])
ax.set_xticks(range(1, 11))
ax.tick_params(labelsize=12, length=5)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_ylabel('-log10 qval', size=12)
# plt.legend()

# plt.show()
plt.savefig('plots/module_enrichment.png', bbox_inches='tight', res=20000, dpi =2000)
    
    
    

    
    
    
    
    
    if len(axis_list) == 0:
        axis_list = list(self.latent.columns)[:3]
        
    ax.set_xlabel(axis_list[0])
        
    if len(axis_list) == 1:
        ax  = plt.subplot(1,1,1)
        ax.scatter(self.latent[axis_list[0]], self.adata.uns['psi'][exon], s=50, c='navy', alpha=0.90, linewidth=0)
        
    elif len(axis_list) == 2:
        ax  = plt.subplot(1,1,1)
        ax.scatter(self.latent[axis_list[0]], self.latent[axis_list[1]], c = self.adata.uns['psi'][exon], 
                   s=50, c=cmap, alpha=0.90, linewidth=0)
        
        ax.set_ylabel(axis_list[1])
        
    elif len(axis_list) == 3:
        
        ax = fig.add_subplot(111, projection='3d')

        ax.grid(False)

        ax.scatter(self.latent[axis_list[0]], self.latent[axis_list[1]], self.latent[axis_list[2]],
                   c=self.adata.uns['psi'][exon], cmap=cmap)
        
        ax.set_ylabel(axis_list[1])
        ax.set_zlabel(axis_list[2])

        
        
    ax.set_title(exon)

        
        
        
        
        
        
        
        
    