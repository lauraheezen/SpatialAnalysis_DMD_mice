# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 18:53:03 2021

@author: trmabdelaal
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import seaborn as sns

adata = scv.read('D2mdx_adata.h5ad')

sc.pp.filter_genes(adata, min_cells=10)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata) 

sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata)  
sc.pl.umap(adata, color='relabelled')

sc.pl.scatter(adata, basis='xy_loc',color='relabelled')
plt.gca().invert_yaxis()

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Linear model
scv.tl.velocity(adata)

# Dynamical model
# scv.tl.recover_dynamics(adata)
# scv.tl.velocity(adata, mode='dynamical')

scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='relabelled',legend_loc='right')

scv.pl.velocity_embedding_stream(adata, basis='xy_loc', color='relabelled',legend_loc='right')
plt.gca().invert_yaxis()

scv.pl.velocity_embedding(adata,basis='xy_loc', color='relabelled')
plt.gca().invert_yaxis()

scv.pl.proportions(adata,groupby='relabelled')

# Proportions per cell
layers = ["spliced", "unspliced", "ambigious"]
layers_keys = [key for key in layers if key in adata.layers.keys()]
counts_layers = [np.sum(adata.layers[key], axis=1) for key in layers_keys]
counts_total = np.sum(counts_layers, 0)
counts_total += counts_total == 0
counts_layers = np.array([counts / counts_total for counts in counts_layers])

fig,ax=plt.subplots()
cax= ax.scatter(adata.obsm['X_xy_loc'][:,0], adata.obsm['X_xy_loc'][:,1],c=counts_layers[1,:]*100,cmap=sns.cm.flare)
ax.invert_yaxis()
#cax.set_clim(0,np.max(counts_layers[1,:])*100)
cax.set_clim(0,35)
cbar = fig.colorbar(cax)
plt.title('Precentage of unspliced per spot')


# RNA velocity Magnitude
Mag_xy = np.linalg.norm(adata.obsm['velocity_xy_loc'],axis=1)

fig,ax=plt.subplots()
cax= ax.scatter(adata.obsm['X_xy_loc'][:,0], adata.obsm['X_xy_loc'][:,1],c=np.log1p(Mag_xy),cmap=sns.cm.flare)
ax.invert_yaxis()
cax.set_clim(0,np.max(np.log1p(Mag_xy)))
cbar = fig.colorbar(cax, ticks=[0,2,4])
plt.title('Magnitude of Spatial RNA velocity vectors (Log_scaled)')

# Top contributing genes driving RNA velocity
Velo_genes = adata.layers['velocity']
Velo_genes = Velo_genes[adata.obs.relabelled=='Inflamed and/or calcified fibers',:]

adata.obs['Selected'] = adata.obs.relabelled=='Inflamed and/or calcified fibers'
scv.pl.scatter(adata,basis='xy_loc', color='Selected')
plt.gca().invert_yaxis()

Velo_genes[np.isnan(Velo_genes)] = 0
Mag = np.linalg.norm(Velo_genes,axis=1)
for i in range(Velo_genes.shape[0]):
    Velo_genes[i,:] = Velo_genes[i,:]/Mag[i]
Velo_genes = np.square(Velo_genes)

Contributing_genes=[]
for i in range(Velo_genes.shape[0]):
    rank = np.argsort(-1*Velo_genes[i,:])[:5] # Top 5
    Contributing_genes.append(np.array(adata.var_names[rank]))
flat_list = pd.Series([item for sublist in Contributing_genes for item in sublist])

flat_list.value_counts() # count list

flat_list.value_counts()/np.sum(adata.obs.relabelled=='Inflamed and/or calcified fibers')*100 # Percentage list
