import os
import numpy as np
from typing import Iterable
import scipy
from scipy.spatial import KDTree
import anndata as ad


## Extract results from complex dictionary
def get_all_values(d):
    if isinstance(d, dict):
        for v in d.values():
            yield from get_all_values(v)
    elif isinstance(d, Iterable) and not isinstance(d, str): # or list, set, ... only
        for v in d:
            yield from get_all_values(v)
    else:
        yield d 

## Subsample data using centroid distance in a latent space (scVI, UMAP, ...)
def cluster_centroid_subsampling(adata, cluster_name, latent_space, n=200):
    ##
    adata.obs[cluster_name] = adata.obs[cluster_name].astype("str")
    ## 
    sub_sampling = {}
    ## Loop through clusters and subsample to "n"
    print(f"Total number of {cluster_name} level hierarchy: {len(adata.obs[cluster_name].unique())}")
    for cluster in adata.obs[cluster_name].unique():
        if adata.obs[cluster_name].value_counts()[cluster] > n:
            ##
            adata_cluster = adata[(adata.obs[cluster_name] == cluster)]
            ##
            adata_centroid = np.median(adata_cluster.obsm[latent_space], axis=0)
            ##
            nnTree = scipy.spatial.KDTree(adata_cluster.obsm[latent_space])
            dist, nn_centroid = nnTree.query(adata_centroid, k=np.minimum(n, adata_cluster.shape[0]), workers=-1) ## This will use all CPUs!
            ##
            sub_sampling[cluster] = adata_cluster.obs.iloc[nn_centroid].index
        else:
            sub_sampling[cluster] = adata[(adata.obs[cluster_name] == cluster)].obs.index
    ##
    adata.strings_to_categoricals()
    ## 
    return sub_sampling


## Runs centroid subsampling and returns indexes of cell to include
def run_cluster_centroid_subsampling(adata, save_adata_filename, cluster_name="Cluster", latent_space="X_UMAP", n=200):
    adata_subsampling = cluster_centroid_subsampling(adata, cluster_name, latent_space, n)
    adata_subsample_index = list(get_all_values(adata_subsampling))
    adata_subsample = adata[adata_subsample_index].copy()
    print(f"Saving new subsampled h5ad file to: {save_adata_filename}")
    adata_subsample.write_h5ad(save_adata_filename)
    return adata_subsample
