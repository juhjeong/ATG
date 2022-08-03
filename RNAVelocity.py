import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

#Import loom files and metadata
sample_one=anndata.read_loom("CONTROL.loom", validate=False)
sample_two=anndata.read_loom("GEM.loom",validate=False)
sample_three=anndata.read_loom("AT.loom",validate=False)
sample_one.var_names_make_unique()
sample_two.var_names_make_unique()
sample_three.var_names_make_unique()
sample_obs=pd.read_csv("cellID_obs.csv")
umap_cord=pd.read_csv("cell_embeddings.csv")
cell_clusters=pd.read_csv("clusters.csv")

#Filter cells of interest and concatenate
cellID_obs_sample_one=sample_obs[sample_obs["x"].str.contains("CONTROL_")]
cellID_obs_sample_two=sample_obs[sample_obs["x"].str.contains("GEM_")]
cellID_obs_sample_three=sample_obs[sample_obs["x"].str.contains("AT_")]
sample_one=sample_one[np.isin(sample_one.obs.index,cellID_obs_sample_one)]
sample_two=sample_two[np.isin(sample_two.obs.index,cellID_obs_sample_two)]
sample_three=sample_three[np.isin(sample_three.obs.index,cellID_obs_sample_three)]
sample_one=sample_one.concatenate(sample_two,sample_three)
sample_one.obs.set_index('Unnamed: 0',inplace=True)
sample_one_index=pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {'Unnamed: 0':'Cell ID'})
umap=umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered=sample_one_index.merge(umap, on = "Cell ID")
cell_clusters=cell_clusters.rename(columns={'Unnamed: 0':'Cell ID'})
cell_clusters=cell_clusters.rename(columns={'x':'seurat_clusters'})
cell_clusters = sample_one_index.merge(cell_clusters,on="Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap']=umap_ordered.values
cell_clusters_ordered=cell_clusters.iloc[:,1:]
sample_one.obs['cell_clusters']=cell_clusters_ordered.values

#Calculate velocity and embedding
scv.pp.filter_and_normalize(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one,mode="stochastic")
scv.tl.velocity_graph(sample_one)
scv.pl.velocity_embedding(sample_one, basis='umap', color="cell_clusters")
ident_colours = ["#F48FEB","#FF514C","#04B4EC","#FF9F25","#a433f5","#6ACB61"]
scv.pl.velocity_embedding_stream(sample_one,basis="umap",color="cell_clusters",palette = ident_colours,size=70,legend_loc='none')

