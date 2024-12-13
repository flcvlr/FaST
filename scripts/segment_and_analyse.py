import numpy
import gc
import seaborn as sns
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import sys
import spateo as st
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
seed1 = 78
seed2 = seed1+10
dataset= sys.argv[1]
my_bin = int(sys.argv[2])

EM_kernel=int(sys.argv[3])
output='seg_k_'+str(EM_kernel)+'_binsize_'+str(my_bin)
if not os.path.exists(dataset+'/'+output):
    os.mkdir(dataset+'/'+output) 

nuc_markers_file = open(sys.path[0]+'/../reference/'+sys.argv[4]+'/'+sys.argv[4]+'.MIR_HG_and_APEX_nuc_markers')
nuc_markers = nuc_markers_file.read().strip().split("\n")

adata = st.io.read_bgi_agg(dataset+'/dge/whole_dataset.txt',gene_agg={'nuclear': nuc_markers})

### nuclear markers

st.cs.segment_densities(adata, 'nuclear', my_bin, k=3, dk=3, distance_threshold=3)
with open(dataset+'/images/nuclear_bins.txt', "w") as txt_file:
    for line in adata.layers['nuclear_bins']:
        line_str = [str(n) for n in line]
        txt_file.write(" ".join(line_str) + "\n")
 
st.cs.score_and_mask_pixels(adata, 'nuclear', k=EM_kernel, method='EM+BP', em_kwargs={'precision': 1e-06 , 'max_iter':15000 , 'seed' : seed1})
st.cs.label_connected_components(adata,'nuclear', min_area=100 ,area_threshold=600, max_area=2000)

st.cs.score_and_mask_pixels(adata, 'nuclear', k=EM_kernel, method='EM+BP',certain_layer='nuclear_labels', scores_layer = 'nuclear_scores_bis', mask_layer='nuclear_mask_bis', em_kwargs={'precision': 1e-06 , 'max_iter':15000 , 'seed' : seed2})
st.cs.label_connected_components(adata,'nuclear_mask_bis',min_area=100, area_threshold=600, max_area=2000, seed_layer='nuclear_labels', out_layer='nuclear_labels_bis')

del adata.layers['nuclear']
del adata.layers['nuclear_bins']
del adata.layers['nuclear_labels']
del adata.layers['nuclear_scores']
del adata.layers['nuclear_scores_bis']
del adata.layers['nuclear_mask']
gc.collect() 

st.cs.segment_densities(adata, 'unspliced', my_bin, k=3, dk=3, distance_threshold=3)
with open(dataset+'/images/unspliced_bins.txt', "w") as txt_file:
    for line in adata.layers['unspliced_bins']:
        line_str = [str(n) for n in line]
        txt_file.write(" ".join(line_str) + "\n")

st.cs.score_and_mask_pixels( adata, 'unspliced', k=EM_kernel, method='EM+BP', certain_layer='nuclear_labels_bis',em_kwargs={'precision': 1e-06, 'max_iter':15000 , 'seed' : seed1})
st.cs.label_connected_components(adata, 'unspliced', min_area=100 ,area_threshold=600, max_area=2000, seed_layer='nuclear_labels_bis')

st.cs.score_and_mask_pixels( adata, 'unspliced', k=EM_kernel, method='EM+BP', certain_layer='unspliced_labels', scores_layer = 'unspliced_scores_bis', mask_layer='unspliced_mask_bis',em_kwargs={'precision': 1e-06, 'max_iter':15000 , 'seed' : seed2})
st.cs.label_connected_components(adata, 'unspliced_mask_bis', min_area=100, area_threshold=600, max_area=2000, seed_layer='unspliced_labels', out_layer='unspliced_labels_bis')

del adata.layers['nuclear_mask_bis']
del adata.layers['unspliced']
del adata.layers['unspliced_bins']
del adata.layers['unspliced_labels']
del adata.layers['unspliced_scores']
del adata.layers['unspliced_scores_bis']
del adata.layers['unspliced_mask']
gc.collect() 

st.cs.segment_densities(adata, 'X', my_bin, k=3, distance_threshold=3, dk=3)
with open(dataset+'/images/X_bins.txt', "w") as txt_file:
    for line in adata.layers['X_bins']:
        line_str = [str(n) for n in line]
        txt_file.write(" ".join(line_str) + "\n")
 

st.cs.score_and_mask_pixels(adata, 'X', k=EM_kernel, method='EM+BP', certain_layer='unspliced_labels_bis',em_kwargs={'precision': 1e-06, 'max_iter':15000 , 'seed' : seed1}) 
st.cs.label_connected_components(adata, 'X', seed_layer='unspliced_labels_bis',min_area=500, area_threshold=2000, max_area=7000)

del adata.layers['unspliced_mask_bis']

unique, counts = numpy.unique(adata.layers['X_labels'], return_counts=True)
cells = dict(zip(unique, counts))
for x in range(0,adata.n_obs):
    for y in range(0,adata.n_vars):
        if (cells[adata.layers['X_labels'][x,y]] < 50) or (cells[adata.layers['X_labels'][x,y]] > 5000):
            adata.layers['X_labels'][x,y]=0
            
st.cs.expand_labels(adata, 'X',distance=5,max_area=5000)

with open(dataset+'/images/final_masks.txt', "w") as txt_file:
    for line in adata.layers['X_labels_expanded']:
        line_str = [str(n) for n in line]
        txt_file.write(" ".join(line_str) + "\n")

cell_adata = st.io.read_bgi(dataset+'/dge/whole_dataset.txt', segmentation_adata=adata, labels_layer='X_labels_expanded')

if sys.argv[4] == 'mouse':
    cell_adata.var["mt"] = cell_adata.var_names.str.startswith("mt-")
   

if sys.argv[4] == 'human':
    cell_adata.var["mt"] = cell_adata.var_names.str.startswith("MT-")
   
    
sc.pp.calculate_qc_metrics(cell_adata, qc_vars=["mt"], inplace=True)

cell_adata.uns['spatial']=''

sc.set_figure_params(facecolor="white", figsize=(8, 8))

with open(dataset+'/'+output+'/log', "w") as txt_file:
    txt_file.write('#cells before filtering: '+str(cell_adata.n_obs)+'\ngenes before filtering:'+ str(cell_adata.n_vars)+'\n')

print(f"#cells before filtering: {cell_adata.n_obs}")
print(f"#genes before filtering: {cell_adata.n_vars}")

fig, axs = plt.subplots(1, 2, figsize=(15, 4))
sns.histplot(cell_adata.obs["total_counts"][cell_adata.obs["total_counts"] < 5000],kde=False, bins=40,ax=axs[0],binrange =(0,2000))
sns.histplot(cell_adata.obs["n_genes_by_counts"][cell_adata.obs["n_genes_by_counts"] < 2000],kde=False, bins=60, ax=axs[1],binrange=(0,2000))
fig.savefig(dataset+'/'+output+'/histograms_RNA_segmented_pre_filtering.png') 

fig, axs = plt.subplots(1, 1, figsize=(8, 4))
sns.histplot(cell_adata.obs["area"][cell_adata.obs["area"]<3000],kde=False, bins=60,ax=axs)
fig.savefig(dataset+'/'+output+'/area_histogram_RNA_segmented_pre_filtering.png') 


sc.pp.filter_cells(cell_adata, min_counts=200)
sc.pp.filter_cells(cell_adata, max_counts=3000)
cell_adata = cell_adata[cell_adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.filter_genes(cell_adata, min_cells=10)


with open(dataset+'/'+output+'/log', "a") as txt_file:
    txt_file.write('#cells after filtering: '+str(cell_adata.n_obs)+'\n#genes after filtering: '+str(cell_adata.n_vars)+'\nmean counts_per_cell: '+str(cell_adata.obs["total_counts"].mean())+'\ntotal UMIs in cells: '+str(cell_adata.obs["total_counts"].sum())+'\n')

print(f"#cells after filtering: {cell_adata.n_obs}")
print(f"#genes after filtering: {cell_adata.n_vars}")
print(f'mean counts_per_cell: {cell_adata.obs["total_counts"].mean()}')
print(f'total UMIs in cells: {cell_adata.obs["total_counts"].sum()}')

fig, axs = plt.subplots(1, 2, figsize=(15, 4))
sns.histplot(cell_adata.obs["total_counts"][cell_adata.obs["total_counts"] < 5000],kde=False, bins=40,ax=axs[0],binrange =(0,2000))
sns.histplot(cell_adata.obs["n_genes_by_counts"][cell_adata.obs["n_genes_by_counts"] < 2000],kde=False, bins=60, ax=axs[1],binrange=(0,2000))
fig.savefig(dataset+'/'+output+'/histograms_RNA_segmented.png') 

fig, axs = plt.subplots(1, 1, figsize=(8, 4))
sns.histplot(cell_adata.obs["area"][cell_adata.obs["area"]<3000],kde=False, bins=60,ax=axs)
fig.savefig(dataset+'/'+output+'/area_histogram_RNA_segmented.png') 

outfile = dataset+'/'+output+'/segmented.h5ad'
cell_adata.write_h5ad(outfile)

sc.settings.figdir =dataset+'/'+output+'/'

sc.pp.normalize_total(cell_adata, inplace=True)
sc.pp.log1p(cell_adata)
sc.pp.highly_variable_genes(cell_adata, flavor="seurat", n_top_genes=2000)
#sc.pp.scale(cell_adata, max_value=10)
sc.pp.pca(cell_adata)
sc.pp.neighbors(cell_adata)
sc.tl.umap(cell_adata)
sc.tl.leiden(cell_adata, key_added="clusters", directed=False, n_iterations=2)
plt.rcParams["figure.figsize"] = (4, 4)
cell_adata.uns['spatial']
sc.pl.umap(cell_adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, show=False, save='umap.png')

plt.rcParams["figure.figsize"] = (8, 8)

spatial_plot=sc.pl.spatial(cell_adata, color=["total_counts", "n_genes_by_counts"], spot_size=30,save = 'spatial_RNA_segmented.png',show=False)

sc.pl.spatial(cell_adata, color=["clusters", "n_genes_by_counts"],spot_size=35,show=False, save='showspatial_w_clustersRNA_segmented.png')

if os.path.exists(dataset+'/dge/whole_dataset_shift.txt'):
    cell_adata_shuf = st.io.read_bgi(dataset+'/dge/whole_dataset_shift.txt', segmentation_adata=adata, labels_layer='X_labels_expanded')
    if sys.argv[4] == 'mouse':
        cell_adata_shuf.var["mt"] = cell_adata_shuf.var_names.str.startswith("mt-")
    if sys.argv[4] == 'human':
        cell_adata_shuf.var["mt"] = cell_adata_shuf.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(cell_adata_shuf, qc_vars=["mt"], inplace=True)
    cell_adata_shuf.uns['spatial']=''
    sc.pp.filter_cells(cell_adata_shuf, min_counts=200)
    sc.pp.filter_cells(cell_adata_shuf, max_counts=5000)
    #cell_adata_shuf = cell_adata_shuf[cell_adata_shuf.obs["pct_counts_mt"] < 20].copy()
    sc.pp.filter_genes(cell_adata_shuf, min_cells=10)
    outfile_shuf = dataset+'/'+output+'/segmented_shuf.h5ad'
    cell_adata_shuf.write_h5ad(outfile_shuf)

