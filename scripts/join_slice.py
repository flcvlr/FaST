import decoupler as dc
import anndata as ad
import sys
import spateo as st

with open(sys.argv[1]+'/tiles_info') as tiles_file:
    line = tiles_file.readline()
    tile_data = line.strip().split(",")
    all_adata = st.read_h5ad(sys.argv[1]+'/dge/K_'+sys.argv[3]+'/'+tile_data[0]+'.h5ad')
    for i in range (0, len(all_adata.obsm['spatial'])):
        all_adata.obsm['spatial'][i][0] += int(tile_data[1])
        all_adata.obsm['spatial'][i][1] += int(tile_data[2])

    for line in tiles_file:
        tile_data = line.strip().split(",")
        adata = st.read_h5ad(sys.argv[1]+'/dge/K_'+sys.argv[3]+'/'+tile_data[0]+'.h5ad')
        for i in range (0, len(adata.obsm['spatial'])):
            adata.obsm['spatial'][i][0] += int(tile_data[1])
            adata.obsm['spatial'][i][1] += int(tile_data[2])
        adata_list = [all_adata,adata]
        all_adata = ad.concat(adata_list, join = 'outer')  
        
outfile = sys.argv[1]+'/RNA_segmented.h5ad'
all_adata.write_h5ad(outfile)

#### analyse dataset

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
outdir = sys.argv[1]+'/images/K_'+sys.argv[3]+'/'
sc.settings.figdir = outdir
sc.set_figure_params(facecolor="white", figsize=(8, 8))

if sys.argv[2] == 'human':
    all_adata.var["mt"] = all_adata.var_names.str.startswith("MT-")

if sys.argv[2] == 'mouse':
    all_adata.var["mt"] = all_adata.var_names.str.startswith("Mt-")

print(f"#cells before filtering: {all_adata.n_obs}")
print(f"#genes before filtering: {all_adata.n_vars}")

sc.pp.calculate_qc_metrics(all_adata, qc_vars=["mt"], inplace=True)
fig, axs = plt.subplots(1, 2, figsize=(15, 4))
sns.histplot(all_adata.obs["total_counts"][all_adata.obs["total_counts"] < 5000],kde=False, bins=40,ax=axs[0],binrange =(0,2000))
sns.histplot(all_adata.obs["n_genes_by_counts"][all_adata.obs["n_genes_by_counts"] < 2000],kde=False, bins=60, ax=axs[1],binrange=(0,2000))
fig.savefig(outdir+'histograms_RNA_segmented_pre_filtering.png') 

fig, axs = plt.subplots(1, 1, figsize=(8, 4))
sns.histplot(all_adata.obs["area"][all_adata.obs["area"]<3000],kde=False, bins=60)
fig.savefig(outdir+'area_histogram_RNA_segmented_pre_filtering.png') 

sc.pp.filter_cells(all_adata, min_counts=200)
sc.pp.filter_cells(all_adata, max_counts=2000)
all_adata = all_adata[all_adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.filter_genes(all_adata, min_cells=10)
print(f"#cells after filtering: {all_adata.n_obs}")
print(f"#genes after filtering: {all_adata.n_vars}")

print(f'mean counts_per_cell: {all_adata.obs["total_counts"].mean()}')
print(f'total UMIs in cells: {all_adata.obs["total_counts"].sum()}')

fig, axs = plt.subplots(1, 2, figsize=(15, 4))
sns.histplot(all_adata.obs["total_counts"][all_adata.obs["total_counts"] < 5000],kde=False, bins=40,ax=axs[0],binrange =(0,2000))
sns.histplot(all_adata.obs["n_genes_by_counts"][all_adata.obs["n_genes_by_counts"] < 2000],kde=False, bins=60, ax=axs[1],binrange=(0,2000))
fig.savefig(outdir+'histograms_RNA_segmented.png') 

fig, axs = plt.subplots(1, 1, figsize=(8, 4))
sns.histplot(all_adata.obs["area"][all_adata.obs["area"]<3000],kde=False, bins=60)
fig.savefig(outdir+'area_histogram_RNA_segmented.png') 

import numpy

def overlap(gene1, gene2,adata, dataset):
    overlap_data=[numpy.matrix.sum(adata[:,gene1].X.todense() > 0), numpy.matrix.sum(adata[:,gene2].X.todense() > 0),numpy.matrix.sum((adata[:,gene1].X.todense() > 0) & (adata[:,gene2].X.todense() > 0))]
    cat=[gene1,gene2,'Double_Positive']
    plt.rcParams['font.size'] = 18
    fig, ax=plt.subplots()
    ax.bar(cat,overlap_data, color=['grey','blue','black'])
    ax.set_title('Number of cells with at least one count in '+dataset)
    perc=str(round((overlap_data[2]/(overlap_data[0]+overlap_data[1]))*100,2))
    ax.annotate('BOTH ='+perc+'%',xy=(2,5*overlap_data[2]), horizontalalignment='center',size='large')
    fig.savefig(outdir+gene1+'_'+gene2+'_overlap.png')

if (sys.argv[2] == 'human'):
    overlap('CD3E','CD19',all_adata, 'FaST')
    overlap('EPCAM','ACTA2',all_adata,'FaST')
    overlap('CCR7','CCL19',all_adata,'FaST')
    overlap('CD3E','KRT6A',all_adata,'FaST')

sc.pp.normalize_total(all_adata, inplace=True)
sc.pp.log1p(all_adata)
sc.pp.highly_variable_genes(all_adata, flavor="seurat", n_top_genes=2000)
#sc.pp.scale(adata, max_value=10)
sc.pp.pca(all_adata)
sc.pp.neighbors(all_adata)
sc.tl.umap(all_adata)
sc.tl.leiden(all_adata, key_added="clusters", directed=False, n_iterations=2)
plt.rcParams["figure.figsize"] = (4, 4)
umap_plot=sc.pl.umap(all_adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save='umap_RNA_segmented.png', show=False)

plt.rcParams["figure.figsize"] = (8, 8)
spatial_plot=sc.pl.spatial(all_adata, color=["total_counts", "n_genes_by_counts"], spot_size=40,save = 'spatial_RNA_segmented.png',show=False)

spatial_clusters = sc.pl.spatial(all_adata, color=["clusters", "n_genes_by_counts"],spot_size=50, save = 'spatial_w_clustersRNA_segmented.png',show=False)

sc.tl.rank_genes_groups(all_adata, "clusters", method="t-test")
heatmap_plot = sc.pl.rank_genes_groups_heatmap(all_adata, n_genes=10, groupby="clusters", save = 'heatmap_RNA_segmented.png',show=False)
#KRT6A_plot = sc.pl.spatial(all_adata, color=["clusters", "KRT6A"],spot_size=40, save = 'KRT6A_RNA_segmented.png',show=False)

all_adata.write_h5ad(outdir+'merged.h5ad')
