from threadpoolctl import threadpool_limits
import seaborn as sns
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import sys
import spateo as st
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

EM_kernel= int(sys.argv[6])

if sys.argv[5] == "apex":
    nuc_markers_file = open(sys.argv[3]+'/reference/'+sys.argv[4]+'/'+sys.argv[4]+'.MIR_HG_and_APEX_nuc_markers')
else:
    nuc_markers_file = open(sys.argv[1]+'/nuclear_markers')
 
    
if sys.argv[7] == 'NO':
    EM_thresh=None
else:
    EM_thresh=float(sys.argv[7])
      
nuc_markers = nuc_markers_file.read().strip().split("\n")
adata = st.io.read_bgi_agg(sys.argv[1]+'/dge/'+sys.argv[2]+'.spateo.txt',gene_agg={'nuclear': nuc_markers} )

### nuclear markers

with threadpool_limits(limits=4):
    st.cs.segment_densities(adata, 'nuclear', 15, k=5, dk=3, distance_threshold=3, background=False)
    st.pl.contours(adata, 'nuclear_bins', scale=0.15)
    st.pl.imshow(adata, 'nuclear_bins', labels=True, save_show_or_return='save', save_kwargs={'path':  sys.argv[1]+'/images/', 'prefix': sys.argv[2]+'_contours', 'dpi': None, 'ext': 'png', 'transparent': False, 'close': True, 'verbose': True})
    st.cs.score_and_mask_pixels(adata, 'nuclear',  threshold=EM_thresh, k=EM_kernel, method='EM+BP',em_kwargs={'max_iter': 15000})
    st.cs.label_connected_components(adata,'nuclear',area_threshold=2000, max_area=4000, min_area=300)
 #   if sys.argv[4] == 'human':
    st.segmentation.utils.filter_cell_labels_by_area(adata,'nuclear',150)
    st.pl.imshow(adata, 'nuclear_labels', labels=True, save_show_or_return='save', save_kwargs={'path':  sys.argv[1]+'/images/', 'prefix': sys.argv[2]+'_nuc_masks', 'dpi': None, 'ext': 'png', 'transparent': False, 'close': True, 'verbose': True})

### unspliced
with threadpool_limits(limits=4):
    st.cs.segment_densities(adata, 'unspliced', 15, k=5, dk=3, distance_threshold=3, background=False)
    st.pl.contours(adata, 'unspliced_bins', scale=0.15)
    st.pl.imshow(adata, 'unspliced_bins', labels=True, save_show_or_return='save', save_kwargs={'path':  sys.argv[1]+'/images/', 'prefix': sys.argv[2]+'_unspliced_contours', 'dpi': None, 'ext': 'png', 'transparent': False, 'close': True, 'verbose': True})
    st.cs.score_and_mask_pixels( adata, 'unspliced',  threshold=EM_thresh, k=EM_kernel, method='EM+BP', certain_layer='nuclear_labels',em_kwargs={'max_iter': 15000})
    st.cs.label_connected_components(adata, 'unspliced',area_threshold=2000, min_area=300, max_area=4000, seed_layer='nuclear_labels')
#    if sys.argv[4] == 'human':
    st.segmentation.utils.filter_cell_labels_by_area(adata,'unspliced',150)
    st.pl.imshow(adata, 'unspliced_labels', labels=True, save_show_or_return='save', save_kwargs={'path':  sys.argv[1]+'/images/', 'prefix': sys.argv[2]+'_unspliced_masks', 'dpi': None, 'ext': 'png', 'transparent': False, 'close': True, 'verbose': True})

### cytoplasm

with threadpool_limits(limits=4):
    st.cs.segment_densities(adata, 'X', 15, k=5, distance_threshold=3, dk=3, background=False)
    st.pl.contours(adata, 'X_bins', scale=0.15)
    st.cs.score_and_mask_pixels(adata, 'X', threshold=EM_thresh, k=EM_kernel+2, method='EM+BP', certain_layer='unspliced_labels',em_kwargs={'max_iter': 15000}) 
    st.cs.label_connected_components(adata, 'X', seed_layer='unspliced_labels',min_area=4000, area_threshold=8000, max_area=20000)
 #   if sys.argv[4] == 'human':
    st.segmentation.utils.filter_cell_labels_by_area(adata,'X',200)
    st.pl.imshow(adata, 'X_labels', labels=True, save_show_or_return='save', save_kwargs={'path':  sys.argv[1]+'/images/', 'prefix': sys.argv[2]+'_final_masks', 'dpi': None, 'ext': 'png', 'transparent': False, 'close': True, 'verbose': True})

cell_adata = st.io.read_bgi(sys.argv[1]+'/dge/'+sys.argv[2]+'.spateo.txt', segmentation_adata=adata, labels_layer='X_labels')
cell_adata=cell_adata[cell_adata.obs['area']> 50]
del cell_adata.obsm['contour']

if sys.argv[4] == 'human':
    cell_adata.var["mt"] = cell_adata.var_names.str.startswith("MT-")

if sys.argv[4] == 'mouse':
    cell_adata.var["mt"] = cell_adata.var_names.str.startswith("Mt-")
    
sc.pp.calculate_qc_metrics(cell_adata, qc_vars=["mt"], inplace=True)

fig, axs = plt.subplots(1, 2, figsize=(15, 4))
sns.histplot(cell_adata.obs["total_counts"][cell_adata.obs["total_counts"] < 5000],kde=False, bins=40,ax=axs[0],binrange =(0,5000))
sns.histplot(cell_adata.obs["n_genes_by_counts"][cell_adata.obs["n_genes_by_counts"] < 2000],kde=False, bins=60, ax=axs[1],binrange=(0,2000))
fig.savefig(sys.argv[1]+'/images/'+sys.argv[2]+'_histograms_RNA_segmented_pre_filtering.png') 

fig, axs = plt.subplots(1, 1, figsize=(8, 4))
sns.histplot(cell_adata.obs["area"][cell_adata.obs["area"]<3000],kde=False, bins=60)
fig.savefig(sys.argv[1]+'/images/'+sys.argv[2]+'_area_histogram_RNA_segmented_pre_filtering.png') 
cat = sys.argv[2]+'_'+cell_adata.obs_names
cell_adata.obs_names = cat

outfile = sys.argv[1]+'/dge/'+sys.argv[2]+'.h5ad'
cell_adata.write_h5ad(outfile)




