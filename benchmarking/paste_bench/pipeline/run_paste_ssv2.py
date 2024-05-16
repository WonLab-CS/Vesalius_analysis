import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scanpy as sc
import squidpy as sq
import paste as pst
import os
import pandas as pd



#-----------------------------------------------------------------------------#
# load data sets - same as with Vesalius
#-----------------------------------------------------------------------------#
input = '/common/wonklab/SSv2/'

seed_counts = os.path.join(input, 'Puck_200115_08.digital_expression.txt.gz')
seed_counts = pd.read_csv(seed_counts, sep='\t', index_col=0)
seed_coord = os.path.join(input, 'Puck_200115_08_bead_locations.csv')
seed_coord = pd.read_csv(seed_coord, skiprows = 1, header=None)
seed_coord.columns = ['barcodes','x','y']

seed = sc.AnnData(seed_counts.T, dtype = np.float32)
seed.var_names_make_unique()
coord_df = seed_coord.loc[:, ['x', 'y']]
seed.obsm["spatial"] = coord_df.to_numpy()
sc.pp.calculate_qc_metrics(seed, inplace=True)
sc.pp.highly_variable_genes(seed, flavor="seurat_v3", n_top_genes=2000)
sc.pp.normalize_total(seed, target_sum=1e4)
sc.pp.log1p(seed)
sc.pp.pca(seed)
sc.pp.neighbors(seed)
sc.tl.umap(seed)
sc.tl.leiden(seed, key_added="clusters", resolution = 0.8)


query_counts = os.path.join(input, 'Puck_190921_21.digital_expression.txt.gz')
query_counts = pd.read_csv(query_counts, sep ='\t', skiprows =3, index_col =0)
query_coord = os.path.join(input, 'Puck_190921_21_bead_locations.csv')
query_coord = pd.read_csv(query_coord, skiprows = 1, header=None)
query_coord.columns = ['barcodes','x','y']

query = sc.AnnData(query_counts.T, dtype = np.float32)
query.var_names_make_unique()
coord_df = query_coord.loc[:, ['x', 'y']]
query.obsm["spatial"] = coord_df.to_numpy()
sc.pp.calculate_qc_metrics(query, inplace=True)
sc.pp.highly_variable_genes(query, flavor="seurat_v3", n_top_genes=2000)
sc.pp.normalize_total(query, target_sum=1e4)
sc.pp.log1p(query)
sc.pp.pca(query)
sc.pp.neighbors(query)
sc.tl.umap(query)
sc.tl.leiden(query, key_added="clusters", resolution = 0.3)

#-----------------------------------------------------------------------------#
# Plots 
#-----------------------------------------------------------------------------#
# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# import matplotlib.colors as mcolors

# custom_palette = ["#E69F00",
#       "#56B4E9",
#       "#009E73",
#       "#F0E442",
#       "#0072B2",
#       "#D55E00",
#       "#CC79A7",
#       "#999999"]


# # Convert hex values to RGB values
# rgb_colors = [mcolors.hex2color(hex_color) for hex_color in custom_palette]

# # Create a dictionary mapping values to colors for the gradient
# n = len(rgb_colors) - 1
# color_dict = {
#     'red': [(i/n, rgb_colors[i][0], rgb_colors[i][0]) for i in range(n+1)],
#     'green': [(i/n, rgb_colors[i][1], rgb_colors[i][1]) for i in range(n+1)],
#     'blue': [(i/n, rgb_colors[i][2], rgb_colors[i][2]) for i in range(n+1)]
# }

# # Create the LinearSegmentedColormap
# categories = list(new_seed.obs["clusters"].cat.categories)
# gradient_cmap = LinearSegmentedColormap('gradient_cmap', segmentdata=color_dict, N=len(categories))
# gradient_colors = gradient_cmap(np.linspace(0, 1, len(categories)))
# sq.pl.spatial_scatter(new_seed,
#                       shape = None,
#                       color = "clusters", 
#                       size = 0.5, 
#                       palette = gradient_colors,
#                       save = '/common/martinp4/output_plots/paste_ssv2_seed.pdf')

# categories = list(query.obs["clusters"].cat.categories)
# cmap = plt.cm.colors.ListedColormap(custom_palette, N = len(categories))
# sq.pl.spatial_scatter(query, shape = None, color = "clusters", size = 0.5,  save = '/common/martinp4/output_plots/paste_ssv2_query.pdf')

#-----------------------------------------------------------------------------#
# pairwise alignment
#-----------------------------------------------------------------------------#
pi12 = pst.pairwise_align(seed, query)
# To visualize the alignment you can stack the slices 
# according to the alignment pi
slices, pis = [seed, query], [pi12]
new_slices = pst.stack_slices_pairwise(slices, pis)

new_seed = new_slices[0]
sc.pp.calculate_qc_metrics(new_seed, inplace=True)
sc.pp.highly_variable_genes(new_seed, flavor="seurat_v3", n_top_genes=2000)
sc.pp.normalize_total(new_seed, target_sum=1e4)
sc.pp.log1p(new_seed)
sc.pp.pca(new_seed)
sc.pp.neighbors(new_seed)
sc.tl.umap(new_seed)
sc.tl.leiden(new_seed, key_added="clusters", resolution = 0.8)

new_query = new_slices[1]
sc.pp.calculate_qc_metrics(new_query, inplace=True)
sc.pp.highly_variable_genes(new_query, flavor="seurat_v3", n_top_genes=2000)
sc.pp.normalize_total(new_query, target_sum=1e4)
sc.pp.log1p(new_query)
sc.pp.pca(new_query)
sc.pp.neighbors(new_query)
sc.tl.umap(new_query)
sc.tl.leiden(new_query, key_added="clusters", resolution = 0.3)

sq.pl.spatial_scatter(new_seed, shape = None, color = "clusters", size = 0.5, save = '/common/martinp4/paste_aligned_ssv2_seed.pdf')
sq.pl.spatial_scatter(new_query, shape = None, color = "clusters", size = 0.5,  save = '/common/martinp4/paste_aligned_ssv2_query.pdf')

