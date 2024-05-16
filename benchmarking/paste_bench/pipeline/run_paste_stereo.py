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
#-----------------------------------------------------------------------------#
# load data
#-----------------------------------------------------------------------------#
input = '/common/wonklab/Stereo_seq/'
output_data = '/common/martinp4/output_data/'
adata = os.path.join(input, 'Mouse_embryo_all_stage.h5ad')
adata = sc.read_h5ad(adata)

seed = adata[adata.obs["timepoint"] == "E16.5"]

query = adata[adata.obs["timepoint"] == "E15.5"]
#-----------------------------------------------------------------------------#
# pairwise alignment
#-----------------------------------------------------------------------------#
print('Start pasting')
pi12 = pst.pairwise_align(seed, query)
# To visualize the alignment you can stack the slices 
# according to the alignment pi
slices, pis = [seed, query], [pi12]
new_slices = pst.stack_slices_pairwise(slices, pis)

new_seed = new_slices[0]
seed_spa = pd.DataFrame(new_seed.obsm['spatial'])
seed_spa = seed_spa.set_axis(new_seed.obs.index.values, axis = 0)
seed_spa = seed_spa.set_axis(['x','y'], axis = 1)
barcodes = pd.DataFrame(seed_spa.index.values)
barcodes = barcodes.set_axis(seed_spa.index.values, axis = 0)
barcodes = barcodes.set_axis(['barcodes'], axis = 1)
export_seed = pd.concat([barcodes, seed_spa, new_seed.obs['annotation']], axis = 1)
export_file = os.join(output_data,"stereo_seed_paste.csv")
export_seed.to_csv(export_file)


new_query = new_slices[1]
query_spa = pd.DataFrame(new_query.obsm['spatial'])
query_spa = query_spa.set_axis(new_query.obs.index.values, axis = 0)
query_spa = query_spa.set_axis(['x','y'], axis = 1)
barcodes = pd.DataFrame(query_spa.index.values)
barcodes = barcodes.set_axis(query_spa.index.values, axis = 0)
barcodes = barcodes.set_axis(['barcodes'], axis = 1)
export_query = pd.concat([barcodes, query_spa, new_query.obs['annotation']], axis = 1)
export_file = os.join(output_data,"stereo_query_paste.csv")
export_query.to_csv(export_file)


# sq.pl.spatial_scatter(new_seed, shape = None, color = "cells", size = 0.5, save = '/common/martinp4/paste_aligned_FISH_seed.pdf')
# sq.pl.spatial_scatter(new_query, shape = None, color = "cells", size = 0.5,  save = '/common/martinp4/paste_aligned_FISH_query.pdf')
