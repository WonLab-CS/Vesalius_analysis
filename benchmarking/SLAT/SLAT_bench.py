#-----------------------------------------------------------------------------#
# Imports
#-----------------------------------------------------------------------------#
import os
import argparse
import re
import torch
import pandas as pd
import scanpy as sc
import numpy as np
import pandas as pd
import scSLAT
from scSLAT.model import Cal_Spatial_Net, load_anndatas, run_SLAT, spatial_match


device = torch.device('cpu')
torch.set_num_threads(1)

#-----------------------------------------------------------------------------#
# Utilities 
#-----------------------------------------------------------------------------#
def get_files(files, id, sim_type, n_samp: int =12):
    counts = [i  for i in files if re.search(sim_type,i) and re.search('gene_counts',i)]
    counts.sort()
    coord = [i  for i in files if re.search(sim_type,i) and re.search('spatial_coordinates',i)]
    coord.sort()
    vec_1 = np.tile(np.array(range(n_samp)), n_samp)
    vec_2 = np.repeat(np.array(range(n_samp)), n_samp)
    file_1 = [counts[vec_1[int(id) - 1]], coord[vec_1[int(id) - 1]]]
    file_2 = [counts[vec_2[int(id) - 1]], coord[vec_2[int(id) - 1]]]
    return file_1, file_2

def load_sim(input_path, file):
    counts = os.path.join(input_path, file[0])
    counts = pd.read_csv(counts, index_col=0).iloc[:, 1:]
    coord = os.path.join(input_path, file[1])
    coord = pd.read_csv(coord, index_col = 0)
    adata = sc.AnnData(counts.T, dtype = np.float32)
    adata.var_names_make_unique()
    coord_df = coord.loc[:, ['x', 'y']]
    adata.obsm["spatial"] = coord_df.to_numpy()
    adata.obs["Territory"] = coord['Territory'].to_numpy()
    adata.obs['cell_labels'] = coord['cell_labels'].to_numpy()
    return adata

def export_match(seed, query, matched, output_path, tag):
    spa = seed[['x','y']]
    spa = spa = spa.set_axis(query.index.values, axis = 0)
    barcodes = pd.DataFrame(query.index.values)
    barcodes = barcodes.set_axis(query.index.values, axis = 0)
    barcodes = barcodes.set_axis(['barcodes'], axis = 1)
    best_match = query['cell_labels'][matched[1]]
    best_match = best_match.set_axis(query.index.values,axis =0)
    export_adata = pd.concat([barcodes, spa, best_match], axis = 1)
    export_file = os.path.join(output_path, tag)
    export_adata.to_csv(export_file)

#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='SLAT benchmarking on synthetic spatial data')
    parser.add_argument('slurm_id', help='Slurm ID from array sub')
    parser.add_argument('type', help = 'type of synthetic data requested')
    args = parser.parse_args()
    output_path = '/common/martinp4/benchmarking_out/SLAT/report'
    files = os.listdir("/common/wonklab/synthetic_spatial")
    files = get_files(files, args.slurm_id, args.type)
    seed = load_sim("/common/wonklab/synthetic_spatial", files[0])
    query = load_sim("/common/wonklab/synthetic_spatial", files[1])
    Cal_Spatial_Net(seed, k_cutoff=10, model='KNN')
    Cal_Spatial_Net(query, k_cutoff=10, model='KNN')
    edges, features = load_anndatas([seed, query], feature='DPCA')
    embd0, embd1, time = run_SLAT(features, edges)
    best, index, distance = spatial_match(features, adatas=[seed,query], reorder=False)
    adata1_df = pd.DataFrame({'index': range(embd0.shape[0]),
                        'x': seed.obsm['spatial'][:,0],
                        'y': seed.obsm['spatial'][:,1],
                        'cell_labels': seed.obs['cell_labels']})
    adata2_df = pd.DataFrame({'index': range(embd1.shape[0]),
                        'x': query.obsm['spatial'][:,0],
                        'y': query.obsm['spatial'][:,1],
                        'cell_labels': query.obs['cell_labels']})
    matching = np.array([range(index.shape[0]), best])
    tag = f'SLAT_aligned_{re.sub(".csv","",files[0][0])}_{re.sub(".csv","",files[1][0])}.csv'
    tag = re.sub('spatial_territories_gene_counts_','',tag)
    export_match(adata1_df,adata2_df, matching, output_path, tag)
    return 0

if __name__ == "__main__":
    main()