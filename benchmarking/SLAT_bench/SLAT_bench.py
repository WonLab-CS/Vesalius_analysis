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
def get_files(files, id, sim_type):
    sim_type = re.escape(sim_type) + '_cells'
    counts = [i  for i in files if re.search(sim_type,i) and re.search('gene_counts',i)]
    counts.sort()
    coord = [i  for i in files if re.search(sim_type,i) and re.search('spatial_coordinates',i)]
    coord.sort()
    n_samp = len(counts)
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
    adata = sc.AnnData(counts.T)
    adata.var_names_make_unique()
    coord_df = coord.loc[:, ['x', 'y']]
    adata.obsm["spatial"] = coord_df.to_numpy()
    adata.obs["Territory"] = coord['Territory'].to_numpy()
    adata.obs['cell_labels'] = coord['cell_labels'].to_numpy()
    adata.obs['interactions'] = coord['interactions'].to_numpy()
    return adata

def export_match(seed, query, best, output_path, tag):
    spa = pd.DataFrame(seed[['x','y']])
    query = query.iloc[best]
    spa = spa.set_axis(query.index.values, axis = 0)
    export_adata = pd.concat([query['barcodes'],spa,query[['cell_labels','interactions']]], axis = 1)
    export_file = os.path.join(output_path, tag)
    export_adata.to_csv(export_file)
    return 0

#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='SLAT benchmarking on synthetic spatial data')
    parser.add_argument('--slurm_id', help='Slurm ID from array sub')
    parser.add_argument('--type', help = 'type of synthetic data requested')
    parser.add_argument('--input', help = 'Location of input data')
    parser.add_argument('--output', help = 'Location of out data')
    args = parser.parse_args()
    output_path = args.output
    files = os.listdir(args.input)
    files = get_files(files, args.slurm_id, args.type)
    if files[0][0] == files[1][0]:
        return 0
    
    seed = load_sim(args.input, files[0])
    query = load_sim(args.input, files[1])
    Cal_Spatial_Net(seed, k_cutoff = 6, model = 'KNN')
    Cal_Spatial_Net(query, k_cutoff = 6, model = 'KNN')
    edges, features = load_anndatas([seed, query], feature='DPCA', check_order=False)
    embd0, embd1, time = run_SLAT(features, edges,epochs = 25,LGCN_layer=5)
    best, index, distance = spatial_match([embd1,embd0], adatas=[query,seed], reorder=False)
    seed = pd.DataFrame({'barcodes': seed.obs.index.values,
                        'x': seed.obsm['spatial'][:,0],
                        'y': seed.obsm['spatial'][:,1],
                        'cell_labels': seed.obs['cell_labels'],
                        'interactions' : seed.obs['interactions']})
    query = pd.DataFrame({'barcodes': query.obs.index.values,
                        'x': query.obsm['spatial'][:,0],
                        'y': query.obsm['spatial'][:,1],
                        'cell_labels': query.obs['cell_labels'],
                        'interactions' : query.obs['interactions']})
    tag = f'SLAT_aligned_{args.type}_{re.sub(".csv","",files[0][0])}_{re.sub(".csv","",files[1][0])}.csv'
    tag = re.sub('gene_counts_','',tag)
    if 'computational_performance' in tag:
        return 0
    else :
        export_match(seed,query, best, output_path, tag)
        return 0

if __name__ == "__main__":
    main()
    