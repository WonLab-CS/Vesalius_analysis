#-----------------------------------------------------------------------------#
# Imports
#-----------------------------------------------------------------------------#
import numpy as np
import scanpy as sc
import squidpy as sq
import paste as pst
import os
import pandas as pd
import argparse
import re
from scipy.spatial import KDTree



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

def export_coord(adata, output_path, tag):
    spa = pd.DataFrame(adata.obsm['spatial'])
    spa = spa.set_axis(adata.obs.index.values, axis = 0)
    spa = spa.set_axis(['x','y'], axis = 1)
    barcodes = pd.DataFrame(spa.index.values)
    barcodes = barcodes.set_axis(spa.index.values, axis = 0)
    barcodes = barcodes.set_axis(['barcodes'], axis = 1)
    export_adata = pd.concat([barcodes, spa, adata.obs[['cell_labels','interactions']]], axis = 1)
    export_file = os.path.join(output_path, tag)
    export_adata.to_csv(export_file)


#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='PASTE benchmarking on synthetic spatial data')
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
    pi12 = pst.pairwise_align(seed, query, numItermax = 10000, use_gpu = False) 
    slices, pis = [seed, query], [pi12]
    new_slices = pst.stack_slices_pairwise(slices, pis)
    new_query = new_slices[1]
    tag = f'PASTE_aligned_{args.type}_{re.sub(".csv","",files[0][0])}_{re.sub(".csv","",files[1][0])}.csv'
    tag = re.sub('gene_counts_','',tag)
    if 'computational_performance' in tag:
        return 0
    else :
        export_coord(new_query, output_path, tag)
        return 0

if __name__ == "__main__":
    main()