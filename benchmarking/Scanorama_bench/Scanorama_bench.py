#-----------------------------------------------------------------------------#
# Imports
#-----------------------------------------------------------------------------#
import numpy as np
import scanpy as sc
import scanorama as sca
import os
import pandas as pd
import argparse
import re
from scipy.spatial.distance import cdist


os.environ['NUMEXPR_MAX_THREADS'] = '1'
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
    adata.obs["cell_labels"] = coord['cell_labels'].to_numpy()
    adata.obs["interactions"] = coord['interactions'].to_numpy()
    sc.pp.log1p(adata)
    return adata



def get_best_match(corrected):
    seed = corrected[0]
    query = corrected[1]
    seed_ls = seed.obsm['X_scanorama']
    query_ls = query.obsm['X_scanorama']
    distances = cdist(seed_ls, query_ls, metric='euclidean')
    ids = np.argmin(distances, axis=0)
    barcodes = pd.DataFrame(query.obs.index.values)
    barcodes = barcodes.set_axis(query.obs.index.values, axis = 0)
    barcodes = barcodes.set_axis(['barcodes'], axis = 1)
    spa = pd.DataFrame(seed.obsm['spatial'][ids])
    spa = spa.set_axis(query.obs.index.values, axis = 0)
    spa = spa.set_axis(['x','y'], axis = 1)
    meta = query.obs
    best_match = pd.concat([barcodes, spa, meta], axis = 1)
    return best_match



def export_coord(best_match, output_path, tag):
    export_file = os.path.join(output_path, tag)
    best_match.to_csv(export_file)



#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='Scanorama benchmarking on synthetic spatial data')
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
    adats = [seed, query]
    sca.integrate_scanpy(adats)
    best_match = get_best_match(adats)
    tag = f'Scanorama_aligned_{args.type}_{re.sub(".csv","",files[0][0])}_{re.sub(".csv","",files[1][0])}.csv'
    tag = re.sub('gene_counts_','',tag)
    if 'computational_performance' in tag:
        return 0
    else :
        export_coord(best_match, output_path, tag)
        return 0

if __name__ == "__main__":
    main()