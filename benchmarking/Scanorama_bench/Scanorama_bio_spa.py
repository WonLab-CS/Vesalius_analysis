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



#-----------------------------------------------------------------------------#
# Utilities 
#-----------------------------------------------------------------------------#
# yes I know this code is stupid. I just don't want to re write the other functions
def get_files(files,data_type, seed, query):
    seed_counts = [i  for i in files if re.search(seed,i) and 
                   re.search('gene_counts',i) and 
                   re.search(data_type,i)]
    seed_coord = [i  for i in files if re.search(seed,i) and 
                  re.search('spatial_coordinates',i) and
                  re.search(data_type,i)]
    query_counts = [i  for i in files if re.search(query,i) and 
                    re.search('gene_counts',i) and
                    re.search(data_type,i)]
    query_coord = [i  for i in files if re.search(query,i) and 
                   re.search('spatial_coordinates',i) and
                   re.search(data_type,i)]
    file_1 = [seed_counts[0], seed_coord[0]]
    file_2 = [query_counts[0], query_coord[0]]
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
    adata.obs['cell_labels'] = coord['cell_labels'].to_numpy()
    adata.obs['interactions'] = coord['interactions'].to_numpy()
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
    parser.add_argument('data_type', help='Data Type')
    parser.add_argument('input', help = 'Location of input data')
    parser.add_argument('output', help = 'Location of out data')
    parser.add_argument('seed_tag', help = 'Reference file tag')
    parser.add_argument('query_tag', help = 'Query file tag')
    args = parser.parse_args()
    output_path = args.output
    files = os.listdir(args.input)
    files = get_files(files,args.data_type, args.seed_tag, args.query_tag)
    seed = load_sim(args.input, files[0])
    query = load_sim(args.input, files[1])
    adats = [seed, query]
    sca.integrate_scanpy(adats)
    corrected = sca.correct_scanpy(adats, return_dimred = True)
    best_match = get_best_match(corrected)
    tag = f'Scanorama_aligned_{args.data_type}_{args.seed_tag}_{args.query_tag}.csv'
    export_coord(best_match, output_path, tag)
    return 0

if __name__ == "__main__":
    main()