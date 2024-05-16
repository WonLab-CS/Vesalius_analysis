#-----------------------------------------------------------------------------#
# Imports
#-----------------------------------------------------------------------------#
import numpy as np
import scanpy as sc
import tangram as tg
import os
import pandas as pd
import argparse
import re


os.environ['NUMEXPR_MAX_THREADS'] = '1'
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
    adata.obs["cell_labels"] = coord['cell_labels'].to_numpy()
    return adata

def  max_idx(df):
    return df.idxmax()

def export_coord(ad_map, query, seed, output_path, tag):
    locs = np.argmax(ad_map.X, axis = 1)
    spa = pd.DataFrame(seed.obsm['spatial'][locs])
    spa = spa.set_axis(query.obs.index.values, axis = 0)
    spa = spa.set_axis(['x','y'], axis = 1)
    barcodes = pd.DataFrame(query.obs.index.values)
    barcodes = barcodes.set_axis(query.obs.index.values, axis = 0)
    barcodes = barcodes.set_axis(['barcodes'], axis = 1)
    best_match = query.obs['cell_labels']
    best_match = best_match.set_axis(query.obs.index.values,axis =0)
    export_adata = pd.concat([barcodes, spa, best_match], axis = 1)
    export_file = os.path.join(output_path, tag)
    export_adata.to_csv(export_file)



#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='Tangram benchmarking on synthetic spatial data')
    parser.add_argument('slurm_id', help='Slurm ID from array sub')
    parser.add_argument('type', help = 'type of synthetic data requested')
    args = parser.parse_args()
    output_path = '/common/martinp4/benchmarking_out/Tangram/report/'
    files = os.listdir("/common/wonklab/synthetic_spatial")
    files = get_files(files, args.slurm_id, args.type)
    seed = load_sim("/common/wonklab/synthetic_spatial", files[0])
    query = load_sim("/common/wonklab/synthetic_spatial", files[1])
    tg.pp_adatas(query, seed, genes=None)
    ad_map = tg.map_cells_to_space(query, seed, density_prior = 'rna_count_based',num_epochs = 1000, mode = "cells", device = 'cpu')
    #tg.project_cell_annotations(ad_map, seed, annotation = 'cell_labels')
    tag = f'Tangram_aligned_{re.sub(".csv","",files[0][0])}_{re.sub(".csv","",files[1][0])}.csv'
    tag = re.sub('spatial_territories_gene_counts_','',tag)
    export_coord(ad_map,query, seed, output_path, tag)
    return 0

if __name__ == "__main__":
    main()