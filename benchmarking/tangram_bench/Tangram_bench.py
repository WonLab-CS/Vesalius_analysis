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
    adata = sc.AnnData(counts.T, dtype = 'float32')
    adata.var_names_make_unique()
    coord_df = coord.loc[:, ['x', 'y']]
    adata.obsm["spatial"] = coord_df.to_numpy()
    adata.obs['cell_labels'] = coord['cell_labels'].to_numpy()
    adata.obs['interactions'] = coord['interactions'].to_numpy()
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000)
    var_genes = adata.var[adata.var['highly_variable']].index.to_numpy()
    return adata, var_genes

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
    interaction = query.obs['interactions']
    interaction = interaction.set_axis(query.obs.index.values,axis =0)
    export_adata = pd.concat([barcodes, spa, best_match, interaction], axis = 1)
    export_file = os.path.join(output_path, tag)
    export_adata.to_csv(export_file)
    return 0



#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='Tangram benchmarking on synthetic spatial data')
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
    seed, seed_genes = load_sim(args.input, files[0])
    query, query_genes = load_sim(args.input, files[1])
    tg.pp_adatas(query, seed, genes = list(set(seed_genes) & set(query_genes)))
    ad_map = tg.map_cells_to_space(query, seed, density_prior = 'rna_count_based',num_epochs = 1000, mode = "cells", device = 'cpu')
    tag = f'Tangram_aligned_{args.type}_{re.sub(".csv","",files[0][0])}_{re.sub(".csv","",files[1][0])}.csv'
    tag = re.sub('gene_counts_','',tag)
    export_coord(ad_map,query, seed, output_path, tag)
    return 0

if __name__ == "__main__":
    main()