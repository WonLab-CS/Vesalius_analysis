#-----------------------------------------------------------------------------#
# Imports
#-----------------------------------------------------------------------------#

import pandas as pd
import argparse
import re
import numpy as np
import csv
import os


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
    counts = pd.read_csv(counts, index_col=0)
    coord = os.path.join(input_path, file[1])
    coord = pd.read_csv(coord, index_col = 0)
    barcodes = pd.DataFrame(coord.index.values)
    barcodes = barcodes.set_axis(coord.index.values, axis = 0)
    barcodes = barcodes.set_axis(['barcodes'], axis = 1)
    labels = pd.concat([barcodes,coord['Territory']], axis = 1)
    coord = pd.concat([barcodes,coord[['x','y']]], axis = 1)
    return counts, labels, coord

def synth_to_sc(counts, labels, data_type, output_path, slurm_id):
    counts = counts.rename(columns = {'genes': 'GENES'})
    labels = labels.rename(columns={'barcodes':'SpotID', 'Territory' : 'CellType'})
    labels['CellType'] = np.repeat('celltype_', labels.shape[0]) + labels['CellType'].astype(str)
    scRNA = os.path.join(output_path, f'scRNA_{str(slurm_id)}_{data_type}.txt')
    scLabels = os.path.join(output_path, f'scLabels_{str(slurm_id)}_{data_type}.txt')
    counts.to_csv(scRNA, index = False, sep = "\t", quoting = csv.QUOTE_NONE)
    labels.to_csv(scLabels, index = False, sep = "\t", quoting = csv.QUOTE_NONE)
    
def synth_to_spa(counts, labels, coord, data_type, output_path, slurm_id):
    counts = counts.rename(columns = {'genes': 'GENES'})
    coord = coord.rename(columns={'barcodes':'SpotID', 'x' : 'col','y' :'row'})
    labels = labels.rename(columns={'barcodes':'SpotID', 'Territory' : 'CellType'})
    labels = labels.rename(columns={'barcodes':'SpotID', 'cell_labels' : 'CellType'}) 
    labels['CellType'] = np.repeat('celltype_', labels.shape[0]) + labels['CellType'].astype(str)
    stRNA = os.path.join(output_path, f'stRNA_{str(slurm_id)}_{data_type}.txt')
    stLabels = os.path.join(output_path, f'stLabels_{str(slurm_id)}_{data_type}.txt')
    stCoord = os.path.join(output_path, f'stCoord_{str(slurm_id)}_{data_type}.txt')
    counts.to_csv(stRNA, index = False, sep = "\t", quoting = csv.QUOTE_NONE)
    coord.to_csv(stCoord, index = False, sep = "\t", quoting = csv.QUOTE_NONE)
    labels.to_csv(stLabels, index = False, sep = "\t", quoting = csv.QUOTE_NONE)
#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='CtyoSpace benchmarking on synthetic spatial data')
    parser.add_argument('slurm_id', help='Slurm ID from array sub')
    parser.add_argument('type', help = 'type of synthetic data requested')
    parser.add_argument('tmp', help='tmp output directory')
    args = parser.parse_args()
    output_path = args.tmp
    files = os.listdir("/common/wonklab/synthetic_spatial")
    files = get_files(files, args.slurm_id, args.type)
    files_used = os.path.join(output_path,'files_used.txt')
    text_file = open(files_used, "w")
    text_file.write(f"Reference:{files[0][0]}\n")
    text_file.write(f"Query:{files[1][0]}\n")
    text_file.close()
    # tmp query data
    tmp = os.path.join("/common/wonklab/synthetic_spatial",files[1][1])
    tmp = pd.read_csv(tmp)
    tmp_loc = os.path.join(output_path,'CytoSpace_tmp_query.csv')
    tmp.to_csv(tmp_loc, index = False, quoting = csv.QUOTE_NONE)
    # Creating single cell like data from synthetic data sets
    counts, labels, coord = load_sim("/common/wonklab/synthetic_spatial", files[1])
    synth_to_sc(counts, labels, args.type, output_path, args.slurm_id)
    counts, labels, coord = load_sim("/common/wonklab/synthetic_spatial", files[0])
    synth_to_spa(counts, labels, coord, args.type, output_path, args.slurm_id)
    return 0

if __name__ == "__main__":
    main()