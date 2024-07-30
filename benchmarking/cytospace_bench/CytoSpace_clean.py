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
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='CytoSpace benchmarking on synthetic spatial data - file clean')
    parser.add_argument('tmp', help='tmp output directory')
    parser.add_argument('out', help='final out')
    args = parser.parse_args()
    input_path = args.tmp
    output_path = args.out
    files_used = os.path.join(input_path,'files_used.txt')
    files_used = open(files_used,'r')
    files_used = files_used.readlines()
    aligned = os.path.join(input_path, 'assigned_locations.csv')
    aligned = pd.read_csv(aligned)
    aligned = aligned.rename(columns={'OriginalCID':'barcodes', 'col' : 'x','row' :'y', 'CellType': 'Territory'})
    aligned = aligned.set_index('barcodes', drop = False).rename_axis(None)
    query = os.path.join(input_path, "CytoSpace_tmp_query.csv")
    query = pd.read_csv(query)
    query = query.set_index('barcodes', drop = False).rename_axis(None)
    new_df = pd.concat([aligned[['barcodes','x','y']],query[['cell_labels']]],axis =1)
    tag = f'CytoSpace_aligned_{re.sub("Reference:","",files_used[0])}_{re.sub("Query:","",files_used[1])}'
    tag = re.sub('\n','',tag)
    tag = re.sub('spatial_territories_gene_counts_','',tag)
    tag = re.sub('.csv','',tag)
    tag = tag + '.csv'
    out_file = os.path.join(output_path,tag)
    new_df.to_csv(out_file,quoting = csv.QUOTE_NONE)
    
    
    

if __name__ == "__main__":
    main()