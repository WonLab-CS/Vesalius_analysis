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

def dup_merge(aligned, query):
  query = query.loc[aligned.index]
  query = query.reset_index(drop = True)
  aligned = aligned.reset_index(drop = True)
  new_df = pd.concat([aligned[['barcodes','x','y']],query[['cell_labels','interactions']]],axis =1)
  return new_df

def main():
    parser = argparse.ArgumentParser(description='CytoSpace benchmarking on synthetic spatial data - file clean')
    parser.add_argument('tmp', help='tmp output directory')
    parser.add_argument('out', help='final out')
    parser.add_argument('data_type', help='data type')
    parser.add_argument('seed_tag', help='seed tag')
    parser.add_argument('query_tag', help='query tag')
    args = parser.parse_args()
    input_path = args.tmp
    output_path = args.out
    if 'trueLab' in str(os.listdir(input_path)):
        lab = f'{args.data_type}trueLab'
    else :
        lab = f'{args.data_type}noLab'
    files_used = os.path.join(input_path,'files_used.txt')
    files_used = open(files_used,'r')
    files_used = files_used.readlines()
    aligned = os.path.join(input_path, 'assigned_locations.csv')
    aligned = pd.read_csv(aligned)
    aligned = aligned.rename(columns={'OriginalCID':'barcodes', 'col' : 'x','row' :'y', 'CellType': 'cell_labels'})
    aligned = aligned.set_index('barcodes', drop = False).rename_axis(None)
    query = os.path.join(input_path, "CytoSpace_tmp_query.csv")
    query = pd.read_csv(query)
    query = query.set_index('barcodes', drop = False).rename_axis(None)
    new_df = dup_merge(aligned,query)
    tag = f'CytoSpace_aligned_{lab}_{args.seed_tag}_{args.query_tag}.csv'
    out_file = os.path.join(output_path,tag)
    new_df.to_csv(out_file,quoting = csv.QUOTE_NONE)
    
    
    

if __name__ == "__main__":
    main()