import numpy as np
import anndata as ad
import pandas as pd
import os
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import argparse
import re
from gpsa import VariationalGPSA
from gpsa import matern12_kernel, rbf_kernel
from gpsa.plotting import callback_twod
import pickle


#-----------------------------------------------------------------------------#
# Configure for cpu 
#-----------------------------------------------------------------------------#
device = torch.device('cpu')
torch.set_num_threads(5)
os.environ["OMP_NUM_THREADS"] = "5" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "5" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "5" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "5" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "5" # export NUMEXPR_NUM_THREADS=1



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
    #adata.obs["Territory"] = coord['Territory'].to_numpy()
    adata.obs['cell_labels'] = coord['cell_labels'].to_numpy()
    adata.obs['interactions'] = coord['interactions'].to_numpy()
    return adata

def export_coord(query,coord,output_path, tag):
    spa = pd.DataFrame(coord)
    spa = spa.set_axis(query.obs.index.values, axis = 0)
    spa = spa.set_axis(['x','y'], axis = 1)
    barcodes = pd.DataFrame(query.obs.index.values)
    barcodes = barcodes.set_axis(query.obs.index.values, axis = 0)
    barcodes = barcodes.set_axis(['barcodes'], axis = 1)
    export_adata = pd.concat([barcodes, spa, query.obs[['cell_labels','interactions']]], axis = 1)
    export_file = os.path.join(output_path, tag)
    export_adata.to_csv(export_file)


# Function taken from in GPSA tutorial and modified for this purpose
def train(model, loss_fn, optimizer ,data_dict, view_idx, Ns, x):
    model.train()
    # Forward pass
    G_means, G_samples, F_latent_samples, F_samples = model.forward(
        {"expression": x}, view_idx = view_idx, Ns = Ns, S=5
    )
    # Compute loss
    loss = loss_fn(data_dict, F_samples)
    # Compute gradients and take optimizer step
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    ls = loss.item()
    return ls, G_means, G_samples, F_latent_samples, F_samples

#-----------------------------------------------------------------------------#
# Main 
#-----------------------------------------------------------------------------#

def main():
    parser = argparse.ArgumentParser(description='GPSA benchmarking on synthetic spatial data')
    parser.add_argument('data_type', help='Slurm ID from array sub')
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
    
    data = {0: seed, 1: query}
    data = ad.concat(data, label = "batch")

    X = data.obsm["spatial"]
    Y = data.X
    view_idx = [np.where(data.obs.batch.values == ii)[0] for ii in range(2)]
    n_samples_list = [len(x) for x in view_idx]

    x = torch.from_numpy(X).float().clone()
    y = torch.from_numpy(Y).float().clone()

    data_dict = {
        "expression": {
            "spatial_coords": x,
            "outputs": y,
            "n_samples_list": n_samples_list,
        }
    }
    
    N_EPOCHS = 50
    model = VariationalGPSA(
        data_dict,
        n_spatial_dims=2,
        m_X_per_view=50,
        m_G=50,
        data_init=True,
        minmax_init=False,
        grid_init=False,
        n_latent_gps={"expression": None},
        mean_function="identity_fixed",
        kernel_func_warp=rbf_kernel,
        kernel_func_data=rbf_kernel,
        fixed_view_idx=0,
    ).to(device)
    view_idx, Ns, _, _ = model.create_view_idx_dict(data_dict)

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
    for t in range(N_EPOCHS):
        loss, G_means, G_samples, F_latent_samples, F_samples = train(model, model.loss_fn, optimizer, data_dict,view_idx, Ns,x)
    
    curr_aligned_coords = G_means["expression"].detach().numpy()
    curr_aligned_coords_slice2 = curr_aligned_coords[view_idx["expression"][1]]
    tag = f'GPSA_aligned_{args.data_type}_{args.seed_tag}_{args.query_tag}.csv'
    export_coord(query, curr_aligned_coords_slice2, output_path, tag)

if __name__ == "__main__":
    main()