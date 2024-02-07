import scanpy as sc
import numpy as np
import scgen as scg
import pandas as pd
import sys
sys.path.append("/work/users/username/projects/cov/integration/analysis")
import helper_VAE as hVAE
import scipy.sparse
import scgen
import matplotlib.pyplot as plt
import seaborn as sns
import harmonypy as hpy
import anndata
import random
import os

base_pr = "/work/users/username/projects/cov/integration/run_proof_principle_228_1/backup/pr/" 
species = 'pr'
object_number = '228_1'
script_number = 'B02'
celltypes = ['Immature_Neutrophils_type1']
for celltype in celltypes:
    model_save_string = base_pr + celltype + '/228_1_B01_' + celltype + '_model.pt' 
    model = scgen.SCGEN.load(model_save_string)

    adata = model.adata
    lat_adata = hVAE.get_latent_representation_object(model,adata)


    lat_adata_hamster = lat_adata[lat_adata.obs['dataset'] == 'hamsterPR']
    lat_adata_human = lat_adata[lat_adata.obs['dataset'] == 'human']
    delta = hVAE.get_delta_in_latent_space(lat_adata_hamster,lat_adata_human)


    ### test set 

    test_adata = sc.read_h5ad(celltype + '/controls_test_pp_228_1_B01_' + celltype +'.h5ad')



    lat_test_adata = hVAE.get_latent_representation_object(model,test_adata)

    lat_hamster_test = lat_test_adata[lat_test_adata.obs['dataset'] == 'hamsterPR']

    lat_hamster_test_shifted = hVAE.shift_adata_in_latent_space(lat_hamster_test,delta)

    lat_hamster_test_shifted.obs['dataset'] = 'humanized hamster'

    lat_test_predicted = hVAE.merge_adata(lat_test_adata,lat_hamster_test_shifted)


    lat_test_predicted = hVAE.merge_adata(lat_test_adata,lat_hamster_test_shifted)


    decoded_test = hVAE.decode_latent_object(model,lat_test_predicted,test_adata)

    gene_expression_df = pd.DataFrame(decoded_test.X,columns = decoded_test.var.index,index = decoded_test.obs.index)
    obs_df = decoded_test.obs[['species','celltype','dataset']]
    obs_df.to_csv(celltype + '/B05_meta_' + celltype + '_' + species + '.csv')
    gene_expression_df.to_csv(celltype + '/B05_gene_expression_' + celltype + '_' + species + '.csv')
