import scanpy as sc
import numpy as np
import scgen as scg
import pandas as pd
import sys
sys.path.append("/work/users/mh823zote/projects/cov/integration/analysis")
import helper_VAE as hVAE
import scipy.sparse
import scgen
import matplotlib.pyplot as plt
import seaborn as sns
import harmonypy as hpy
import anndata
import random
import os

base = "/work/users/mh823zote/projects/cov/integration" 
species = 'pr'
object_number = '228_1'
script_number = 'B01'
#celltypes = ['B_Cells','Classical_Monocytes','NK_Cells','Non_Classical_Monocytes']
celltypes = ['Non_Classical_Monocytes']
for celltype in celltypes:
    hVAE.make_folder_if_not_existing(celltype)
    

    path_counts = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_PR.h5ad'
    path_obs = base + '/data/05_newPrepro/04_moreOrthologues/v228_1_PR_obs.csv'
    path_var = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_PR_var.csv'


    _,hamster = hVAE.process_v228_1_proof_principle(path_counts,path_obs,path_var)
    celltype_string = 'uniform_name_overview'
    hamster_ct = hVAE.filter_adata_obs(hamster,col_name=celltype_string,val=celltype)
    
    sc.pp.highly_variable_genes(hamster_ct,layer = 'log_counts')
    hamster_ct.var['highly_variable']
    hvg_hamster_ct = list(hamster_ct.var['highly_variable'][hamster_ct.var['highly_variable'] == True].index)
    hamster_ct = hamster_ct[:,hvg_hamster_ct]

    #hamster_ct

    binary_hamster = list(hamster_ct.obs['hamster'] == 'Ha1')
    hamster_set_variable = ['test' if value else 'train' for value in binary_hamster]
    hamster_ct.obs['set'] = hamster_set_variable

    hamster_ct.obs = hamster_ct.obs[['uniform_name_overview','timepoint','hamster','set']]

    hamster_ct.obs[['set','hamster']].value_counts()

    hamster_ct_train = hVAE.filter_adata_obs(hamster_ct,col_name='set',val='train')

    hamster_ct_test = hVAE.filter_adata_obs(hamster_ct,col_name='set',val='test')

    hamster_ct_train.write_h5ad(celltype + '/H0_01_train_' + celltype + '_.h5ad')

    hamster_ct_test.write_h5ad(celltype + '/H0_01_test_' + celltype + '_.h5ad')
    
    # Train VAE

import os

batch_key = 'timepoint'
labels_key = 'uniform_name_overview'
model_save_string = os.getcwd()

model_save_string = model_save_string + '/' + celltype + '/H0_01_model_' + celltype + '.pt'

hamster_ct_train = hamster_ct_train.copy()
model = hVAE.train_VAE(hamster_ct_train,
                       batch_key=batch_key,
                       labels_key= labels_key,
                       n_latent = 10,
                       save_string = model_save_string)
    
    