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

base = "/work/users/username/projects/cov/integration" 
species = 'pr'
celltype = 'Non_Classical_Monocytes'
save_location = base +'/run_VAE_228_1/results/' + species + '/' + celltype

path_counts_PR = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_PR.h5ad'
path_obs_PR = base + '/data/05_newPrepro/04_moreOrthologues/v228_1_PR_obs.csv'
path_var_PR = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_PR_var.csv'

human,hamster = hVAE.process_v228_1(path_counts_PR,path_obs_PR,path_var_PR)


human_ct = hVAE.filter_adata_obs(human,col_name='uniform_name_overview',val=celltype)
hamster_ct = hVAE.filter_adata_obs(hamster,col_name='uniform_name_overview',val=celltype)
for setting_str in ['B4']:
    gene_set,n_latent,do_batch_removal_human = hVAE.choose_setting(setting_str)
    for num in range(1,6):
        experiment_name = setting_str + '_run_' + str(num)
        hVAE.create_folder_structure(save_location,experiment_name)
        hVAE.run_VAE_celltype(hamster_ct,
                 human_ct,
                 species = species,
                 celltype = celltype,
                 save_input=False,
                 save_location = save_location,
                 experiment_name = experiment_name,
                 gene_set = gene_set,
                 n_latent = n_latent,
                 make_umaps = False,
                 batch_key = 'compare_umap',
                 labels_key = 'uniform_name_overview')
