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
species = 'ma'
celltype = 'Immature_Neutrophils'
save_location = base +'/run_VAE_228_1/results/' + species + '/' + celltype + '_acurate'

path_counts = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_MA.h5ad'
path_obs = base + '/data/05_newPrepro/04_moreOrthologues/v228_1_MA_obs.csv'
path_obs_v228_2 = base + '/data/05_newPrepro/04_moreOrthologues/h228_2_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.txt'
path_var = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_MA_var.csv'


human,hamster = hVAE.process_v228_1_Neutrophil(path_counts,path_obs,path_var)
path_h270_1 = base + '/data/05_newPrepro/04_moreOrthologues/h270_1_clusterlabels_in_h228_3_seurat.txt'
df_h270 = pd.read_csv(path_h270_1,index_col = 0)
human.obs = human.obs.join(df_h270[['uniform_name_overview3']])
hamster.obs = hamster.obs.join(df_h270[['uniform_name_overview3']])
celltype_string = 'uniform_name_overview3'
#subtype 1
human_ct = hVAE.filter_adata_obs(human,col_name=celltype_string,val='Immature Neutrophils 1')
hamster_ct = hVAE.filter_adata_obs(hamster,col_name=celltype_string,val='Immature Neutrophils 1')

for setting_str in ['B4']:
    gene_set,n_latent,do_batch_removal_human = hVAE.choose_setting(setting_str)
    for num in range(1,6):
        experiment_name = setting_str + 'type1_h270_run_' + str(num)
        hVAE.create_folder_structure(save_location,experiment_name)
        hVAE.run_VAE_celltype(hamster_ct,
                 human_ct,
                 species = species,
                 celltype = 'Immature Neutrophils 1',
                 save_input=False,
                 save_location = save_location,
                 experiment_name = experiment_name,
                 gene_set = gene_set,
                 n_latent = n_latent,
                 make_umaps = False,
                 batch_key = 'compare_umap',
                 labels_key = celltype_string)
