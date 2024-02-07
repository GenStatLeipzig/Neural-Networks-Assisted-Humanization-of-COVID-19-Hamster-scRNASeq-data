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

base = "/work/users/username/projects/cov/integration" 
species = 'ma'
object_number = '228_1'
script_number = 'B01'
celltype = 'Immature_Neutrophils'
hVAE.make_folder_if_not_existing(celltype + '_type1')

path_counts = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_MA.h5ad'
path_obs = base + '/data/05_newPrepro/04_moreOrthologues/v228_1_MA_obs.csv'
path_var = base + '/data/05_newPrepro/04_moreOrthologues/v227_2_MA_var.csv'
path_obs_v228_2 = base + '/data/05_newPrepro/04_moreOrthologues/h228_2_qc_integrated_rpca_cellanno_FILTEREDqc_FILTEREDpureCells.txt'


human,hamster = hVAE.process_v228_1_Neutrophil_proof_principle(path_counts,path_obs,path_var)
path_h270_1 = base + '/data/05_newPrepro/04_moreOrthologues/h270_1_clusterlabels_in_h228_3_seurat.txt'
df_h270 = pd.read_csv(path_h270_1,index_col = 0)
human.obs = human.obs.join(df_h270[['uniform_name_overview3']])
hamster.obs = hamster.obs.join(df_h270[['uniform_name_overview3']])
celltype_string = 'uniform_name_overview3'


#subtype 1
human_ct = hVAE.filter_adata_obs(human,col_name=celltype_string,val='Immature Neutrophils 1')
hamster_ct = hVAE.filter_adata_obs(hamster,col_name=celltype_string,val='Immature Neutrophils 1')


human.obs['donor'].value_counts()

sc.pp.highly_variable_genes(hamster_ct,layer = 'log_counts')
sc.pp.highly_variable_genes(human_ct,layer = 'log_counts')

hamster_ct.var['highly_variable']
hvg_hamster_ct = list(hamster_ct.var['highly_variable'][hamster_ct.var['highly_variable'] == True].index)

human_ct.var['highly_variable']
hvg_human_ct = list(human_ct.var['highly_variable'][human_ct.var['highly_variable'] == True].index)

hvg_hamster_ct = set(hvg_hamster_ct)
hvg_human_ct = set(hvg_human_ct)
union_hvg = list(hvg_hamster_ct.union(hvg_human_ct))

hamster_ct = hamster_ct[:,union_hvg]

human_ct = human_ct[:,union_hvg]

hamster_ct_0 = hVAE.filter_adata_obs(hamster_ct,col_name='timepoint',val = 'd0')

human_ct_0 = hVAE.filter_adata_obs(human_ct,col_name='who_per_sample',val = '0')

binary_hamster = list(hamster_ct_0.obs['hamster'] == 'Ha1')
hamster_set_variable = ['test' if value else 'train' for value in binary_hamster]
hamster_ct_0.obs['set'] = hamster_set_variable


binary_human = list(human_ct_0.obs['donor'] == 'BN-31')
human_set_variable = ['test' if value else 'train' for value in binary_human]
human_ct_0.obs['set'] = human_set_variable

controls = hVAE.merge_adata(hamster_ct_0,human_ct_0)
controls.obs['celltype'] = celltype





## train test split 


controls_train = controls[controls.obs['set'] == 'train']
controls_test = controls[controls.obs['set'] == 'test']

controls_train.write_h5ad(celltype + '_type1' + '/controls_train_pp_' + hVAE.make_identifier(object_number=object_number,
                     script_number=script_number,
                     celltype = celltype + '_type1') + '.h5ad') 
controls_test.write_h5ad(celltype + '_type1' + '/controls_test_pp_' + hVAE.make_identifier(object_number=object_number,
                     script_number=script_number,
                     celltype = celltype + '_type1') + '.h5ad') 


batch_key = 'dataset'
labels_key = 'celltype'

model_save_string = celltype + '_type1' + '/' + hVAE.make_identifier(object_number=object_number,
                     script_number=script_number,
                     celltype = celltype + '_type1') + '_model.pt'


controls_train = controls_train.copy()
model = hVAE.train_VAE(controls_train,batch_key=batch_key,labels_key= labels_key,n_latent = 10,save_string = model_save_string)
