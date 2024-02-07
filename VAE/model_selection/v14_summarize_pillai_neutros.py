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

pre = 'v13'
save_folder = '/work/users/username/projects/cov/integration/run_VAE_228_1/results/model_selection/pillais_trace/'
base = "/work/users/username/projects/cov/integration" 
celltypes = ['CD4+_T_Cells', 'Classical_Monocytes', 'B_Cells','CD8+_T_Cells','NK_Cells','Non_Classical_Monocytes','Neutrophils_acurate','Immature_Neutrophils_acurate']
#celltypes = ['Neutrophils_acurate','Immature_Neutrophils_acurate']
all_species = ['ma','pr']

df_pre = pd.read_csv(save_folder + 'v13' + '_' + 'ma' + '_' + 'Classical_Monocytes' + '.csv',index_col = 0)
df = pd.DataFrame(columns = df_pre.columns)
species = 'ma'
for celltype in celltypes:
    df_temp = pd.DataFrame(pd.read_csv(save_folder + pre + '_' + species + '_' + celltype + '.csv',index_col = 0).loc['Value']).T
    df_temp = df_temp.rename(index={'Value':celltype})
    df = pd.concat([df,df_temp])
df.to_csv(save_folder + 'v14' + '_' + species + '_' + 'all_pillai.csv')


df_pre = pd.read_csv(save_folder + 'v13' + '_' + 'ma' + '_' + 'Classical_Monocytes' + '.csv',index_col = 0)
df = pd.DataFrame(columns = df_pre.columns)
species = 'pr'
for celltype in celltypes:
    df_temp = pd.DataFrame(pd.read_csv(save_folder + pre + '_' + species + '_' + celltype + '.csv',index_col = 0).loc['Value']).T
    df_temp = df_temp.rename(index={'Value':celltype})
    df = pd.concat([df,df_temp])
df.to_csv(save_folder + 'v14' + '_' + species + '_' + 'all_pillai.csv')
