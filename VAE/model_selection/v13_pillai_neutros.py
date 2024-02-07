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

celltypes = ['Neutrophils_acurate','Immature_Neutrophils_acurate']
#celltypes = ['CD4+_T_Cells', 'Classical_Monocytes', 'B_Cells','CD8+_T_Cells', 'NK_Cells', 'Non_Classical_Monocytes']
all_species = ['ma','pr']

base_results = '/work/users/username/projects/cov/integration/run_VAE_228_1/results/'

for species in all_species:
    for celltype in celltypes:
        pillais_overview = pd.DataFrame(index = ['Value', 'Num DF', 'Den DF', 'F Value', 'Pr > F'])
        settings = ['B4']
        runs = [1,2,3,4,5]
        for setting in settings:
            for run in runs:
                p_df = hVAE.v227_2_get_pillais_trace(base_results,species,celltype,setting,run)
                pillais_overview = pd.merge(p_df , pillais_overview, left_index=True, right_index=True)
        pillais_overview.to_csv(save_folder + pre + '_' + species + '_' + celltype + '.csv')
