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

def prepare_supplement(path,experiment,hamster_species = 'ma'):
    model = scgen.SCGEN.load(path  + experiment + '/model/model.pt')

    latent_adata = hVAE.get_latent_representation_object(model,model.adata)

    latent_adata_hamster = hVAE.leave_out_adata_obs(latent_adata,'cells','Whole_blood')
    latent_adata_human = hVAE.filter_adata_obs(latent_adata,'cells','Whole_blood')

    if hamster_species == 'ma':
        latent_adata_hamster = hVAE.control_infection_annotation_ma(latent_adata_hamster)
    if hamster_species == 'pr':     
            latent_adata_hamster = hVAE.control_infection_annotation_pr(latent_adata_hamster)
    latent_adata_human = hVAE.control_infection_annotation_human(latent_adata_human)
    lat_merged = hVAE.merge_adata(latent_adata_human,latent_adata_hamster)

    lat_hamster_control = hVAE.filter_adata_obs(latent_adata_hamster,'general_type','hamster_control')
    lat_hamster_infection = hVAE.filter_adata_obs(latent_adata_hamster,'general_type','hamster_infection')
    lat_human_control = hVAE.filter_adata_obs(latent_adata_human,'general_type','human_control')
    lat_human_infection = hVAE.filter_adata_obs(latent_adata_human,'general_type','human_infection')
    delta = hVAE.get_delta_in_latent_space(lat_hamster_control,lat_human_control)
    lat_adata_hamster_shifted = hVAE.shift_adata_in_latent_space(latent_adata_hamster,delta)    

    if hamster_species == 'ma':
        lat_adata_hamster_shifted = hVAE.control_infection_annotation_ma_shifted(lat_adata_hamster_shifted)
    if hamster_species == 'pr':
        lat_adata_hamster_shifted = hVAE.control_infection_annotation_pr_shifted(lat_adata_hamster_shifted)

    lat_adata_hamster_shifted.obs['compare_umap'] = lat_adata_hamster_shifted.obs['compare_umap'].apply(lambda x: x + '_shifted')

    lat_all_overview = hVAE.merge_adata(lat_merged,lat_adata_hamster_shifted)
    return lat_all_overview

color_dict = {}

color_dict['d2_shifted'] = "#1d4497"
color_dict['d3_shifted'] = "#5773c0"
color_dict['d5_shifted'] = "#7da7ea"
color_dict['e14_shifted'] = "#8cc8bc"

color_dict['hd_D2_shifted'] = "#3a507f"
color_dict['hd_D3_shifted'] = "#369cc9"
color_dict['ld_D2_shifted'] = "#31c7ba"
color_dict['ld_D3_shifted'] = "#1c9d7c"

color_dict['3'] = "#f7c267"
color_dict['4'] = "#d39a2d"
color_dict['5'] = "#9b332b"
color_dict['7'] = "#591c19"

run_dict_ma = {}
run_dict_ma['Non_Classical_Monocytes'] = 'B4_run_2'
run_dict_ma['NK_Cells'] = 'B4_run_2'
run_dict_ma['Neutrophils_acurate'] = 'B4_run_3'
run_dict_ma['B_Cells'] = 'B4_run_5'
run_dict_ma['CD4+_T_Cells'] = 'B4_run_5'
run_dict_ma['CD8+_T_Cells'] = 'B4_run_4'
run_dict_ma['Classical_Monocytes'] = 'B4_run_2'
run_dict_ma['Immature_Neutrophils_acurate'] = 'B4type1_h270_run_4'

run_dict_pr = {}
run_dict_pr['Non_Classical_Monocytes'] = 'B4_run_3'
run_dict_pr['NK_Cells'] = 'B4_run_5'
run_dict_pr['Neutrophils_acurate'] = 'B4_run_3'
run_dict_pr['B_Cells'] = 'B4_run_4'
run_dict_pr['CD4+_T_Cells'] = 'B4_run_5'
run_dict_pr['CD8+_T_Cells'] = 'B4_run_1'
run_dict_pr['Classical_Monocytes'] = 'B4_run_3'
run_dict_pr['Immature_Neutrophils_acurate'] = 'B4type1_h270_run_5'

celltypes = ['Non_Classical_Monocytes', 'NK_Cells', 'Neutrophils_acurate', 'B_Cells', 'CD4+_T_Cells', 'CD8+_T_Cells','Classical_Monocytes','Immature_Neutrophils_acurate'] 


#ma
for celltype in celltypes:
    path = '/work/users/username/projects/cov/integration/run_VAE_228_1/backup_results/results_v228_1/ma/' + celltype
    hamster_species = 'ma'
    experiment = '/'  + str(run_dict_ma[celltype])
    lat_all_overview  = prepare_supplement(path,experiment,hamster_species)
    hamster_degrees = ['d2_shifted','d3_shifted','d5_shifted','e14_shifted']
    fig,axes = plt.subplots(2, 4, figsize=(20, 8))
    current_axes = -1
    for hamster_degree in hamster_degrees:
        current_axes += 1
        comp_degrees = ['3','4','5','7']
        comp_degrees.append(hamster_degree)
        lat_comp = lat_all_overview[lat_all_overview.obs['compare_umap'].isin(comp_degrees)]
        hVAE.prepare_umap(lat_comp)
        palette_demuth = [color_dict['3'], color_dict['4'], color_dict['5'],color_dict['7'],color_dict[hamster_degree]]
        sc.pl.umap(lat_comp,color = 'compare_umap', palette=palette_demuth,ax = axes[0,current_axes],show= False)

    comp_degrees = ['3','4','5','7']
    current_axes = -1
    for comp_degree in comp_degrees:
        current_axes += 1
        hamster_degrees = ['d2_shifted','d3_shifted','d5_shifted','e14_shifted']
        hamster_degrees.append(comp_degree)
        lat_comp = lat_all_overview[lat_all_overview.obs['compare_umap'].isin(hamster_degrees)]
        hVAE.prepare_umap(lat_comp)
        palette_demuth = [color_dict[comp_degree], color_dict['d2_shifted'], color_dict['d3_shifted'],color_dict['d5_shifted'],color_dict['e14_shifted']]
        sc.pl.umap(lat_comp,color = 'compare_umap', palette=palette_demuth,ax = axes[1,current_axes],show= False)
    plt.tight_layout()
    plt.savefig('ma/s1_' + + str(celltype) + '_supp.pdf')
  
    
#pr
for celltype in celltypes:
    path = '/work/users/username/projects/cov/integration/run_VAE_228_1/backup_results/results_v228_1/pr/' + celltype
    hamster_species = 'pr'
    experiment = '/'  + str(run_dict_pr[celltype])
    lat_all_overview  = prepare_supplement(path,experiment,hamster_species = 'pr')
    hamster_degrees = ['hd_D2_shifted','hd_D3_shifted','ld_D2_shifted','ld_D3_shifted']
    fig,axes = plt.subplots(2, 4, figsize=(20, 8))
    current_axes = -1
    for hamster_degree in hamster_degrees:
        current_axes += 1
        comp_degrees = ['3','4','5','7']
        comp_degrees.append(hamster_degree)
        lat_comp = lat_all_overview[lat_all_overview.obs['compare_umap'].isin(comp_degrees)]
        hVAE.prepare_umap(lat_comp)
        palette_demuth = [color_dict['3'], color_dict['4'], color_dict['5'],color_dict['7'],color_dict[hamster_degree]]
        sc.pl.umap(lat_comp,color = 'compare_umap', palette=palette_demuth,ax = axes[0,current_axes],show= False)

    comp_degrees = ['3','4','5','7']
    current_axes = -1
    for comp_degree in comp_degrees:
        current_axes += 1
        hamster_degrees = ['hd_D2_shifted','hd_D3_shifted','ld_D2_shifted','ld_D3_shifted']
        hamster_degrees.append(comp_degree)
        lat_comp = lat_all_overview[lat_all_overview.obs['compare_umap'].isin(hamster_degrees)]
        hVAE.prepare_umap(lat_comp)
        palette_demuth = [color_dict[comp_degree], color_dict['hd_D2_shifted'], color_dict['hd_D3_shifted'],color_dict['ld_D2_shifted'],color_dict['ld_D3_shifted']]
        sc.pl.umap(lat_comp,color = 'compare_umap', palette=palette_demuth,ax = axes[1,current_axes],show= False)
    plt.tight_layout()
    plt.savefig('pr/s1_' + str(celltype) + '_supp.pdf')

