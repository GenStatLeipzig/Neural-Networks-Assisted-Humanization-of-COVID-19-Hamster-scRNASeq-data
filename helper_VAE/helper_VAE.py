import numpy as np
import pandas as pd
import scanpy as sc
import scgen
from scipy import stats
import torch
import anndata as ad
from anndata import AnnData
import random
import matplotlib.pyplot as plt
import scvi
import os
import seaborn as sns
import harmonypy as hpy
import scipy.stats as scistats
from statsmodels.multivariate.manova import MANOVA

def filter_adata_obs(adata,col_name,val):
    '''
    adata    :   anndata object
    col_name :   string, column name
    val      :   string or float, value in column to filter for
    '''
    return adata[adata.obs[col_name] == val]

def leave_out_adata_obs(adata,col_name,val):
    '''
    Leave out specific value
    adata    :   anndata object
    col_name :   string, column name
    val      :   string or float, value in column to filter for
    '''
    return adata[~(adata.obs[col_name] == val)]

def merge_adata(adata1,adata2):
    '''
    adata1 : anndata object
    adata2 : anndata object
    '''
    return adata1.concatenate(adata2)

def show_unique(adata,col_name):
    print(pd.unique(adata.obs[col_name]))
    

def prepare_umap(adata):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
def get_latent_representation_object(model,adata):
    latent_X = model.get_latent_representation(adata=adata)
    latent_adata = sc.AnnData(X=latent_X, obs=adata.obs.copy())
    return latent_adata

def get_latent_shift(model,adata,delta):
    pred_latent_X = model.get_latent_representation(adata=adata) + delta
    latent_pred_adata = sc.AnnData(X=pred_latent_X, obs=adata.obs.copy())
    return latent_pred_adata

def decode(model,latent_adata,input_adata):
    predicted_cells = (
                model.module.generative(torch.Tensor(latent_adata.X))["px"].cpu().detach().numpy()
    )
    predicted_adata = AnnData(
                X=predicted_cells,
                obs=input_adata.obs.copy(),
                var=input_adata.var.copy(),
                obsm=input_adata.obsm.copy(),
            )
    return predicted_adata

def make_string_column_total(adata,column):
    adata.obs[column] = make_string_column(adata.obs[column])

def make_string_column(column):
    new_column = []
    for j in range(len(column)):
        new_column.append(str(column[j]))
    return new_column

def calculate_pseudotime_median(adata,repetitions,column_name,column_value):
    sc.pp.neighbors(adata)
    sc.tl.diffmap(adata)
    degrees = sorted(pd.unique(adata.obs[column_name]))
    random_keys = random.sample(list(np.flatnonzero(adata.obs[column_name] == column_value)),repetitions)
    
    df_dict = {}
    for j in range(repetitions):
        pseudotime_dict = {}
        adata_temp = adata.copy()
        root_loc = random_keys[j]
        adata_temp.uns['iroot'] = root_loc
        sc.tl.dpt(adata_temp,copy=False)
        pseudo_list = []
        for DEG in degrees:
            pt_DEG = filter_adata_obs(adata_temp,column_name,DEG).obs['dpt_pseudotime']
            pt_rem_inf = list(filter(lambda val: val !=  np.inf, pt_DEG))
            pseudo_list.append(np.median(pt_rem_inf))
        pseudotime_dict[column_value] = pseudo_list
        df_dict[j] = pd.DataFrame(pseudotime_dict,index = degrees)
    df_final = df_dict[0].copy()/repetitions
    for q in range(1,repetitions):
        df_final += df_dict[q]/repetitions
    return df_final


def calculate_pseudotime(adata,repetitions,column_name,column_value):
    sc.pp.neighbors(adata)
    sc.tl.diffmap(adata)
    degrees = sorted(pd.unique(adata.obs[column_name]))
    random_keys = random.sample(list(np.flatnonzero(adata.obs[column_name] == column_value)),repetitions)
    
    df_dict = {}
    for j in range(repetitions):
        pseudotime_dict = {}
        adata_temp = adata.copy()
        root_loc = random_keys[j]
        adata_temp.uns['iroot'] = root_loc
        sc.tl.dpt(adata_temp,copy=False)
        pseudo_list = []
        for DEG in degrees:
            pt_DEG = filter_adata_obs(adata_temp,column_name,DEG).obs['dpt_pseudotime']
            pt_rem_inf = list(filter(lambda val: val !=  np.inf, pt_DEG))
            pseudo_list.append(np.mean(pt_rem_inf))
        pseudotime_dict[column_value] = pseudo_list
        df_dict[j] = pd.DataFrame(pseudotime_dict,index = degrees)
    df_final = df_dict[0].copy()/repetitions
    for q in range(1,repetitions):
        df_final += df_dict[q]/repetitions
    return df_final

def make_compare_column_schulte(adata):
    compare_column = []
    for q in range(len(adata.obs)):
        compare_column.append(str(adata.obs['who_per_sample'][q]))
    adata.obs['compare_umap'] =  compare_column 
    
def make_compare_column_yapeng(adata):
    compare_column = []
    for q in range(len(adata.obs)):
        compare_column.append(str(adata.obs['Who.Ordinal.Scale'][q]))
    adata.obs['compare_umap'] =  compare_column 
    

def make_compare_column(adata):
    compare_column = []
    for q in range(len(adata.obs)):
        if adata.obs['species'][q] == 'humanSchulte':
            compare_column.append(str(adata.obs['who_per_sample'][q]))
        else:
            compare_column.append(str(adata.obs['timepoint'][q]))
    adata.obs['compare_umap'] =  compare_column 
    
def make_compare_column_new(adata):
    compare_column = []
    for q in range(len(adata.obs)):
        if adata.obs['species'][q] == 'humanSchulte':
            compare_column.append(str(adata.obs['who_per_sample'][q]))
        else:
            compare_column.append(str(adata.obs['timepoint_new'][q]))
    adata.obs['compare_umap'] =  compare_column 
    
def get_delta(adata_base,adata_dest,model):
    control_base = np.mean(model.get_latent_representation(adata_base), axis=0)
    control_dest = np.mean(model.get_latent_representation(adata_dest), axis=0)
    delta = control_dest - control_base
    return delta

def pseudo_plot(df):
    plt.bar(['0', '3', '4', '5', '6', '7'],df.iloc[0:6][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
def pseudo_plot_cohort2(df):
    plt.bar(['0', '3', '4', '5', '7'],df.iloc[0:5][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
def pseudo_plot_hamster(df):
    plt.bar(['d0', 'd2', 'd3', 'd5', 'e14'],df.iloc[0:5][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
def pseudo_plot_hamster_infected(df):
    plt.bar(['d2', 'd3', 'd5', 'e14'],df.iloc[0:4][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
def pseudo_plot_hamsterPR_infected(df):
    plt.bar(['hd_d2', 'hd_d3', 'ld_d2', 'ld_d3'],df.iloc[0:4][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
def pseudo_plot_hamsterPR(df):
    plt.bar(['d0', 'hd_d2', 'hd_d3', 'ld_d2', 'ld_d3'],df.iloc[0:5][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
def Intersection(list1, list2):
    return set(list1).intersection(list2)
def filter_for_common_genes(adata,common_genes):
    gene_pos = []
    for gene in common_genes:
        pos = list(adata.var.index).index(gene)
        gene_pos.append(pos)
    return gene_pos
import scipy.sparse

def add_column(base_adata,integrated_adata,column_name):
    tags = list(base_adata.obs.index)
    column_list = []
    for tag in tags:
        try:
            column_list.append(integrated_adata.obs.loc[tag][column_name])
        except:
            column_list.append('not available')
    base_adata.obs[column_name] = column_list
    return base_adata

def do_uniform_name_accurate_schulte_dict():
    schulte_uniform_accurate_dict = {}
    schulte_uniform_accurate_dict[0] =  'Classical_Monocytes'
    schulte_uniform_accurate_dict[1] =  'HLA-DR+CD83+_Monocytes'
    schulte_uniform_accurate_dict[2] =  'CD163+_Monocytes'
    schulte_uniform_accurate_dict[3] = 'HLA-DR-S100A+_Monocytes'
    schulte_uniform_accurate_dict[4] = 'Non_Classical_Monocytes'
    schulte_uniform_accurate_dict[5] = 'Neutrophils'
    schulte_uniform_accurate_dict[6] = 'Immature_Neutrophils'
    schulte_uniform_accurate_dict[7] = 'mDCs'
    schulte_uniform_accurate_dict[8] = 'pDCs'
    schulte_uniform_accurate_dict[9] = 'CD4+_T_Cells'
    schulte_uniform_accurate_dict[10] = 'CD4+_T_Cells'
    schulte_uniform_accurate_dict[11] = 'CD4+_T_Cells'
    schulte_uniform_accurate_dict[12] = 'CD8+_T_Cells'
    schulte_uniform_accurate_dict[13] = 'CD8+_T_Cells'
    schulte_uniform_accurate_dict[14] = 'CD8+_T_Cells'
    schulte_uniform_accurate_dict[15] = 'NK_Cells'
    schulte_uniform_accurate_dict[16] = 'B_Cells'
    schulte_uniform_accurate_dict[17] = 'B_Cells'
    schulte_uniform_accurate_dict[18] = 'B_Cells'
    schulte_uniform_accurate_dict[19] = 'Plasmablasts'
    schulte_uniform_accurate_dict[20] = 'Megakaryocyte'
    schulte_uniform_accurate_dict[19] = 'Plasmablasts'
    schulte_uniform_accurate_dict[20] = 'Megakaryocyte'
    schulte_uniform_accurate_dict[21] = 'undefined'
    schulte_uniform_accurate_dict[22] = 'undefined'
    return schulte_uniform_accurate_dict


def do_uniform_name_overview_schulte_dict():
    schulte_uniform_dict = {}
    schulte_uniform_dict[0] =  'Classical_Monocytes'
    schulte_uniform_dict[1] =  'Classical_Monocytes'
    schulte_uniform_dict[2] =  'Classical_Monocytes'
    schulte_uniform_dict[3] = 'Classical_Monocytes'
    schulte_uniform_dict[4] = 'Non_Classical_Monocytes'
    schulte_uniform_dict[5] = 'Neutrophils'
    schulte_uniform_dict[6] = 'Neutrophils'
    schulte_uniform_dict[7] = 'DC'
    schulte_uniform_dict[8] = 'DC'
    schulte_uniform_dict[9] = 'T_Cells'
    schulte_uniform_dict[10] = 'T_Cells'
    schulte_uniform_dict[11] = 'T_Cells'
    schulte_uniform_dict[12] = 'T_Cells'
    schulte_uniform_dict[13] = 'T_Cells'
    schulte_uniform_dict[14] = 'T_Cells'
    schulte_uniform_dict[15] = 'NK_Cells'
    schulte_uniform_dict[16] = 'B_Cells'
    schulte_uniform_dict[17] = 'B_Cells'
    schulte_uniform_dict[18] = 'B_Cells'
    schulte_uniform_dict[19] = 'Plasmablasts'
    schulte_uniform_dict[20] = 'Megakaryocyte'
    schulte_uniform_dict[19] = 'Plasmablasts'
    schulte_uniform_dict[20] = 'Megakaryocyte'
    schulte_uniform_dict[21] = 'undefined'
    schulte_uniform_dict[22] = 'undefined'
    return schulte_uniform_dict


def do_uniform_name_accurate_hamsterMA_dict():
    accurate_dict = {}
    accurate_dict['Immature neutrophil'] =  'Immature_Neutrophils'
    accurate_dict['Neutrophil'] =  'Neutrophils'
    accurate_dict['B'] =  'B_cells'
    accurate_dict['mixed2'] = 'undefined'
    accurate_dict['T'] = 'T_cells'
    accurate_dict['Activated T'] = 'Activated_T_cells'
    accurate_dict['NK'] = 'NK_Cells'
    accurate_dict['Classical monocyte'] = 'Classical_Monocytes'
    accurate_dict['mDC'] = 'mDCs'
    accurate_dict['Platelet'] = 'Platelet'
    accurate_dict['mixed1'] = 'undefined'
    accurate_dict['Non-classical monocyte'] = 'Non_Classical_Monocytes'
    accurate_dict['mixed3'] = 'undefined'
    accurate_dict['mixed4'] = 'undefined'
    accurate_dict['mixed5'] = 'undefined'
    accurate_dict['mixed6'] = 'undefined'
    return accurate_dict

def do_uniform_name_overview_hamsterMA_dict():
    overview_dict = {}
    overview_dict['Immature neutrophil'] =  'Neutrophils'
    overview_dict['Neutrophil'] =  'Neutrophils'
    overview_dict['B'] =  'B_cells'
    overview_dict['mixed2'] = 'undefined'
    overview_dict['T'] = 'T_cells'
    overview_dict['Activated T'] = 'T_cells'
    overview_dict['NK'] = 'NK_Cells'
    overview_dict['Classical monocyte'] = 'Classical_Monocytes'
    overview_dict['mDC'] = 'DC'
    overview_dict['Platelet'] = 'Platelet'
    overview_dict['mixed1'] = 'undefined'
    overview_dict['Non-classical monocyte'] = 'Non_Classical_Monocytes'
    overview_dict['mixed3'] = 'undefined'
    overview_dict['mixed4'] = 'undefined'
    overview_dict['mixed5'] = 'undefined'
    overview_dict['mixed6'] = 'undefined'
    return overview_dict

def do_uniform_name_accurate_hamsterDwarf_dict():
    accurate_dict = {}
    accurate_dict['Immature neutrophil'] =  'Immature_Neutrophils'
    accurate_dict['Neutrophil'] =  'Neutrophils'
    accurate_dict['mixed1'] = 'undefined'
    accurate_dict['B'] =  'B_cells'
    accurate_dict['Platelet'] = 'Platelet'
    accurate_dict['T'] = 'T_cells'
    accurate_dict['Non-classical monocyte'] = 'Non_Classical_Monocytes'
    accurate_dict['Classical monocyte'] = 'Classical_Monocytes'
    accurate_dict['mDC'] = 'mDCs'
    accurate_dict['NK'] = 'NK_Cells'
    accurate_dict['mixed2'] = 'undefined'
    accurate_dict['mixed3'] = 'undefined'
    accurate_dict['mixed4'] = 'undefined'
    accurate_dict['mixed5'] = 'undefined'
    accurate_dict['mixed6'] = 'undefined'
    accurate_dict['mixed7'] = 'undefined'
    accurate_dict['mixed8'] = 'undefined'
    return accurate_dict


def do_uniform_name_overview_hamsterDwarf_dict():
    overview_dict = {}
    overview_dict['Immature neutrophil'] =  'Neutrophils'
    overview_dict['Neutrophil'] =  'Neutrophils'
    overview_dict['mixed1'] = 'undefined'
    overview_dict['B'] =  'B_cells'
    overview_dict['Platelet'] = 'Platelet'
    overview_dict['T'] = 'T_cells'
    overview_dict['Non-classical monocyte'] = 'Non_Classical_Monocytes'
    overview_dict['Classical monocyte'] = 'Classical_Monocytes'
    overview_dict['mDC'] = 'DC'
    overview_dict['NK'] = 'NK_Cells'
    overview_dict['mixed2'] = 'undefined'
    overview_dict['mixed3'] = 'undefined'
    overview_dict['mixed4'] = 'undefined'
    overview_dict['mixed5'] = 'undefined'
    overview_dict['mixed6'] = 'undefined'
    overview_dict['mixed7'] = 'undefined'
    overview_dict['mixed8'] = 'undefined'
    return overview_dict


def do_uniform_name_accurate_Yapeng_dict():
    accurate_dict = {}
    accurate_dict[0] =  'B_cells'
    accurate_dict[1] =  'DC'
    accurate_dict[2] = 'HSC'
    accurate_dict[3] =  'Megakaryocyte'
    accurate_dict[4] = 'Classical_Monocytes'
    accurate_dict[5] = 'Non_Classical_Monocytes'
    accurate_dict[6] =  'NK_Cells'
    accurate_dict[7] =  'CD4+_T_Cells'
    accurate_dict[8] =  'CD8+_T_Cells'
    return accurate_dict
    
def do_uniform_name_overview_Yapeng_dict():
    overview_dict = {}
    overview_dict[0] =  'B_cells'
    overview_dict[1] =  'DC'
    overview_dict[2] = 'HSC'
    overview_dict[3] =  'Megakaryocyte'
    overview_dict[4] = 'Classical_Monocytes'
    overview_dict[5] = 'Non_Classical_Monocytes'
    overview_dict[6] =  'NK_Cells'
    overview_dict[7] =  'T_Cells'
    overview_dict[8] =  'T_Cells'
    return overview_dict




def do_uniform_name_accurate_schulteBD_dict():
    schulte_uniform_accurate_dict = {}
    schulte_uniform_accurate_dict[0] =  'Neutrophils'
    schulte_uniform_accurate_dict[1] =  'Neutrophils'
    schulte_uniform_accurate_dict[2] =  'Neutrophils'
    schulte_uniform_accurate_dict[3] = 'Neutrophils'
    schulte_uniform_accurate_dict[4] = 'Immature_Neutrophils'
    schulte_uniform_accurate_dict[5] = 'Immature_Neutrophils'
    schulte_uniform_accurate_dict[6] = 'Eosinophils'
    schulte_uniform_accurate_dict[7] = 'Classical_Monocytes'
    schulte_uniform_accurate_dict[8] = 'Classical_Monocytes'
    schulte_uniform_accurate_dict[9] = 'Classical_Monocytes'
    schulte_uniform_accurate_dict[10] = 'Non_Classical_Monocytes'
    schulte_uniform_accurate_dict[11] = 'CD8+_T_Cells'
    schulte_uniform_accurate_dict[12] = 'CD4+_T_Cells'
    schulte_uniform_accurate_dict[13] = 'CD4+_T_Cells'
    schulte_uniform_accurate_dict[14] = 'CD4+_T_Cells'
    schulte_uniform_accurate_dict[15] = 'Proliferating_Cells'
    schulte_uniform_accurate_dict[16] = 'NK_Cells'
    schulte_uniform_accurate_dict[17] = 'B_Cells'
    schulte_uniform_accurate_dict[18] = 'B_Cells'
    schulte_uniform_accurate_dict[19] = 'Plasmablasts'
    schulte_uniform_accurate_dict[20] = 'Megakaryocyte'
    schulte_uniform_accurate_dict[21] = 'mDCs'
    schulte_uniform_accurate_dict[22] = 'pDCs'
    schulte_uniform_accurate_dict[23] = 'HSC'
    schulte_uniform_accurate_dict[24] = 'undefined'
    return schulte_uniform_accurate_dict

def do_uniform_name_overview_schulteBD_dict():
    schulte_uniform_overview_dict = {}
    schulte_uniform_overview_dict[0] =  'Neutrophils'
    schulte_uniform_overview_dict[1] =  'Neutrophils'
    schulte_uniform_overview_dict[2] =  'Neutrophils'
    schulte_uniform_overview_dict[3] = 'Neutrophils'
    schulte_uniform_overview_dict[4] = 'Neutrophils'
    schulte_uniform_overview_dict[5] = 'Neutrophils'
    schulte_uniform_overview_dict[6] = 'Eosinophils'
    schulte_uniform_overview_dict[7] = 'Classical_Monocytes'
    schulte_uniform_overview_dict[8] = 'Classical_Monocytes'
    schulte_uniform_overview_dict[9] = 'Classical_Monocytes'
    schulte_uniform_overview_dict[10] = 'Non_Classical_Monocytes'
    schulte_uniform_overview_dict[11] = 'T_Cells'
    schulte_uniform_overview_dict[12] = 'T_Cells'
    schulte_uniform_overview_dict[13] = 'T_Cells'
    schulte_uniform_overview_dict[14] = 'T_Cells'
    schulte_uniform_overview_dict[15] = 'Proliferating_Cells'
    schulte_uniform_overview_dict[16] = 'NK_Cells'
    schulte_uniform_overview_dict[17] = 'B_Cells'
    schulte_uniform_overview_dict[18] = 'B_Cells'
    schulte_uniform_overview_dict[19] = 'Plasmablasts'
    schulte_uniform_overview_dict[20] = 'Megakaryocyte'
    schulte_uniform_overview_dict[21] = 'DC'
    schulte_uniform_overview_dict[22] = 'DC'
    schulte_uniform_overview_dict[23] = 'HSC'
    schulte_uniform_overview_dict[24] = 'undefined'
    return schulte_uniform_overview_dict











def add_uniform_name_accurate(adata,translate_dict,orig_column,new_column = 'uniform_name_acurate'):
    uniform_accurate = []
    for ID in adata.obs[orig_column]:
        uniform_accurate.append(translate_dict[ID])
    adata.obs[new_column] = uniform_accurate 
    return adata

def get_delta_in_latent_space(adata_lat_base,adata_lat_dest):
    lat_mean_base = np.mean(adata_lat_base.X,axis = 0)
    lat_mean_dest = np.mean(adata_lat_dest.X,axis = 0)
    delta = lat_mean_dest - lat_mean_base
    return delta

def shift_adata_in_latent_space(adata_base,delta):
    shifted_adata = adata_base.copy()
    shifted_adata.X += delta
    return shifted_adata

def train_VAE(adata,batch_key,labels_key,n_latent,save_string,save = True):
    scvi.data.setup_anndata(adata,batch_key=batch_key,labels_key=labels_key)
    model = scgen.SCGEN(adata,n_latent = n_latent)
    model.train(max_epochs=100,batch_size=32,early_stopping=True,early_stopping_patience=25)
    if save:
        model.save(save_string, overwrite=True,save_anndata=True)
    return model
    

def make_compare_column_schulte(adata):
    compare_column = []
    for q in range(len(adata.obs)):
        compare_column.append(str(adata.obs['who_per_sample'][q]))
    adata.obs['compare_umap'] =  compare_column 
    
def get_min_max_plotting(df):
    min_val = min(df.values)[0]
    max_val = max(df.values)[0]
    vis_shift = 0.1*(max_val - min_val)
    return min_val-vis_shift,max_val+vis_shift


def pseudo_plot_3_to_7(df):
    plt.bar(['3', '4', '5', '6', '7'],df.iloc[0:5][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
    
def pseudo_plot_3_to_7_coh_2(df):
    plt.bar(['3', '4', '5', '7'],df.iloc[0:4][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
    
def basic_preprocessing(adata,min_genes=200,target_sum = 1e4):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
def prepare_hamsterMA_coh1(hamsterMA,schulte_cohort1):
    translate_hamsterMA_acc = do_uniform_name_accurate_hamsterMA_dict()
    translate_hamsterMA_overview = do_uniform_name_overview_hamsterMA_dict()
    add_uniform_name_accurate(adata = hamsterMA,
                                   translate_dict = translate_hamsterMA_acc,
                                   orig_column='celltype',
                                   new_column = 'uniform_name_acurate')

    add_uniform_name_accurate(adata = hamsterMA,
                                   translate_dict = translate_hamsterMA_overview,
                                   orig_column='celltype',
                                   new_column = 'uniform_name_overview')
    
    translate_schulte_cohort1_acc = do_uniform_name_accurate_schulte_dict()
    translate_schulte_cohort1_overview = do_uniform_name_overview_schulte_dict()
    add_uniform_name_accurate(adata = schulte_cohort1,
                                   translate_dict = translate_schulte_cohort1_acc,
                                   orig_column='id.celltype',
                                   new_column = 'uniform_name_acurate')

    add_uniform_name_accurate(adata = schulte_cohort1,
                                   translate_dict = translate_schulte_cohort1_overview,
                                   orig_column='id.celltype',
                                   new_column = 'uniform_name_overview')
    
    hamsterMA_genes = list(hamsterMA.var.index)
    human_cohort1_genes = list(schulte_cohort1.var.index)
    common_genes = list(Intersection(hamsterMA_genes,human_cohort1_genes))
    hamsterMA = hamsterMA[:,filter_for_common_genes(hamsterMA,common_genes)]
    hamsterMA.__dict__['_raw'].__dict__['_var'] = hamsterMA.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    schulte_cohort1 = schulte_cohort1[:,filter_for_common_genes(schulte_cohort1,common_genes)]
    return hamsterMA,schulte_cohort1


def make_yapeng_float(column):
    new_col = []
    for q in column:
        if q == '1 or 2':
            new_col.append(1.5)
        else:
            new_col.append(float(q))
    return new_col

def pseudo_plot_1_or2_to_7(df):
    plt.bar(['1 or 2', '3', '4', '5','6', '7'],df.iloc[0:6][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
    
def adapt_yapeng_compare_umap_2(adata):
    new_col = []
    for deg in adata.obs['compare_umap']:
        if ((deg == '1') or (deg == '2')):
            new_col.append('1 or 2')
        else:
            new_col.append(deg)
    adata.obs['compare_umap_2'] = new_col   
    
def cleanup(adata):
    try:
        del adata.obsm['X_diffmap']
    except:
        pass
    try:
        del adata.obsm['X_umap']
    except:
        pass
    try:
        del adata.obsm['X_pca']
    except:
        pass

def adapt_timepoint_hamsterPR(adata):
    tp_new = []
    for j in range(len(adata.obs['timepoint'])):
        if adata.obs['timepoint'][j] != 'hd_D2':
            tp_new.append(adata.obs['timepoint'][j])
        else:
            if adata.obs['hd_D2_cluster_1'][j] == 1:
                tp_new.append(adata.obs['timepoint'][j])
            else:
                tp_new.append('hd_D2_like_hd_D3')
    adata.obs['timepoint_new'] = tp_new
    
def pseudo_plot_hamsterPR_new(df):
    plt.bar(['d0', 'hd_d2','hd_D2_like_hd_D3','hd_d3', 'ld_d2', 'ld_d3'],df.iloc[0:6][df.columns.values[0]].values)
    plt.title(df.columns.values[0])
    
def make_compare_column_hamsterPR(adata):
    compare_column = []
    for q in range(len(adata.obs)):
        compare_column.append(str(adata.obs['timepoint_new'][q]))
    adata.obs['compare_umap'] =  compare_column 

def add_closeness_score_human(deg):
    dpts = deg.values[:,0]
    min_val = min(dpts)
    max_val = max(dpts)
    vis_shift = 0.1*(max_val - min_val)
    closeness = (dpts - min_val) + vis_shift
    closeness_score = 1-closeness/closeness.sum(0)
    deg['closeness_score'] = closeness_score
    
def add_closeness_score_PR(df_PR):
    dpts = df_PR.values[:,1]
    min_val = min(dpts)
    max_val = max(dpts)
    vis_shift = 0.1*(max_val - min_val)
    closeness = (dpts - min_val) + vis_shift
    closeness_score = 1-closeness/closeness.sum(0)
    df_PR['closeness_score'] = closeness_score
    
def make_vs_rest_column(adata,column,cov_degree):
    new_col = []
    for cov_deg in class_monos.obs['Who.Ordinal.Scale']:
        if cov_deg == cov_degree:
            new_col.append(cov_degree)
        else:
            new_col.append('rest')
    adata.obs[cov_degree + '_vs_rest'] = new_col
    
def add_rank(pt_array):
    sort = pt_array.argsort()
    ranks = np.arange(len(pt_array))[sort.argsort()] + 1
    return ranks

def PKF_score_PR(experiment_location):
    hd_D2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D2.csv')[0:4]
    hd_D3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D3.csv')[0:4]
    ld_D2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_ld_D2.csv')[0:4]
    ld_D3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_ld_D3.csv')[0:4]
    hd_D2['ranks'] = add_rank(hd_D2['hd_D2 shifted'])
    hd_D3['ranks'] = add_rank(hd_D3['hd_D3 shifted'])
    ld_D2['ranks'] = add_rank(ld_D2['ld_D2 shifted'])
    ld_D3['ranks'] = add_rank(ld_D3['ld_D3 shifted'])
    top_rank_hd_D2 = float(hd_D2[hd_D2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_hd_D3 = float(hd_D3[hd_D3['ranks'] == 1]['cov_degree'].values[0])
    top_rank_ld_D2 = float(hd_D2[ld_D2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_ld_D3 = float(hd_D3[ld_D3['ranks'] == 1]['cov_degree'].values[0])
    pkf_score = 0
    if top_rank_hd_D2 > 3:
        pkf_score += .1

    if top_rank_hd_D3 > 3:
        pkf_score += .1

    if top_rank_hd_D2 > 4:
        pkf_score += .1

    if top_rank_hd_D3 > 4:
        pkf_score += .1

    if top_rank_ld_D2 < 5:
        pkf_score += .1

    if top_rank_ld_D3 < 5:
        pkf_score += .1

    if max(top_rank_ld_D2,top_rank_ld_D3) <=  min(top_rank_hd_D2,top_rank_hd_D3):
        pkf_score += .1

    if max(top_rank_ld_D2,top_rank_ld_D3) <  min(top_rank_hd_D2,top_rank_hd_D3):
        pkf_score += .1

    if min(top_rank_ld_D2,top_rank_ld_D3) <  min(top_rank_hd_D2,top_rank_hd_D3):
        pkf_score += .1

    if (np.mean([top_rank_hd_D2,top_rank_hd_D3]) > np.mean([top_rank_ld_D2,top_rank_ld_D3])):
        pkf_score += .1
    return np.round(pkf_score,2)

def PKF_score_MA(experiment_location):
    d2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_d2.csv')[0:4]
    d3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_d3.csv')[0:4]
    d5 = pd.read_csv(experiment_location + '/figures/dpt_df/df_d5.csv')[0:4]
    e14 = pd.read_csv(experiment_location + '/figures/dpt_df/df_e14.csv')[0:4]
    d2['ranks'] = add_rank(d2['d2 shifted'])
    d3['ranks'] = add_rank(d3['d3 shifted'])
    d5['ranks'] = add_rank(d5['d5 shifted'])
    e14['ranks'] = add_rank(e14['e14 shifted'])
    top_rank_d2 = float(d2[d2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_d3 = float(d3[d3['ranks'] == 1]['cov_degree'].values[0])
    top_rank_d5 = float(d5[d5['ranks'] == 1]['cov_degree'].values[0])
    top_rank_e14 = float(e14[e14['ranks'] == 1]['cov_degree'].values[0])
    pkf_score = 0
    if top_rank_d2 > 3:
        pkf_score += .1
    if top_rank_d2 > 4:
        pkf_score += .1
    if top_rank_d3 > 3:
        pkf_score += .1
    if top_rank_d3 > 4:
        pkf_score += .1
    if top_rank_e14 < 5:
        pkf_score += .1
    if top_rank_e14 < 4:
        pkf_score += .1
    if max(top_rank_d2,top_rank_d3) > top_rank_e14:
        pkf_score += .1
    if min(top_rank_d2,top_rank_d3) > top_rank_e14:
        pkf_score += .1
    return np.round(pkf_score/0.8,2)

def add_pkf_score(experiment_location,hamster_type):
    #experiment_location: path to experiment result folder e.g '~/archive_experiment_folders/experiment_1_2'
    #hamster_type: str, 'MA' or 'PR'
    if hamster_type == 'PR':
        pkf_score = PKF_score_PR(experiment_location)
    if hamster_type == 'MA':
        pkf_score = PKF_score_MA(experiment_location)
    out = pd.DataFrame()
    out['pkf_score'] = [pkf_score]
    out.to_csv(experiment_location + '/description/pkf_score.csv')
    
def add_pkf_score_relative(experiment_location,hamster_type,class_mono_PR = False):
    #experiment_location: path to experiment result folder e.g '~/archive_experiment_folders/experiment_1_2'
    #hamster_type: str, 'ma' or 'pr'
    if hamster_type == 'pr':
        #if class_mono_PR == True:
        #    pkf_score = PKF_score_PR_relative_divided(experiment_location)
        #else:
        pkf_score = PKF_score_PR_relative(experiment_location)
    if hamster_type == 'ma':
        pkf_score = PKF_score_MA_relative(experiment_location)
    out = pd.DataFrame()
    out['pkf_score'] = [pkf_score]
    out.to_csv(experiment_location + '/description/pkf_score_relative.csv')
    
def PKF_score_PR_relative(experiment_location):
    hd_D2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D2.csv')[0:4]
    hd_D3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D3.csv')[0:4]
    ld_D2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_ld_D2.csv')[0:4]
    ld_D3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_ld_D3.csv')[0:4]
    hd_D2['ranks'] = add_rank(hd_D2['hd_D2 shifted'])
    hd_D3['ranks'] = add_rank(hd_D3['hd_D3 shifted'])
    ld_D2['ranks'] = add_rank(ld_D2['ld_D2 shifted'])
    ld_D3['ranks'] = add_rank(ld_D3['ld_D3 shifted'])
    top_rank_hd_D2 = float(hd_D2[hd_D2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_hd_D3 = float(hd_D3[hd_D3['ranks'] == 1]['cov_degree'].values[0])
    top_rank_ld_D2 = float(ld_D2[ld_D2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_ld_D3 = float(hd_D3[ld_D3['ranks'] == 1]['cov_degree'].values[0])
    pkf_score = 0
  

    if max(top_rank_ld_D2,top_rank_ld_D3) <=  min(top_rank_hd_D2,top_rank_hd_D3):
        pkf_score += .1

    if max(top_rank_ld_D2,top_rank_ld_D3) <  min(top_rank_hd_D2,top_rank_hd_D3):
        pkf_score += .1

    if min(top_rank_ld_D2,top_rank_ld_D3) <  min(top_rank_hd_D2,top_rank_hd_D3):
        pkf_score += .1

    if (np.mean([top_rank_hd_D2,top_rank_hd_D3]) > np.mean([top_rank_ld_D2,top_rank_ld_D3])):
        pkf_score += .1
    
    if (np.mean([top_rank_hd_D2,top_rank_hd_D3]) > np.mean([top_rank_ld_D2,top_rank_ld_D3])):
        pkf_score += .1
    
    if not (top_rank_hd_D2 == top_rank_hd_D3 == top_rank_ld_D2 == top_rank_ld_D3):
        pkf_score += .1
    
    
    
    return np.round(pkf_score/0.6,2)


def PKF_score_MA_relative(experiment_location):
    d2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_d2.csv')[0:4]
    d3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_d3.csv')[0:4]
    d5 = pd.read_csv(experiment_location + '/figures/dpt_df/df_d5.csv')[0:4]
    e14 = pd.read_csv(experiment_location + '/figures/dpt_df/df_e14.csv')[0:4]
    d2['ranks'] = add_rank(d2['d2 shifted'])
    d3['ranks'] = add_rank(d3['d3 shifted'])
    d5['ranks'] = add_rank(d5['d5 shifted'])
    e14['ranks'] = add_rank(e14['e14 shifted'])
    top_rank_d2 = float(d2[d2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_d3 = float(d3[d3['ranks'] == 1]['cov_degree'].values[0])
    top_rank_d5 = float(d5[d5['ranks'] == 1]['cov_degree'].values[0])
    top_rank_e14 = float(e14[e14['ranks'] == 1]['cov_degree'].values[0])
    pkf_score = 0
    if max(top_rank_d2,top_rank_d3) > top_rank_e14:
        pkf_score += .25
    if min(top_rank_d2,top_rank_d3) > top_rank_e14:
        pkf_score += .25
    if not (top_rank_d2 == top_rank_d3 == top_rank_e14):
        pkf_score += .25
    if not (top_rank_d2 == top_rank_d3):
        pkf_score += .25
    return pkf_score

def PKF_score_PR_relative_divided(experiment_location):
    hd_D2_outlier = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D2_outlier.csv')[0:4]
    hd_D2_standard = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D2_standard.csv')[0:4]
    hd_D3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_hd_D3.csv')[0:4]
    ld_D2 = pd.read_csv(experiment_location + '/figures/dpt_df/df_ld_D2.csv')[0:4]
    ld_D3 = pd.read_csv(experiment_location + '/figures/dpt_df/df_ld_D3.csv')[0:4]
    hd_D2_outlier['ranks'] = add_rank(hd_D2_outlier['hd_D2_outlier shifted'])
    hd_D2_standard['ranks'] = add_rank(hd_D2_standard['hd_D2_standard shifted'])
    hd_D3['ranks'] = add_rank(hd_D3['hd_D3 shifted'])
    ld_D2['ranks'] = add_rank(ld_D2['ld_D2 shifted'])
    ld_D3['ranks'] = add_rank(ld_D3['ld_D3 shifted'])
    top_rank_hd_D2_outlier = float(hd_D2_outlier[hd_D2_outlier['ranks'] == 1]['cov_degree'].values[0])
    top_rank_hd_D2_standard = float(hd_D2_standard[hd_D2_standard['ranks'] == 1]['cov_degree'].values[0])
    top_rank_hd_D3 = float(hd_D3[hd_D3['ranks'] == 1]['cov_degree'].values[0])
    top_rank_ld_D2 = float(ld_D2[ld_D2['ranks'] == 1]['cov_degree'].values[0])
    top_rank_ld_D3 = float(hd_D3[ld_D3['ranks'] == 1]['cov_degree'].values[0])
    pkf_score = 0
  

    if max(top_rank_ld_D2,top_rank_ld_D3) <=  min(top_rank_hd_D2_outlier,top_rank_hd_D2_standard,top_rank_hd_D3):
        pkf_score += .1

    if max(top_rank_ld_D2,top_rank_ld_D3) <  min(top_rank_hd_D2_outlier,top_rank_hd_D2_standard,top_rank_hd_D3):
        pkf_score += .1

    if min(top_rank_ld_D2,top_rank_ld_D3) <  min(top_rank_hd_D2_outlier,top_rank_hd_D2_standard,top_rank_hd_D3):
        pkf_score += .1

    if (np.mean([top_rank_hd_D2_outlier,top_rank_hd_D2_standard,top_rank_hd_D3]) > np.mean([top_rank_ld_D2,top_rank_ld_D3])):
        pkf_score += .1
    
    if (np.mean([top_rank_hd_D2_outlier,top_rank_hd_D2_standard,top_rank_hd_D3]) > np.mean([top_rank_ld_D2,top_rank_ld_D3])):
        pkf_score += .1
    
    if not (top_rank_hd_D2_outlier == top_rank_hd_D2_standard == top_rank_hd_D3 == top_rank_ld_D2 == top_rank_ld_D3):
        pkf_score += .1
    
    
    
    return np.round(pkf_score/0.6,2)


def make_translate_dict_plot():
    translate_dict_deg = {}
    translate_dict_deg[0] = '3'
    translate_dict_deg[1] = '4'
    translate_dict_deg[2] = '5'
    translate_dict_deg[3] = '7'
    return translate_dict_deg

def annotate_dwarf(hamsterPR,T_sub):
    translate_hamsterPR_acc = do_uniform_name_accurate_hamsterDwarf_dict()
    translate_hamsterPR_overview = do_uniform_name_overview_hamsterDwarf_dict()
    add_uniform_name_accurate(adata = hamsterPR,
                               translate_dict = translate_hamsterPR_acc,
                               orig_column='celltype',
                               new_column = 'uniform_name_acurate')

    add_uniform_name_accurate(adata = hamsterPR,
                                   translate_dict = translate_hamsterPR_overview,
                                   orig_column='celltype',
                                   new_column = 'uniform_name_overview')
    hamsterPR = subset_T_cells_PR(hamsterPR,T_sub)
    return hamsterPR

def annotate_golden(hamsterMA,T_sub):
    translate_hamsterMA_acc = do_uniform_name_accurate_hamsterMA_dict()
    translate_hamsterMA_overview = do_uniform_name_overview_hamsterMA_dict()
    add_uniform_name_accurate(adata = hamsterMA,
                               translate_dict = translate_hamsterMA_acc,
                               orig_column='celltype',
                               new_column = 'uniform_name_acurate')

    add_uniform_name_accurate(adata = hamsterMA,
                                   translate_dict = translate_hamsterMA_overview,
                                   orig_column='celltype',
                                   new_column = 'uniform_name_overview')
    hamsterMA = subset_T_cells_MA(hamsterMA,T_sub)
    return hamsterMA
    
def subset_T_cells_MA(hamsterMA,T_sub):
    T_MA = filter_adata_obs(T_sub,'orig.ident','ma')
    new_cell_label = []
    for q in range(len(hamsterMA.obs)):
        ct = hamsterMA.obs['uniform_name_acurate'][q]
        if ct == 'T_cells':
            ind = hamsterMA.obs.index[q]
            try:
                new_label = T_MA[T_MA.obs.index == ind].obs['T_cell_Hamster_acurate'].values[0]
                if new_label == 'Neutrophils':
                    new_cell_label.append('T_neutro')
                else:
                    new_cell_label.append(new_label)
            except:
                new_cell_label.append('undefined')
        else:
            new_cell_label.append(ct)
    hamsterMA.obs['uniform_name_acurate'] = new_cell_label
    return hamsterMA

def subset_T_cells_PR(hamsterPR,T_sub):
    T_PR = filter_adata_obs(T_sub,'orig.ident','pr')
    new_cell_label = []
    for q in range(len(hamsterPR.obs)):
        ct = hamsterPR.obs['uniform_name_acurate'][q]
        if ct == 'T_cells':
            ind = hamsterPR.obs.index[q]
            try:
                new_label = T_PR[T_PR.obs.index == ind].obs['T_cell_Hamster_acurate'].values[0]
                if new_label == 'Neutrophils':
                    new_cell_label.append('T_neutro')
                else:
                    new_cell_label.append(new_label)
            except:
                new_cell_label.append('undefined')
        else:
            new_cell_label.append(ct)
    hamsterPR.obs['uniform_name_acurate'] = new_cell_label
    return hamsterPR

def annotate_cohort2(schulte_cohort2):
    cohort2_WB = filter_adata_obs(schulte_cohort2,'cells','Whole_blood')
    translate_schulte_cohort2_acc = do_uniform_name_accurate_schulteBD_dict()
    translate_schulte_cohort2_overview = do_uniform_name_overview_schulteBD_dict()
    add_uniform_name_accurate(adata = cohort2_WB,
                                   translate_dict = translate_schulte_cohort2_acc,
                                   orig_column='cluster_labels_res.0.8',
                                   new_column = 'uniform_name_acurate')

    add_uniform_name_accurate(adata = cohort2_WB,
                                   translate_dict = translate_schulte_cohort2_overview,
                                   orig_column='cluster_labels_res.0.8',
                                   new_column = 'uniform_name_overview')
    return cohort2_WB

def golden_human_preprocess(hamsterMA,cohort2_WB):
    sc.pp.filter_cells(cohort2_WB, min_genes=200)
    sc.pp.filter_genes(cohort2_WB, min_cells=20)
    sc.pp.filter_cells(hamsterMA, min_genes=200)
    sc.pp.filter_genes(hamsterMA, min_cells=20)
    
    hamsterMA_genes = list(hamsterMA.var.index)
    human_cohort2_genes = list(cohort2_WB.var.index)
    common_genes = list(Intersection(hamsterMA_genes,human_cohort2_genes))
    cohort2_WB = cohort2_WB[:,common_genes]
    cohort2_WB.__dict__['_raw'].__dict__['_var'] = cohort2_WB.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    hamsterMA = hamsterMA[:,common_genes]
    hamsterMA.__dict__['_raw'].__dict__['_var'] = hamsterMA.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    
    sc.pp.normalize_total(cohort2_WB, target_sum=1e4,exclude_highly_expressed = True)
    sc.pp.log1p(cohort2_WB)
    
    sc.pp.normalize_total(hamsterMA, target_sum=1e4,exclude_highly_expressed = True)
    sc.pp.log1p(hamsterMA)
    return hamsterMA,cohort2_WB

def dwarf_human_preprocess(hamsterPR,cohort2_WB):
    sc.pp.filter_cells(cohort2_WB, min_genes=200)
    sc.pp.filter_genes(cohort2_WB, min_cells=20)
    
    sc.pp.filter_cells(hamsterPR, min_genes=200)
    sc.pp.filter_genes(hamsterPR, min_cells=20)
    
    hamsterPR_genes = list(hamsterPR.var.index)
    human_cohort2_genes = list(cohort2_WB.var.index)
    common_genes = list(Intersection(hamsterPR_genes,human_cohort2_genes))
    cohort2_WB = cohort2_WB[:,common_genes]
    cohort2_WB.__dict__['_raw'].__dict__['_var'] = cohort2_WB.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    hamsterPR = hamsterPR[:,common_genes]
    hamsterPR.__dict__['_raw'].__dict__['_var'] = hamsterPR.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

    sc.pp.normalize_total(cohort2_WB, target_sum=1e4,exclude_highly_expressed = True)
    sc.pp.log1p(cohort2_WB)

    sc.pp.normalize_total(hamsterPR, target_sum=1e4,exclude_highly_expressed = True)
    sc.pp.log1p(hamsterPR)
    return hamsterPR,cohort2_WB

def preprocess(hamster_rawfile,human_rawfile,T_sub,
               species = 'ma',
               celltype = 'Classical_Monocytes',
               create_folder = True,
               save_input = True,
               save_location = None,
               experiment_name = None,
              do_batch_removal_human = False,
              do_batch_removal_human_full_set = False,
              mito_filtering = False,
              plot_qc_metics = False):
    #species 'ma' for golden hamster, 'pr' for dwarf hamster
    if (do_batch_removal_human & do_batch_removal_human_full_set):
        raise Exception("Decide on batch removal") 
    if species == 'ma':
        hamster = annotate_golden(hamster_rawfile,T_sub)
        human = annotate_cohort2(human_rawfile)
        if mito_filtering:
            hamster = remove_mito(hamster,plot_qc_metics = plot_qc_metics)
            human = remove_mito(human,plot_qc_metics = plot_qc_metics)
        hamster,human = golden_human_preprocess(hamster,human)
        
    if species == 'pr':
        hamster = annotate_dwarf(hamster_rawfile,T_sub)
        human = annotate_cohort2(human_rawfile)
        if mito_filtering:
            hamster = remove_mito(hamster,plot_qc_metics = plot_qc_metics)
            human = remove_mito(human,plot_qc_metics = plot_qc_metics)
        hamster,human = dwarf_human_preprocess(hamster,human)
        
    if do_batch_removal_human_full_set:
        sc.pp.combat(human,'orig.ident')
        
    if ((species == 'ma') & (celltype == 'Classical_Monocytes')):
        hamster = filter_adata_obs(hamster,'uniform_name_acurate','Classical_Monocytes')
        hamster = leave_out_adata_obs(hamster,'seurat_clusters',15)
    else:
        hamster = filter_adata_obs(hamster,'uniform_name_acurate',celltype)
        
    if celltype == 'B_cells':
        human = filter_adata_obs(human,'uniform_name_acurate','B_Cells')
    else:
        human = filter_adata_obs(human,'uniform_name_acurate',celltype)
    if create_folder:
        if (save_location != None) & (experiment_name != None):
            create_folder_structure(save_location,experiment_name)
        else:
            print('No save location or experiment name given!')
    if do_batch_removal_human:
        sc.pp.combat(human,'orig.ident')
    if save_input:
        if save_location != None:
            hamster.write_h5ad(save_location + '/' + str(experiment_name) + '/data/hamster.h5ad')
            human.write_h5ad(save_location + '/' + str(experiment_name) + '/data/human.h5ad')
    return hamster,human
    
def create_folder_structure(location,experiment_name):
    os.mkdir(location + '/' + str(experiment_name))
    os.mkdir(location + '/' + str(experiment_name) + '/model')
    os.mkdir(location + '/' + str(experiment_name) + '/data')
    os.mkdir(location + '/' + str(experiment_name) + '/figures')
    os.mkdir(location + '/' + str(experiment_name) + '/figures/dpt_df')
    os.mkdir(location + '/' + str(experiment_name) + '/figures/train_figures')
    os.mkdir(location + '/' + str(experiment_name) + '/figures/dpt_figures')
    os.mkdir(location + '/' + str(experiment_name) + '/figures/LISI_figures')
    os.mkdir(location + '/' + str(experiment_name) + '/figures/heatmap_figures')
    os.mkdir(location + '/' + str(experiment_name) + '/figures/umaps')
    os.mkdir(location + '/' + str(experiment_name) + '/scripts')
    os.mkdir(location + '/' + str(experiment_name) + '/description')
    
    
def train_autoencoder_model_PR(hamster,human,
                            celltype = 'Classical_Monocytes',
                            save_input = True,
                            save_location = None,
                            experiment_name = None,
                           gene_set = 'full',
                           make_umaps = False,
                           batch_key = 'compare_umap',
                           labels_key = 'uniform_name_acurate',
                           n_latent = 5):
    #gene_set = 'full' or 'union_hv' or 'intersection_hv'
    loc_dpt_df = save_location + '/' + str(experiment_name) + '/figures/dpt_df'
    loc_dpt_figures = save_location + '/' + str(experiment_name) + '/figures/dpt_figures'
    loc_heatmap = save_location + '/' + str(experiment_name) + '/figures/heatmap_figures'
    loc_umaps = save_location + '/' + str(experiment_name) + '/figures/umaps'
    model_save_string = save_location + '/' + str(experiment_name) + '/model/model.pt'
    
    make_compare_column(hamster)
    make_compare_column_schulte(human)
    
    if gene_set == 'union_hv':
        sc.pp.highly_variable_genes(hamster,layer = 'log_counts')
        sc.pp.highly_variable_genes(human,layer = 'log_counts')
        hvg_hamster = list(hamster.var['highly_variable'][hamster.var['highly_variable'] == True].index)
        hvg_human = list(human.var['highly_variable'][human.var['highly_variable'] == True].index)
        hvg = list(set(hvg_hamster).union(hvg_human))
        hamster = hamster[:,hvg]
        human = human[:,hvg]
    
    if gene_set == 'intersection_hv':
        sc.pp.highly_variable_genes(hamster,layer = 'log_counts')
        sc.pp.highly_variable_genes(human,layer = 'log_counts')
        hvg_hamster = list(hamster.var['highly_variable'][hamster.var['highly_variable'] == True].index)
        hvg_human = list(human.var['highly_variable'][human.var['highly_variable'] == True].index)
        hvg = list(set(hvg_hamster).intersection(hvg_human))
        hamster = hamster[:,hvg]
        human = human[:,hvg]
        
    if make_umaps:
        prepare_umap(hamster)
        sc.pl.umap(hamster,color = 'compare_umap',save =  loc_umaps + '/hamster_GE_space.pdf')

        prepare_umap(human)
        sc.pl.umap(human,color = 'compare_umap',save = loc_umaps + '/human_GE_space.pdf')

    merged = merge_adata(human,hamster)

    if make_umaps:
        prepare_umap(merged)
        sc.pl.umap(merged,color = 'compare_umap',save = loc_umaps + '/merged_GE_space.pdf')
    try:
        del merged.obsm['X_diffmap']
    except:
        pass
    try:
        del merged.obsm['X_umap']
    except:
        pass
    try:
        del merged.obsm['X_pca']
    except:
        pass
    
    model = train_VAE(merged,batch_key,labels_key,n_latent,save_string = model_save_string)
    data = model.adata.copy()
    lat_data = get_latent_representation_object(model,data)
    
    if make_umaps:
        prepare_umap(lat_data)
        sc.pl.umap(lat_data,color = 'compare_umap',save = loc_umaps + '/merged_lat_space.pdf')
        
    lat_hamster = filter_adata_obs(lat_data,'dataset','hamsterPR')
    lat_human = leave_out_adata_obs(lat_data,'dataset','hamsterPR')

    lat_hamster_d0 = filter_adata_obs(lat_hamster,'compare_umap','D0')
    lat_hamster_infectious = leave_out_adata_obs(lat_hamster,'compare_umap','D0')
    lat_human_0 = filter_adata_obs(lat_human,'compare_umap','0')
    lat_human_infectious = leave_out_adata_obs(lat_human,'compare_umap','0')

    delta = get_delta_in_latent_space(lat_hamster_d0,lat_human_0) 


    lat_hamster_hd_D2 = filter_adata_obs(lat_hamster,'timepoint','hd_D2')
    lat_hamster_hd_D2_shifted_control = shift_adata_in_latent_space(lat_hamster_hd_D2,delta)
    lat_hamster_hd_D2_shifted_control.obs['compare_umap'] = 'hd_D2 shifted'
    hd_D2_shift_human_infectious = merge_adata(lat_hamster_hd_D2_shifted_control,lat_human_infectious)
    df_hd_D2 = calculate_pseudotime(hd_D2_shift_human_infectious,10,'compare_umap','hd_D2 shifted')
    df_hd_D2.to_csv(loc_dpt_df + '/df_hd_D2.csv',index_label = 'cov_degree')

    lat_hamster_hd_D3 = filter_adata_obs(lat_hamster,'timepoint','hd_D3')
    lat_hamster_hd_D3_shifted_control = shift_adata_in_latent_space(lat_hamster_hd_D3,delta)
    lat_hamster_hd_D3_shifted_control.obs['compare_umap'] = 'hd_D3 shifted'
    hd_D3_shift_human_infectious = merge_adata(lat_hamster_hd_D3_shifted_control,lat_human_infectious)
    df_hd_D3 = calculate_pseudotime(hd_D3_shift_human_infectious,10,'compare_umap','hd_D3 shifted')
    df_hd_D3.to_csv(loc_dpt_df + '/df_hd_D3.csv',index_label = 'cov_degree')

    lat_hamster_ld_D3 = filter_adata_obs(lat_hamster,'timepoint','ld_D3')
    lat_hamster_ld_D3_shifted_control = shift_adata_in_latent_space(lat_hamster_ld_D3,delta)
    lat_hamster_ld_D3_shifted_control.obs['compare_umap'] = 'ld_D3 shifted'
    ld_D3_shift_human_infectious = merge_adata(lat_hamster_ld_D3_shifted_control,lat_human_infectious)
    df_ld_D3 = calculate_pseudotime(ld_D3_shift_human_infectious,10,'compare_umap','ld_D3 shifted')
    df_ld_D3.to_csv(loc_dpt_df + '/df_ld_D3.csv',index_label = 'cov_degree')

    lat_hamster_ld_D2 = filter_adata_obs(lat_hamster,'timepoint','ld_D2')
    lat_hamster_ld_D2_shifted_control = shift_adata_in_latent_space(lat_hamster_ld_D2,delta)
    lat_hamster_ld_D2_shifted_control.obs['compare_umap'] = 'ld_D2 shifted'
    ld_D2_shift_human_infectious = merge_adata(lat_hamster_ld_D2_shifted_control,lat_human_infectious)
    df_ld_D2 = calculate_pseudotime(ld_D2_shift_human_infectious,10,'compare_umap','ld_D2 shifted')
    df_ld_D2.to_csv(loc_dpt_df + '/df_ld_D2.csv',index_label = 'cov_degree')

    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_hd_D2[0:4])
    plt.ylim(get_min_max_plotting(df_hd_D2[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_hd_D2.pdf')
    plt.close()
   

    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_hd_D3[0:4])
    plt.ylim(get_min_max_plotting(df_hd_D3[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_hd_D3.pdf')
    plt.close()


    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_ld_D2[0:4])
    plt.ylim(get_min_max_plotting(df_ld_D2[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_ld_D2.pdf')
    plt.close()


    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_ld_D3[0:4])
    plt.ylim(get_min_max_plotting(df_ld_D3[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_ld_D3.pdf')
    plt.close()
    
    df_hd_D2 = pd.read_csv(loc_dpt_df + '/df_hd_D2.csv')[0:4]
    add_closeness_score_PR(df_hd_D2)
    df_hd_D2['day'] = 'hd D2'

    df_hd_D3 = pd.read_csv(loc_dpt_df + '/df_hd_D3.csv')[0:4]
    add_closeness_score_PR(df_hd_D3)
    df_hd_D3['day'] = 'hd D3'

    df_ld_D2 = pd.read_csv(loc_dpt_df + '/df_ld_D2.csv')[0:4]
    add_closeness_score_PR(df_ld_D2)
    df_ld_D2['day'] = 'ld D2'

    df_ld_D3 = pd.read_csv(loc_dpt_df + '/df_ld_D3.csv')[0:4]
    add_closeness_score_PR(df_ld_D3)
    df_ld_D3['day'] = 'ld D3'


    df_hamsterPR = df_hd_D2.append(df_hd_D3).append(df_ld_D2).append(df_ld_D3)
    df_hamsterPR_pivot = df_hamsterPR.pivot("day", "cov_degree", "closeness_score").astype(float)

    plt.figure()
    sns.heatmap(df_hamsterPR_pivot,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('WHO degree',fontsize = 12)
    plt.ylabel('dose / timepoint',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_heatmap + '/heatmap.pdf')
    plt.close()    
    
    
def train_autoencoder_model_MA(hamster,human,celltype = 'Classical_Monocytes',
                            save_input = True,
                            save_location = None,
                            experiment_name = None,
                           gene_set = 'full',
                           make_umaps = False,
                           batch_key = 'compare_umap',
                           labels_key = 'uniform_name_acurate',
                           n_latent = 5):
    #gene_set = 'full' or 'union_hv' or 'intersection_hv'
    loc_dpt_df = save_location + '/' + str(experiment_name) + '/figures/dpt_df'
    loc_dpt_figures = save_location + '/' + str(experiment_name) + '/figures/dpt_figures'
    loc_heatmap = save_location + '/' + str(experiment_name) + '/figures/heatmap_figures'
    loc_umaps = save_location + '/' + str(experiment_name) + '/figures/umaps'
    model_save_string = save_location + '/' + str(experiment_name) + '/model/model.pt'
    
    make_compare_column(hamster)
    make_compare_column_schulte(human)
    
    if gene_set == 'union_hv':
        sc.pp.highly_variable_genes(hamster,layer = 'log_counts')
        sc.pp.highly_variable_genes(human,layer = 'log_counts')
        hvg_hamster = list(hamster.var['highly_variable'][hamster.var['highly_variable'] == True].index)
        hvg_human = list(human.var['highly_variable'][human.var['highly_variable'] == True].index)
        hvg = list(set(hvg_hamster).union(hvg_human))
        hamster = hamster[:,hvg]
        human = human[:,hvg]
    
    if gene_set == 'intersection_hv':
        sc.pp.highly_variable_genes(hamster,layer = 'log_counts')
        sc.pp.highly_variable_genes(human,layer = 'log_counts')
        hvg_hamster = list(hamster.var['highly_variable'][hamster.var['highly_variable'] == True].index)
        hvg_human = list(human.var['highly_variable'][human.var['highly_variable'] == True].index)
        hvg = list(set(hvg_hamster).intersection(hvg_human))
        hamster = hamster[:,hvg]
        human = human[:,hvg]
        
    if make_umaps:
        prepare_umap(hamster)
        sc.pl.umap(hamster,color = 'compare_umap',save =  loc_umaps + '/hamster_GE_space.pdf')

        prepare_umap(human)
        sc.pl.umap(human,color = 'compare_umap',save = loc_umaps + '/human_GE_space.pdf')

    merged = merge_adata(human,hamster)

    if make_umaps:
        prepare_umap(merged)
        sc.pl.umap(merged,color = 'compare_umap',save = loc_umaps + '/merged_GE_space.pdf')
    try:
        del merged.obsm['X_diffmap']
    except:
        pass
    try:
        del merged.obsm['X_umap']
    except:
        pass
    try:
        del merged.obsm['X_pca']
    except:
        pass
    
    model = train_VAE(merged,batch_key,labels_key,n_latent,save_string = model_save_string)
    data = model.adata.copy()
    lat_data = get_latent_representation_object(model,data)
    
    if make_umaps:
        prepare_umap(lat_data)
        sc.pl.umap(lat_data,color = 'compare_umap',save = loc_umaps + '/merged_lat_space.pdf')
        
    lat_hamster = filter_adata_obs(lat_data,'dataset','hamsterMA')
    lat_human = leave_out_adata_obs(lat_data,'dataset','hamsterMA')

    lat_hamster_d0 = filter_adata_obs(lat_hamster,'compare_umap','d0')
    lat_hamster_infectious = leave_out_adata_obs(lat_hamster,'compare_umap','d0')
    lat_human_0 = filter_adata_obs(lat_human,'compare_umap','0')
    lat_human_infectious = leave_out_adata_obs(lat_human,'compare_umap','0')

    delta = get_delta_in_latent_space(lat_hamster_d0,lat_human_0) 


    lat_hamster_d2 = filter_adata_obs(lat_hamster,'timepoint','d2')
    lat_hamster_d2_shifted_control = shift_adata_in_latent_space(lat_hamster_d2,delta)
    lat_hamster_d2_shifted_control.obs['compare_umap'] = 'd2 shifted'
    d2_shift_human_infectious = merge_adata(lat_hamster_d2_shifted_control,lat_human_infectious)
    df_d2 = calculate_pseudotime(d2_shift_human_infectious,10,'compare_umap','d2 shifted')
    df_d2.to_csv(loc_dpt_df + '/df_d2.csv',index_label = 'cov_degree')

    lat_hamster_d3 = filter_adata_obs(lat_hamster,'timepoint','d3')
    lat_hamster_d3_shifted_control = shift_adata_in_latent_space(lat_hamster_d3,delta)
    lat_hamster_d3_shifted_control.obs['compare_umap'] = 'd3 shifted'
    d3_shift_human_infectious = merge_adata(lat_hamster_d3_shifted_control,lat_human_infectious)
    df_d3 = calculate_pseudotime(d3_shift_human_infectious,10,'compare_umap','d3 shifted')
    df_d3.to_csv(loc_dpt_df + '/df_d3.csv',index_label = 'cov_degree')

    lat_hamster_d5 = filter_adata_obs(lat_hamster,'timepoint','d5')
    lat_hamster_d5_shifted_control = shift_adata_in_latent_space(lat_hamster_d5,delta)
    lat_hamster_d5_shifted_control.obs['compare_umap'] = 'd5 shifted'
    d5_shift_human_infectious = merge_adata(lat_hamster_d5_shifted_control,lat_human_infectious)
    df_d5 = calculate_pseudotime(d5_shift_human_infectious,10,'compare_umap','d5 shifted')
    df_d5.to_csv(loc_dpt_df + '/df_d5.csv',index_label = 'cov_degree')

    lat_hamster_e14 = filter_adata_obs(lat_hamster,'timepoint','e14')
    lat_hamster_e14_shifted_control = shift_adata_in_latent_space(lat_hamster_e14,delta)
    lat_hamster_e14_shifted_control.obs['compare_umap'] = 'e14 shifted'
    e14_shift_human_infectious = merge_adata(lat_hamster_e14_shifted_control,lat_human_infectious)
    df_e14 = calculate_pseudotime(e14_shift_human_infectious,10,'compare_umap','e14 shifted')
    df_e14.to_csv(loc_dpt_df + '/df_e14.csv',index_label = 'cov_degree')

    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_d2[0:4])
    plt.ylim(get_min_max_plotting(df_d2[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_d2.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_d3[0:4])
    plt.ylim(get_min_max_plotting(df_d3[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_d3.pdf')
    plt.close()


    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_d5[0:4])
    plt.ylim(get_min_max_plotting(df_d5[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_d5.pdf')
    plt.close()


    plt.figure()
    pseudo_plot_3_to_7_coh_2(df_e14[0:4])
    plt.ylim(get_min_max_plotting(df_e14[0:4]))
    plt.savefig(loc_dpt_figures + '/dpt_e14.pdf')
    plt.close()
    
    df_d2 = pd.read_csv(loc_dpt_df + '/df_d2.csv')[0:4]
    add_closeness_score_PR(df_d2)
    df_d2['day'] = 'd2'

    df_d3 = pd.read_csv(loc_dpt_df + '/df_d3.csv')[0:4]
    add_closeness_score_PR(df_d3)
    df_d3['day'] = 'd3'

    df_d5 = pd.read_csv(loc_dpt_df + '/df_d5.csv')[0:4]
    add_closeness_score_PR(df_d5)
    df_d5['day'] = 'd5'

    df_e14 = pd.read_csv(loc_dpt_df + '/df_e14.csv')[0:4]
    add_closeness_score_PR(df_e14)
    df_e14['day'] = 'e14'


    df_hamster = df_d2.append(df_d3).append(df_d5).append(df_e14)
    df_hamster_pivot = df_hamster.pivot("day", "cov_degree", "closeness_score").astype(float)

    plt.figure()
    sns.heatmap(df_hamster_pivot,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('WHO degree',fontsize = 12)
    plt.ylabel('dose / timepoint',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_heatmap + '/heatmap.pdf')
    plt.close()
    
def reverse_mappingMA(celltype = 'Classical_Monocytes',save_location = None,experiment_name = None):
    model_save_string = save_location + '/' + str(experiment_name) + '/model/model.pt'
    loc_reverse_df = save_location + '/' + str(experiment_name) + '/figures/dpt_df'
    loc_reverse_figures = save_location + '/' + str(experiment_name) + '/figures/dpt_figures'
    loc_reverse_heatmap = save_location + '/' + str(experiment_name) + '/figures/heatmap_figures'
    model = scgen.SCGEN.load(model_save_string)
    data = model.adata.copy()
    lat_data = get_latent_representation_object(model,data)
    lat_hamster = filter_adata_obs(lat_data,'dataset','hamsterMA')
    lat_human = leave_out_adata_obs(lat_data,'dataset','hamsterMA')

    lat_hamster_d0 = filter_adata_obs(lat_hamster,'compare_umap','d0')
    lat_hamster_infectious = leave_out_adata_obs(lat_hamster,'compare_umap','d0')
    lat_human_0 = filter_adata_obs(lat_human,'compare_umap','0')
    lat_human_infectious = leave_out_adata_obs(lat_human,'compare_umap','0')

    delta = get_delta_in_latent_space(lat_hamster_d0,lat_human_0)
    lat_hamster_infectious = leave_out_adata_obs(lat_hamster,'compare_umap','d0')
    lat_hamster_infected_shifted = shift_adata_in_latent_space(lat_hamster_infectious,delta)

    lat_human_3 = filter_adata_obs(lat_human,'who_per_sample','3')
    human_3_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_3)
    df_3 = calculate_pseudotime(human_3_infected_shifted,10,'compare_umap','3')[1:]

    lat_human_4 = filter_adata_obs(lat_human,'who_per_sample','4')
    human_4_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_4)
    df_4 = calculate_pseudotime(human_4_infected_shifted,10,'compare_umap','4')[1:]

    lat_human_5 = filter_adata_obs(lat_human,'who_per_sample','5')
    human_5_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_5)
    df_5 = calculate_pseudotime(human_5_infected_shifted,10,'compare_umap','5')[1:]

    lat_human_7 = filter_adata_obs(lat_human,'who_per_sample','7')
    human_7_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_7)
    df_7 = calculate_pseudotime(human_7_infected_shifted,10,'compare_umap','7')[1:]

    df_3.to_csv(loc_reverse_df + '/df_3.csv',index_label = 'timepoint')
    df_4.to_csv(loc_reverse_df + '/df_4.csv',index_label = 'timepoint')
    df_5.to_csv(loc_reverse_df + '/df_5.csv',index_label = 'timepoint')
    df_7.to_csv(loc_reverse_df + '/df_7.csv',index_label = 'timepoint')

    plt.figure()
    pseudo_plot_hamster_infected(df_3)
    plt.ylim(get_min_max_plotting(df_3))
    plt.savefig(loc_reverse_figures + '/dpt_3.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_hamster_infected(df_4)
    plt.ylim(get_min_max_plotting(df_4))
    plt.savefig(loc_reverse_figures + '/dpt_4.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_hamster_infected(df_5)
    plt.ylim(get_min_max_plotting(df_5))
    plt.savefig(loc_reverse_figures + '/dpt_5.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_hamster_infected(df_7)
    plt.ylim(get_min_max_plotting(df_7))
    plt.savefig(loc_reverse_figures + '/dpt_7.pdf')
    plt.close()

    df_3 = pd.read_csv(loc_reverse_df + '/df_3.csv')
    add_closeness_score_PR(df_3)
    df_3['WHO'] = '3'

    df_4 = pd.read_csv(loc_reverse_df + '/df_4.csv')
    add_closeness_score_PR(df_4)
    df_4['WHO'] = '4'

    df_5 = pd.read_csv(loc_reverse_df + '/df_5.csv')
    add_closeness_score_PR(df_5)
    df_5['WHO'] = '5'

    df_7 = pd.read_csv(loc_reverse_df + '/df_7.csv')
    add_closeness_score_PR(df_7)
    df_7['WHO'] = '7'

    df_all = df_3.append(df_4).append(df_5).append(df_7)
    df_all_pivot = df_all.pivot("WHO", "timepoint", "closeness_score").astype(float)

    plt.figure()
    sns.heatmap(df_all_pivot,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('timepoint',fontsize = 12)
    plt.ylabel('WHO degree',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_reverse_heatmap + '/reverse_heatmap.pdf')
    plt.close()  

    df_all = df_3.append(df_4).append(df_5).append(df_7)
    df_all_pivot = df_all.pivot("timepoint","WHO", "closeness_score").astype(float)

    plt.figure()
    sns.heatmap(df_all_pivot,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('timepoint',fontsize = 12)
    plt.ylabel('WHO degree',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_reverse_heatmap + '/reverse_heatmap_columnwise.pdf')
    plt.close()  
    
def reverse_mappingPR(celltype = 'Classical_Monocytes',save_location = None,experiment_name = None):
    model_save_string = save_location + '/' + str(experiment_name) + '/model/model.pt'
    loc_reverse_df = save_location + '/' + str(experiment_name) + '/figures/dpt_df'
    loc_reverse_figures = save_location + '/' + str(experiment_name) + '/figures/dpt_figures'
    loc_reverse_heatmap = save_location + '/' + str(experiment_name) + '/figures/heatmap_figures'
    model = scgen.SCGEN.load(model_save_string)
    data = model.adata.copy()
    lat_data = get_latent_representation_object(model,data)
    lat_hamster = filter_adata_obs(lat_data,'dataset','hamsterPR')
    lat_human = leave_out_adata_obs(lat_data,'dataset','hamsterPR')

    lat_hamster_d0 = filter_adata_obs(lat_hamster,'compare_umap','D0')
    lat_hamster_infectious = leave_out_adata_obs(lat_hamster,'compare_umap','D0')
    lat_human_0 = filter_adata_obs(lat_human,'compare_umap','0')
    lat_human_infectious = leave_out_adata_obs(lat_human,'compare_umap','0')

    delta = get_delta_in_latent_space(lat_hamster_d0,lat_human_0) 
    lat_hamster_infectious = leave_out_adata_obs(lat_hamster,'compare_umap','D0')
    lat_hamster_infected_shifted = shift_adata_in_latent_space(lat_hamster_infectious,delta)

    lat_human_3 = filter_adata_obs(lat_human,'who_per_sample','3')
    human_3_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_3)
    df_3 = calculate_pseudotime(human_3_infected_shifted,10,'compare_umap','3')[1:]

    lat_human_4 = filter_adata_obs(lat_human,'who_per_sample','4')
    human_4_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_4)
    df_4 = calculate_pseudotime(human_4_infected_shifted,10,'compare_umap','4')[1:]

    lat_human_5 = filter_adata_obs(lat_human,'who_per_sample','5')
    human_5_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_5)
    df_5 = calculate_pseudotime(human_5_infected_shifted,10,'compare_umap','5')[1:]

    lat_human_7 = filter_adata_obs(lat_human,'who_per_sample','7')
    human_7_infected_shifted = merge_adata(lat_hamster_infected_shifted,lat_human_7)
    df_7 = calculate_pseudotime(human_7_infected_shifted,10,'compare_umap','7')[1:]

    df_3.to_csv(loc_reverse_df + '/df_3.csv',index_label = 'timepoint')
    df_4.to_csv(loc_reverse_df + '/df_4.csv',index_label = 'timepoint')
    df_5.to_csv(loc_reverse_df + '/df_5.csv',index_label = 'timepoint')
    df_7.to_csv(loc_reverse_df + '/df_7.csv',index_label = 'timepoint')

    plt.figure()
    pseudo_plot_hamsterPR_infected(df_3)
    plt.ylim(get_min_max_plotting(df_3))
    plt.savefig(loc_reverse_figures + '/dpt_3.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_hamsterPR_infected(df_4)
    plt.ylim(get_min_max_plotting(df_4))
    plt.savefig(loc_reverse_figures + '/dpt_4.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_hamsterPR_infected(df_5)
    plt.ylim(get_min_max_plotting(df_5))
    plt.savefig(loc_reverse_figures + '/dpt_5.pdf')
    plt.close()

    plt.figure()
    pseudo_plot_hamsterPR_infected(df_7)
    plt.ylim(get_min_max_plotting(df_7))
    plt.savefig(loc_reverse_figures + '/dpt_7.pdf')
    plt.close()

    df_3 = pd.read_csv(loc_reverse_df + '/df_3.csv')
    add_closeness_score_PR(df_3)
    df_3['WHO'] = '3'

    df_4 = pd.read_csv(loc_reverse_df + '/df_4.csv')
    add_closeness_score_PR(df_4)
    df_4['WHO'] = '4'

    df_5 = pd.read_csv(loc_reverse_df + '/df_5.csv')
    add_closeness_score_PR(df_5)
    df_5['WHO'] = '5'

    df_7 = pd.read_csv(loc_reverse_df + '/df_7.csv')
    add_closeness_score_PR(df_7)
    df_7['WHO'] = '7'

    df_all = df_3.append(df_4).append(df_5).append(df_7)
    df_all_pivot = df_all.pivot("WHO", "timepoint", "closeness_score").astype(float)

    plt.figure()
    sns.heatmap(df_all_pivot,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('timepoint',fontsize = 12)
    plt.ylabel('WHO degree',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_reverse_heatmap + '/reverse_heatmap.pdf')
    plt.close()  

    df_all = df_3.append(df_4).append(df_5).append(df_7)
    df_all_pivot = df_all.pivot("timepoint","WHO", "closeness_score").astype(float)

    plt.figure()
    sns.heatmap(df_all_pivot,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('timepoint',fontsize = 12)
    plt.ylabel('WHO degree',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_reverse_heatmap + '/reverse_heatmap_columnwise.pdf')
    plt.close()  
    
    
def undirected_mappingMA(celltype = 'Classical_Monocytes',save_location = None,experiment_name = None):
    model_save_string = save_location + '/' + str(experiment_name) + '/model/model.pt'
    loc_reverse_df = save_location + '/' + str(experiment_name) + '/figures/dpt_df'
    loc_reverse_figures = save_location + '/' + str(experiment_name) + '/figures/dpt_figures'
    loc_reverse_heatmap = save_location + '/' + str(experiment_name) + '/figures/heatmap_figures'
    df_3 = pd.read_csv(loc_reverse_df + '/df_3.csv')
    add_closeness_score_PR(df_3)
    df_3['WHO'] = '3'

    df_4 = pd.read_csv(loc_reverse_df + '/df_4.csv')
    add_closeness_score_PR(df_4)
    df_4['WHO'] = '4'

    df_5 = pd.read_csv(loc_reverse_df + '/df_5.csv')
    add_closeness_score_PR(df_5)
    df_5['WHO'] = '5'

    df_7 = pd.read_csv(loc_reverse_df + '/df_7.csv')
    add_closeness_score_PR(df_7)
    df_7['WHO'] = '7'

    df_all_reverse = df_3.append(df_4).append(df_5).append(df_7)
    df_all_pivot_reverse = df_all_reverse.pivot("timepoint","WHO", "closeness_score").astype(float)

    df_d2 = pd.read_csv(loc_reverse_df + '/df_d2.csv')[0:4]
    add_closeness_score_PR(df_d2)
    df_d2['day'] = 'd2'

    df_d3 = pd.read_csv(loc_reverse_df + '/df_d3.csv')[0:4]
    add_closeness_score_PR(df_d3)
    df_d3['day'] = 'd3'

    df_d5 = pd.read_csv(loc_reverse_df + '/df_d5.csv')[0:4]
    add_closeness_score_PR(df_d5)
    df_d5['day'] = 'd5'

    df_e14 = pd.read_csv(loc_reverse_df + '/df_e14.csv')[0:4]
    add_closeness_score_PR(df_e14)
    df_e14['day'] = 'e14'


    df_hamster = df_d2.append(df_d3).append(df_d5).append(df_e14)
    df_hamster_pivot = df_hamster.pivot("day", "cov_degree", "closeness_score").astype(float)

    df_pivot_undirected = df_hamster_pivot*df_all_pivot_reverse

    plt.figure()
    sns.heatmap(df_pivot_undirected,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('timepoint',fontsize = 12)
    plt.ylabel('WHO degree',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_reverse_heatmap + '/undirected_heatmap.pdf')
    plt.close()  
    
    
def undirected_mappingPR(celltype = 'Classical_Monocytes',save_location = None,experiment_name = None):
    model_save_string = save_location + '/' + str(experiment_name) + '/model/model.pt'
    loc_reverse_df = save_location + '/' + str(experiment_name) + '/figures/dpt_df'
    loc_reverse_figures = save_location + '/' + str(experiment_name) + '/figures/dpt_figures'
    loc_reverse_heatmap = save_location + '/' + str(experiment_name) + '/figures/heatmap_figures'
    df_3 = pd.read_csv(loc_reverse_df + '/df_3.csv')
    add_closeness_score_PR(df_3)
    df_3['WHO'] = '3'

    df_4 = pd.read_csv(loc_reverse_df + '/df_4.csv')
    add_closeness_score_PR(df_4)
    df_4['WHO'] = '4'

    df_5 = pd.read_csv(loc_reverse_df + '/df_5.csv')
    add_closeness_score_PR(df_5)
    df_5['WHO'] = '5'

    df_7 = pd.read_csv(loc_reverse_df + '/df_7.csv')
    add_closeness_score_PR(df_7)
    df_7['WHO'] = '7'

    df_all_reverse = df_3.append(df_4).append(df_5).append(df_7)
    df_all_pivot_reverse = df_all_reverse.pivot("timepoint","WHO", "closeness_score").astype(float)

    df_hd_D2 = pd.read_csv(loc_reverse_df + '/df_hd_D2.csv')[0:4]
    add_closeness_score_PR(df_hd_D2)
    df_hd_D2['day'] = 'hd D2'
    df_hd_D2['timepoint'] = 'hd_D2'

    df_hd_D3 = pd.read_csv(loc_reverse_df + '/df_hd_D3.csv')[0:4]
    add_closeness_score_PR(df_hd_D3)
    df_hd_D3['day'] = 'hd D3'
    df_hd_D3['timepoint'] = 'hd_D3'

    df_ld_D2 = pd.read_csv(loc_reverse_df + '/df_ld_D2.csv')[0:4]
    add_closeness_score_PR(df_ld_D2)
    df_ld_D2['day'] = 'ld D2'
    df_ld_D2['timepoint'] = 'ld_D2'

    df_ld_D3 = pd.read_csv(loc_reverse_df + '/df_ld_D3.csv')[0:4]
    add_closeness_score_PR(df_ld_D3)
    df_ld_D3['day'] = 'ld D3'
    df_ld_D3['timepoint'] = 'ld_D3'


    df_hamsterPR = df_hd_D2.append(df_hd_D3).append(df_ld_D2).append(df_ld_D3)
    df_hamsterPR_pivot = df_hamsterPR.pivot("timepoint", "cov_degree", "closeness_score").astype(float)

    df_pivot_undirected = df_hamsterPR_pivot*df_all_pivot_reverse

    plt.figure()
    sns.heatmap(df_pivot_undirected,cmap = "Reds",annot=True)
    plt.yticks(rotation = 0)
    plt.xlabel('timepoint',fontsize = 12)
    plt.ylabel('WHO degree',fontsize = 12)
    plt.title(celltype,fontsize = 14)
    plt.savefig(loc_reverse_heatmap + '/undirected_heatmap.pdf')
    plt.close()  
    
    
def full_analysis(hamster_rawfile,
                 human_rawfile,
                 T_sub_file, 
                 species = 'ma',
                 celltype = 'Classical_Monocytes',
                 save_input=True,
                 save_location = None,
                 experiment_name = None,
                 gene_set = 'full',
                 n_latent = 5,
                 make_umaps = False,
                 batch_key = 'compare_umap',
                 labels_key = 'uniform_name_acurate',
                 do_batch_removal_human_full_set = False,
                 do_batch_removal_human = False,
                 mito_filtering = False,
              plot_qc_metics = False):
    #hamster_rawfile: AnnData, hamster raw input file
    #human_rawfile:  AnnData, human raw input file
    #T_sub_file:  AnnData, raw input file for T cell subclustering
    #species: str, 'ma' or 'pr' --------- 'ma' for syrian goldhamster and 'pr' roborovski dwarf hamster
    #celltype: str, category in adata.obs['uniform_name_acurate']
    #save_input: True or False - whether to save preprocessed data
    #save_location: str path to folder where to save experiment
    #eperiment_name: str, name of experiment, folder gets created in save_location directory
    #gene_set: 'full' or 'union_hv' or 'intersection_hv' - possible gene set inputs
    #make_umaps: True or False, whether to generate umaps
    #batch_key: column ins adata.obs considered as batch for scGEN 
    #labels_key: column ins adata.obs considered as labels_key
    #do_batch_removal_human: whether tp perform batch removal for 'orig.ident' column in adata.obs
    hamster,human = preprocess(hamster_rawfile = hamster_rawfile,
                               human_rawfile = human_rawfile,
                               T_sub = T_sub_file,
                               species = species,
                               celltype =  celltype,
                               create_folder = True,
                               save_input = save_input,
                               save_location = save_location,
                               experiment_name = experiment_name,
                               do_batch_removal_human_full_set = do_batch_removal_human_full_set,
                              do_batch_removal_human = do_batch_removal_human,
                               mito_filtering = mito_filtering,
                               plot_qc_metics = plot_qc_metics)
    if species == 'ma':
        train_autoencoder_model_MA(hamster,
                                   human,
                                   celltype = celltype,
                                   save_input = save_input,
                                   save_location = save_location,
                                   experiment_name = experiment_name,
                                   gene_set = gene_set,
                                   make_umaps = make_umaps,
                                   batch_key = batch_key,
                                   labels_key = labels_key,
                                   n_latent = n_latent)
        reverse_mappingMA(celltype = celltype,
                          save_location = save_location,
                          experiment_name = experiment_name)
        undirected_mappingMA(celltype = celltype,
                             save_location = save_location,
                             experiment_name = experiment_name)
        
    if species == 'pr':
        train_autoencoder_model_PR(hamster,
                                   human,
                                   celltype = celltype,
                                   save_input = save_input,
                                   save_location = save_location,
                                   experiment_name = experiment_name,
                                   gene_set = gene_set,
                                   make_umaps = make_umaps,
                                   batch_key = batch_key,
                                   labels_key = labels_key,
                                   n_latent = n_latent)
        reverse_mappingPR(celltype = celltype,
                          save_location = save_location,
                          experiment_name = experiment_name)
        undirected_mappingPR(celltype = celltype,
                             save_location = save_location,
                             experiment_name = experiment_name)
        
    add_pkf_score_relative(save_location + '/' + str(experiment_name),species)
    write_input_to_csv(species = species,
                   celltype = celltype,
                   save_location = save_location,
                   experiment_name = experiment_name,
                   gene_set = gene_set,
                   n_latent = n_latent,
                   batch_key = batch_key,
                   labels_key = labels_key)
    
    
def show_pkf_score(path_to_folder):
    files = os.listdir(path_to_folder)
    for file in files:
        pkf_score_loc = path_to_folder + '/' + file + '/description/pkf_score_relative.csv'
        try:
            pkf_score = pd.read_csv(pkf_score_loc)['pkf_score'][0]
            print(file + ' : ' + str(pkf_score))
        except:
            pass
        
def write_input_to_csv(species,celltype,save_location,experiment_name,gene_set,n_latent,batch_key,labels_key):
    categories = ['species','celltype','save_location','experiment_name','gene_set','n_latent','batch_key','labels_key']
    values = [species,celltype,save_location,experiment_name,gene_set,n_latent,batch_key,labels_key]
    param_df = pd.DataFrame([values],columns = categories)
    param_df.to_csv(save_location + '/' + str(experiment_name) + '/description/parameter_df.csv',index = False)
    
    
def remove_mito(adata,plot_qc_metics = True, save_qc_metric = None,mad_threshold = 5,mt_pct_threshold =  8):
    adata.var["mt" ] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
    if plot_qc_metics:
        #sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
        if save_qc_metric is not None:
            p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
        else:
            p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
        # sc.pl.violin(adata, 'total_counts')
        p2 = sc.pl.violin(adata, "pct_counts_mt")
        p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", mad_threshold)
    | is_outlier(adata, "log1p_n_genes_by_counts", mad_threshold)
    | is_outlier(adata, "pct_counts_in_top_20_genes", mad_threshold))
    adata.obs.outlier.value_counts()
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > mt_pct_threshold
)
    adata.obs.mt_outlier.value_counts()
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * M.mad()) | (
        np.median(M) + nmads * M.mad() < M
    )
    return outlier


def setting_A1():
    gene_set = 'intersection_hv'
    n_latent = 5
    do_batch_removal_human = True
    return gene_set,n_latent,do_batch_removal_human

def setting_A2():
    gene_set = 'intersection_hv'
    n_latent = 10
    do_batch_removal_human = True
    return gene_set,n_latent,do_batch_removal_human

def setting_A3():
    gene_set = 'intersection_hv'
    n_latent = 5
    do_batch_removal_human = False
    return gene_set,n_latent,do_batch_removal_human

def setting_A4():
    gene_set = 'intersection_hv'
    n_latent = 10
    do_batch_removal_human = False
    return gene_set,n_latent,do_batch_removal_human

def setting_B1():
    gene_set = 'union_hv'
    n_latent = 5
    do_batch_removal_human = True
    return gene_set,n_latent,do_batch_removal_human

def setting_B2():
    gene_set = 'union_hv'
    n_latent = 10
    do_batch_removal_human = True
    return gene_set,n_latent,do_batch_removal_human

def setting_B3():
    gene_set = 'union_hv'
    n_latent = 5
    do_batch_removal_human = False
    return gene_set,n_latent,do_batch_removal_human

def setting_B4():
    gene_set = 'union_hv'
    n_latent = 10
    do_batch_removal_human = False
    return gene_set,n_latent,do_batch_removal_human

def setting_A5():
    gene_set = 'intersection_hv'
    n_latent = 5
    do_batch_removal_human = False
    do_batch_removal_human_full_set = True
    return gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set

def setting_A6():
    gene_set = 'intersection_hv'
    n_latent = 10
    do_batch_removal_human = False
    do_batch_removal_human_full_set = True
    return gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set

def setting_B5():
    gene_set = 'union_hv'
    n_latent = 5
    do_batch_removal_human = False
    do_batch_removal_human_full_set = True
    return gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set

def setting_B6():
    gene_set = 'union_hv'
    n_latent = 10
    do_batch_removal_human = False
    do_batch_removal_human_full_set = True
    return gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set

def choose_setting(setting_str):
    #setting_str: 'A1','A2','A3','A4','B1','B2','B3','B4'
    if setting_str == 'A1':
        gene_set,n_latent,do_batch_removal_human = setting_A1()
    if setting_str == 'A2':
        gene_set,n_latent,do_batch_removal_human = setting_A2()
    if setting_str == 'A3':
        gene_set,n_latent,do_batch_removal_human = setting_A3()
    if setting_str == 'A4':
        gene_set,n_latent,do_batch_removal_human = setting_A4()
    if setting_str == 'B1':
        gene_set,n_latent,do_batch_removal_human = setting_B1()
    if setting_str == 'B2':
        gene_set,n_latent,do_batch_removal_human = setting_B2()
    if setting_str == 'B3':
        gene_set,n_latent,do_batch_removal_human = setting_B3()
    if setting_str == 'B4':
        gene_set,n_latent,do_batch_removal_human = setting_B4()
    return gene_set,n_latent,do_batch_removal_human

def choose_setting_out(setting_str):
    #setting_str: 'A5','A6','B5','B6'
    if setting_str == 'A5':
        gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set = setting_A5()
    if setting_str == 'A6':
        gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set = setting_A6()
    if setting_str == 'B5':
        gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set = setting_B5()
    if setting_str == 'B6':
        gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set = setting_B6()
    return gene_set,n_latent,do_batch_removal_human,do_batch_removal_human_full_set

def LISI_preprocess_human(human_rawfile,
                          celltype = 'Classical_Monocytes',
                          mito_filtering = True,
                          plot_qc_metics = False,
                          do_batch_removal_human = False):
    
    human = annotate_cohort2(human_rawfile)
    sc.pp.filter_cells(human, min_genes=200)
    sc.pp.filter_genes(human, min_cells=20)
    if mito_filtering:
        human = remove_mito(human,plot_qc_metics = plot_qc_metics)
    sc.pp.normalize_total(human, target_sum=1e4,exclude_highly_expressed = True)
    sc.pp.log1p(human)
    if do_batch_removal_human:
        sc.pp.combat(human,'orig.ident')
    human = filter_adata_obs(human,'uniform_name_acurate',celltype)
    #if do_batch_removal_human:
    #    sc.pp.combat(human,'orig.ident')
    return human

def LISI_perform_human(preprocessed_file,celltype,save_string):
    sc.pp.highly_variable_genes(preprocessed_file,layer = 'log_counts')
    hvg_human = list(preprocessed_file.var['highly_variable'][preprocessed_file.var['highly_variable'] == True].index)
    human = preprocessed_file[:,hvg_human]
    human.obs['LISI_score'] = hpy.lisi.compute_lisi(human.X,human.obs.filter(items=['who_per_sample']),['who_per_sample'])
    human.obs[['who_per_sample','LISI_score']].to_csv('df_' + save_string + '.csv',index = True)
    plt.figure()
    sns.boxplot(data = human.obs,x = "LISI_score",y="who_per_sample",
                order=["0","3","4","5","7"])
    plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig('box_' + save_string + '.pdf',bbox_inches='tight')
    plt.close()
    return human

def LISI_perform_hamster(preprocessed_file,celltype,save_string):
    sc.pp.highly_variable_genes(preprocessed_file,layer = 'log_counts')
    hvg_hamster = list(preprocessed_file.var['highly_variable'][preprocessed_file.var['highly_variable'] == True].index)
    hamster = preprocessed_file[:,hvg_hamster]
    hamster.obs['LISI_score'] = hpy.lisi.compute_lisi(hamster.X,hamster.obs.filter(items=['timepoint']),['timepoint'])
    hamster.obs[['timepoint','LISI_score']].to_csv('df_' + save_string + '.csv',index = True)
    plt.figure()
    sns.boxplot(data = hamster.obs,x = "LISI_score",y="timepoint")
    plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig('box_' + save_string + '.pdf',bbox_inches='tight')
    plt.close()
    return hamster

def LISI_perform_human_latent(human,celltype,save_string):
    human.obs['LISI_score'] = hpy.lisi.compute_lisi(human.X,human.obs.filter(items=['who_per_sample']),['who_per_sample'])
    human.obs[['who_per_sample','LISI_score']].to_csv('df_' + save_string + '_latent.csv',index = True)
    plt.figure()
    sns.boxplot(data = human.obs,x = "LISI_score",y="who_per_sample",
                order=["0","3","4","5","7"])
    plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig('box_' + save_string + '_latent.pdf',bbox_inches='tight')
    plt.close()
    return human

def LISI_umap_human(adata,save_string):
    prepare_umap(adata)
    sc.pl.umap(adata,color = ['who_per_sample','LISI_score'],save = save_string + '_umap.pdf')
    
def LISI_umap_hamster(adata,save_string):
    prepare_umap(adata)
    sc.pl.umap(adata,color = ['timepoint','LISI_score'],save = save_string + '_umap.pdf')
    
def make_pkf_score_df(path_to_folder,species,celltype):
    files = sorted(os.listdir(path_to_folder))
    res_df = pd.DataFrame(columns = ['setting','run','pkf_score'])
    for file in files:
        try:
            setting = file[0:2]
            run = file[-1:]
            pkf_score_loc = path_to_folder + '/' + file + '/description/pkf_score_relative.csv'
            pkf_score = pd.read_csv(pkf_score_loc)['pkf_score'][0]
            res_df.loc[file] = [setting,run,pkf_score] 
        except:
            pass
    res_df['celltype'] = celltype
    res_df['species'] = species
    return res_df

def get_LISI_latent_human(celltype,hamster_species,setting,run,save_string):
    base_exp = '/work/users/mh823zote/projects/cov/integration/run_VAE/results/' + hamster_species + '/grid/' + celltype + '/' + setting + '_run_' + str(run)
    model = scgen.SCGEN.load(base_exp + '/model/model.pt')
    adata = model.adata.copy()
    adata = filter_adata_obs(adata,'cells','Whole_blood')
    lat_adata = get_latent_representation_object(model,adata)
    LISI_perform_human_latent(lat_adata,celltype,save_string)
    return lat_adata

def get_LISI_latent_hamster(celltype,hamster_species,setting,run,save_string):
    base_exp = '/work/users/mh823zote/projects/cov/integration/run_VAE/results/' + hamster_species + '/grid/' + celltype + '/' + setting + '_run_' + str(run)
    model = scgen.SCGEN.load(base_exp + '/model/model.pt')
    adata = model.adata.copy()
    adata = leave_out_adata_obs(adata,'cells','Whole_blood')
    lat_adata = get_latent_representation_object(model,adata)
    LISI_perform_hamster_latent(lat_adata,celltype,save_string)
    return lat_adata

def LISI_perform_hamster_latent(hamster,celltype,save_string):
    hamster.obs['LISI_score'] = hpy.lisi.compute_lisi(hamster.X,hamster.obs.filter(items=['timepoint']),['timepoint'])
    hamster.obs[['timepoint','LISI_score']].to_csv('df_' + save_string + '_latent.csv',index = True)
    plt.figure()
    sns.boxplot(data = hamster.obs,x = "LISI_score",y="timepoint")
    plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig('box_' + save_string + '_latent.pdf',bbox_inches='tight')
    plt.close()
    return hamster

def LISI_preprocess_hamster(hamster_rawfile,T_sub,species,
                          celltype = 'Classical_Monocytes',
                          mito_filtering = True,
                          plot_qc_metics = False):
    if species == 'ma':
        hamster = annotate_golden(hamster_rawfile,T_sub)
    if species == 'pr':
        hamster = annotate_dwarf(hamster_rawfile,T_sub)
    sc.pp.filter_cells(hamster, min_genes=200)
    sc.pp.filter_genes(hamster, min_cells=20)
    if mito_filtering:
        hamster = remove_mito(hamster,plot_qc_metics = plot_qc_metics)
    sc.pp.normalize_total(hamster, target_sum=1e4,exclude_highly_expressed = True)
    sc.pp.log1p(hamster)
    hamster = filter_adata_obs(hamster,'uniform_name_acurate',celltype)
    if ((species == 'ma') & (celltype == 'Classical_Monocytes')):
        hamster = leave_out_adata_obs(hamster,'seurat_clusters',15)
    return hamster

def analyse_Autoencoder_latent_space(celltype,hamster_species,setting,run):
    #outputs in current folder
    save_string_human = setting + '_run_' + str(run) + '_' +  hamster_species + '_human'
    save_string_hamster = setting + '_run_' + str(run) + '_' +  hamster_species + '_hamster'
    
    #human
    lat_adata_human = get_LISI_latent_human(celltype,hamster_species,setting,run,save_string_human)
    LISI_umap_human(lat_adata_human,save_string_human)
    
    #hamster
    lat_adata_hamster = get_LISI_latent_hamster(celltype,hamster_species,setting,run,save_string_hamster)
    LISI_umap_hamster(lat_adata_hamster,save_string_hamster)

def return_analysis_input(triple):
    # e.g triple = ['ma','Classical_Monocytes','A5_run_3']
    hamster_species = triple[0]
    celltype = triple[1]
    setting = triple[2][0:2]
    run = int(triple[2][-1:])
    return celltype,hamster_species,setting,run


def control_infection_annotation_pr(adata):
    infection_string = []
    for q in range(len(adata)):
        if adata.obs['timepoint'][q] == 'D0':
            infection_string.append('hamster_control')
        else:
            infection_string.append('hamster_infection')
    adata.obs['general_type'] = infection_string
    return adata
    
def control_infection_annotation_ma(adata):
    infection_string = []
    for q in range(len(adata)):
        if adata.obs['timepoint'][q] == 'd0':
            infection_string.append('hamster_control')
        else:
            infection_string.append('hamster_infection')
    adata.obs['general_type'] = infection_string
    return adata
    
def control_infection_annotation_human(adata):
    infection_string = []
    for q in range(len(adata)):
        if adata.obs['who_per_sample'][q] == '0':
            infection_string.append('human_control')
        else:
            infection_string.append('human_infection')
    adata.obs['general_type'] = infection_string
    return adata

def control_infection_annotation_ma_shifted(adata):
    infection_string = []
    for q in range(len(adata)):
        if adata.obs['timepoint'][q] == 'd0':
            infection_string.append('hamster_control_shifted')
        else:
            infection_string.append('hamster_infection_shifted')
    adata.obs['general_type'] = infection_string
    return adata

def control_infection_annotation_pr_shifted(adata):
    infection_string = []
    for q in range(len(adata)):
        if adata.obs['timepoint'][q] == 'D0':
            infection_string.append('hamster_control_shifted')
        else:
            infection_string.append('hamster_infection_shifted')
    adata.obs['general_type'] = infection_string
    return adata

def show_latent_shift(celltype,hamster_species,setting,run):
    save_base = hamster_species + '_' + celltype + '_' + setting + '_run_' + str(run)
    save_string_human = setting + '_run_' + str(run) + '_' +  hamster_species + '_human'
    save_string_hamster = setting + '_run_' + str(run) + '_' +  hamster_species + '_hamster'
    lat_adata_human = get_LISI_latent_human(celltype,hamster_species,setting,run,save_string_human)
    lat_adata_hamster = get_LISI_latent_hamster(celltype,hamster_species,setting,run,save_string_hamster)
    
    if hamster_species == 'ma':
        lat_adata_hamster = control_infection_annotation_ma(lat_adata_hamster)
    if hamster_species == 'pr':     
        lat_adata_hamster = control_infection_annotation_pr(lat_adata_hamster)
    lat_adata_human = control_infection_annotation_human(lat_adata_human)
    lat_merged = merge_adata(lat_adata_human,lat_adata_hamster)
    
    lat_hamster_control = filter_adata_obs(lat_adata_hamster,'general_type','hamster_control')
    lat_hamster_infection = filter_adata_obs(lat_adata_hamster,'general_type','hamster_infection')
    lat_human_control = filter_adata_obs(lat_adata_human,'general_type','human_control')
    lat_human_infection = filter_adata_obs(lat_adata_human,'general_type','human_infection')
    
    delta = get_delta_in_latent_space(lat_hamster_control,lat_human_control)
    
    lat_adata_hamster_shifted = shift_adata_in_latent_space(lat_adata_hamster,delta)
    
    if hamster_species == 'ma':
        lat_adata_hamster_shifted = control_infection_annotation_ma_shifted(lat_adata_hamster_shifted)
    if hamster_species == 'pr':
        lat_adata_hamster_shifted = control_infection_annotation_pr_shifted(lat_adata_hamster_shifted)
    

    prepare_umap(lat_merged)
    sc.pl.umap(lat_merged,color = 'general_type',save = save_base + '_lat.pdf')
    
    lat_all_overview = merge_adata(lat_merged,lat_adata_hamster_shifted)
    prepare_umap(lat_all_overview)
    sc.pl.umap(lat_all_overview,color = 'general_type',save = save_base + '_lat_all_shift.pdf')
    
    lat_only_shift = merge_adata(lat_adata_human,lat_adata_hamster_shifted)
    prepare_umap(lat_only_shift)
    sc.pl.umap(lat_only_shift,color = 'general_type',save = save_base + '_lat_only_shift.pdf')
    
def save_heatmap_csv(triple,base_path,save_path_base):
    #triple experiment info ['pr', 'NK_Cells', 'A5_run_4']
    #base_path: script location
    #save_path_base: save location 
    hamster_species = triple[0]
    ct = triple[1]
    nr_run = triple[2]
    loc_reverse_df = base_path + '/run_VAE/results/' + hamster_species + '/grid/' + ct + '/' + nr_run + '/figures/dpt_df'
    
    if hamster_species == 'ma':
        df_3 = pd.read_csv(loc_reverse_df + '/df_3.csv')
        add_closeness_score_PR(df_3)
        df_3['WHO'] = '3'

        df_4 = pd.read_csv(loc_reverse_df + '/df_4.csv')
        add_closeness_score_PR(df_4)
        df_4['WHO'] = '4'

        df_5 = pd.read_csv(loc_reverse_df + '/df_5.csv')
        add_closeness_score_PR(df_5)
        df_5['WHO'] = '5'

        df_7 = pd.read_csv(loc_reverse_df + '/df_7.csv')
        add_closeness_score_PR(df_7)
        df_7['WHO'] = '7'

        df_all_reverse = df_3.append(df_4).append(df_5).append(df_7)
        df_all_pivot_reverse = df_all_reverse.pivot("timepoint","WHO", "closeness_score").astype(float)

        df_d2 = pd.read_csv(loc_reverse_df + '/df_d2.csv')[0:4]
        add_closeness_score_PR(df_d2)
        df_d2['day'] = 'd2'

        df_d3 = pd.read_csv(loc_reverse_df + '/df_d3.csv')[0:4]
        add_closeness_score_PR(df_d3)
        df_d3['day'] = 'd3'

        df_d5 = pd.read_csv(loc_reverse_df + '/df_d5.csv')[0:4]
        add_closeness_score_PR(df_d5)
        df_d5['day'] = 'd5'

        df_e14 = pd.read_csv(loc_reverse_df + '/df_e14.csv')[0:4]
        add_closeness_score_PR(df_e14)
        df_e14['day'] = 'e14'


        df_hamster = df_d2.append(df_d3).append(df_d5).append(df_e14)
        df_hamster_pivot = df_hamster.pivot("day", "cov_degree", "closeness_score").astype(float)

        df_pivot_undirected = df_hamster_pivot*df_all_pivot_reverse
        
        df_hamster_pivot.to_csv(save_path_base + '/' + hamster_species + '/' + ct + '/df_pivot_hamster_on_human.csv')
        df_all_pivot_reverse.to_csv(save_path_base + '/' + hamster_species + '/' + ct + '/df_pivot_human_on_hamster.csv')
        df_pivot_undirected.to_csv(save_path_base + '/' + hamster_species + '/' + ct + '/df_pivot_undirected.csv')

    if hamster_species == 'pr':
        df_3 = pd.read_csv(loc_reverse_df + '/df_3.csv')
        add_closeness_score_PR(df_3)
        df_3['WHO'] = '3'

        df_4 = pd.read_csv(loc_reverse_df + '/df_4.csv')
        add_closeness_score_PR(df_4)
        df_4['WHO'] = '4'

        df_5 = pd.read_csv(loc_reverse_df + '/df_5.csv')
        add_closeness_score_PR(df_5)
        df_5['WHO'] = '5'

        df_7 = pd.read_csv(loc_reverse_df + '/df_7.csv')
        add_closeness_score_PR(df_7)
        df_7['WHO'] = '7'

        df_all_reverse = df_3.append(df_4).append(df_5).append(df_7)
        df_all_pivot_reverse = df_all_reverse.pivot("timepoint","WHO", "closeness_score").astype(float)

        df_hd_D2 = pd.read_csv(loc_reverse_df + '/df_hd_D2.csv')[0:4]
        add_closeness_score_PR(df_hd_D2)
        df_hd_D2['day'] = 'hd D2'
        df_hd_D2['timepoint'] = 'hd_D2'

        df_hd_D3 = pd.read_csv(loc_reverse_df + '/df_hd_D3.csv')[0:4]
        add_closeness_score_PR(df_hd_D3)
        df_hd_D3['day'] = 'hd D3'
        df_hd_D3['timepoint'] = 'hd_D3'

        df_ld_D2 = pd.read_csv(loc_reverse_df + '/df_ld_D2.csv')[0:4]
        add_closeness_score_PR(df_ld_D2)
        df_ld_D2['day'] = 'ld D2'
        df_ld_D2['timepoint'] = 'ld_D2'

        df_ld_D3 = pd.read_csv(loc_reverse_df + '/df_ld_D3.csv')[0:4]
        add_closeness_score_PR(df_ld_D3)
        df_ld_D3['day'] = 'ld D3'
        df_ld_D3['timepoint'] = 'ld_D3'


        df_hamsterPR = df_hd_D2.append(df_hd_D3).append(df_ld_D2).append(df_ld_D3)
        df_hamsterPR_pivot = df_hamsterPR.pivot("timepoint", "cov_degree", "closeness_score").astype(float)

        df_pivot_undirected = df_hamsterPR_pivot*df_all_pivot_reverse

        df_hamsterPR_pivot.to_csv(save_path_base + '/' + hamster_species + '/' + ct + '/df_pivot_hamster_on_human.csv')
        df_all_pivot_reverse.to_csv(save_path_base + '/' + hamster_species + '/' + ct + '/df_pivot_human_on_hamster.csv')
        df_pivot_undirected.to_csv(save_path_base + '/' + hamster_species + '/' + ct + '/df_pivot_undirected.csv')
     
    
def print_filter_results(adata_before,adata_after):
    print('cells')
    print('-------------------------------------------------------------')
    print('cell number before : ' + str(len(adata_before)))
    print('cell number after  : ' + str(len(adata_after)))
    print('cells filtered     : ' + str(len(adata_before)-len(adata_after)) + ' (' + str(np.round(100*(len(adata_before)-len(adata_after))/len(adata_before),4)) + '%)')
    print(' ')
    print('genes')
    print('-------------------------------------------------------------')
    print('genes before       : ' + str(len(adata_before.var)))
    print('genes after        : ' + str(len(adata_after.var)))
    print('genes filtered     : ' + str(len(adata_before.var)-len(adata_after.var)) + ' (' + str(np.round(100*(len(adata_before.var)-len(adata_after.var))/len(adata_before.var),4)) + '%)')
    print(' ')
    

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * M.mad()) | (
        np.median(M) + nmads * M.mad() < M
    )
    return outlier

def outlier_MAD_analysis(adata,obs_string,nmads,print_overview = True):
    outlier = is_outlier(adata, obs_string, nmads)
    adata.obs['out_' + obs_string] = outlier
    min_non_outlier = np.round(np.min(adata[adata.obs['out_' + obs_string] == False].obs[obs_string]),4)
    max_non_outlier = np.round(np.max(adata[adata.obs['out_' + obs_string] == False].obs[obs_string]),4)
    min_outlier = np.round(np.min(adata[adata.obs['out_' + obs_string] == True].obs[obs_string]),4)
    max_outlier = np.round(np.max(adata[adata.obs['out_' + obs_string] == True].obs[obs_string]),4)
    nr_outliers = len(adata[adata.obs['out_' + obs_string] == True])
    nr_non_outliers = len(adata[adata.obs['out_' + obs_string] == False])
    fraction_outliers = np.round(nr_outliers / (nr_outliers + nr_non_outliers),6)
    if print_overview == True:
        print('range non outlier : [' + str(min_non_outlier) +  ' ; ' + str(max_non_outlier) + ']')
        if sum(outlier) != 0:
            if (min_outlier >= max_non_outlier) | (max_outlier <= min_non_outlier):
                print('range outlier     : [' + str(min_outlier) +  ' ; ' + str(max_outlier) + ']')
            else:
                lower_max = np.max(adata[adata.obs[obs_string] <= min_non_outlier].obs[obs_string])
                upper_min = np.min(adata[adata.obs[obs_string] >= max_non_outlier].obs[obs_string])
                print('range outlier     : [' + str(min_outlier) +  ' ; ' + str(lower_max) + '] and ['  + str(upper_min) +  ' ; ' + str(max_outlier) + ']')
            print('number outliers   : ' + str(nr_outliers))
            print('fraction outliers : ' + str(100*fraction_outliers) + '%')
        else:
            print('range outlier     :' + ' no outliers')
    return adata

def get_normalization_UMAP(df_UMAP):
    #UMAP euclidean distances
    df_UMAP_val = df_UMAP.values
    np.min(df_UMAP_val)
    min_val = np.min(df_UMAP_val)
    max_val = np.max(df_UMAP_val)
    shift_UMAP = 0.1*(max_val - min_val)
    norm_sim =  1-(df_UMAP.values - min_val + (0.5)*shift_UMAP) / ((max_val - min_val) + shift_UMAP)
    return pd.DataFrame(norm_sim)

def MSE_UMAP_scGEN_df(df_UMAP_norm,df_scGEN):
    MSE = ((df_scGEN - df_UMAP_norm)**2).values.mean()
    return MSE

def matrix_metric_UMAP_VAE(UMAP_euclidean,VAE_similarity):
    UMAP_similarity = get_normalization_UMAP(UMAP_euclidean)
    return MSE_UMAP_scGEN_df(UMAP_similarity,VAE_similarity)

def matrix_metric_UMAP_VAE(UMAP_euclidean,VAE_similarity,plot=False):
    UMAP_similarity = get_normalization_UMAP(UMAP_euclidean)
    
    #mse
    mse = MSE_UMAP_scGEN_df(UMAP_similarity,VAE_similarity)
    
    #spearman r correlation
    flat_UMAP = UMAP_similarity.to_numpy().flatten()
    flat_VAE = VAE_similarity.to_numpy().flatten()
    spearman = scistats.spearmanr(a=flat_UMAP,b=flat_VAE)
    if plot == True:
        plt.plot(flat_UMAP,flat_VAE,'x')
        plt.xlabel('UMAP similarity')
        plt.ylabel('VAE similarity')
        plt.ylim([0,1])
        plt.xlim([0,1])
        plt.plot([0,1],[0,1])
    return mse,spearman

def print_matrix_sim_results(mse,spearman):
    print('MSE                  : ' + str(mse))
    print('spearman correlation : '  + str(spearman[0]))
    print('p-value              : '  + str(spearman[1]))
    
def log_normalize(adata,target_sum = 1e4):
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
def run_VAE_celltype(hamster_ct_file,
                 human_ct_file,
                 species = 'ma',
                 celltype = 'Classical_Monocytes',
                 save_input=False,
                 save_location = None,
                 experiment_name = None,
                 gene_set = 'full',
                 n_latent = 5,
                 make_umaps = False,
                 batch_key = 'compare_umap',
                 labels_key = 'uniform_name_acurate',
                 do_LISI = True,
                 do_LISI_umap = True,
                    latent_shift_vis = True):
    #hamster_rawfile: AnnData, hamster raw input file
    #human_rawfile:  AnnData, human raw input file
    #T_sub_file:  AnnData, raw input file for T cell subclustering
    #species: str, 'ma' or 'pr' --------- 'ma' for syrian goldhamster and 'pr' roborovski dwarf hamster
    #celltype: str, category in adata.obs['uniform_name_acurate']
    #save_input: True or False - whether to save preprocessed data
    #save_location: str path to folder where to save experiment
    #eperiment_name: str, name of experiment, folder gets created in save_location directory
    #gene_set: 'full' or 'union_hv' or 'intersection_hv' - possible gene set inputs
    #make_umaps: True or False, whether to generate umaps
    #batch_key: column ins adata.obs considered as batch for scGEN 
    #labels_key: column ins adata.obs considered as labels_key
    #do_batch_removal_human: whether tp perform batch removal for 'orig.ident' column in adata.obs
    
    #preprocessing??
    hamster = hamster_ct_file.copy()
    human = human_ct_file.copy()
    
    if species == 'ma':
        train_autoencoder_model_MA(hamster,
                                   human,
                                   celltype = celltype,
                                   save_input = save_input,
                                   save_location = save_location,
                                   experiment_name = experiment_name,
                                   gene_set = gene_set,
                                   make_umaps = make_umaps,
                                   batch_key = batch_key,
                                   labels_key = labels_key,
                                   n_latent = n_latent)
        reverse_mappingMA(celltype = celltype,
                          save_location = save_location,
                          experiment_name = experiment_name)
        undirected_mappingMA(celltype = celltype,
                             save_location = save_location,
                             experiment_name = experiment_name)
        
    if species == 'pr':
        train_autoencoder_model_PR(hamster,
                                   human,
                                   celltype = celltype,
                                   save_input = save_input,
                                   save_location = save_location,
                                   experiment_name = experiment_name,
                                   gene_set = gene_set,
                                   make_umaps = make_umaps,
                                   batch_key = batch_key,
                                   labels_key = labels_key,
                                   n_latent = n_latent)
        reverse_mappingPR(celltype = celltype,
                          save_location = save_location,
                          experiment_name = experiment_name)
        undirected_mappingPR(celltype = celltype,
                             save_location = save_location,
                             experiment_name = experiment_name)
        
    add_pkf_score_relative(save_location + '/' + str(experiment_name),species)
    write_input_to_csv(species = species,
                   celltype = celltype,
                   save_location = save_location,
                   experiment_name = experiment_name,
                   gene_set = gene_set,
                   n_latent = n_latent,
                   batch_key = batch_key,
                   labels_key = labels_key)
    save_train_history(save_location,experiment_name)
    if do_LISI:
        LISI_perform_human_celltype(human,save_location,experiment_name,celltype = celltype,species = species,do_LISI_umap = do_LISI_umap)
        LISI_perform_hamster_celltype(hamster,save_location,experiment_name,celltype = celltype,species = species,do_LISI_umap = do_LISI_umap)
        LISI_latent_human_celltype(save_location,experiment_name,celltype = celltype,species= species,do_LISI_umap = do_LISI_umap)
        LISI_latent_hamster_celltype(save_location,experiment_name,celltype = celltype,species= species,do_LISI_umap = do_LISI_umap)
    if latent_shift_vis == True:
        latent_shift_analysis_celltype(save_location,experiment_name,celltype = celltype,species= species)
    
def save_train_history(save_location,experiment_name):
    model = scgen.SCGEN.load(save_location + '/' + experiment_name + '/model/model.pt')
    # loss
    #-----------------------------------------------------------------------------------------
    plt.figure()
    plt.title('loss')
    plt.plot(model.history['train_loss_epoch'],label = 'train')
    plt.plot(model.history['validation_loss'],label = 'val')
    plt.xlabel('epoch')
    plt.legend()
    plt.savefig(save_location + '/' + experiment_name + '/figures/train_figures/loss.pdf')
    plt.close()
    #-----------------------------------------------------------------------------------------
    plt.title('reconstruction loss')
    plt.plot(model.history['reconstruction_loss_train'],label = 'train')
    plt.plot(model.history['reconstruction_loss_validation'],label = 'val')
    plt.xlabel('epoch')
    plt.legend()
    plt.savefig(save_location + '/' + experiment_name + '/figures/train_figures/reconstruction_loss.pdf')
    plt.close()
    #-----------------------------------------------------------------------------------------
    plt.plot(model.history['elbo_train'], label = 'elbo')
    plt.plot(model.history['kl_local_train'],label = 'kl_local')
    plt.xlabel('epoch')
    plt.legend()
    plt.savefig(save_location + '/' + experiment_name + '/figures/train_figures/elbo.pdf')
    plt.close()
    
    o = pd.DataFrame()
    o['train_loss_epoch'] = model.history['train_loss_epoch']
    o['validation_loss'] = model.history['validation_loss']
    o['reconstruction_loss_train'] = model.history['reconstruction_loss_train']
    o['reconstruction_loss_validation'] = model.history['reconstruction_loss_validation']
    o['elbo_train'] = model.history['elbo_train']
    o['kl_local_train'] = model.history['kl_local_train']
    
    o.to_csv(save_location + '/' + experiment_name + '/description/train_scores.csv')
    
def LISI_perform_human_celltype(preprocessed_file,save_location,experiment_name,celltype = None,species= None,do_LISI_umap = True):
    sc.pp.highly_variable_genes(preprocessed_file,layer = 'log_counts')
    hvg_human = list(preprocessed_file.var['highly_variable'][preprocessed_file.var['highly_variable'] == True].index)
    human = preprocessed_file[:,hvg_human]
    human.obs['LISI_score'] = hpy.lisi.compute_lisi(human.X,human.obs.filter(items=['who_per_sample']),['who_per_sample'])
    human.obs[['who_per_sample','LISI_score']].to_csv(save_location + '/' + experiment_name + '/description/LISI_scores_human_GE.csv',index = True)
    plt.figure()
    sns.boxplot(data = human.obs,x = "LISI_score",y="who_per_sample",
                order=["0","3","4","5","7"])
    #plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig(save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_boxplot_human_GE.pdf',bbox_inches='tight')
    plt.close()
    
    if do_LISI_umap:
        prepare_umap(human)
        sc.pl.umap(human,color = ['who_per_sample','LISI_score'],save = celltype + '_' + species + '_' + experiment_name + '.pdf')
        src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + '.pdf'
        dst_path = save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_umap_human_GE.pdf'
        os.rename(src_path,dst_path)
    
    
def LISI_perform_human_latent_celltype(human,save_location,experiment_name,celltype = None,species= None,do_LISI_umap = True):
    human.obs['LISI_score'] = hpy.lisi.compute_lisi(human.X,human.obs.filter(items=['who_per_sample']),['who_per_sample'])
    human.obs[['who_per_sample','LISI_score']].to_csv(save_location + '/' + experiment_name + '/description/LISI_scores_human_latent.csv',index = True)
    plt.figure()
    sns.boxplot(data = human.obs,x = "LISI_score",y="who_per_sample",
                order=["0","3","4","5","7"])
    #plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig(save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_boxplot_human_latent.pdf',bbox_inches='tight')
    plt.close()
    
    if do_LISI_umap:
        prepare_umap(human)
        sc.pl.umap(human,color = ['who_per_sample','LISI_score'],save = celltype + '_' + species + '_' + experiment_name + '.pdf')
        src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + '.pdf'
        dst_path = save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_umap_human_latent.pdf'
        os.rename(src_path,dst_path)
        
def LISI_latent_human_celltype(save_location,experiment_name,celltype = None,species= None,do_LISI_umap = True):
    model = scgen.SCGEN.load(save_location + '/' + experiment_name + '/model/model.pt')
    adata = model.adata.copy()
    adata = filter_adata_obs(adata,'cells','Whole_blood')
    lat_adata = get_latent_representation_object(model,adata)
    LISI_perform_human_latent_celltype(lat_adata,save_location,experiment_name,celltype = celltype,species= species,do_LISI_umap = do_LISI_umap)
    
    
def LISI_perform_hamster_celltype(preprocessed_file,save_location,experiment_name,celltype = None,species= None,do_LISI_umap = True):
    sc.pp.highly_variable_genes(preprocessed_file,layer = 'log_counts')
    hvg_hamster = list(preprocessed_file.var['highly_variable'][preprocessed_file.var['highly_variable'] == True].index)
    hamster = preprocessed_file[:,hvg_hamster]
    hamster.obs['LISI_score'] = hpy.lisi.compute_lisi(hamster.X,hamster.obs.filter(items=['timepoint']),['timepoint'])
    hamster.obs[['timepoint','LISI_score']].to_csv(save_location + '/' + experiment_name + '/description/LISI_scores_hamster_GE.csv',index = True)
    plt.figure()
    sns.boxplot(data = hamster.obs,x = "LISI_score",y="timepoint")
    #plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig(save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_boxplot_hamster_GE.pdf',bbox_inches='tight')
    plt.close()
    if do_LISI_umap:
        prepare_umap(hamster)
        sc.pl.umap(hamster,color = ['timepoint','LISI_score'],save = celltype + '_' + species + '_' + experiment_name + 'hamster.pdf')
        src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + 'hamster.pdf'
        dst_path = save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_umap_hamster_GE.pdf'
        os.rename(src_path,dst_path)
        
def LISI_perform_hamster_latent_celltype(hamster,save_location,experiment_name,celltype = None,species= None,do_LISI_umap = True):
    hamster.obs['LISI_score'] = hpy.lisi.compute_lisi(hamster.X,hamster.obs.filter(items=['timepoint']),['timepoint'])
    hamster.obs[['timepoint','LISI_score']].to_csv(save_location + '/' + experiment_name + '/description/LISI_scores_hamster_latent.csv',index = True)
    plt.figure()
    sns.boxplot(data = hamster.obs,x = "LISI_score",y="timepoint")
    #plt.title(celltype)
    plt.xlim([0.95,4.5])
    plt.savefig(save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_boxplot_hamster_latent.pdf',bbox_inches='tight')
    plt.close()
    if do_LISI_umap:
        prepare_umap(hamster)
        sc.pl.umap(hamster,color = ['timepoint','LISI_score'],save = celltype + '_' + species + '_' + experiment_name + 'hamster_latent.pdf')
        src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + 'hamster_latent.pdf'
        dst_path = save_location + '/' + experiment_name + '/figures/LISI_figures/LISI_umap_hamster_latent.pdf'
        os.rename(src_path,dst_path)
        
def LISI_latent_hamster_celltype(save_location,experiment_name,celltype = None,species= None,do_LISI_umap = True):
    model = scgen.SCGEN.load(save_location + '/' + experiment_name + '/model/model.pt')
    adata = model.adata.copy()
    adata = leave_out_adata_obs(adata,'cells','Whole_blood')
    lat_adata = get_latent_representation_object(model,adata)
    LISI_perform_hamster_latent_celltype(lat_adata,save_location,experiment_name,celltype = celltype,species= species,do_LISI_umap = do_LISI_umap)
    
def latent_shift_analysis_celltype(save_location,experiment_name,celltype = None,species= None):
    model = scgen.SCGEN.load(save_location + '/' + experiment_name + '/model/model.pt')
    hamster_species = species

    latent_adata = get_latent_representation_object(model,model.adata)

    latent_adata_hamster = leave_out_adata_obs(latent_adata,'cells','Whole_blood')
    latent_adata_human = filter_adata_obs(latent_adata,'cells','Whole_blood')

    if hamster_species == 'ma':
        latent_adata_hamster = control_infection_annotation_ma(latent_adata_hamster)
    if hamster_species == 'pr':     
            latent_adata_hamster = control_infection_annotation_pr(latent_adata_hamster)
    latent_adata_human = control_infection_annotation_human(latent_adata_human)
    lat_merged = merge_adata(latent_adata_human,latent_adata_hamster)

    lat_hamster_control = filter_adata_obs(latent_adata_hamster,'general_type','hamster_control')
    lat_hamster_infection = filter_adata_obs(latent_adata_hamster,'general_type','hamster_infection')
    lat_human_control = filter_adata_obs(latent_adata_human,'general_type','human_control')
    lat_human_infection = filter_adata_obs(latent_adata_human,'general_type','human_infection')
    delta = get_delta_in_latent_space(lat_hamster_control,lat_human_control)
    lat_adata_hamster_shifted = shift_adata_in_latent_space(latent_adata_hamster,delta)    

    if hamster_species == 'ma':
        lat_adata_hamster_shifted = control_infection_annotation_ma_shifted(lat_adata_hamster_shifted)
    if hamster_species == 'pr':
        lat_adata_hamster_shifted = control_infection_annotation_pr_shifted(lat_adata_hamster_shifted)

    prepare_umap(lat_merged)
    sc.pl.umap(lat_merged,color = 'general_type',save = celltype + '_' + species + '_' + experiment_name + 'joint_latent.pdf')
    src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + 'joint_latent.pdf'
    dst_path = save_location + '/' + experiment_name + '/figures/umaps/joint_latent.pdf'
    os.rename(src_path,dst_path)

    lat_all_overview = merge_adata(lat_merged,lat_adata_hamster_shifted)
    prepare_umap(lat_all_overview)
    sc.pl.umap(lat_all_overview,color = 'general_type',save = celltype + '_' + species + '_' + experiment_name + 'joint_latent_shifted_all.pdf')
    src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + 'joint_latent_shifted_all.pdf'
    dst_path = save_location + '/' + experiment_name + '/figures/umaps/joint_latent_shifted_all.pdf'
    os.rename(src_path,dst_path)

    lat_only_shift = merge_adata(latent_adata_human,lat_adata_hamster_shifted)
    prepare_umap(lat_only_shift)
    sc.pl.umap(lat_only_shift,color = 'general_type',save = celltype + '_' + species + '_' + experiment_name + 'joint_latent_only_shifted.pdf')
    src_path = 'figures/umap' + celltype + '_' + species + '_' + experiment_name + 'joint_latent_only_shifted.pdf'
    dst_path = save_location + '/' + experiment_name + '/figures/umaps/joint_latent_only_shifted.pdf'
    os.rename(src_path,dst_path)
    
    
def spearman_matrix_UMAP_VAE_no_norm(UMAP_euclidean,VAE_similarity,plot=False):
    #spearman r correlation
    flat_UMAP = UMAP_euclidean.to_numpy().flatten()
    flat_VAE = VAE_similarity.to_numpy().flatten()
    spearman = scistats.spearmanr(a=flat_UMAP,b=flat_VAE)
    if plot == True:
        plt.plot(flat_UMAP,flat_VAE,'x')
        plt.xlabel('UMAP Euclidean distance')
        plt.ylabel('VAE similarity')
        #plt.ylim([0,1])
        #plt.xlim([0,1])
        #plt.plot([0,1],[0,1])
    return spearman

def build_adata_scratch(expression_matrix,meta_obs = None, meta_var = None):
    adata = ad.AnnData(expression_matrix)
    if meta_obs is not None:
        adata.obs = meta_obs
    if meta_var is not None:
        adata.var = meta_var
    return adata
    
def unnamed_as_index(df):
    df = df.rename(columns={'Unnamed: 0': 'index'})
    df.set_index('index', inplace=True)
    return df


def process_v223_2(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)
    
    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)
    
    adata = sc.read_h5ad(path_counts)
    
    adata.obs =  obs_df
    adata.var =  var_df
    adata.obs['dataset'] = adata.obs['species'].copy()
    
    adata.obs = adata.obs[['cells','species','uniform_name_overview','who_per_sample','timepoint','dataset']]
    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')
    adata.layers['log_counts'] = adata.X.copy()
    
    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})

    hamster_ma = filter_adata_obs(hamster,col_name='species',val='ma')
    hamster_pr = filter_adata_obs(hamster,col_name='species',val='pr')
    
    sc.pp.log1p(hamster_ma,layer = 'log_counts')
    sc.pp.log1p(hamster_pr, layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    
    return human, hamster_ma, hamster_pr


def process_v227_2(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)

    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)

    adata = sc.read_h5ad(path_counts)

    adata.obs =  obs_df
    adata.var =  var_df

    adata.obs['dataset'] = adata.obs['species'].copy()

    adata.obs = adata.obs[['cells','species','uniform_name_overview','who_per_sample','timepoint','dataset']]

    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')

    adata.layers['log_counts'] = adata.X.copy()

    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})

    sc.pp.log1p(hamster,layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    return human,hamster

def v227_2_get_pillais_trace(base_results,species,celltype,setting,run):
    model_path = base_results + species + '/' + celltype + '/' + setting + '_run_' + str(run) + '/' + 'model/model.pt'
    model = scgen.SCGEN.load(model_path)
    data = model.adata.copy()
    lat_data = get_latent_representation_object(model,data)
    lat_data.X.shape
    c,g = lat_data.X.shape
    column_names = [f'dimension_{i+1}' for i in range(g)]
    lat_df = pd.DataFrame(data=lat_data.X, columns=column_names)
    lat_df['severity_group'] = lat_data.obs['compare_umap'].values
    predictors = ' + '.join([f'dimension_{i}' for i in range(1, g+1)])
    formula = f'{predictors} ~ severity_group'
    
    #fit MANOVA
    fit = MANOVA.from_formula(formula, data=lat_df)
    results = fit.mv_test()
    
    #get pillais trace results
    pillai_trace = results.results['severity_group']['stat'].loc['Pillai\'s trace']
    pillai_trace_df = pd.DataFrame(pillai_trace)
    pillai_trace_df = pillai_trace_df.rename(columns={'Pillai\'s trace' : setting + '_run_' + str(run)})
    return pillai_trace_df
    
def v227_2_pkf_comparison(base_results,species,celltype,setting,run):
    pkf_path = base_results + species + '/' + celltype + '/' + setting + '_run_' + str(run) + '/' + 'description/pkf_score_relative.csv'
    pkf_df = pd.read_csv(pkf_path,index_col = 0)
    pkf_df = pkf_df.rename(columns={'pkf_score' : setting + '_run_' + str(run)})
    return pkf_df

def add_Neutrophil_resolution_227_2(adata):
    column_human = 'cluster_labels_res.0.8'
    human_neutrophil_names = ['Neutrophils_1','Neutrophils_2','Neutrophils_3','Neutrophils_4']
    human_immature_neutrophil_names = ['Immature Neutrophils_1','Immature Neutrophils_2']

    column_hamster = 'celltype'
    hamster_neutrophil_names = ['Neutrophil']
    hamster_immature_neutrophil_names = ['Immature neutrophil']

    uniform_name_acurate = []
    for j in range(len(adata)):
        if adata.obs['species'][j] == 'human':
            if adata.obs[column_human][j] in human_neutrophil_names:
                uniform_name_acurate.append('Neutrophils')
            if adata.obs[column_human][j] in human_immature_neutrophil_names:
                uniform_name_acurate.append('Immature_Neutrophils')
            if (adata.obs[column_human][j] not in human_neutrophil_names) & (adata.obs[column_human][j] not in human_immature_neutrophil_names):
                uniform_name_acurate.append('no_Neutrophils')
        if adata.obs['species'][j] != 'human':
            if adata.obs[column_hamster][j] in hamster_neutrophil_names:
                uniform_name_acurate.append('Neutrophils')
            if adata.obs[column_hamster][j] in hamster_immature_neutrophil_names:
                uniform_name_acurate.append('Immature_Neutrophils')
            if (adata.obs[column_hamster][j] not in hamster_neutrophil_names) & (adata.obs[column_hamster][j] not in hamster_immature_neutrophil_names):
                uniform_name_acurate.append('no_Neutrophils')
    adata.obs['uniform_name_acurate_N'] = uniform_name_acurate
    return adata

def process_v227_2_Neutrophil(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)

    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)

    adata = sc.read_h5ad(path_counts)

    adata.obs =  obs_df
    adata.var =  var_df

    adata.obs['dataset'] = adata.obs['species'].copy()

    adata.obs = adata.obs[['cells','species','cluster_labels_res.0.8','celltype','uniform_name_overview','who_per_sample','timepoint','dataset']]

    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')
    adata = add_Neutrophil_resolution_227_2(adata)

    adata.layers['log_counts'] = adata.X.copy()

    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})

    sc.pp.log1p(hamster,layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    return human,hamster

def process_v228_1(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)
    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)
    adata = sc.read_h5ad(path_counts)
    adata.obs =  obs_df
    adata.var =  var_df
    #-----------------------------------------------------------------------
    adata = filter_adata_obs(adata,
                      col_name='uniform_name_overview_keep',val=True)
    #-----------------------------------------------------------------------
    adata.obs['dataset'] = adata.obs['species'].copy()

    adata.obs = adata.obs[['cells','species','uniform_name_overview','who_per_sample','timepoint','dataset']]

    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')

    adata.layers['log_counts'] = adata.X.copy()

    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})


    sc.pp.log1p(hamster,layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    return human,hamster

def process_v228_1_Neutrophil(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)

    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)

    adata = sc.read_h5ad(path_counts)

    adata.obs =  obs_df
    adata.var =  var_df
    #-----------------------------------------------------------------------
    adata = filter_adata_obs(adata,
                      col_name='uniform_name_overview_keep',val=True)
    #-----------------------------------------------------------------------

    adata.obs['dataset'] = adata.obs['species'].copy()

    adata.obs = adata.obs[['cells','species','cluster_labels_res.0.8','celltype','uniform_name_overview','who_per_sample','timepoint','dataset']]

    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')
    adata = add_Neutrophil_resolution_227_2(adata)

    adata.layers['log_counts'] = adata.X.copy()

    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})

    sc.pp.log1p(hamster,layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    return human,hamster

def process_v228_1_proof_principle(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)
    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)
    adata = sc.read_h5ad(path_counts)
    adata.obs =  obs_df
    adata.var =  var_df
    #-----------------------------------------------------------------------
    adata = filter_adata_obs(adata,
                      col_name='uniform_name_overview_keep',val=True)
    #-----------------------------------------------------------------------
    adata.obs['dataset'] = adata.obs['species'].copy()

    adata.obs = adata.obs[['cells','donor','species','uniform_name_overview','who_per_sample','timepoint','dataset']]

    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')

    adata.layers['log_counts'] = adata.X.copy()

    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})


    sc.pp.log1p(hamster,layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    return human,hamster

def decode_latent_object(model,latent_object,GE_object):
    #GE object only for gene names
    decoded_cells = (
                    model.module.generative(torch.Tensor(latent_object.X))["px"].cpu().detach().numpy()
        )
    decoded_adata = AnnData(
                    X=decoded_cells,
                    obs=latent_object.obs.copy(),
                    var = pd.DataFrame(index=GE_object.var.index)
                )
    return decoded_adata

def make_identifier(object_number,script_number,celltype = None):
    if celltype is not None:
        return object_number + '_' + script_number + '_' + celltype
    else:
        return object_number + '_' + script_number
    
def make_folder_if_not_existing(folder_path):
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
        print('folder created!')
    else:
        print('folder exists!')
    
def process_v228_1_Neutrophil_proof_principle(path_counts,path_obs,path_var):
    obs_df = pd.read_csv(path_obs,low_memory=False)
    obs_df.set_index('rn', inplace=True)

    var_df = pd.read_csv(path_var,low_memory=False)
    var_df.set_index('rn', inplace=True)

    adata = sc.read_h5ad(path_counts)

    adata.obs =  obs_df
    adata.var =  var_df
    #-----------------------------------------------------------------------
    adata = filter_adata_obs(adata,
                      col_name='uniform_name_overview_keep',val=True)
    #-----------------------------------------------------------------------

    adata.obs['dataset'] = adata.obs['species'].copy()

    adata.obs = adata.obs[['cells','species','cluster_labels_res.0.8','celltype','uniform_name_overview','who_per_sample','timepoint','dataset','donor','hamster']]

    adata.obs['cells'] = adata.obs['cells'].fillna('hamster_WB')
    adata = add_Neutrophil_resolution_227_2(adata)

    adata.layers['log_counts'] = adata.X.copy()

    human = filter_adata_obs(adata,col_name='cells',val='Whole_blood')
    hamster = leave_out_adata_obs(adata,col_name='cells',val='Whole_blood')

    human.obs = human.obs.drop(['timepoint'], axis = 1)
    human.obs['who_per_sample'] = human.obs['who_per_sample'].astype(int).astype(str)

    hamster.obs = hamster.obs.drop(['who_per_sample'], axis = 1)
    hamster.obs['species'] = hamster.obs['species'].replace({'hamsterMA': 'ma', 'hamsterPR': 'pr'})

    sc.pp.log1p(hamster,layer = 'log_counts')
    sc.pp.log1p(human, layer = 'log_counts')
    return human,hamster
    
    
    


