import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

save_folder_pillai = '/work/users/username/projects/cov/integration/run_VAE_228_1/results/model_selection/pillais_trace/'


pillai_ma = pd.read_csv(save_folder_pillai + 'v14' + '_' + 'ma' + '_' + 'all_pillai.csv',index_col = 0)
pillai_pr = pd.read_csv(save_folder_pillai + 'v14' + '_' + 'pr' + '_' + 'all_pillai.csv',index_col = 0)

ma_B4 = pillai_ma.filter(regex='^B4')
pr_B4 = pillai_pr.filter(regex='^B4')

best_model_pr_B4 = pr_B4.idxmax(axis=1)
best_model_ma_B4 = ma_B4.idxmax(axis=1)

best_models_B4 = pd.DataFrame({'best_model_ma': best_model_ma_B4,'best_model_pr': best_model_pr_B4})
best_models_B4.to_csv(save_folder_pillai + 'v15_best_models_B4.csv')

best_models_B4 = pd.DataFrame({'best_model_ma': best_model_ma_B4,'best_model_pr': best_model_pr_B4})
best_models_B4.to_csv(save_folder_pillai + 'v15_best_models_B4.csv')
