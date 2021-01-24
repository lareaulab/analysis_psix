import importlib
import numpy as np
import pandas as pd
import os
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from scipy import stats as st
import seaborn as sns

import numpy.random as r

import numpy as np
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

from scipy.stats import combine_pvalues
# %run -i '../../utils/Kruskal_Wallis_test_functions.py'

from tqdm import tqdm

############

from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import euclidean_distances


def get_distance_matrix(pca, k=10):
    nbrs = NearestNeighbors(n_neighbors=k).fit(pca[['PC_1', 'PC_2']])
    distances, indices = nbrs.kneighbors(pca[['PC_1', 'PC_2']])
    
    cells = list(pca.index)
    
    W = pd.DataFrame(np.zeros((len(cells), len(cells))))
    W.columns = cells
    W.index = cells
    
    for i in tqdm(range(len(cells))):
        cell_i = cells[i]
        sigma = np.max(distances[i])
        for j in range(len(distances[i])):
            cell_j = cells[indices[i][j]]
            d = distances[i][j]
            w = np.exp(-(d**2)/(sigma**2))        
            W.loc[cell_i, cell_j] = w
    
    return W


def get_signature_matrix(PSI_tab):
    return (PSI_tab - PSI_tab.mean())/PSI_tab.std()


def make_mock_C_scores(norm_PSI, Ws, exon_list, total_cells, mock=100000):
    exon_out_list = []
    C_scores = []
    for i in tqdm(range(mock)):
        mock_run = True
        while mock_run:
            exon = r.choice(exon_list, 1)[0]
#             print(exon)
            scramble_cells = r.choice(norm_PSI.columns, total_cells, replace=False)
            mock_PSI = pd.DataFrame(norm_PSI.loc[exon, scramble_cells]).T

#             print(mock_PSI.shape)
#             print(norm_PSI.shape)

            mock_PSI.columns = norm_PSI.columns
            mock_df = mock_PSI.loc[exon]
#             print(type(mock_df))
#             print(mock_df)
            mock_score = get_C(mock_df, Ws)
            
#             print(mock_score)

            if mock_score >= 0:
                C_scores.append(mock_score)
                exon_out_list.append('mock_'+exon+'_'+str(i))
                mock_run = False
    return exon_out_list, C_scores
                    
def get_C(exon_score, W):
    exon_score = exon_score.dropna()
    obs_cells = exon_score.index
    x = (exon_score.values.reshape(-1, 1) - exon_score.values.reshape(1, -1))
    w = W.loc[obs_cells, obs_cells]
    num = (len(obs_cells)-1)*((w*(x**2)).sum().sum())
    den = (2*w.sum().sum())*np.sum((exon_score - exon_score.mean())**2)
    C = num/den
    score = 1 - C
    return score
    
    
##############################


# Ahora si el bueno

def get_mock_dict(PSI_tab, norm_PSI, Ws, mock=200):

    total_cells = len(PSI_tab.columns)

    exons_05_10 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.4) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.45)] 
    exons_10_20 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.3) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.40)] 
    exons_20_30 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.2) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.30)] 
    exons_30_40 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))>0.1) & (np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.20)] 
    exons_40_50 = PSI_tab.index[(np.abs(0.5 - PSI_tab.mean(axis = 1))<=0.1)] 

    exons_obs_10_20 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.9) & (PSI_tab.isna().mean(axis=1) > 0.8)]
    exons_obs_20_30 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.8) & (PSI_tab.isna().mean(axis=1) > 0.7)]
    exons_obs_30_40 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.7) & (PSI_tab.isna().mean(axis=1) > 0.6)]
    exons_obs_40_50 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.6) & (PSI_tab.isna().mean(axis=1) > 0.5)]
    exons_obs_50_60 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.5) & (PSI_tab.isna().mean(axis=1) > 0.4)]
    exons_obs_60_70 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.4) & (PSI_tab.isna().mean(axis=1) > 0.3)]
    exons_obs_70_80 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.3) & (PSI_tab.isna().mean(axis=1) > 0.2)]
    exons_obs_80_90 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.2) & (PSI_tab.isna().mean(axis=1) > 0.1)]
    exons_obs_90_100 = PSI_tab.index[(PSI_tab.isna().mean(axis=1) <= 0.1)]

    list1 = [exons_05_10, exons_10_20, exons_20_30, exons_30_40, exons_40_50]
    list2 = [exons_obs_10_20, exons_obs_20_30, exons_obs_30_40, exons_obs_40_50,
             exons_obs_50_60, exons_obs_60_70, exons_obs_70_80, exons_obs_80_90, exons_obs_90_100]


    exon_out_list = []
    C_score_list = []

    for lista_1 in list1:
        for lista_2 in list2:
            combination = lista_1 & lista_2
            if len(combination) > 0:
                exon_out, C_scores = make_mock_C_scores(norm_PSI, Ws, lista_1&lista_2, total_cells, mock=mock)
                exon_out_list.append(exon_out)
                C_score_list.append(C_scores)


    psi_key = ['psi_05_10', 'psi_10_20', 'psi_20_30', 'psi_30_40', 'psi_40_50']
    obs_key = ['obs_10_20', 'obs_20_30', 'obs_30_40', 'obs_40_50', 
               'obs_50_60', 'obs_60_70', 'obs_70_80', 'obs_80_90', 'obs_90_100']

    counter = 0
    mock_dict = {}
    for pk in psi_key:
        obs_dict = {}
        for ok in obs_key:

            #fit_alpha, fit_loc, fit_beta = st.gamma.fit(C_score_list[counter])
            #random_data = st.gamma.rvs(fit_alpha, loc=fit_loc, scale=fit_beta, size=1000000)
            random_data = C_score_list[counter]
            obs_dict.update({ok:random_data})
            counter += 1
        mock_dict.update({pk:obs_dict})

    return mock_dict



#######################


def get_C_score_pval_gamma(PSI_tab, norm_PSI, Ws, exon_list, total_cells, mock_dict):
    
    
    exon_out_list = []
    C_list = []
    p_list = []
    
    for exon in tqdm(exon_list):
        psi_mean = PSI_tab.loc[exon].mean()
        obs_mean = PSI_tab.loc[exon].isna().mean()
        
        
        
        exon_df = norm_PSI.loc[exon] # to make things faster
        exon_score = get_C(exon_df, Ws)
        if exon_score >= 0:
            C_list.append(exon_score)
            exon_out_list.append(exon)
            
            if (np.abs(0.5 - psi_mean) > 0.4) and (np.abs(0.5 - psi_mean) <= 0.45):
                pk = 'psi_05_10'
            elif (np.abs(0.5 - psi_mean) > 0.3) and (np.abs(0.5 - psi_mean) <= 0.4):
                pk = 'psi_10_20'
            elif (np.abs(0.5 - psi_mean) > 0.2) and (np.abs(0.5 - psi_mean) <= 0.3):
                pk = 'psi_20_30'
            elif (np.abs(0.5 - psi_mean) > 0.1) and (np.abs(0.5 - psi_mean) <= 0.2):
                pk = 'psi_30_40'
            elif (np.abs(0.5 - psi_mean) <= 0.1):
                pk = 'psi_40_50'
             
            if (obs_mean <= 0.9) and (obs_mean > 0.8):
                ok = 'obs_10_20'
            elif (obs_mean <= 0.8) and (obs_mean > 0.7):
                ok = 'obs_20_30'
            elif (obs_mean <= 0.7) and (obs_mean > 0.6):
                ok = 'obs_30_40'
            elif (obs_mean <= 0.6) and (obs_mean > 0.5):
                ok = 'obs_40_50'
            elif (obs_mean <= 0.5) and (obs_mean > 0.4):
                ok = 'obs_50_60'
            elif (obs_mean <= 0.4) and (obs_mean > 0.3):
                ok = 'obs_60_70'
            elif (obs_mean <= 0.3) and (obs_mean > 0.2):
                ok = 'obs_70_80'
            elif (obs_mean <= 0.2) and (obs_mean > 0.1):
                ok = 'obs_80_90'
            elif (obs_mean <= 0.1):
                ok = 'obs_90_100'
            else:
                ok = 'obs_10_20'
                
            
            random_data = mock_dict[pk][ok]
            
            x = np.sum(random_data > exon_score)
            n = len(random_data)
            pv = (x+1)/(n+1)
            
            p_list.append(pv)
            
    pval_df = pd.DataFrame()
    pval_df['C_score'] = C_list
    pval_df['pval'] = p_list
    pval_df.index = exon_out_list
    return pval_df

     
tiklova_PSI = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/tiklova/skipped_exons_psi.tab',
                         sep='\t', index_col=0) 

tiklova_mrna = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation//tiklova/mrna_counts.tab', sep='\t', index_col=0)

tiklova_mrna_event = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/tiklova/mrna_per_event.tab',
                                  sep='\t', index_col=0)

tiklova_pca = pd.read_csv('/mnt/lareaulab/cfbuenabadn/data_sc_regulation/data_autocorrelation/tiklova/rd_pc2.tab', sep='\t', index_col=0)

discard = [x for x in tiklova_mrna.index if ((x[:3] in ['mt-', 'Gm0', 'Gm1', 'Gm2', 'Gm3', 'Gm4', 'Gm5', 
                                    'Gm6', 'Gm7', 'Gm8', 'Gm9', 'Mir']) or (x[-3:] in ['Rik', '-ps']) or (x [-4:-1] == '-ps'))]

good_genes = [x for x in tiklova_mrna.index if x not in discard]
good_exons = [x for x in tiklova_PSI.index if x.split('_')[0] in good_genes]

tiklova_PSI = tiklova_PSI.loc[good_exons]

print(tiklova_PSI.shape)
print(len(good_exons))

tiklova_exons = tiklova_PSI.index[np.abs(0.5 - tiklova_PSI.mean(axis=1)) <= 0.45] & tiklova_PSI.index[tiklova_PSI.isna().mean(axis=1) < 0.90]

k = 100 #int(round(np.sqrt(len(tiklova_PSI.columns))))

W_tiklova = get_distance_matrix(tiklova_pca, k=k)
tiklova_norm_PSI = get_signature_matrix(tiklova_PSI)

# observed_exons_1 = tiklova_PSI.index[tiklova_PSI[tiklova_pca.index[tiklova_pca.AC==0]].isna().mean(axis=1) <= (1-0.5)]
# observed_exons_2 = tiklova_PSI.index[tiklova_PSI[tiklova_pca.index[tiklova_pca.AC==1]].isna().mean(axis=1) <= (1-0.5)]
# observed_exons_3 = tiklova_PSI.index[tiklova_PSI[tiklova_pca.index[tiklova_pca.AC==2]].isna().mean(axis=1) <= (1-0.5)]
# observed_exons_4 = tiklova_PSI.index[tiklova_PSI[tiklova_pca.index[tiklova_pca.AC==3]].isna().mean(axis=1) <= (1-0.5)]
# observed_exons_5 = tiklova_PSI.index[tiklova_PSI[tiklova_pca.index[tiklova_pca.AC==4]].isna().mean(axis=1) <= (1-0.5)]

# test_exons = []
# for exon in good_exons:
#     exon_counts = 0
#     exon_counts += (exon in observed_exons_1)
#     exon_counts += (exon in observed_exons_2)
#     exon_counts += (exon in observed_exons_3)
#     exon_counts += (exon in observed_exons_4)
#     exon_counts += (exon in observed_exons_5)
#     if exon_counts >= 3:
#         test_exons.append(exon)



tiklova_mock_dict = get_mock_dict(tiklova_PSI.loc[good_exons], tiklova_norm_PSI.loc[good_exons], W_tiklova, mock=10000)

pgamma_tiklova = get_C_score_pval_gamma(tiklova_PSI.loc[good_exons], tiklova_norm_PSI.loc[good_exons], 
                                W_tiklova, good_exons, len(tiklova_PSI.columns), tiklova_mock_dict)

pgamma_tiklova.to_csv('tiklova_GearyC_k100.tab', sep='\t', header=True, index=True)

    
    
