from __future__ import division
import os, sys, json
import pandas as pd
import numpy as np

import  matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms


import phylo_tools as pt

import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.stats.api as sms
import statsmodels.formula.api as smf
from statsmodels.compat import lzip

#from sklearn.metrics import pairwise_distances
#from skbio.stats.ordination import pcoa

import parse_file
import timecourse_utils
import mutation_spectrum_utils


np.random.seed(123456789)

treatments=pt.treatments
replicates = pt.replicates



df = pd.read_csv(pt.get_path() + '/data/competition_2018-1-9-count.txt', sep = '\t')
rows_to_keep = []
for index, row in df.iterrows():
    wt = row['WT']
    spo0a = row['spoA']
    if (wt == 'TMTC') or (spo0a == 'TMTC') :
        continue
    wt_spo0a = int(wt) + int(spo0a)
    if (wt_spo0a < 30):
        continue
    rows_to_keep.append(index)


def cfus_ml(column, conc):
    return column * (10 ** (int(conc) * -1))


df_keep = df.iloc[rows_to_keep]
df_keep.WT = df_keep.WT.astype(int)
df_keep.spoA = df_keep.spoA.astype(int)

df_keep['WT_cfus_ml'] = df_keep.apply(lambda x: cfus_ml(x.WT, x.Concentration), axis=1)
df_keep['spoA_cfus_ml'] = df_keep.apply(lambda x: cfus_ml(x.spoA, x.Concentration), axis=1)

df_keep = df_keep.drop(['Concentration', 'WT', 'spoA', 'Rep'], 1)
df_keep = df_keep.groupby(['Day','Flask'], as_index=False).mean()

flask_1 = df_keep.loc[df_keep['Flask'] == 1]
flask_2 = df_keep.loc[df_keep['Flask'] == 2]
flask_3 = df_keep.loc[df_keep['Flask'] == 3]

# mean initia CFU counts from notebook

relative_fitness_1 = np.log((flask_1['spoA_cfus_ml'].values / flask_1['WT_cfus_ml'].values) * ( 0.51/0.49))
relative_fitness_2 = np.log((flask_2['spoA_cfus_ml'].values / flask_2['WT_cfus_ml'].values) * ( 0.48/0.52))
relative_fitness_3 = np.log((flask_3['spoA_cfus_ml'].values / flask_3['WT_cfus_ml'].values) * ( 0.54/0.46))

relative_fitness_per_time_1 = np.log((flask_1['spoA_cfus_ml'].values / flask_1['WT_cfus_ml'].values) * ( 0.51/0.49)) /  flask_1['Day'].values
relative_fitness_per_time_2 = np.log((flask_2['spoA_cfus_ml'].values / flask_2['WT_cfus_ml'].values) * ( 0.48/0.52)) /  flask_2['Day'].values
relative_fitness_per_time_3 = np.log((flask_3['spoA_cfus_ml'].values / flask_3['WT_cfus_ml'].values) * ( 0.54/0.46)) /  flask_3['Day'].values


zipped_relative = list(zip(list(relative_fitness_1), list(relative_fitness_2), list(relative_fitness_3)))
relative_mean = (relative_fitness_1 + relative_fitness_2 + relative_fitness_3) / 3
relative_se_list = []
for i , item in enumerate(zipped_relative):
    relative_se_list.append(2*np.std(np.asarray(item)) / np.sqrt(len(item)))

zipped_relative_time = list(zip(list(relative_fitness_per_time_1), list(relative_fitness_per_time_2), list(relative_fitness_per_time_3)))
relative_time_mean = (relative_fitness_per_time_1 + relative_fitness_per_time_2 + relative_fitness_per_time_3) / 3
relative_time_se_list = []
for i , item in enumerate(zipped_relative_time):
    relative_time_se_list.append(2*np.std(np.asarray(item)) / np.sqrt(len(item)))




fig = plt.figure(figsize = (10, 9))

ax_relative = plt.subplot2grid((2, 1), (0, 0), colspan=1)
ax_time_relative = plt.subplot2grid((2, 1), (1, 0), colspan=1)

ax_relative.axhline(y=0, color='grey', linestyle='--', lw = 3)
ax_relative.errorbar(flask_1['Day'].values, relative_mean, relative_se_list, linestyle='-', marker='o', lw = 3)
ax_relative.set_ylim(-5, 5)
ax_relative.set_ylabel(pt.latex_dict['S'] + '\nrelative fitness at ' + r'$t$' + ', ' + r'$X(t)$'  , fontsize = 16)
ax_relative.set_xscale('log', base=10)

ax_relative.text(-0.1, 1.07, pt.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_relative.transAxes)


ax_time_relative.axhline(y=0, color='grey', linestyle='--', lw = 3)
ax_time_relative.errorbar(flask_1['Day'].values, relative_time_mean, relative_time_se_list, linestyle='-', marker='o', lw = 3)
ax_time_relative.set_ylim(-0.35, 0.75)
ax_time_relative.set_ylabel(pt.latex_dict['S'] + '\nrelative fitness, ' + r'$\Delta X$'  , fontsize = 16)
ax_time_relative.set_xscale('log', base=10)

ax_time_relative.text(-0.1, 1.07, pt.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_time_relative.transAxes)


ax_time_relative.set_xlabel('Days, ' + r'$t$', fontsize = 20)

fig_name = pt.get_path() + '/figs/fitness_spo0a.jpg'
fig.savefig(fig_name, format='jpg', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
