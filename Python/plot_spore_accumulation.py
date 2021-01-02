from __future__ import division
import os, sys, json
import pandas as pd
import numpy as np

import  matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import matplotlib.transforms as transforms
from matplotlib import cm

import phylo_tools as pt

import scipy.stats as stats

import parse_file
import timecourse_utils
import mutation_spectrum_utils


np.random.seed(123456789)

treatments=pt.treatments
replicates = pt.replicates


color_range =  np.linspace(0.0, 1.0, 10)
rgb_blue = cm.get_cmap('Blues')( color_range )
rgb_red = cm.get_cmap('Reds')( color_range )




path_IN = pt.get_path() + '/data/spore_assay/Sporulation_170912_long.txt'
IN = pd.read_csv(path_IN, sep = '\t')
IN = IN.loc[IN['Time_hours'] <= 400]
#d100
IN_0B1_100 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 100)]
IN_2B1_100 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 100)]
IN_mean_0B1_100 = IN_0B1_100['Vegetative_percent'].groupby(IN_0B1_100['Time_hours']).mean().reset_index()
IN_mean_2B1_100 = IN_2B1_100['Vegetative_percent'].groupby(IN_2B1_100['Time_hours']).mean().reset_index()
IN_std_0B1_100 = IN_0B1_100['Vegetative_percent'].groupby(IN_0B1_100['Time_hours']).std().reset_index()
IN_std_2B1_100 = IN_2B1_100['Vegetative_percent'].groupby(IN_2B1_100['Time_hours']).std().reset_index()
# Day 500
IN_0B1_500 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 500)]
IN_2B1_500 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 500)]
IN_mean_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).mean().reset_index()
IN_mean_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).mean().reset_index()
IN_std_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).std().reset_index()
IN_std_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).std().reset_index()

fig = plt.figure(figsize=(6,3))

#plt.scatter(IN_mean_0B1.Time_hours.values, IN_mean_0B1.Vegetative_percent.values, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
plt.plot(IN_mean_0B1_100.Time_hours.values, 1.001- IN_mean_0B1_100.Vegetative_percent.values, \
    'b-',  c=rgb_blue[5])
plt.plot(IN_mean_2B1_100.Time_hours.values, 1.001- IN_mean_2B1_100.Vegetative_percent.values, \
    'b-',  c =rgb_red[5])


plt.errorbar(IN_mean_0B1_100.Time_hours.values, 1.001- IN_mean_0B1_100.Vegetative_percent.values, \
    IN_std_0B1_100.Vegetative_percent.values,  linestyle='None', marker='o', c=rgb_blue[5], elinewidth=1.5, label="1-day WT, day 100",)
plt.errorbar(IN_mean_2B1_100.Time_hours.values, 1.001- IN_mean_2B1_100.Vegetative_percent.values, \
    IN_std_2B1_100.Vegetative_percent.values, linestyle='None', marker='o', c=rgb_red[5], elinewidth=1.5, label="100-day WT, day 100",)


plt.plot(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
    'b-',  c=rgb_blue[9])
plt.plot(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
    'b-',  c = rgb_red[9])

plt.errorbar(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
    IN_std_0B1_500.Vegetative_percent.values,  linestyle='None', marker='o', c=rgb_blue[9], elinewidth=1.5, label="1-day WT, day 500",)
plt.errorbar(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
    IN_std_2B1_500.Vegetative_percent.values, linestyle='None', marker='o', c = rgb_red[9], elinewidth=1.5, label="100-day WT, day 500",)

#plt.title('Bacillus sporulation', fontsize = 24)
plt.xlabel('Time (hours)', fontsize = 14)
plt.ylabel('Percent spores', fontsize = 14)
plt.ylim(0.0008, 1.1)
plt.xlim(-5, 400)
plt.yscale('log')
plt.legend(numpoints=1, prop={'size':8},  loc='lower right', frameon=False)

plt.title(pt.latex_dict['B'], fontsize=14)

fig_name = pt.get_path() + '/figs/spore_assay.pdf'
fig.savefig(fig_name, format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
