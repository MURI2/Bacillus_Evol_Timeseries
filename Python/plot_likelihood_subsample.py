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


ntot_subsample=50
subsamples=10000
# ntot_subsample minimum number of mutations

G_subsample_dict = {}

G_all_mutations_dict = {}

for taxon in ['B', 'S']:

    for treatment in treatments:

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        populations = [treatment+taxon + replicate for replicate in replicates ]

        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)

        G_subsample_list = []
        for i in range(subsamples):

            G_subsample = mutation_spectrum_utils.calculate_subsampled_total_parallelism(gene_parallelism_statistics, ntot_subsample=ntot_subsample)

            G_subsample_list.append(G_subsample)

        G_subsample_list.sort()

        G_subsample_dict[treatment+taxon] = G_subsample_list


        G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

        G_all_mutations_dict[treatment+taxon] = G


fig = plt.figure(figsize = (8, 8))

taxon_xaxis_dict = {'B': -0.1, 'S':0.1}

for treatment in treatments:

    for taxon in ['B', 'S']:

        G_subsample_mean = np.mean(G_subsample_dict[treatment+taxon])
        G_subsample_025 = G_subsample_dict[treatment+taxon][ int( 0.025 * subsamples)  ]
        G_subsample_975 = G_subsample_dict[treatment+taxon][ int( 0.975 * subsamples)  ]

        #xerr1 = [ [z_lclb_mpd_null_mean - lclb_mpd_025, z_lcpl_mpd_null_mean - lcpl_mpd_025, z_hclb_mpd_null_mean - hclb_mpd_025, z_hcpl_mpd_null_mean - hcpl_mpd_025 ] ,
        #        [lclb_mpd_975 - z_lclb_mpd_null_mean, lcpl_mpd_975 - z_lcpl_mpd_null_mean, hclb_mpd_975 - z_hclb_mpd_null_mean, hcpl_mpd_975 -z_hcpl_mpd_null_mean ]]

        plt.errorbar(int(treatment) + taxon_xaxis_dict[taxon], G_subsample_mean, yerr = [ [G_subsample_mean-G_subsample_025], [ G_subsample_975-G_subsample_mean]], \
                fmt = 'o', alpha = 1, barsabove = True, marker = 's', \
                mfc = 'white', mec = 'white', lw=3.5, c = 'k', zorder=1, ms=17)

        plt.scatter(int(treatment) + taxon_xaxis_dict[taxon], G_subsample_mean, marker='s', s = 250, \
            linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), alpha=1, zorder=2)


        plt.scatter(int(treatment) + taxon_xaxis_dict[taxon], G_all_mutations_dict[treatment+taxon], marker=pt.plot_species_marker(taxon), s = 250, \
            linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), alpha=1, zorder=2)

plt.xlabel("Transfer time (days)", fontsize = 20)

plt.xticks((0,1,2), ('1', '10', '100'), fontsize=14  )
plt.rc('ytick',labelsize=12)

plt.ylim([1.2,6.2])


plt.ylabel("Net increase in log-likelihood, " r'$\Delta \ell$' , fontsize = 20)



legend_elements = [Line2D([0], [0], color = 'none', marker='o', label=pt.latex_dict['B'],
                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], marker='o', color='none', label=pt.latex_dict['S'],
                    markerfacecolor='w', markersize=13, markeredgewidth=2),
                Line2D([0], [0], color = 'none', marker='s', label=pt.latex_dict['B'] + ', sub-sampled',
                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], marker='s', color='none', label=pt.latex_dict['S'] + ', sub-sampled',
                    markerfacecolor='w', markersize=13, markeredgewidth=2)]
# Create the figure
plt.legend(handles=legend_elements, loc='upper left')

fig.subplots_adjust() #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/G_score_subsample.jpg"
fig.savefig(fig_name, format='jpg', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
