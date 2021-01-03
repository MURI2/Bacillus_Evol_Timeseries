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



sys.stderr.write("Loading mutation data...\n")

mutation_trajectories = {}
fixed_mutation_trajectories = {}
delta_mutation_trajectories = {}
#transit_times = {}
taxa = ['B', 'S']

for treatment in treatments:
    for taxon in taxa:
        for replicate in replicates:

            population = treatment + taxon + replicate
            sys.stderr.write("Processing %s...\t" % population)
            times, Ms, fixed_Ms = parse_file.get_mutation_fixation_trajectories(population)
            fixed_mutation_trajectories[population] = (times, fixed_Ms)
            mutation_trajectories[population] = (times,np.log10(Ms))
            delta_mutation_trajectories[population] = (times[1:], np.log10(Ms[1:]/Ms[:-1] ))

            sys.stderr.write("analyzed %d mutations!\n" % len(Ms))



fig = plt.figure(figsize = (10, 10))

column_count = 0

for treatment in treatments:

    ax_t_vs_M = plt.subplot2grid((3, 3), (0, column_count), colspan=1)

    ax_t_vs_delta_M = plt.subplot2grid((3, 3), (1, column_count), colspan=1)

    ax_M_vs_F = plt.subplot2grid((3, 3), (2, column_count), colspan=1)

    ax_t_vs_M.text(-0.1, 1.07, pt.sub_plot_labels[column_count], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_t_vs_M.transAxes)
    ax_t_vs_delta_M.text(-0.1, 1.07, pt.sub_plot_labels[column_count+3], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_t_vs_delta_M.transAxes)
    ax_M_vs_F.text(-0.1, 1.07, pt.sub_plot_labels[column_count+6], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_M_vs_F.transAxes)

    #ax_M_vs_F.plot([0,300],[0,300],'--',linewidth=1,color='k', zorder=1)

    for taxon_i, taxon in enumerate(taxa):

        treatment_taxon_populations = []

        for replicate in replicates:

            population = treatment + taxon + replicate

            Mts,Ms = mutation_trajectories[population]
            fixed_Mts, fixed_Ms = fixed_mutation_trajectories[population]
            deltaMts, deltaMs = delta_mutation_trajectories[population]

            ax_t_vs_M.plot(Mts, 10**Ms, 'o-',color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            ax_t_vs_M.set_yscale('log', base=10)
            ax_t_vs_M.tick_params(axis='x', labelsize=8)

            # back transform to format plot axes
            ax_t_vs_delta_M.plot(deltaMts, 10**deltaMs, color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon))
            ax_t_vs_delta_M.set_yscale('log', base=10)

            if isinstance(fixed_Ms, np.floating):
                fixed_Ms = [fixed_Ms]*len(fixed_Mts)
            ax_M_vs_F.plot(fixed_Mts, fixed_Ms, 'o-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=1, markersize=7,linewidth=3, markeredgewidth=1.5, zorder=1)
            #ax_M_vs_F.set_xlabel('Days, ' + r'$t$', fontsize = 12)

            treatment_taxon_populations.append(population)

        avg_Mts, avg_Ms = timecourse_utils.average_trajectories([mutation_trajectories[population] for population in treatment_taxon_populations])

        avg_deltaMts, avg_deltaMs = timecourse_utils.average_trajectories([delta_mutation_trajectories[population] for population in treatment_taxon_populations])

        if taxon == 'B':
            ls = '--'
        else:
            ls = ':'

        ax_t_vs_delta_M.axhline(y=1, c='grey', linestyle=':', lw=3, zorder=1)
        ax_t_vs_M.plot(avg_Mts, 10**avg_Ms, ls,color='k', marker=" ", alpha=1, linewidth=4, zorder=2)
        ax_t_vs_delta_M.plot(avg_deltaMts, 10**avg_deltaMs, ls,color='k', marker=" ", alpha=1, linewidth=4, zorder=2)


        if (taxon_i == 0) and (column_count==0):
            legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= r'$\overline{M}_{WT} (t)$'),
                               Line2D([0], [0], ls=':', color='k', lw=1.5, label= r'$\overline{M}_{\Delta \mathrm{spo0A}} (t)$')]
            ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=12)

    ax_t_vs_M.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

    if treatment == '2':
        ax_M_vs_F.yaxis.set_major_locator(MaxNLocator(integer=True))

    if column_count == 0:

        ax_t_vs_M.set_ylabel('Mutations, ' + r'$M(t)$', fontsize = 15)
        ax_M_vs_F.set_ylabel('Fixed mutations', fontsize = 15)
        ax_t_vs_delta_M.set_ylabel('Change in mutations,\n' + r'$M(t)/M(t-1)$', fontsize = 15)

        legend_elements = [Line2D([0], [0], color = 'none', marker='o', label=pt.latex_dict['B'],
                            markerfacecolor='k', markersize=10),
                        Line2D([0], [0], marker='o', color='none', label=pt.latex_dict['S'],
                            markerfacecolor='w', markersize=10, markeredgewidth=2)]

        ax_M_vs_F.legend(handles=legend_elements, loc='upper left', fontsize=8)


    ax_t_vs_M.set_ylim([0.008 , 300])
    ax_t_vs_delta_M.set_ylim([0.1 , 70])
    #ax_M_vs_F.set_ylim([-0.1 , 44])

    column_count += 1

fig.text(0.53, 0.02, 'Days, ' + r'$t$', ha='center', fontsize=28)


fig_name = pt.get_path() + '/figs/mutation_accumulation_B_S.jpg'
fig.savefig(fig_name, format='jpg',  bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
