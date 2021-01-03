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

from scipy.special import gammaln

from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

from sklearn.decomposition import PCA


np.random.seed(123456789)

treatments=pt.treatments
replicates = pt.replicates


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']


#legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= r'$\overline{M}_{WT} (t)$'),
#                   Line2D([0], [0], ls=':', color='k', lw=1.5, label= r'$\overline{M}_{\Delta \mathrm{spo0A}} (t)$')]
#ax_t_vs_M.legend(handles=legend_elements, loc='lower right', fontsize=12)



















def allele_survival():
    frequencies = np.linspace(0.1,0.9,201)
    df = 0.05
    fstar = 0.5

    #tstar = 20025

    null_num_in_bin = np.zeros_like(frequencies)
    null_avg_f = np.zeros_like(frequencies)

    origin_fixation_times = {}

    taxa = ['B', 'S']

    #fig = plt.figure()
    fig = plt.figure(figsize = (12, 6))

    for treatment in ['0']:
        for taxon_idx, taxon in enumerate(taxa):
            ax_i = plt.subplot2grid((1, 2), (0, taxon_idx), colspan=1)
            for replicate in ['1', '2', '3', '4', '5']:
                population = treatment + taxon + replicate

                if population in pt.populations_to_ignore:
                    continue

                num_runs = []

                num_in_bin = np.zeros_like(null_num_in_bin)
                avg_f = np.zeros_like(null_avg_f)

                origin_fixation_times[population] = ([],[],[])

                sys.stderr.write("Processing fixation probabilities for %s...\n" % population)

                # load mutations
                mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
                population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
                state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

                num_runs = []

                for mutation_idx in range(0,len(mutations)):

                    location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]

                    state_Ls = state_trajectories[mutation_idx]

                    good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

                    freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

                    masked_times = times[good_idxs]
                    masked_freqs = freqs[good_idxs]
                    masked_depths = depths[good_idxs]
                    masked_state_Ls = state_Ls[good_idxs]

                    # Estimate appearance and fixation times
                    if masked_state_Ls[-1] in parse_file.well_mixed_fixed_states:
                        t0,tf,dt = timecourse_utils.calculate_appearance_fixation_time_from_hmm(masked_times, masked_freqs, masked_state_Ls)
                        origin_fixation_times[population][0].append(t0)
                        origin_fixation_times[population][1].append(tf)
                        origin_fixation_times[population][2].append(dt)

                    # Now split the trajectory into independent polymorphic runs
                    # each of which contains a single final state (unless end of experiment)
                    independent_runs = timecourse_utils.split_well_mixed_hmm(masked_times,masked_freqs, masked_state_Ls)
                    num_runs.append(len(independent_runs))

                    for run_idxs in independent_runs:

                        if len(run_idxs)<2:
                            # need at least one polymorphic state and one final state
                            continue
                        # initial time
                        t = masked_times[run_idxs[0]]

                        # get final state
                        final_state = masked_state_Ls[run_idxs[-1]]

                        # get frequency of parent clade during run
                        #if final_state == parse_file.well_mixed_hmm_states['F'] or final_state==parse_file.well_mixed_hmm_states['P']:
                        #    parent_freqs = masked_freqs
                        #else:
                        #parent_freqs = masked_freqs
                        #elif final_state == parse_file.clade_hmm_states['F'] or final_state == parse_file.well_mixed_hmm_states['P']:
                        #    parent_freqs = masked_fmajors
                        #elif final_state == parse_file.clade_hmm_states['F'] or final_state == parse_file.clade_hmm_states['Pm']:
                        #    parent_freqs = masked_fmajors
                        #else:
                        #    parent_freqs = masked_fextincts


                        # don't neet to renormalize the freqs because we have no population structure

                        #renormalized_freqs = np.clip(masked_freqs[run_idxs]/parent_freqs[run_idxs],0,1)

                        # get fixation weight
                        if final_state in parse_file.well_mixed_fixed_states:
                            fixation_weight = 1.0
                        elif final_state in parse_file.well_mixed_extinct_states:
                            fixation_weight = 0.0
                        else:
                            fixation_weight = masked_freqs[-1] > fstar

                        individual_bin_weights = np.exp(-np.power((masked_freqs[:,None]-frequencies[None,:])/df,2))
                        individual_bin_weights *= (individual_bin_weights>0.01)

                        bin_weights = individual_bin_weights.sum(axis=0)
                        fixation_weights = (individual_bin_weights*fixation_weight).sum(axis=0)

                        #num_in_bin += bin_weights
                        num_in_bin += bin_weights
                        avg_f += fixation_weights

                        null_num_in_bin += bin_weights
                        null_avg_f += fixation_weights


                avg_f = avg_f/(num_in_bin+(num_in_bin<1))

                ax_i.plot(frequencies[(num_in_bin>=1)], avg_f[(num_in_bin>=1)], '.-', color= pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), fillstyle=pt.plot_species_fillstyle(taxon), alpha=0.9, linewidth=2.5, markersize = 8, zorder=2)
                ax_i.plot([0, 1], [0, 1], '--', c = 'black', alpha=0.9, markersize = 10, zorder=1)

                ax_i.set_xlim([0,1])
                ax_i.set_ylim([0,1])
                ax_i.set_xlabel('Allele frequency', fontsize = 14)
                ax_i.set_ylabel(r'$\mathrm{Pr}\left [ \mathrm{survival} \right ]$', fontsize = 14)

                ax_i.spines['top'].set_visible(False)

                line, = ax_i.plot([0,1],[1,1],'k:',linewidth=5)
                line.set_dashes((0.5, 0.5))

                ax_i.set_title( latex_dict[taxon], fontweight='bold', fontsize=17)

                #fig.suptitle(latex_dict[strain], fontsize=28, fontweight='bold')

            if (taxon_idx == 0):
                legend_elements = [Line2D([0], [0], ls='--', color='k', lw=1.5, label= 'Quasi-neutral'),
                                   Line2D([0], [0], ls=':', color='k', lw=1.5, label= 'Hitchhiking' )]
                ax_i.legend(handles=legend_elements, loc='upper left', fontsize=12)

    fig_name = pt.get_path() + '/figs/freq_vs_prob_fixation.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()





def likelihood_plot(ntot_subsample=50, subsamples=10000):
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



    legend_elements = [Line2D([0], [0], color = 'none', marker='o', label=latex_dict['B'],
                        markerfacecolor='k', markersize=13),
                    Line2D([0], [0], marker='o', color='none', label=latex_dict['S'],
                        markerfacecolor='w', markersize=13, markeredgewidth=2),
                    Line2D([0], [0], color = 'none', marker='s', label=latex_dict['B'] + ', sub-sampled',
                        markerfacecolor='k', markersize=13),
                    Line2D([0], [0], marker='s', color='none', label=latex_dict['S'] + ', sub-sampled',
                        markerfacecolor='w', markersize=13, markeredgewidth=2)]
    # Create the figure
    plt.legend(handles=legend_elements, loc='upper left')

    fig.subplots_adjust() #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/G_score_subsample.pdf"
    fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
