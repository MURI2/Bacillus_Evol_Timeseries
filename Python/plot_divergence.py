from __future__ import division
import os, sys, pickle, random
import numpy as np

import  matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import ColorConverter

import scipy.stats as stats
import statsmodels.api as sm

import parse_file
import timecourse_utils
import mutation_spectrum_utils
import phylo_tools as pt

#import get_random_matrix

import phik

np.random.seed(123456789)


permutations_divergence = 10000

treatment_pairs = [['0','1'],['0','2'],['1','2']]

significant_multiplicity_dict = {}
significant_n_mut_dict = {}
gene_size_dict = {}
gene_mean_size_dict = {}
for taxon in pt.taxa:
    significant_multiplicity_dict[taxon] = {}
    significant_n_mut_dict[taxon] = {}
    gene_size_dict[taxon] = {}

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % ('0'+taxon)))
    Ltot = 0
    for gene_name in sorted(convergence_matrix.keys()):
        Lmin=0
        L = max([convergence_matrix[gene_name]['length'],Lmin])
        Ltot += L
    Lavg = Ltot*1.0/len(convergence_matrix.keys())

    gene_mean_size_dict[taxon] = Lavg

    for treatment_idx, treatment in enumerate(pt.treatments):

        significant_multiplicity_taxon_path = pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
        if os.path.exists(significant_multiplicity_taxon_path) == False:
            continue
        significant_multiplicity_taxon = open(significant_multiplicity_taxon_path, "r")
        for i, line in enumerate( significant_multiplicity_taxon ):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            gene_size_dict[taxon][items[0]] = float(items[-5])
            if items[0] not in significant_multiplicity_dict[taxon]:
                significant_multiplicity_dict[taxon][items[0]] = {}

            if items[0] not in significant_n_mut_dict[taxon]:
                significant_n_mut_dict[taxon][items[0]] = {}

            significant_multiplicity_dict[taxon][items[0]][treatment] = float(items[-2])
            significant_n_mut_dict[taxon][items[0]][treatment] = float(items[-4])



def calculate_divergence_correlations_between_taxa():

    sys.stdout.write("Starting divergence tests...\n")

    divergence_dict = {}
    for treatment_idx, treatment in enumerate(pt.treatments):
        all_genes = set(significant_n_mut_dict['B'].keys()) & significant_n_mut_dict['S'].keys()
        result = []
        for gene in all_genes:
            if (treatment in significant_n_mut_dict['B'][gene]) and (treatment in significant_n_mut_dict['S'][gene]):
                result.append((significant_n_mut_dict['B'][gene][treatment], significant_n_mut_dict['S'][gene][treatment], gene ))

        n_x = [int(x[0]) for x in result]
        n_y = [int(x[1]) for x in result]
        gene_names = [x[2] for x in result]

        gene_sizes_taxon = [gene_size_dict['B'][gene_i] for gene_i in gene_names]
        gene_sizes_taxon = np.asarray(gene_sizes_taxon)
        taxon_Lmean = gene_mean_size_dict['B']

        n_matrix = np.asarray([n_x, n_y])
        mult_matrix = n_matrix * (taxon_Lmean / gene_sizes_taxon)
        rel_mult_matrix = mult_matrix/mult_matrix.sum(axis=1)[:,None]
        pearsons_corr = np.corrcoef(rel_mult_matrix[0,:], rel_mult_matrix[1,:])[1,0]
        pearsons_corr_squared = pearsons_corr**2
        pearsons_corr_null = []
        pearsons_corr_squared_null = []
        for k in range(permutations_divergence):

            if (k % 2000 == 0) and (k>0):

                sys.stdout.write("%d iterations\n" % (k))

            n_matrix_random = phik.simulation.sim_2d_data_patefield(n_matrix)
            mult_matrix_random = n_matrix_random * (taxon_Lmean / gene_sizes_taxon)
            rel_mult_matrix_random = mult_matrix_random/mult_matrix_random.sum(axis=1)[:,None]
            pearsons_corr_random = np.corrcoef(rel_mult_matrix_random[0,:], rel_mult_matrix_random[1,:])[1,0]
            pearsons_corr_squared_random = pearsons_corr_random**2

            pearsons_corr_null.append(pearsons_corr_random)
            pearsons_corr_squared_null.append(pearsons_corr_squared_random)

        pearsons_corr_null = np.asarray(pearsons_corr_null)
        pearsons_corr_squared_null = np.asarray(pearsons_corr_squared_null)

        pearsons_corr_null_abs = np.absolute(pearsons_corr_null)
        pearsons_corr_squared_null_abs = np.absolute(pearsons_corr_squared_null)

        Z_corr_squared = (pearsons_corr_squared - np.mean(pearsons_corr_squared_null)) / np.std(pearsons_corr_squared_null)
        Z_corr = (pearsons_corr - np.mean(pearsons_corr_null)) / np.std(pearsons_corr_null)

        P_corr_squared = (len(pearsons_corr_squared_null_abs[pearsons_corr_squared_null_abs < np.absolute(pearsons_corr_squared)]) +1) / (permutations_divergence+1)
        P_corr = (len(pearsons_corr_null_abs[pearsons_corr_null_abs < np.absolute(pearsons_corr)]) +1) / (permutations_divergence+1)


        divergence_dict[treatment] = {}
        divergence_dict[treatment]['pearsons_corr_squared'] = pearsons_corr_squared
        divergence_dict[treatment]['P_value_corr_squared'] = P_corr_squared
        divergence_dict[treatment]['Z_corr_squared'] = Z_corr_squared

        divergence_dict[treatment]['pearsons_corr'] = pearsons_corr
        divergence_dict[treatment]['P_value_corr'] = P_corr
        divergence_dict[treatment]['Z_corr'] = Z_corr

        sys.stdout.write("%d-day, WT vs. spo0A: rho=%f, P=%f, Z=%f\n" % (10**int(treatment), pearsons_corr, P_corr, Z_corr))

    sys.stdout.write("Dumping pickle......\n")
    with open(pt.get_path()+'/data/divergence_pearsons_between_taxa.pickle', 'wb') as handle:
        pickle.dump(divergence_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stdout.write("Done!\n")






def calculate_divergence_correlations_between_treatments():

    sys.stdout.write("Starting divergence tests...\n")

    divergence_dict = {}

    for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):

        treatment_pair_set = (treatment_pair[0], treatment_pair[1])

        divergence_dict[treatment_pair_set] = {}

        for taxon in pt.taxa:

            #result = [(x[treatment_pair[0]],x[treatment_pair[1]]) for x in significant_multiplicity_dict[taxon].values() if (treatment_pair[0] in x) and (treatment_pair[1] in x)]
            #result = [(x[treatment_pair[0]],x[treatment_pair[1]], x) for x in significant_n_mut_dict[taxon].values() if (treatment_pair[0] in x) and (treatment_pair[1] in x)]
            result = [(dicts[treatment_pair[0]],dicts[treatment_pair[1]], keys) for keys, dicts in significant_n_mut_dict[taxon].items() if (treatment_pair[0] in dicts) and (treatment_pair[1] in dicts)]

            n_x = [int(x[0]) for x in result]
            n_y = [int(x[1]) for x in result]
            gene_names = [x[2] for x in result]

            gene_sizes_taxon_treatment_pair = [gene_size_dict[taxon][gene_i] for gene_i in gene_names]
            gene_sizes_taxon_treatment_pair = np.asarray(gene_sizes_taxon_treatment_pair)
            taxon_Lmean = gene_mean_size_dict[taxon]

            n_matrix = np.asarray([n_x, n_y])
            mult_matrix = n_matrix * (taxon_Lmean / gene_sizes_taxon_treatment_pair)
            rel_mult_matrix = mult_matrix/mult_matrix.sum(axis=1)[:,None]
            pearsons_corr = np.corrcoef(rel_mult_matrix[0,:], rel_mult_matrix[1,:])[1,0]
            pearsons_corr_squared = pearsons_corr**2

            pearsons_corr_null = []
            pearsons_corr_squared_null = []
            for k in range(permutations_divergence):

                if (k % 2000 == 0) and (k>0):

                    sys.stdout.write("%d iterations\n" % (k))

                n_matrix_random = phik.simulation.sim_2d_data_patefield(n_matrix)
                mult_matrix_random = n_matrix_random * (taxon_Lmean / gene_sizes_taxon_treatment_pair)
                rel_mult_matrix_random = mult_matrix_random/mult_matrix_random.sum(axis=1)[:,None]
                pearsons_corr_random = np.corrcoef(rel_mult_matrix_random[0,:], rel_mult_matrix_random[1,:])[1,0]
                pearsons_corr_squared_random = pearsons_corr_random**2

                pearsons_corr_null.append(pearsons_corr_random)
                pearsons_corr_squared_null.append(pearsons_corr_squared_random)

            pearsons_corr_null = np.asarray(pearsons_corr_null)
            pearsons_corr_squared_null = np.asarray(pearsons_corr_squared_null)

            pearsons_corr_null_abs = np.absolute(pearsons_corr_null)
            pearsons_corr_squared_null_abs = np.absolute(pearsons_corr_squared_null)

            Z_corr_squared = (pearsons_corr_squared - np.mean(pearsons_corr_squared_null)) / np.std(pearsons_corr_squared_null)
            Z_corr = (pearsons_corr - np.mean(pearsons_corr_null)) / np.std(pearsons_corr_null)

            P_corr_squared = (len(pearsons_corr_squared_null_abs[pearsons_corr_squared_null_abs < np.absolute(pearsons_corr_squared)]) +1) / (permutations_divergence+1)
            P_corr = (len(pearsons_corr_null_abs[pearsons_corr_null_abs < np.absolute(pearsons_corr)]) +1) / (permutations_divergence+1)

            divergence_dict[treatment_pair_set][taxon] = {}
            divergence_dict[treatment_pair_set][taxon]['pearsons_corr_squared'] = pearsons_corr_squared
            divergence_dict[treatment_pair_set][taxon]['P_value_corr_squared'] = P_corr_squared
            divergence_dict[treatment_pair_set][taxon]['Z_corr_squared'] = Z_corr_squared

            divergence_dict[treatment_pair_set][taxon]['pearsons_corr'] = pearsons_corr
            divergence_dict[treatment_pair_set][taxon]['P_value_corr'] = P_corr
            divergence_dict[treatment_pair_set][taxon]['Z_corr'] = Z_corr

            sys.stdout.write("%d vs %d-day, %s: rho=%f, P=%f, Z=%f\n" % (10**int(treatment_pair[0]), 10**int(treatment_pair[1]), taxon, pearsons_corr, P_corr, Z_corr))

    sys.stdout.write("Dumping pickle......\n")
    with open(pt.get_path()+'/data/divergence_pearsons_between_treatments.pickle', 'wb') as handle:
        pickle.dump(divergence_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.stdout.write("Done!\n")


#calculate_divergence_correlations_between_taxa()
#calculate_divergence_correlations_between_treatments()

with open(pt.get_path()+'/data/divergence_pearsons_between_taxa.pickle', 'rb') as handle:
    divergence_dict_between_taxa = pickle.load(handle)

with open(pt.get_path()+'/data/divergence_pearsons_between_treatments.pickle', 'rb') as handle:
    divergence_dict_between_treatments = pickle.load(handle)




gs = gridspec.GridSpec(nrows=2, ncols=1)

fig = plt.figure(figsize = (10, 13))
ax_between_taxa = fig.add_subplot(gs[0, 0])
ax_between_treatments = fig.add_subplot(gs[1, 0])

ax_between_taxa.text(-0.1, 1.07, pt.sub_plot_labels[0], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_between_taxa.transAxes)
ax_between_treatments.text(-0.1, 1.07, pt.sub_plot_labels[1], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_between_treatments.transAxes)



for treatment_idx, treatment in enumerate(divergence_dict_between_taxa.keys()):

    Z_corr = divergence_dict_between_taxa[treatment]['Z_corr']

    marker_style = dict(color=pt.get_colors(treatment),
                        markerfacecoloralt='white',
                        markerfacecolor=pt.get_colors(treatment),
                        mew=3)

    ax_between_taxa.plot(treatment, Z_corr, markersize = 30, marker = 'o',  \
        linewidth=0.4,  alpha=1, fillstyle='left', zorder=2 , **marker_style)



ax_between_taxa.set_xlim([-0.5,2.5])
ax_between_taxa.set_ylim([-33,5])


ax_between_taxa.axhline( y=0, color='k', lw=3, linestyle=':', alpha = 1, zorder=1)
ax_between_taxa.text(0.125, 0.91, 'Convergence', fontsize=15, fontweight='bold', ha='center', va='center', transform=ax_between_taxa.transAxes)
ax_between_taxa.text(0.115, 0.83, 'Divergence', fontsize=15 , fontweight='bold', ha='center', va='center', transform=ax_between_taxa.transAxes)
ax_between_taxa.set_ylabel("Standardized correlation, "+ r'$Z_{\rho}$' , fontsize = 16)

ax_between_taxa.set_xticks([0, 1, 2])
ax_between_taxa.set_xticklabels(['1-day', '10-days', '100-days'], fontweight='bold', fontsize=18 )



ax_between_taxa.set_title("B. subtilis " + r'$\mathbf{WT}$' + " vs. B. subtilis "  + r'$\mathbf{\Delta spo0A}$', style='italic', fontsize=16, fontweight='bold')




count = 0
for taxon in pt.taxa:
    for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):

        Z_corr = divergence_dict_between_treatments[tuple(treatment_pair)][taxon]['Z_corr']


        marker_style = dict(color='k', marker='o',
                markerfacecoloralt=pt.get_colors(treatment_pair[1]),
                markerfacecolor=pt.get_colors(treatment_pair[0]),
                mew=2)

        ax_between_treatments.plot(count, Z_corr, markersize = 28,   \
            linewidth=2,  alpha=1, zorder=3, fillstyle='left', **marker_style)


        count+=1





ax_between_treatments.set_xticks([0, 1, 2, 3, 4, 5])

ax_between_treatments.set_xticklabels(['1-day vs.\n10-days', '1-day vs.\n100-days', '1-day vs.\n100-days', '1-day vs.\n10-days', '1-day vs.\n100-days', '1-day vs.\n100-days'], fontweight='bold', fontsize=13 )




ax_between_treatments.axhline( y=0, color='k', lw=3, linestyle=':', alpha = 1, zorder=1)
ax_between_treatments.axvline( x=2.5, color='k', lw=3, linestyle='-', alpha = 1, zorder=1)

ax_between_treatments.text(0.125, 0.86, 'Convergence', fontsize=15, fontweight='bold', ha='center', va='center', transform=ax_between_treatments.transAxes)
ax_between_treatments.text(0.115, 0.77, 'Divergence', fontsize=15 , fontweight='bold', ha='center', va='center', transform=ax_between_treatments.transAxes)
ax_between_treatments.set_xlim([-0.5,5.5])
ax_between_treatments.set_ylim([-18,4])

ax_between_treatments.set_ylabel("Standardized correlation, "+ r'$Z_{\rho}$' , fontsize = 16)


ax_between_treatments.text(0.25, -0.155, "B. subtilis " + r'$\mathbf{WT}$', style='italic', fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_between_treatments.transAxes)
ax_between_treatments.text(0.75, -0.155, "B. subtilis "  + r'$\mathbf{\Delta spo0A}$', style='italic', fontsize=16, fontweight='bold', ha='center', va='center', transform=ax_between_treatments.transAxes)


fig.subplots_adjust(hspace=0.15,wspace=0.2) #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/divergence.pdf"
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



sys.stderr.write("Done with figure!\n")
