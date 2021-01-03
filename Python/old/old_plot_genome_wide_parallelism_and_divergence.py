from __future__ import division
import os, sys
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




np.random.seed(123456789)

subsamples=10000




def calculate_parallelism_statistics_partition(taxon, treatment, fmax_partition=0.5):

    convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

    populations = [treatment+taxon + replicate for replicate in pt.replicates ]

    significant_genes = []

    significant_multiplicity_taxon_path = pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
    #if os.path.exists(significant_multiplicity_taxon_path) == False:
    #    continue
    significant_multiplicity_taxon = open(significant_multiplicity_taxon_path, "r")
    for i, line in enumerate( significant_multiplicity_taxon ):
        if i == 0:
            continue
        line = line.strip()
        items = line.split(",")
        if items[0] not in significant_genes:
            significant_genes.append(items[0])

    # Now calculate gene counts
    #Ltot = 0
    #Ngenes = 0
    ntot_less = 0
    ntot_greater = 0
    fmax_true_false = []
    positions = [0]
    n_greater_all = []
    n_less_all = []
    for gene_name in sorted(convergence_matrix.keys()):

        if gene_name not in significant_genes:
            continue

        #L = max([convergence_matrix[gene_name]['length'],Lmin])
        n_less = 0
        n_greater = 0
        #num_pops = 0

        for population in populations:

            # filter by cutoff for maximum allele Frequency
            convergence_matrix_mutations_population_filtered_greater = [k for k in convergence_matrix[gene_name]['mutations'][population] if (k[-1] >= fmax_partition ) ]
            convergence_matrix_mutations_population_filtered_less = [k for k in convergence_matrix[gene_name]['mutations'][population] if (k[-1] < fmax_partition ) ]

            new_muts_greater = len(convergence_matrix_mutations_population_filtered_greater)
            new_muts_less = len(convergence_matrix_mutations_population_filtered_less)

            #fmax_greater =

            if (new_muts_greater > 0) and (new_muts_less > 0):

                n_greater += new_muts_greater
                n_less += new_muts_less

                #num_pops += 1
                #n += new_muts
                #for t,l,f,f_max in convergence_matrix_mutations_population_filtered_greater:
                #    times.append(t)
                #    # get maximum allele frequency

        if (n_greater == 0) or (new_muts_less == 0):
            continue

        fmax_true_false.extend([True] * n_greater)
        fmax_true_false.extend([False] * n_less)

        ntot_less += n_less
        ntot_greater += n_greater


        n_greater_all.append(n_greater)
        n_less_all.append(n_less)


        positions.append(ntot_less+ntot_greater)


    n_less_all = np.asarray(n_less_all)
    n_greater_all = np.asarray(n_greater_all)
    n_all = n_less_all + n_greater_all

    positions = np.asarray(positions)
    fmax_true_false = np.asarray(fmax_true_false)


    ntot = ntot_less + ntot_greater

    likelihood_partition = sum((n_less_all*np.log((n_less_all*ntot)/(ntot_less*n_all) )) + (n_greater_all*np.log((n_greater_all*ntot)/(ntot_greater*n_all) )))

    print("likelihood", fmax_partition, likelihood_partition)

    null_likelihood_partition = []

    for i in range(1000):

        fmax_true_false_permute = np.random.permutation(fmax_true_false)

        #print(fmax_true_false_permute)

        n_less_permute_all = []
        n_greater_permute_all = []

        for gene_idx in range(len(positions)-1):

            gene_fmax_permute = fmax_true_false_permute[positions[gene_idx]:positions[gene_idx+1]]

            n_greater_permute = len(gene_fmax_permute[gene_fmax_permute==True])
            n_less_permute = len(gene_fmax_permute[gene_fmax_permute!=True])

            if (n_greater_permute > 0 ) and (n_less_permute > 0):

                n_greater_permute_all.append(n_greater_permute)
                n_less_permute_all.append(n_less_permute)

        n_less_permute_all = np.asarray(n_less_permute_all)
        n_greater_permute_all = np.asarray(n_greater_permute_all)

        n_permute_all = n_less_permute_all + n_greater_permute_all

        ntot_less_permute = len(n_less_permute_all)
        ntot_greater_permute = len(n_greater_permute_all)

        ntot_permute = ntot_less_permute + ntot_greater_permute


        likelihood_partition_permute = sum((n_less_permute_all*np.log((n_less_permute_all*ntot_permute)/(ntot_less_permute*n_permute_all) )) + (n_greater_permute_all*np.log((n_greater_permute_all*ntot_permute)/(ntot_greater_permute*n_permute_all) )))

        null_likelihood_partition.append(likelihood_partition_permute)

    null_likelihood_partition = np.asarray(null_likelihood_partition)

    #print((likelihood_partition - np.mean(null_likelihood_partition) ) / np.std(null_likelihood_partition) )

    #print(fmax_true_false)





def likelihood_subsample(taxon, treatment, ntot_subsample=50, fmax_cutoff=0.8, fmin_cutoff=0.0, subsamples=10000):
    # ntot_subsample minimum number of mutations

    # Load convergence matrix
    convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

    populations = [treatment+taxon + replicate for replicate in pt.replicates ]

    gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations, fmax_min=fmax_cutoff)

    G_subsample_list = []
    for i in range(subsamples):

        G_subsample = mutation_spectrum_utils.calculate_subsampled_total_parallelism(gene_parallelism_statistics, ntot_subsample=ntot_subsample)

        G_subsample_list.append(G_subsample)

    G_subsample_list.sort()

    G_CIs_dict = {}

    G_subsample_mean = np.mean(G_subsample_list)
    G_subsample_025 = G_subsample_list[ int( 0.025 * subsamples)  ]
    G_subsample_975 = G_subsample_list[ int( 0.975 * subsamples)  ]

    G_CIs_dict['G_mean'] = G_subsample_mean
    G_CIs_dict['G_025'] = G_subsample_025
    G_CIs_dict['G_975'] = G_subsample_975

    return G_CIs_dict



def calculate_likelihood_ratio_fmax(taxon, treatment, ntot_subsample=50, fmax_partition=0.8, subsamples=10000):

    convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

    populations = [treatment+taxon + replicate for replicate in pt.replicates ]

    gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations, fmax_min=fmax_cutoff)

    G_subsample_list = []




#fmax_cutoffs = np.asarray([0,0.2,0.4,0.6])
#fmax_cutoffs = np.asarray([0,0.1,0.2,0.3,0.4,0.5])
fmax_cutoffs = np.asarray([0.05, 0.1, 0.15,0.2,0.25,0.3])

G_dict_all = {}
taxa=['B','C','D','F','J','P']
treatments = ['0','1']
#for taxon in taxa:

#    sys.stdout.write("Sub-sampling taxon: %s\n" % (taxon))

#    for treatment in treatments:

#        if treatment+taxon in pt.treatment_taxa_to_ignore:
#            continue

#        for fmax_cutoff in fmax_cutoffs:
#
#            calculate_parallelism_statistics_partition(taxon, treatment, fmax_partition=fmax_cutoff)



ntotal_dict = {}
for taxon in taxa:

    sys.stdout.write("Sub-sampling taxon: %s\n" % (taxon))

    G_dict_all[taxon] = {}
    if taxon == 'J':
        ntotal = 50
    else:
        # calculate ntot for all frequency cutoffs
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % ('1'+taxon)))
        populations = ['1'+taxon + replicate for replicate in pt.replicates ]
        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,fmax_min=max(fmax_cutoffs))
        ntotal = 0
        for gene_i, gene_parallelism_statistics_i in gene_parallelism_statistics.items():
            ntotal+=gene_parallelism_statistics_i['observed']
    ntotal_dict[taxon] = ntotal
    for treatment in treatments:
        if treatment+taxon in pt.treatment_taxa_to_ignore:
            continue

        G_dict_all[taxon][treatment] = {}

        for fmax_cutoff in fmax_cutoffs:

            fmax_cutoff_dict = likelihood_subsample(taxon, treatment, ntot_subsample=ntotal,fmax_cutoff=fmax_cutoff, subsamples=subsamples)

            G_dict_all[taxon][treatment][fmax_cutoff] = fmax_cutoff_dict



fig = plt.figure(figsize = (9, 12))
gs = gridspec.GridSpec(nrows=4, ncols=3)
ax_count=0
for taxon_list_idx, taxon_list in enumerate([['B','C','J'],['D','F','P']]):
    for taxon_idx, taxon in enumerate(taxon_list):
        if taxon == '':
            continue
        ax = fig.add_subplot(gs[taxon_list_idx, taxon_idx])
        #ax.set_xlim([-0.05,max(fmax_cutoffs)+0.05])
        ax.set_xlim([-0.03,max(fmax_cutoffs)+0.05])
        ax.text(-0.2, 1.07, pt.sub_plot_labels[ax_count], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(pt.latex_genus_bold_dict[taxon] + ' ('+ r'$n_{total}=$' + str(ntotal_dict[taxon]) + ')' , fontsize=12)

        ax_count+=1

        for treatment in treatments:
            if treatment+taxon in pt.treatment_taxa_to_ignore:
                continue

            delta_l_list = []
            delta_025 = []
            delta_975 = []

            for fmax_cutoff in fmax_cutoffs:
                delta_l_list.append(G_dict_all[taxon][treatment][fmax_cutoff]['G_mean'])
                delta_025.append(G_dict_all[taxon][treatment][fmax_cutoff]['G_025'])
                delta_975.append(G_dict_all[taxon][treatment][fmax_cutoff]['G_975'])

            delta_l_list = np.asarray(delta_l_list)
            delta_025 = np.asarray(delta_025)
            delta_975 = np.asarray(delta_975)

            ax.errorbar(fmax_cutoffs, delta_l_list, yerr = [ delta_l_list-delta_025,  delta_975-delta_l_list] , \
                    fmt = 'o', alpha = 1, barsabove = True, marker = pt.plot_species_marker(taxon), \
                    mfc = 'white', mec = 'white', lw=2, c = 'k', zorder=1, ms=17)

            ax.scatter(fmax_cutoffs, delta_l_list, marker=pt.plot_species_marker(taxon), s = 150, \
                linewidth=3, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), alpha=1, zorder=2)


# now do divergence

significant_multiplicity_dict = {}

for taxon in pt.taxa:
    significant_multiplicity_dict[taxon] = {}
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
            if items[0] not in significant_multiplicity_dict[taxon]:
                significant_multiplicity_dict[taxon][items[0]] = {}
            significant_multiplicity_dict[taxon][items[0]][treatment] = float(items[-2])



fig.text(0.5, 0.485, "Maximum allele frequency (" + r'$f_{max}$' + ") cutoff", ha='center', va='center', fontsize=16)
fig.text(0.05, 0.7 ,"Net increase in log-likelihood, " r'$\Delta \ell$' , ha='center', va='center', rotation='vertical', fontsize=16)



treatment_pairs = [['0','1'],['0','2'],['1','2']]
#fig = plt.figure(figsize = (9, 6))
ax_divergence = fig.add_subplot(gs[2:4, 0:3])
ax_divergence.axhline( y=0, color='k', lw=2.5, linestyle='--', alpha = 1, zorder=1)
ax_count_divergence=0
slopes_dict = {}

record_strs = [",".join(['treatment_pair', 'taxon', 'tree_name', 'slope', 'slope_standard_error'])]



for treatment_pair_idx, treatment_pair in enumerate(treatment_pairs):

    treatment_pair_set = (treatment_pair[0], treatment_pair[1])

    slopes_dict[treatment_pair_set] = []

    # get mean colors
    ccv = ColorConverter()

    color_1 = np.array(ccv.to_rgb( pt.get_colors( treatment_pair[0] ) ))
    color_2 = np.array(ccv.to_rgb( pt.get_colors( treatment_pair[1] ) ))

    mix_color = 0.7 * (color_1 + color_2)
    mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)

    if (treatment_pair[0] == '0') and (treatment_pair[1] == '1'):
        #mix_color = pt.lighten_color(mix_color, amount=2.8)
        mix_color = 'gold'

    ax_count_divergence_treatment_pair = []
    treatment_pair_slopes = []

    for taxon in pt.taxa:

        result = [(x[treatment_pair[0]],x[treatment_pair[1]]) for x in significant_multiplicity_dict[taxon].values() if (treatment_pair[0] in x) and (treatment_pair[1] in x)]
        if len(result) ==0:
            continue

        mult_x = np.log10([x[0] for x in result])
        mult_y = np.log10([x[1] for x in result])

        mult_x = sm.add_constant(mult_x)
        mod = sm.OLS(mult_y, mult_x)
        res = mod.fit()
        slope = res.params[1]
        CI_025 = res.conf_int(0.05)[1][0]
        CI_975 = res.conf_int(0.05)[1][1]

        slope_std_error = res.bse[1]

        slopes_dict[treatment_pair_set].append(slope)

        def flip_slope_and_CIs(slope, CI_025, CI_975,null=1):
            new_slope = abs(slope-null)
            delta_CI_025 = abs(CI_025-slope)
            delta_CI_975 = abs(CI_975-slope)

            if slope < null:
                new_CI_025 = new_slope - delta_CI_975
                new_CI_975 = new_slope + delta_CI_025
            else:
                new_CI_025 = new_slope - delta_CI_025
                new_CI_975 = new_slope + delta_CI_975

            return new_slope,new_CI_025,new_CI_975


        new_slope,new_CI_025,new_CI_975 = flip_slope_and_CIs(slope, CI_025, CI_975)


        ax_divergence.errorbar(ax_count_divergence, new_slope, yerr = [ [new_slope-new_CI_025], [new_CI_975-new_slope]], \
                fmt = 'o', alpha = 1, barsabove = True, marker = pt.plot_species_marker(taxon), \
                mfc = 'white', mec = 'white', lw=3, c = 'k', zorder=2, ms=17)

        ax_divergence.scatter(ax_count_divergence, new_slope, marker=pt.plot_species_marker(taxon), s = 250, \
            linewidth=2, facecolors=mix_color, edgecolors='k', alpha=1, zorder=3)


        ax_count_divergence_treatment_pair.append(ax_count_divergence)
        treatment_pair_slopes.append(new_slope)

        ax_count_divergence+=1

        record_str = ",".join(['%s_%s' % treatment_pair_set,  str(taxon), pt.tree_name_dict[taxon], str(slope), str(slope_std_error)])
        record_strs.append(record_str)


    if treatment_pair_set == ('0','2'):
        xmin_mean = (min(ax_count_divergence_treatment_pair)+0.1)/16
        xmax_mean = (max(ax_count_divergence_treatment_pair)+1-0.1)/16

    else:
        xmin_mean = (min(ax_count_divergence_treatment_pair))/16
        xmax_mean = (max(ax_count_divergence_treatment_pair)+1)/16

    ax_divergence.axhline(y=np.mean(treatment_pair_slopes), xmin=xmin_mean, xmax=xmax_mean, color=mix_color, lw=3, linestyle='--', alpha = 1, zorder=1)

    if treatment_pair_idx < 2:

        ax_divergence.axvline( x=ax_count_divergence-0.5, color='k', lw=2, linestyle='-', alpha = 1, zorder=2)


#ax_divergence.set_xticks([], [])

#ax_divergence.set_xticklabels([2,7.5,13], ['1-day vs. 10-days', '1-day vs. 100-days', '10-days vs. 100-days'], fontsize=13)
ax_divergence.set_xticklabels( [], fontsize=12)

ax_divergence.text(0.15, -0.04, '1-day vs. 10-days', fontsize=15,  ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.5, -0.04, '1-day vs. 100-days', fontsize=15,  ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.84, -0.04, '10-days vs. 100-days', fontsize=15,  ha='center', va='center', transform=ax_divergence.transAxes)

#fig.text(0.2, 0.01, '1-day vs. 10-days', fontsize=14,  ha='center', va='center')
#fig.text(0.5, 0.01, '1-day vs. 100-days', fontsize=14,  ha='center', va='center')
#fig.text(0.6, 0.01, '10-day vs. 100-days', fontsize=14, ha='center', va='center')


ax_divergence.tick_params(axis='x', labelsize=14, length = 0)


ax_divergence.set_ylabel("Degree of divergent evolution, " + r'$\left | \beta_{1} - 1 \right |$' , fontsize = 16)

ax_divergence.text(-0.05, 1.07, pt.sub_plot_labels[ax_count], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_divergence.transAxes)


ax_divergence.text(0.79, 0.93, '$t_{ \mathrm{1\, vs \, 10, \, 1\, vs \, 100} } = %.3f , $' % (0.0367789), fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)
#ax_divergence.text(0.92, 0.93, r'$ P \nless 0.05 $', fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.925, 0.93, '$P = %.3f$' % (0.974), fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)


ax_divergence.text(0.80, 0.88, '$t_{ \mathrm{1\, vs \, 10, \, 10\, vs \, 100} } = %.2f , $' % (-0.334594), fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)
#ax_divergence.text(0.94, 0.88, r'$P \nless 0.05 $', fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.945, 0.88, '$P = %.3f$' % (0.769762), fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)


ax_divergence.text(0.80, 0.83, '$t_{ \mathrm{1\, vs \, 100, \, 10\, vs \, 100} } = %.2f , $' % (-0.37743), fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)
#ax_divergence.text(0.944, 0.83, r'$P \nless 0.05 $', fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)
ax_divergence.text(0.95, 0.83, '$P = %.3f$' % (0.742142), fontsize=9,  ha='center', va='center', transform=ax_divergence.transAxes)



fig.subplots_adjust(hspace=0.35,wspace=0.3) #hspace=0.3, wspace=0.5
fig_name = pt.get_path() + "/figs/G_score_vs_fmax_subsample_divergence.pdf"
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



sys.stderr.write("Done with figure!\n")

sys.stderr.write("Writing intermediate file for phylogenetic t-test...\n")
file = open(pt.get_path()+'/data/divergence_slopes.csv',"w")
record_str = "\n".join(record_strs)
file.write(record_str)
file.close()
sys.stderr.write("Done!\n")
