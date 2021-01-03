import os, sys
import parse_file
import phylo_tools as pt
import timecourse_utils
import numpy

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

import statsmodels.stats.multitest as multitest
import scipy.stats as stats

treatments=pt.treatments
replicates = pt.replicates
taxa = ['B','C','D','F','J','P']

sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']

nonsynonymous_types = set(['missense','nonsense'])
synonymous_types = set(['synonymous'])

non_appeared = {}
non_fixed = {}

syn_appeared = {}
syn_fixed = {}

targeted_Lsyn = {}
targeted_Lnon = {}
targeted_fixed_Lsyn = {}
targeted_fixed_Lnon = {}

taxon_Lsyn_dict = {}
taxon_Lnon_dict = {}

populations = []

for taxon in taxa:
    Lsyn, Lnon, substitution_specific_synonymous_fraction = parse_file.calculate_synonymous_nonsynonymous_target_sizes(taxon)
    taxon_Lsyn_dict[taxon] = Lsyn
    taxon_Lnon_dict[taxon] = Lnon
    for treatment in treatments:
        for replicate in replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            populations.append(population)

            non_appeared[population] = 1
            non_fixed[population] = 1

            syn_appeared[population] = 1
            syn_fixed[population] = 1

            targeted_Lsyn[population] = 1
            targeted_Lnon[population] = 1
            targeted_fixed_Lsyn[population] = 1
            targeted_fixed_Lnon[population] = 1

            mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
            population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
            state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

            num_processed_mutations = 0

            for mutation_idx in range(0,len(mutations)):

                #location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
                location, gene_name, allele, var_type, codon, position_in_codon, AAs_count,  test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]

                state_Ls = state_trajectories[mutation_idx]

                good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

                freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

                masked_times = times[good_idxs]
                masked_freqs = freqs[good_idxs]
                masked_state_Ls = state_Ls[good_idxs]

                fixed_weight = timecourse_utils.calculate_fixed_weight(masked_state_Ls[-1],masked_freqs[-1])

                if var_type in nonsynonymous_types or var_type in synonymous_types:
                    targeted_Lnon[population] += (1-substitution_specific_synonymous_fraction[allele])
                    targeted_fixed_Lnon[population] += fixed_weight*(1-substitution_specific_synonymous_fraction[allele])
                    targeted_Lsyn[population] += substitution_specific_synonymous_fraction[allele]
                    targeted_fixed_Lsyn[population] += fixed_weight*substitution_specific_synonymous_fraction[allele]

                if var_type in nonsynonymous_types:
                    non_appeared[population]+=1
                    non_fixed[population]+=fixed_weight
                    num_processed_mutations+=1

                elif var_type in synonymous_types:
                    syn_appeared[population]+=1
                    syn_fixed[population]+=fixed_weight
                    num_processed_mutations+=1



#total_non_appeared = sum([non_appeared[population] for population in populations])
#total_non_fixed = sum([non_fixed[population] for population in populations])
#total_syn_appeared = sum([syn_appeared[population] for population in populations])
#total_syn_fixed = sum([syn_fixed[population] for population in populations])

#dnds_appeared = total_non_appeared/total_syn_appeared*Lsyn/Lnon
#dnds_fixed = total_non_fixed/total_syn_fixed*Lsyn/Lnon

#total_targeted_Lnon = sum(targeted_Lnon.values())
#total_targeted_Lsyn = sum(targeted_Lsyn.values())
#targeted_dnds_appeared = total_non_appeared/total_syn_appeared*total_targeted_Lsyn/total_targeted_Lnon
#total_targeted_Lnon = sum(targeted_fixed_Lnon.values())
#total_targeted_Lsyn = sum(targeted_fixed_Lsyn.values())
#targeted_dnds_fixed = total_non_fixed/total_syn_fixed*total_targeted_Lsyn/total_targeted_Lnon


#population_dnds_appeared = numpy.array([non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*Lsyn/Lnon for population in populations])
#population_dnds_fixed = numpy.array([non_fixed[population]/(syn_fixed[population]+(syn_fixed[population]==0))*Lsyn/Lnon for population in populations])
#targeted_population_dnds_appeared = numpy.array([non_appeared[population]/syn_appeared[population]*targeted_Lsyn[population]/targeted_Lnon[population] for population in populations])
#targeted_population_dnds_fixed = numpy.array([non_fixed[population]/(syn_fixed[population]+(syn_fixed[population]==0))*targeted_fixed_Lsyn[population]/targeted_fixed_Lnon[population] for population in populations])


#for population in populations:
#    sys.stdout.write("%s: %d non, %d syn, dN/dS = %g, (dN/dS)* = %g\n" % (population, non_appeared[population], syn_appeared[population], non_appeared[population]/syn_appeared[population]*Lsyn/Lnon, non_appeared[population]/syn_appeared[population]*targeted_Lsyn[population]/targeted_Lnon[population]))
#    #sys.stdout.write("%s fixed: %d non, %d syn\n" % (population, non_fixed[population], syn_fixed[population]))





#fig, ax = plt.subplots(figsize=(4,6))
fig = plt.figure(figsize = (10, 5))
gs = gridspec.GridSpec(nrows=2, ncols=3)
#fig = plt.figure(figsize = (12, 6))
anova_pvalues = []
anova_F = []
for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):
    for taxon_idx, taxon in enumerate(taxon_list):
        dnds_samples = []
        for treatment in treatments:

            populations_plot = [ treatment+taxon+replicate for replicate in replicates if treatment+taxon+replicate not in pt.populations_to_ignore ]
            taxon_treatment_dnds_appeared = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*taxon_Lsyn_dict[taxon]/taxon_Lnon_dict[taxon] for population in populations_plot]
            if len(taxon_treatment_dnds_appeared) < 2:
                continue
            dnds_samples.append(taxon_treatment_dnds_appeared)

        if taxon == 'J':
            fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1])

        else:
            fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1], dnds_samples[2])

        sys.stdout.write("%s: dN/dS one-way ANOVA: F = %g, P = %g\n" % (taxon, fvalue, pvalue))

        anova_pvalues.append(pvalue)
        anova_F.append(fvalue)


reject, pvals_corrected, alphacSidak, alphacBonf = multitest.multipletests(anova_pvalues, alpha=0.05, method='fdr_bh')


ax_count = 0

for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):
    for taxon_idx, taxon in enumerate(taxon_list):
        ax = fig.add_subplot(gs[taxon_list_idx, taxon_idx])
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold')

        dnds_samples = []
        for treatment in treatments:

            populations_plot = [ treatment+taxon+replicate for replicate in replicates if treatment+taxon+replicate not in pt.populations_to_ignore ]
            taxon_treatment_dnds_appeared = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*taxon_Lsyn_dict[taxon]/taxon_Lnon_dict[taxon] for population in populations_plot]
            if len(taxon_treatment_dnds_appeared) < 2:
                continue
            ax.scatter( [int(treatment)] * len(taxon_treatment_dnds_appeared), taxon_treatment_dnds_appeared,  marker=pt.plot_species_marker(taxon),  linewidth=2, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), s=100, zorder=2, alpha=0.8)
            if len(taxon_treatment_dnds_appeared) > 2:
                ax.errorbar(int(treatment),numpy.mean(taxon_treatment_dnds_appeared), yerr= 2*numpy.std(taxon_treatment_dnds_appeared) / numpy.sqrt(len(taxon_treatment_dnds_appeared)), linestyle='-', c = 'k', marker=pt.plot_species_marker(taxon), lw = 2.5,  zorder=3)
            #dnds_treatment.append(taxon_treatment_dnds_appeared)

            dnds_samples.append(taxon_treatment_dnds_appeared)

        ax.text(-0.1, 1.07, sub_plot_labels[ax_count], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)
        ax.text(0.7, 0.9, r'$F=$'+ str(round( anova_F[ax_count],3) ), fontsize=10, ha='center', va='center', transform=ax.transAxes)
        ax.text(0.7, 0.8, r'$P_{BH}=$'+ str(round(pvals_corrected[ax_count], 3)) , fontsize=10, ha='center', va='center', transform=ax.transAxes)

        ax_count+=1

        if taxon == 'J':
            ax.set_xticks([0,2])
            ax.set_xticklabels( ['1','100'] )
            ax.set_xlim([-0.3, 2.3])

            #fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1])

        else:
            ax.set_xticks([0,1,2])
            ax.set_xticklabels( ['1','10','100'] )
            ax.set_xlim([-0.3, 2.3])

            #fvalue, pvalue = stats.f_oneway(dnds_samples[0], dnds_samples[1], dnds_samples[2])



legend_elements = [Line2D([0], [0], color = 'none', marker=pt.plot_species_marker('B'), label=pt.latex_genus_dict['B'],
                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], color = 'none', marker=pt.plot_species_marker('C'), label=pt.latex_genus_dict['C'],
                                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], color = 'none', marker=pt.plot_species_marker('D'), label=pt.latex_genus_dict['D'],
                                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], color = 'none', marker=pt.plot_species_marker('F'), label=pt.latex_genus_dict['F'],
                                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], color = 'none', marker=pt.plot_species_marker('J'), label=pt.latex_genus_dict['J'],
                                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], color = 'none', marker=pt.plot_species_marker('P'), label=pt.latex_genus_dict['P'],
                                    markerfacecolor='k', markersize=13)]

# Create the figure
#ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

fig.text(0.5, -0.01, "Transfer time", ha='center', va='center', fontsize=18)
fig.text(0.05, 0.5, 'dN/dS', ha='center', va='center', rotation='vertical', fontsize=18)


fig.subplots_adjust(hspace=0.3, wspace=0.5)
fig_name = pt.get_path() + '/figs/dn_ds.pdf'
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
