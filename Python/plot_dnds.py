import os, sys
import parse_file
import phylo_tools as pt
import timecourse_utils
import numpy

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import scipy.stats as stats



Lsyn, Lnon, substitution_specific_synonymous_fraction = parse_file.calculate_synonymous_nonsynonymous_target_sizes('B')

treatments=pt.treatments
replicates = pt.replicates
taxa = ['B', 'S']

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

populations = []

for taxon in taxa:
    for treatment in treatments:

        nonsynonymous_fmax_all = []
        synonymous_fmax_all = []

        for replicate in replicates:

            population = treatment + taxon + replicate
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

                location, gene_name, allele, var_type, codon, position_in_codon, AAs_count, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
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

                    nonsynonymous_fmax_all.append(max(freqs))

                elif var_type in synonymous_types:
                    syn_appeared[population]+=1
                    syn_fixed[population]+=fixed_weight
                    num_processed_mutations+=1

                    synonymous_fmax_all.append(max(freqs))


total_non_appeared = sum([non_appeared[population] for population in populations])
total_non_fixed = sum([non_fixed[population] for population in populations])
total_syn_appeared = sum([syn_appeared[population] for population in populations])
total_syn_fixed = sum([syn_fixed[population] for population in populations])

dnds_appeared = total_non_appeared/total_syn_appeared*Lsyn/Lnon
dnds_fixed = total_non_fixed/total_syn_fixed*Lsyn/Lnon

total_targeted_Lnon = sum(targeted_Lnon.values())
total_targeted_Lsyn = sum(targeted_Lsyn.values())
targeted_dnds_appeared = total_non_appeared/total_syn_appeared*total_targeted_Lsyn/total_targeted_Lnon
total_targeted_Lnon = sum(targeted_fixed_Lnon.values())
total_targeted_Lsyn = sum(targeted_fixed_Lsyn.values())
targeted_dnds_fixed = total_non_fixed/total_syn_fixed*total_targeted_Lsyn/total_targeted_Lnon


population_dnds_appeared = numpy.array([non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*Lsyn/Lnon for population in populations])

population_dnds_fixed = numpy.array([non_fixed[population]/(syn_fixed[population]+(syn_fixed[population]==0))*Lsyn/Lnon for population in populations])

targeted_population_dnds_appeared = numpy.array([non_appeared[population]/syn_appeared[population]*targeted_Lsyn[population]/targeted_Lnon[population] for population in populations])

targeted_population_dnds_fixed = numpy.array([non_fixed[population]/(syn_fixed[population]+(syn_fixed[population]==0))*targeted_fixed_Lsyn[population]/targeted_fixed_Lnon[population] for population in populations])


#for population in populations:
#    sys.stdout.write("%s: %d non, %d syn, dN/dS = %g, (dN/dS)* = %g\n" % (population, non_appeared[population], syn_appeared[population], non_appeared[population]/syn_appeared[population]*Lsyn/Lnon, non_appeared[population]/syn_appeared[population]*targeted_Lsyn[population]/targeted_Lnon[population]))



jitter_shift = [-0.1, 0.1]

fig, ax = plt.subplots(figsize=(4,4))
#fig = plt.figure(figsize = (12, 6))
for treatment in ['0', '1', '2']:
    #ax_i = plt.subplot2grid((1, 2), (0,taxon_idx), colspan=1)
    #ax_i.set_title( latex_dict[taxon], fontsize=17)
    #ax_i.text(-0.1, 1.07, sub_plot_labels[treatment_idx], fontsize=22, fontweight='bold', ha='center', va='center', transform=ax_i.transAxes)
    dnds_treatment = []
    for taxon_idx, taxon in enumerate(taxa):

        populations_plot = [ treatment+taxon + replicate for replicate in replicates ]

        taxon_treatment_dnds_appeared = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*Lsyn/Lnon for population in populations_plot]

        ax.scatter( [int(treatment) + jitter_shift[taxon_idx] ] * len(taxon_treatment_dnds_appeared), taxon_treatment_dnds_appeared , marker=pt.plot_species_marker(taxon),  linewidth=2, facecolors=pt.get_scatter_facecolor(taxon, treatment), edgecolors=pt.get_colors(treatment), s=120, zorder=2)
        ax.errorbar(int(treatment)+ jitter_shift[taxon_idx], numpy.mean(taxon_treatment_dnds_appeared), yerr= 2*numpy.std(taxon_treatment_dnds_appeared) / numpy.sqrt(len(taxon_treatment_dnds_appeared)), linestyle='-', c = 'k', marker=pt.plot_species_marker(taxon), lw = 2.5)
        dnds_treatment.append(taxon_treatment_dnds_appeared)

    t, p = stats.ttest_ind(dnds_treatment[0], dnds_treatment[1], equal_var=False)
    #print(t, p)


# to t test between treatments for each strain

for treatment in ['0', '1', '2']:

    populations_B = [ treatment+'B' + replicate for replicate in replicates ]
    populations_S = [ treatment+'S' + replicate for replicate in replicates ]

    taxon_treatment_dnds_appeared_B = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*Lsyn/Lnon for population in populations_B]
    taxon_treatment_dnds_appeared_S = [non_appeared[population]/(syn_appeared[population]+(syn_appeared[population]==0))*Lsyn/Lnon for population in populations_S]

    t, p = stats.ttest_ind(taxon_treatment_dnds_appeared_B, taxon_treatment_dnds_appeared_S, equal_var=True)
    print(treatment, t, p)






#ax.set_xlim([-0.25, 1.25])
ax.set_ylim([0.35, 1.1])
ax.axhline(y=1, color='k', linestyle=':', lw=2.5, alpha = 0.8, zorder=1)
#ax_i.set_xscale('log', basex=10)
#ax.set_xlabel("Transfer time (days)", fontsize = 18)
ax.set_ylabel("pN/pS", fontsize = 18)

#x_ticks = [0,1]
ax.set_xticks([0, 1, 2])
ax.set_xticklabels( ['1-day', '10-days', '100-days'], fontweight='bold', fontsize=12)



legend_elements = [Line2D([0], [0], color = 'none', marker='o', label=pt.latex_dict['B'],
                    markerfacecolor='k', markersize=13),
                Line2D([0], [0], marker='o', color='none', label=pt.latex_dict['S'],
                    markerfacecolor='w', markersize=13, markeredgewidth=2)]
# Create the figure
ax.legend(handles=legend_elements, loc='upper right')


fig.subplots_adjust(hspace=0.3, wspace=0.5)
fig_name = pt.get_path() + '/figs/plot_dn_ds.jpg'
fig.savefig(fig_name, format='jpg', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
