from __future__ import division
import os, sys
import numpy as np

import  matplotlib.pyplot as plt
from matplotlib.colors import ColorConverter
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

import scipy.stats as stats

import phylo_tools as pt
import parse_file
import timecourse_utils
import mutation_spectrum_utils


def plot_within_taxon_paralleliism(taxon, slope_null=1):

    fig = plt.figure(figsize = (12, 8))

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    # to get the common gene names for each ID

    ax_multiplicity = plt.subplot2grid((2, 3), (0, 0), colspan=1)
    ax_mult_freq = plt.subplot2grid((2, 3), (0, 1), colspan=1)
    ax_venn = plt.subplot2grid((2, 3), (0, 2), colspan=1)

    ax_multiplicity.set_xscale('log', base=10)
    ax_multiplicity.set_yscale('log', base=10)
    ax_multiplicity.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    ax_multiplicity.set_ylabel('Fraction mutations ' + r'$\geq m$', fontsize=14)
    ax_multiplicity.text(-0.1, 1.07, pt.sub_plot_labels[0], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)


    ax_multiplicity.set_ylim([0.001, 1.1])
    ax_multiplicity.set_xlim([0.07, 130])


    ax_mult_freq.set_xscale('log', base=10)
    ax_mult_freq.set_yscale('log', base=10)
    ax_mult_freq.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    ax_mult_freq.set_ylabel('Mean maximum allele frequency, ' + r'$\overline{f}_{max}$', fontsize=11)
    ax_mult_freq.text(-0.1, 1.07, pt.sub_plot_labels[1], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_mult_freq.transAxes)

    ax_venn.axis('off')
    ax_venn.text(-0.1, 1.07, pt.sub_plot_labels[2], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_venn.transAxes)


    alpha_treatment_dict = {'0':0.5, '1':0.5, '2':0.8}


    significant_multiplicity_dict = {}

    significant_multiplicity_values_dict = {}

    multiplicity_dict = {}

    g_score_p_label_dict = {}

    all_mults = []
    all_freqs = []

    treatments_in_taxon = []

    label_y_axes = [0.3,0.2,0.1]

    for treatment_idx, treatment in enumerate(pt.treatments):


        significan_multiplicity_taxon_path = pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
        if os.path.exists(significan_multiplicity_taxon_path) == False:
            continue
        treatments_in_taxon.append(treatment)
        significan_multiplicity_taxon = open(significan_multiplicity_taxon_path, "r")

        significan_multiplicity_list = []
        for i, line in enumerate( significan_multiplicity_taxon ):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            significan_multiplicity_list.append(items[0])

            if items[0] not in significant_multiplicity_values_dict:
                significant_multiplicity_values_dict[items[0]] = {}
                significant_multiplicity_values_dict[items[0]][treatment] = float(items[-2])
            else:
                significant_multiplicity_values_dict[items[0]][treatment] = float(items[-2])


        significant_multiplicity_dict[treatment] = significan_multiplicity_list

        populations = [treatment+taxon + replicate for replicate in pt.replicates ]

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))
        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)
        #print(gene_parallelism_statistics)
        G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

        sys.stdout.write("Total parallelism for %s = %g (p=%g)\n" % (treatment+taxon, G,pvalue))

        predictors = []
        responses = []

        gene_hits = []
        gene_predictors = []
        mean_gene_freqs = []

        Ls = []

        ax_mult_freqs_x = []
        ax_mult_freqs_y = []

        for gene_name in convergence_matrix.keys():

            convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

            Ls.append(convergence_matrix[gene_name]['length'])
            m = gene_parallelism_statistics[gene_name]['multiplicity']


            if gene_name not in multiplicity_dict:
                multiplicity_dict[gene_name] = {}
                multiplicity_dict[gene_name][treatment] = m
            else:
                multiplicity_dict[gene_name][treatment] = m

            n = 0
            nfixed = 0
            freqs = []
            nf_max = 0

            for population in populations:
                for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:
                    fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                    predictors.append(m)
                    responses.append(fixed_weight)

                    n+=1
                    nfixed+=fixed_weight

                    # get freqs for regression
                    #if L == parse_file.POLYMORPHIC:
                    #freqs.append(f_max)
                    nf_max+=timecourse_utils.calculate_fixed_weight(L,f_max)

            if n > 0.5:
                gene_hits.append(n)
                gene_predictors.append(m)
                #mean_gene_freqs.append(np.mean(freqs))

                if nf_max > 0:
                    ax_mult_freqs_x.append(m)
                    ax_mult_freqs_y.append( nf_max / n )

        Ls = np.asarray(Ls)
        ntot = len(predictors)
        mavg = ntot*1.0/len(Ls)

        predictors, responses = (np.array(x) for x in zip(*sorted(zip(predictors, responses), key=lambda pair: (pair[0]))))

        gene_hits, gene_predictors = (np.array(x) for x in zip(*sorted(zip(gene_hits, gene_predictors), key=lambda pair: (pair[0]))))

        rescaled_predictors = np.exp(np.fabs(np.log(predictors/mavg)))

        null_survival_function = mutation_spectrum_utils.NullMultiplicitySurvivalFunction.from_parallelism_statistics(gene_parallelism_statistics)

        # default base is 10
        theory_ms = np.logspace(-2,2,100)
        theory_survivals = null_survival_function(theory_ms)
        theory_survivals /= theory_survivals[0]

        sys.stderr.write("Done!\n")

        ax_multiplicity.plot(theory_ms, theory_survivals, lw=3, color=pt.get_colors(treatment), alpha=0.8, ls=':',  zorder=1)

        ax_multiplicity.plot(predictors, (len(predictors)-np.arange(0,len(predictors)))*1.0/len(predictors), lw=3, color=pt.get_colors(treatment),alpha=0.8, ls='--', label= str(int(10**int(treatment))) + '-day', drawstyle='steps', zorder=2)

        #ax_multiplicity.text(0.2, 0.3, g_score_p_label_dict['0'], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)
        #ax_multiplicity.text(0.2, 0.2, g_score_p_label_dict['1'], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)
        #ax_multiplicity.text(0.2, 0.1, g_score_p_label_dict['2'], fontsize=25, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)

        if pvalue < 0.001:
            pretty_pvalue = r'$\ll 0.001$'
        else:
            pretty_pvalue = '=' + str(round( pvalue, 4 ))

        g_score_p_label = r'$\Delta \ell_{{{}}}=$'.format(str(10**int(treatment))) + str(round(G, 3)) + ', ' + r'$P$' + pretty_pvalue

        text_color = pt.lighten_color(pt.get_colors(treatment), amount=1.3)

        ax_multiplicity.text(0.26, label_y_axes[treatment_idx], g_score_p_label, fontsize=7, ha='center', va='center',  color='k', transform=ax_multiplicity.transAxes)

        ax_mult_freq.scatter(ax_mult_freqs_x, ax_mult_freqs_y, color=pt.get_colors(treatment), edgecolors=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), alpha=alpha_treatment_dict[treatment])

        all_mults.extend(ax_mult_freqs_x)
        all_freqs.extend(ax_mult_freqs_y)

        #slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(ax_mult_freqs_x), np.log10(ax_mult_freqs_y))
        #print(slope, p_value)

    # make treatment pairs
    treatments_in_taxon.sort(key=float)

    for i in range(0, len(treatments_in_taxon)):

        for j in range(i+1, len(treatments_in_taxon)):

            ax_mult_i_j = plt.subplot2grid((2, 3), (1, i+j-1), colspan=1)
            ax_mult_i_j.set_xscale('log', base=10)
            ax_mult_i_j.set_yscale('log', base=10)
            ax_mult_i_j.set_xlabel( str( 10** int( treatments_in_taxon[i] ) ) + '-day gene multiplicity, ' + r'$m$', fontsize=14)
            ax_mult_i_j.set_ylabel(str( 10** int( treatments_in_taxon[j] ) ) + '-day gene multiplicity, ' + r'$m$', fontsize=14)
            ax_mult_i_j.plot([   0.05, 200  ], [   0.05, 200   ], lw = 3, c='grey', ls = '--', zorder=1 )
            ax_mult_i_j.set_xlim([   0.05, 200  ])
            ax_mult_i_j.set_ylim([   0.05, 200  ])

            ax_mult_i_j.text(-0.1, 1.07, pt.sub_plot_labels[2+i+j], fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_mult_i_j.transAxes)

            multiplicity_pair = [(multiplicity_dict[gene_name][treatments_in_taxon[i]], multiplicity_dict[gene_name][treatments_in_taxon[j]]) for gene_name in sorted(multiplicity_dict) if (multiplicity_dict[gene_name][treatments_in_taxon[i]] > 0) and (multiplicity_dict[gene_name][treatments_in_taxon[j]] > 0) ]
            significant_multiplicity_pair = [(significant_multiplicity_values_dict[gene_name][treatments_in_taxon[i]], significant_multiplicity_values_dict[gene_name][treatments_in_taxon[j]]) for gene_name in sorted(significant_multiplicity_values_dict) if (treatments_in_taxon[i] in significant_multiplicity_values_dict[gene_name]) and (treatments_in_taxon[j] in significant_multiplicity_values_dict[gene_name]) ]

            # get mean colors
            ccv = ColorConverter()

            color_1 = np.array(ccv.to_rgb( pt.get_colors( treatments_in_taxon[i] ) ))
            color_2 = np.array(ccv.to_rgb( pt.get_colors( treatments_in_taxon[j] ) ))

            mix_color = 0.7 * (color_1 + color_2)
            mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)

            if (treatments_in_taxon[i] == '0') and (treatments_in_taxon[j] == '1'):
                #mix_color = pt.lighten_color(mix_color, amount=2.8)
                mix_color = 'gold'

            mult_i = [x[0] for x in multiplicity_pair]
            mult_j = [x[1] for x in multiplicity_pair]

            ax_mult_i_j.scatter(mult_i, mult_j, marker=pt.plot_species_marker(taxon), facecolors=mix_color, edgecolors='none', alpha=0.8,s=90, zorder=2)

            mult_significant_i = [x[0] for x in significant_multiplicity_pair]
            mult_significant_j = [x[1] for x in significant_multiplicity_pair]
            ax_mult_i_j.scatter(mult_significant_i, mult_significant_j, marker=pt.plot_species_marker(taxon), facecolors=mix_color, edgecolors='k', lw=1.5, alpha=0.7,s=90, zorder=3)

            #slope_mult, intercept_mult, r_value_mult, p_value_mult, std_err_mult = stats.linregress(np.log10(mult_significant_i), np.log10(mult_significant_j))

            mult_ij = mult_significant_i+mult_significant_j + mult_i + mult_j

            ax_mult_i_j.set_xlim([min(mult_ij)*0.5, max(mult_ij)*1.5])
            ax_mult_i_j.set_ylim([min(mult_ij)*0.5, max(mult_ij)*1.5])

            # null slope of 1
            #ratio = (slope_mult - slope_null) / std_err_mult
            #p_value_mult_new_null = stats.t.sf(np.abs(ratio), len(mult_significant_j)-2)*2

            #if p_value_mult_new_null < 0.05:
            #    x_log10_fit_range =  np.linspace(np.log10(min(mult_i) * 0.5), np.log10(max(mult_i) * 1.5), 10000)

            #    y_fit_range = 10 ** (slope_mult*x_log10_fit_range + intercept_mult)
            #    ax_mult_i_j.plot(10**x_log10_fit_range, y_fit_range, c='k', lw=3, linestyle='--', zorder=4)

            #ax_mult_i_j.text(0.05, 0.9, r'$\beta_{1}=$'+str(round(slope_mult,3)), fontsize=12, transform=ax_mult_i_j.transAxes)
            #ax_mult_i_j.text(0.05, 0.82, r'$r^{2}=$'+str(round(r_value_mult**2,3)), fontsize=12, transform=ax_mult_i_j.transAxes)
            #ax_mult_i_j.text(0.05, 0.74, pt.get_p_value_latex(p_value_mult_new_null), fontsize=12, transform=ax_mult_i_j.transAxes)

    #if taxon == 'F':
    #    subset_tuple = (len( significant_multiplicity_dict['0']), \
    #                    len( significant_multiplicity_dict['1']), \
    #                    len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['1'])))

    #    venn = venn2(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('1')))
    #    c = venn2_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')

    subset_tuple = (len( significant_multiplicity_dict['0']), \
                    len( significant_multiplicity_dict['1']), \
                    len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['1'])), \
                    len(significant_multiplicity_dict['2']), \
                    len(set(significant_multiplicity_dict['0']) & set(significant_multiplicity_dict['2'])), \
                    len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])),  \
                    len(set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['1']) & set(significant_multiplicity_dict['2'])))

    venn = venn3(subsets = subset_tuple, ax=ax_venn, set_labels=('', '', ''), set_colors=(pt.get_colors('0'), pt.get_colors('1'), pt.get_colors('2')))
    c = venn3_circles(subsets=subset_tuple, ax=ax_venn, linestyle='dashed')

    ax_mult_freq.set_xlim([min(all_mults)*0.5, max(all_mults)*1.5])
    ax_mult_freq.set_ylim([min(all_freqs)*0.5, max(all_freqs)*1.5])

    fig.suptitle(pt.latex_dict[taxon], fontsize=30)

    fig.subplots_adjust(wspace=0.3) #hspace=0.3, wspace=0.5
    fig_name = pt.get_path() + "/figs/multiplicity_%s.jpg" % taxon
    fig.savefig(fig_name, format='jpg', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


for taxon in ['B','S']:

    plot_within_taxon_paralleliism(taxon)
