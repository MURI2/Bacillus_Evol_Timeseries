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


taxa = ['B', 'S']

parallelism_axes = {}


fig = plt.figure(figsize = (12, 12))

gene_data = parse_file.parse_gene_list('B')

gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
# to get the common gene names for each ID

for treatment_idx, treatment in enumerate(pt.treatments):

    ax_multiplicity = plt.subplot2grid((3, 3), (0, treatment_idx), colspan=1)
    ax_regression = plt.subplot2grid((3, 3), (1, treatment_idx), colspan=1)
    ax_venn = plt.subplot2grid((3, 3), (2, treatment_idx), colspan=1)

    ax_multiplicity.set_xlabel('Gene multiplicity, ' + r'$m$', fontsize=14)
    ax_multiplicity.set_ylabel('Fraction mutations ' + r'$\geq m$', fontsize=14)
    ax_multiplicity.set_title( str(10**int(treatment))+ '-day transfers', fontsize=17)

    ax_multiplicity.set_xscale('log', base=10)
    ax_multiplicity.set_yscale('log', base=10)

    ax_multiplicity.set_ylim([0.001, 1.1])
    ax_multiplicity.set_xlim([0.07, 130])
    ax_multiplicity.text(-0.1, 1.07, pt.sub_plot_labels[treatment_idx], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_multiplicity.transAxes)


    ax_regression.set_xlabel('Gene multiplicity, ' + r'$m$' + '\n' + pt.latex_dict['B']   , fontsize=14)
    ax_regression.set_ylabel('Gene multiplicity, ' + r'$m$' + '\n' + pt.latex_dict['S'] , fontsize=14)

    ax_regression.set_xscale('log', base=10)
    ax_regression.set_yscale('log', base=10)
    ax_regression.text(-0.1, 1.07, pt.sub_plot_labels[3+treatment_idx], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_regression.transAxes)


    ax_venn.axis('off')
    ax_venn.text(-0.1, 1.07, pt.sub_plot_labels[6+treatment_idx], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_venn.transAxes)


    mult_taxa_dict = {}

    for taxon in taxa:

        if taxon == 'B':
            multiplicity_label = r'$\mathrm{wt}$'
        else:
            multiplicity_label = r'$\Delta \mathrm{spo0A}$'

        populations = [treatment+taxon + replicate for replicate in pt.replicates ]

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix,populations,Lmin=100)

        G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

        sys.stdout.write("Total parallelism for %s = %g (p=%g)\n" % (treatment+taxon, G,pvalue))

        predictors = []
        responses = []

        gene_hits = []
        gene_predictors = []

        Ls = []

        for gene_name in convergence_matrix.keys():

            # get multiplicities for regression
            if gene_name not in mult_taxa_dict:
                mult_taxa_dict[gene_name] = {}
                mult_taxa_dict[gene_name][taxon] = gene_parallelism_statistics[gene_name]['multiplicity']
            else:
                mult_taxa_dict[gene_name][taxon] = gene_parallelism_statistics[gene_name]['multiplicity']

            convergence_matrix[gene_name]['length'] < 50 and convergence_matrix[gene_name]['length']

            Ls.append(convergence_matrix[gene_name]['length'])
            m = gene_parallelism_statistics[gene_name]['multiplicity']

            n = 0
            nfixed = 0

            for population in populations:
                for t,L,f,f_max in convergence_matrix[gene_name]['mutations'][population]:
                    fixed_weight = timecourse_utils.calculate_fixed_weight(L,f)

                    predictors.append(m)
                    responses.append(fixed_weight)

                    n+=1
                    nfixed+=fixed_weight

            if n > 0.5:
                gene_hits.append(n)
                gene_predictors.append(m)

        Ls = np.asarray(Ls)
        ntot = len(predictors)
        mavg = ntot*1.0/len(Ls)

        predictors, responses = (np.array(x) for x in zip(*sorted(zip(predictors, responses), key=lambda pair: (pair[0]))))

        gene_hits, gene_predictors = (np.array(x) for x in zip(*sorted(zip(gene_hits, gene_predictors), key=lambda pair: (pair[0]))))

        #rescaled_predictors = np.exp(np.fabs(np.log(predictors/mavg)))

        #logit_mod = sm.Logit(responses, sm.add_constant(rescaled_predictors))
        #logit_res = logit_mod.fit()

        #sys.stdout.write("Logistic regression for %ss:\n" % (population))
        #sys.stdout.write("%s\n" % str(logit_res.summary()))
        #sys.stdout.write("Params:\n")
        #sys.stdout.write("%s\n" % str(logit_res.params))

        sys.stdout.write("Avg mutation multiplicity=%g, Avg fixed mutation multiplicity=%g\n" % (predictors.sum()/len(responses), (predictors*responses).sum()/responses.sum()))
        sys.stderr.write("Calculating null distribution...\n")
        null_survival_function = mutation_spectrum_utils.NullMultiplicitySurvivalFunction.from_parallelism_statistics(gene_parallelism_statistics)

        # default base is 10
        theory_ms = np.logspace(-2,2,100)
        theory_survivals = null_survival_function(theory_ms)
        theory_survivals /= theory_survivals[0]

        sys.stderr.write("Done!\n")

        # step function
        ax_multiplicity.plot(predictors, (len(predictors)-np.arange(0,len(predictors)))*1.0/len(predictors), lw=3, color=pt.get_colors(treatment),alpha=0.8, ls=pt.get_taxon_ls(taxon), label='Observed ' + multiplicity_label, drawstyle='steps', zorder=2)

        ax_multiplicity.plot(theory_ms, theory_survivals, lw=3, color='grey',alpha=0.8, ls=pt.get_taxon_ls(taxon), label= 'Null ' +  multiplicity_label, zorder=1)

    if treatment_idx == 0:
        ax_multiplicity.legend( loc='lower left', fontsize=8)

    for gene_name, gene_dict in mult_taxa_dict.items():
        if 'B' not in gene_dict:
            mult_taxa_dict[gene_name]['B'] = 0
        if 'S' not in gene_dict:
            mult_taxa_dict[gene_name]['S'] = 0

    # then get venn diagram
    # import significant genes
    parallel_genes_B = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+'B'), "r")
    parallel_genes_S = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+'S'), "r")

    gene_significant_multiplicity_dict = {}


    parallel_genes_B_list = []
    for i, line in enumerate(parallel_genes_B):
        if i == 0:
            continue
        line = line.strip()
        items = line.split(",")
        parallel_genes_B_list.append(items[0])

        if items[0] not in gene_significant_multiplicity_dict:
            gene_significant_multiplicity_dict[items[0]] = {}
        gene_significant_multiplicity_dict[items[0]]['B'] = float(items[6].strip())

    parallel_genes_S_list = []
    for i, line in enumerate(parallel_genes_S):
        if i == 0:
            continue
        line = line.strip()
        items = line.split(",")
        parallel_genes_S_list.append(items[0])

        if items[0] not in gene_significant_multiplicity_dict:
            gene_significant_multiplicity_dict[items[0]] = {}
        gene_significant_multiplicity_dict[items[0]]['S'] = float(items[6].strip())


    mult_B_S = [(mult_taxa_dict[gene_name]['B'], mult_taxa_dict[gene_name]['S']) for gene_name in sorted(mult_taxa_dict) if (mult_taxa_dict[gene_name]['B'] > 0) and (mult_taxa_dict[gene_name]['S'] > 0) ]
    mult_significant_B_S = [(gene_significant_multiplicity_dict[gene_name]['B'], gene_significant_multiplicity_dict[gene_name]['S']) for gene_name in sorted(gene_significant_multiplicity_dict) if ('B' in gene_significant_multiplicity_dict[gene_name]) and ('S' in gene_significant_multiplicity_dict[gene_name]) ]

    mult_B = [x[0] for x in mult_B_S]
    mult_S = [x[1] for x in mult_B_S]

    mult_significant_B = [x[0] for x in mult_significant_B_S]
    mult_significant_S = [x[1] for x in mult_significant_B_S]

    mult_all = mult_B  + mult_S + mult_significant_B + mult_significant_S

    ax_regression.plot([ 0.0001, 1000  ], [   0.0001, 1000  ], lw = 3, c='k', ls = '--', zorder=1 )


    ax_regression.set_xlim([min(mult_all)*0.5, max(mult_all)*1.5])
    ax_regression.set_ylim([min(mult_all)*0.5, max(mult_all)*1.5])


    ax_regression.scatter(mult_B, mult_S, color=pt.get_colors(treatment),alpha=1,s=90, zorder=2)

    ax_regression.scatter(mult_significant_B, mult_significant_S, linewidth=3, facecolors=pt.get_colors(treatment), edgecolors='k', alpha=1, s=90, zorder=3)

    #ax_regression.set_xlim([  0.05, 200 ])
    #ax_regression.set_ylim([  0.05, 200 ])


    #if len(mult_significant_B) >= 3:
    #    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(mult_significant_B), np.log10(mult_significant_S))

        # hypothetical slope of 1
    #    ratio = (slope - 1) / std_err
    #    pval = stats.t.sf(np.abs(ratio), len(mult_significant_B)-2)*2
        # two sided or one sided?
    #    sys.stderr.write("Treatment %s-day slope t-test = %f, p = %f, df=%d\n" %  (str(10**int(treatment)), round(ratio, 3),  round(pval, 3),  len(mult_significant_B)-2  ) )

    #    ax_regression.text(0.1, 0.9, r'$\beta_{1} = $' + str(round(slope, 3)), fontsize=11, transform=ax_regression.transAxes)
    #    ax_regression.text(0.1, 0.8, r'$P \nless 0.05$', fontsize=11, transform=ax_regression.transAxes)

    venn = venn2(subsets = (len(parallel_genes_B_list), len(parallel_genes_S_list), len(set(parallel_genes_B_list) & set(parallel_genes_S_list))), ax=ax_venn, set_labels=('', ''), set_colors=(pt.get_colors(treatment), pt.get_colors(treatment)))
    c = venn2_circles(subsets=(len(parallel_genes_B_list), len(parallel_genes_S_list), len(set(parallel_genes_B_list) & set(parallel_genes_S_list))), ax=ax_venn, linestyle='dashed')
    #set_colors=(pt.get_colors(treatment), pt.get_colors(treatment)),

    c[0].set_ls('--')
    c[1].set_ls(':')
    c[0].set_lw(5)
    c[1].set_lw(5)
    c[0].set_edgecolor(pt.get_colors(treatment))
    c[1].set_edgecolor(pt.get_colors(treatment))
    #c[0].set_radius(len(parallel_genes_B_list) / 80 )
    #c[1].set_radius(len(parallel_genes_S_list) / 80)

fig.subplots_adjust(hspace=0.3, wspace=0.5)
fig_name = pt.get_path() + '/figs/multiplicity_B_vs_S.pdf'
fig.savefig(fig_name, format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
