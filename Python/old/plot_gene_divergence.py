from __future__ import division
import os, sys
import numpy
from math import log10,log,exp

import  matplotlib.pyplot as plt
#from matplotlib.lines import Line2D
#import matplotlib.gridspec as gridspec

import parse_file
import timecourse_utils
import stats_utils
import mutation_spectrum_utils
import phylo_tools as pt

numpy.random.seed(123456789)

FDR = 0.05
nmin = 3

fmax_divider = 0.5

if len(sys.argv) > 1:
    level = sys.argv[1]
else:
    level = 'gene'

significant_gene_dict_all_taxa = {}

for taxon in ['B']:

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    # to get the common gene names for each ID
    #gene_name_dict = dict(zip(gene_names, genes ))
    #protein_id_dict = dict(zip(gene_names, protein_ids ))

    for treatment in ['0','1']:

        if treatment+taxon in pt.treatment_taxa_to_ignore:
            sys.stderr.write("Skipping %s, too few surviving replicates ...\n" % (treatment+taxon))
            continue

        populations = [treatment+taxon + replicate for replicate in pt.replicates ]

        sys.stderr.write("Analyzing %s level parallelism for %s...\n" % (level,treatment+taxon))

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        # Calculate basic parallellism statistics
        gene_parallelism_statistics_minor = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix, populations,fmax_max=0.5)
        gene_parallelism_statistics_major = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix, populations,fmax_min=0.5)

        # Do same thing for multiplicity statistic
        pooled_multiplicities_minor = numpy.array([gene_parallelism_statistics_minor[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics_minor.keys()])
        pooled_multiplicities_minor.sort()
        pooled_multiplicities_major = numpy.array([gene_parallelism_statistics_major[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics_major.keys()])
        pooled_multiplicities_major.sort()


        gene_logpvalues_minor = mutation_spectrum_utils.calculate_parallelism_logpvalues(gene_parallelism_statistics_minor)
        pooled_pvalues_minor = []
        for gene_name in gene_logpvalues_minor.keys():
            if gene_parallelism_statistics_minor[gene_name]['observed']>=nmin:
                pooled_pvalues_minor.append( gene_logpvalues_minor[gene_name] )
        pooled_pvalues_minor = numpy.asarray(pooled_pvalues_minor)
        if len(pooled_pvalues_minor) == 0:
            continue
        pooled_pvalues_minor.sort()

        # same thing for major freqs
        gene_logpvalues_major = mutation_spectrum_utils.calculate_parallelism_logpvalues(gene_parallelism_statistics_major)
        pooled_pvalues_major = []
        for gene_name in gene_logpvalues_major.keys():
            if gene_parallelism_statistics_major[gene_name]['observed']>=nmin:
                print(gene_parallelism_statistics_major[gene_name])
                pooled_pvalues_major.append( gene_logpvalues_major[gene_name] )
        pooled_pvalues_major = numpy.asarray(pooled_pvalues_major)
        if len(pooled_pvalues_major) == 0:
            continue
        pooled_pvalues_major.sort()

        print(pooled_pvalues_major)


        null_pvalue_survival_minor = mutation_spectrum_utils.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics_minor, nmin=nmin)
        null_pvalue_survival_major = mutation_spectrum_utils.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics_major, nmin=nmin)

        observed_ps_minor, observed_pvalue_survival_minor = stats_utils.calculate_unnormalized_survival_from_vector(pooled_pvalues_minor, min_x=-4)
        observed_ps_major, observed_pvalue_survival_major = stats_utils.calculate_unnormalized_survival_from_vector(pooled_pvalues_major, min_x=-8)

        # Pvalue version
        threshold_idx_minor = numpy.nonzero((null_pvalue_survival_minor(observed_ps_minor)*1.0/observed_pvalue_survival_minor)<FDR)[0][0]
        pstar_minor = observed_ps_minor[threshold_idx_minor] # lowest value where this is true
        num_significant_minor = observed_pvalue_survival_minor[threshold_idx_minor]

        # Pvalue version for minor
        threshold_idx_major = numpy.nonzero((null_pvalue_survival_major(observed_ps_major)*1.0/observed_pvalue_survival_major)<FDR)[0][0]
        pstar_major = observed_ps_major[threshold_idx_major] # lowest value where this is true
        num_significant_major = observed_pvalue_survival_major[threshold_idx_major]

        sys.stdout.write("Found %d significant %ss for fmax >= %f (p* = %g)\n" % (num_significant_minor, level, fmax_divider, exp(-pstar_minor)))
        sys.stdout.write("Found %d significant %ss for fmax < %f (p* = %g)\n" % (num_significant_major, level, fmax_divider, exp(-pstar_major)))


        gene_parallelism_statistics_significant_minor = {}

        for gene_name in sorted(gene_parallelism_statistics_minor, key=lambda x: gene_parallelism_statistics_minor.get(x)['observed'],reverse=True):

            if gene_logpvalues_minor[gene_name] >= pstar_minor and gene_parallelism_statistics_minor[gene_name]['observed']>=nmin:

                gene_parallelism_statistics_significant_minor[gene_name] = {}
                gene_parallelism_statistics_significant_minor[gene_name]['length'] = gene_parallelism_statistics_minor[gene_name]['length']
                gene_parallelism_statistics_significant_minor[gene_name]['observed'] = gene_parallelism_statistics_minor[gene_name]['observed']
                gene_parallelism_statistics_significant_minor[gene_name]['expected'] = gene_parallelism_statistics_minor[gene_name]['expected']
                gene_parallelism_statistics_significant_minor[gene_name]['multiplicity'] = gene_parallelism_statistics_minor[gene_name]['multiplicity']


        gene_parallelism_statistics_significant_major = {}

        for gene_name in sorted(gene_parallelism_statistics_major, key=lambda x: gene_parallelism_statistics_major.get(x)['observed'],reverse=True):


            if treatment == '1':
                if gene_parallelism_statistics_major[gene_name]['observed'] >0:
                    print(gene_parallelism_statistics_major[gene_name])

            if gene_logpvalues_major[gene_name] >= pstar_major and gene_parallelism_statistics_major[gene_name]['observed']>=nmin:

                gene_parallelism_statistics_significant_major[gene_name] = {}
                gene_parallelism_statistics_significant_major[gene_name]['length'] = gene_parallelism_statistics_major[gene_name]['length']
                gene_parallelism_statistics_significant_major[gene_name]['observed'] = gene_parallelism_statistics_major[gene_name]['observed']
                gene_parallelism_statistics_significant_major[gene_name]['expected'] = gene_parallelism_statistics_major[gene_name]['expected']
                gene_parallelism_statistics_significant_major[gene_name]['multiplicity'] = gene_parallelism_statistics_major[gene_name]['multiplicity']

        significant_gene_dict_all_taxa[treatment+taxon] = {}
        significant_gene_dict_all_taxa[treatment+taxon]['minor'] = gene_parallelism_statistics_significant_minor
        significant_gene_dict_all_taxa[treatment+taxon]['major'] = gene_parallelism_statistics_significant_major


print(significant_gene_dict_all_taxa['0C']['major'])


for freq_cutoff in ['minor', 'major']:


    dict_freq_1 = significant_gene_dict_all_taxa['0C'][freq_cutoff]
    dict_freq_10 = significant_gene_dict_all_taxa['1C'][freq_cutoff]

    #print(freq_cutoff, set(dict_freq_10.keys()))
