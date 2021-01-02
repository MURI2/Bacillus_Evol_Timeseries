import numpy
import sys
import parse_file
import mutation_spectrum_utils
from math import log10,log,exp
import stats_utils

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import phylo_tools as pt

FDR = 0.05
nmin = 3
fmax_divider = 0.5

if len(sys.argv) > 1:
    level = sys.argv[1]
else:
    level = 'gene'

########
#
# Set up parallelism results all taxa
#
########

# first get data, then make plot


taxa = [ 'B', 'S']
treatments=pt.treatments
#treatments = ['2']
replicates = pt.replicates


for taxon in taxa:

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    # to get the common gene names for each ID
    gene_name_dict = dict(zip(gene_names, genes ))
    protein_id_dict = dict(zip(gene_names, protein_ids ))

    for treatment in treatments:

        populations = [treatment+taxon + replicate for replicate in replicates ]

        sys.stderr.write("Analyzing %s level parallelism for %s...\n" % (level,treatment+taxon))

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        # Calculate basic parallellism statistics
        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix, populations)

        # Calculate G score for entire gene (G=n*g)
        gene_G_scores = mutation_spectrum_utils.calculate_G_scores(gene_parallelism_statistics)
        pooled_G_scores = numpy.asarray(list(gene_G_scores.values()))

        pooled_G_scores.sort()

        null_G_survival = mutation_spectrum_utils.NullGeneGSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics)

        observed_Gs, observed_G_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_G_scores)

        # Do same thing for multiplicity statistic
        pooled_multiplicities = numpy.array([gene_parallelism_statistics[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics.keys()])
        pooled_multiplicities.sort()

        null_multiplicity_survival = mutation_spectrum_utils.NullGeneMultiplicitySurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )

        observed_ms, observed_multiplicity_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_multiplicities)

        # Do same thing for num hits
        pooled_hits = numpy.array([gene_parallelism_statistics[gene_name]['observed'] for gene_name in gene_parallelism_statistics.keys()])
        pooled_hits.sort()

        null_uniform_hit_survival = mutation_spectrum_utils.NullUniformGeneHitSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )

        null_hit_survival = mutation_spectrum_utils.NullGeneHitSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )

        observed_ns, observed_hit_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_hits)

        # Do same thing for pvalues
        gene_qvalues, gene_pvalues = mutation_spectrum_utils.calculate_parallelism_qvalues(gene_parallelism_statistics)

        gene_logpvalues = mutation_spectrum_utils.calculate_parallelism_logpvalues(gene_parallelism_statistics)
        pooled_pvalues = []
        for gene_name in gene_logpvalues.keys():
            if gene_parallelism_statistics[gene_name]['observed']>=nmin:
                pooled_pvalues.append( gene_logpvalues[gene_name] )
        pooled_pvalues = numpy.asarray(pooled_pvalues)
        if len(pooled_pvalues) == 0:
            continue
        pooled_pvalues.sort()

        null_pvalue_survival = mutation_spectrum_utils.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics, nmin=nmin)

        observed_ps, observed_pvalue_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_pvalues, min_x=-4)

        # Print counts of gene hits, just for reference
        for n in range(0,21):
            sys.stdout.write('%d %d-hit %ss\n' % ((numpy.fabs(pooled_hits-n)<0.5).sum(), n, level))

        # Pvalue version
        threshold_idx = numpy.nonzero((null_pvalue_survival(observed_ps)*1.0/observed_pvalue_survival)<FDR)[0][0]
        pstar = observed_ps[threshold_idx] # lowest value where this is true
        num_significant = observed_pvalue_survival[threshold_idx]

        sys.stdout.write("Found %d significant %ss (p* = %g)\n" % (num_significant, level, exp(-pstar)))


        ntot = 0
        nsignificant = 0
        Ltot = 0
        Lsignificant = 0

        nonsignificant_genes = []

        output_file = open(pt.get_path() +'/data/timecourse_final/' +  ("parallel_%ss_%s.txt" % (level, treatment+taxon)) ,"w")
        # print header
        output_file.write(", ".join(["Locus tag", "RefSeq protein ID", "Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))


        output_file_not_significant = open(pt.get_path() +'/data/timecourse_final/' +  ("parallel_not_significant_%ss_%s.txt" % (level, treatment+taxon)) ,"w")
        # print header
        output_file_not_significant.write(", ".join(["Locus tag", "RefSeq protein ID", "Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))


        sys.stdout.write("-log p^* = %g\n" % pstar)
        #sys.stderr.write("Nonsignificant genes:\n")

        for gene_name in sorted(gene_parallelism_statistics, key=lambda x: gene_parallelism_statistics.get(x)['observed'],reverse=True):

            ntot += gene_parallelism_statistics[gene_name]['observed']
            Ltot += gene_parallelism_statistics[gene_name]['length']

            if gene_logpvalues[gene_name] >= pstar and gene_parallelism_statistics[gene_name]['observed']>=nmin:
            #if gene_G_scores[gene_name]>=Gstar:
            #if gene_qvalues[gene_name]<FDR:

                nsignificant += gene_parallelism_statistics[gene_name]['observed']
                Lsignificant += gene_parallelism_statistics[gene_name]['length']

                output_file.write("\n")
                output_file.write("%s, %s, %s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, protein_id_dict[gene_name], gene_name_dict[gene_name], gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], gene_logpvalues[gene_name]))

            else:

                if gene_parallelism_statistics[gene_name]['observed']>=nmin:
                        #sys.stderr.write("%s, %s,  %s, %0.1f, %d, %0.2f, %0.2f, %g\n" % (gene_name, protein_id_dict[gene_name], gene_name_dict[gene_name], gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], gene_logpvalues[gene_name]))

                        output_file_not_significant.write("\n")
                        output_file_not_significant.write("%s, %s, %s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, protein_id_dict[gene_name], gene_name_dict[gene_name], gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], gene_logpvalues[gene_name]))

                nonsignificant_genes.append(gene_name)

        output_file.close()
        output_file_not_significant.close()

        observed_G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics, nonsignificant_genes)

        sys.stdout.write("Significant %ss:\n" % level)
        sys.stdout.write("n = %g (%0.2f of total), L = %g (%0.2f of total)\n" % (nsignificant, nsignificant*1.0/ntot, Lsignificant, Lsignificant*1.0/Ltot))
        sys.stdout.write("Remaining total parallelism = %g (p=%g)\n" % (observed_G, pvalue))








for taxon in taxa:

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    # to get the common gene names for each ID
    gene_name_dict = dict(zip(gene_names, genes ))
    protein_id_dict = dict(zip(gene_names, protein_ids ))

    for treatment in treatments:


        populations = [treatment+taxon + replicate for replicate in replicates ]

        sys.stderr.write("Analyzing %s level parallelism for %s...\n" % (level,treatment+taxon))

        # Load convergence matrix
        convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+taxon)))

        # Calculate basic parallellism statistics
        gene_parallelism_statistics = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix, populations,fmax_min=0.5)

        ## Calculate G score for entire gene (G=n*g)
        #gene_G_scores = mutation_spectrum_utils.calculate_G_scores(gene_parallelism_statistics)
        #pooled_G_scores = numpy.asarray(list(gene_G_scores.values()))

        #pooled_G_scores.sort()

        #null_G_survival = mutation_spectrum_utils.NullGeneGSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics)

        #observed_Gs, observed_G_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_G_scores)

        # Do same thing for multiplicity statistic
        pooled_multiplicities = numpy.array([gene_parallelism_statistics[gene_name]['multiplicity'] for gene_name in gene_parallelism_statistics.keys()])
        pooled_multiplicities.sort()

        #null_multiplicity_survival = mutation_spectrum_utils.NullGeneMultiplicitySurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )

        #observed_ms, observed_multiplicity_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_multiplicities)

        # Do same thing for num hits
        pooled_hits = numpy.array([gene_parallelism_statistics[gene_name]['observed'] for gene_name in gene_parallelism_statistics.keys()])
        pooled_hits.sort()

        #null_uniform_hit_survival = mutation_spectrum_utils.NullUniformGeneHitSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )

        #null_hit_survival = mutation_spectrum_utils.NullGeneHitSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics )

        #observed_ns, observed_hit_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_hits)

        # Do same thing for pvalues
        gene_qvalues, gene_pvalues = mutation_spectrum_utils.calculate_parallelism_qvalues(gene_parallelism_statistics)

        gene_logpvalues = mutation_spectrum_utils.calculate_parallelism_logpvalues(gene_parallelism_statistics)
        pooled_pvalues = []
        for gene_name in gene_logpvalues.keys():
            if gene_parallelism_statistics[gene_name]['observed']>=nmin:
                pooled_pvalues.append( gene_logpvalues[gene_name] )
        pooled_pvalues = numpy.asarray(pooled_pvalues)
        if len(pooled_pvalues) == 0:
            continue
        pooled_pvalues.sort()

        null_pvalue_survival = mutation_spectrum_utils.NullGeneLogpSurvivalFunction.from_parallelism_statistics( gene_parallelism_statistics, nmin=nmin)

        observed_ps, observed_pvalue_survival = stats_utils.calculate_unnormalized_survival_from_vector(pooled_pvalues, min_x=-4)

        # Pvalue version
        threshold_idx = numpy.nonzero((null_pvalue_survival(observed_ps)*1.0/observed_pvalue_survival)<FDR)[0][0]
        pstar = observed_ps[threshold_idx] # lowest value where this is true
        num_significant = observed_pvalue_survival[threshold_idx]

        sys.stdout.write("Found %d significant %ss (p* = %g)\n" % (num_significant, level, exp(-pstar)))


        ntot = 0
        nsignificant = 0
        Ltot = 0
        Lsignificant = 0

        nonsignificant_genes = []

        output_file = open(pt.get_path() +'/data/timecourse_final/' +  ("parallel_%ss_%s_fmax_min.txt" % (level, treatment+taxon)) ,"w")
        # print header
        output_file.write(", ".join(["Locus tag", "RefSeq protein ID", "Gene", "Length", "Observed", "Expected", "Multiplicity", "-log10(P)"]))


        for gene_name in sorted(gene_parallelism_statistics, key=lambda x: gene_parallelism_statistics.get(x)['observed'],reverse=True):

            ntot += gene_parallelism_statistics[gene_name]['observed']
            Ltot += gene_parallelism_statistics[gene_name]['length']

            if gene_logpvalues[gene_name] >= pstar and gene_parallelism_statistics[gene_name]['observed']>=nmin:
            #if gene_G_scores[gene_name]>=Gstar:
            #if gene_qvalues[gene_name]<FDR:

                nsignificant += gene_parallelism_statistics[gene_name]['observed']
                Lsignificant += gene_parallelism_statistics[gene_name]['length']

                output_file.write("\n")
                output_file.write("%s, %s, %s, %0.1f, %d, %0.2f, %0.2f, %g" % (gene_name, protein_id_dict[gene_name], gene_name_dict[gene_name], gene_parallelism_statistics[gene_name]['length'],  gene_parallelism_statistics[gene_name]['observed'], gene_parallelism_statistics[gene_name]['expected'], gene_parallelism_statistics[gene_name]['multiplicity'], gene_logpvalues[gene_name]))


        output_file.close()
