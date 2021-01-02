from __future__ import division
import os, sys
import phylo_tools as pt
import parse_file
import timecourse_utils
import mutation_spectrum_utils
import numpy as np
from numpy.random import multinomial


treatments=pt.treatments
replicates = pt.replicates

G_subsample_dict = {}

G_all_mutations_dict = {}

ntot_subsample = 200
num_bootstraps = 10000
iter = 1000

for treatment in ['1']:

    # Load convergence matrix
    convergence_matrix_B = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+'B')))
    convergence_matrix_S = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % (treatment+'S')))

    populations_B = [treatment+'B' + replicate for replicate in replicates ]
    populations_S = [treatment+'S' + replicate for replicate in replicates ]

    gene_parallelism_statistics_B = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix_B,populations_B,Lmin=100)
    gene_parallelism_statistics_S = mutation_spectrum_utils.calculate_parallelism_statistics(convergence_matrix_S,populations_S,Lmin=100)

    #G_subsample_list = []
    #for i in range(subsamples):

    #    G_subsample = mutation_spectrum_utils.calculate_subsampled_total_parallelism(gene_parallelism_statistics, ntot_subsample=ntot_subsample)

    #    G_subsample_list.append(G_subsample)

    #G_subsample_list.sort()

    #G_subsample_dict[treatment+taxon] = G_subsample_list

    #G, pvalue = mutation_spectrum_utils.calculate_total_parallelism(gene_parallelism_statistics)

    #G_all_mutations_dict[treatment+taxon] = G

    Ls = []
    ns_B = []
    ns_S = []

    for gene_name in gene_parallelism_statistics_B.keys():

        Ls.append( gene_parallelism_statistics_B[gene_name]['length'] )
        ns_B.append( gene_parallelism_statistics_B[gene_name]['observed'] )
        ns_S.append( gene_parallelism_statistics_S[gene_name]['observed'] )


    Ls = np.asarray(Ls)
    ns_B = np.asarray(ns_B)
    ns_S = np.asarray(ns_S)

    Ltot = Ls.sum()
    ntot_B = ns_B.sum()
    ntot_S = ns_S.sum()

    ps = Ls*1.0/Ltot

    observed_G_list = []

    #for i in range(iter):

    #ns_subsample_B = multinomial(ntot_subsample, ns_B*1.0/ntot_B)
    #ns_subsample_S = multinomial(ntot_subsample, ns_S*1.0/ntot_S)

    #ns_subsample_B = multinomial(ntot_subsample, ns_B*1.0/ntot_B)
    #ns_subsample_S = multinomial(ntot_subsample, ns_S*1.0/ntot_S)

    ns_subsamle_merged = np.asarray([ns_B,  ns_S ]).T
    # remove zeros
    ns_subsamle_merged = ns_subsamle_merged[np.all(ns_subsamle_merged != 0, axis=1)]
    print(ns_subsamle_merged)
    ns_subsamle_merged_observed = np.sum(ns_subsamle_merged, axis=1)

    gs = ns_subsamle_merged_observed * np.abs(np.log((ns_subsamle_merged[:,0] / ns_subsamle_merged[:,0].sum()  ) / (ns_subsamle_merged[:,1] / ns_subsamle_merged[:,1].sum() )  ))

    observed_G = gs.sum() / ns_subsamle_merged_observed.sum()

    bootstrapped_Gs = []
    for bootstrap_idx in range(0,num_bootstraps):
        bootstrapped_ns_B = multinomial(ntot_B,ps)
        bootstrapped_ns_S = multinomial(ntot_S,ps)

        bootstrapped_ns_subsamle_merged = np.asarray([bootstrapped_ns_B,  bootstrapped_ns_S ]).T
        bootstrapped_ns_subsamle_merged = bootstrapped_ns_subsamle_merged[np.all(bootstrapped_ns_subsamle_merged != 0, axis=1)]
        bootstrapped_ns_subsamle_merged_observed = np.sum(bootstrapped_ns_subsamle_merged, axis=1)

        #print(bootstrapped_ns_subsamle_merged[:,0].sum(), bootstrapped_ns_subsamle_merged[:,1].sum())

        bootstrapped_gs = bootstrapped_ns_subsamle_merged_observed * np.abs(np.log((bootstrapped_ns_subsamle_merged[:,0] / bootstrapped_ns_subsamle_merged[:,0].sum()) / (bootstrapped_ns_subsamle_merged[:,1] / bootstrapped_ns_subsamle_merged[:,1].sum())))
        bootstrapped_G = bootstrapped_gs.sum()/bootstrapped_ns_subsamle_merged_observed.sum()

        bootstrapped_Gs.append(bootstrapped_G)

    bootstrapped_Gs = np.array(bootstrapped_Gs)

    pvalue = ((bootstrapped_Gs>=observed_G).sum()+1.0)/(len(bootstrapped_Gs)+1.0)
    print(observed_G, pvalue)


    print(bootstrapped_Gs)
