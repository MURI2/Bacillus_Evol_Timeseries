from __future__ import division
import timecourse_utils
import parse_file
from math import log,exp
import numpy
from scipy.special import gammaln as loggamma
import sys
from scipy.optimize import fmin
from scipy.interpolate import interp1d
from scipy.stats import nbinom

import phylo_tools as pt

import uniform_prior_well_mixed_hmm as well_mixed_model


taxon = sys.argv[1]
treatments = ['0', '1', '2']
reps = ['1', '2', '3', '4', '5']


min_coverage = 5

samples_to_remove_dict = pt.samples_to_remove

for treatment in treatments:
    for replicate in reps:
        population = treatment + taxon + replicate

        if population in pt.populations_to_ignore:
            continue

        if population in pt.samples_to_remove:
            times_to_ignore = pt.samples_to_remove[population]
        else:
            times_to_ignore=None

        mutations = []
        good_mutations = []

        sys.stderr.write("Processing %s\n" % population)

        mutation_data, depth_tuple = parse_file.parse_annotated_timecourse(population)

        #times = mutation_data[0][10]
        times = mutation_data[0][9]

        if times_to_ignore != None:
            times_to_ignore_idx=[numpy.where(times == t)[0] for t in times_to_ignore][0]
        else:
            times_to_ignore_idx=None
        #times = np.delete(times, times_to_ignore_idx)
        #alts = np.delete(alts, times_to_ignore_idx)
        #depths = np.delete(depths, times_to_ignore_idx)

        for mutation_idx in range(0,len(mutation_data)):

            location, gene_name, allele, var_type, codon, position_in_codon, AA_count, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutation_data[mutation_idx]

            good_idxs, masked_alts, masked_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue, min_coverage)
            if times_to_ignore_idx != None:
                times = numpy.delete(times, times_to_ignore_idx)
                good_idxs = numpy.delete(good_idxs, times_to_ignore_idx)
                masked_alts = numpy.delete(masked_alts, times_to_ignore_idx)
                masked_depths = numpy.delete(masked_depths, times_to_ignore_idx)

            max_freq = (masked_alts*1.0/(masked_depths+(masked_depths<1))).max()

            num_good_timepoints = (masked_depths>0.5).sum()

            if num_good_timepoints > 5:
                mutations.append((masked_alts*1.0, masked_depths*1.0))
                if var_type!='sv' and var_type!='indel' and masked_depths[-1]>1:
                #if masked_depths[-1]>1:
                    good_mutations.append((masked_alts*1.0, masked_depths*1.0))
            else:
                mutations.append((numpy.zeros_like(times)*1.0, numpy.zeros_like(times)))


        A = []
        D = []
        for i in range(0,len(mutations)):
            Ai,Di = mutations[i]
            A.append(Ai)
            D.append(Di)

        A = numpy.array(A)
        D = numpy.array(D)

        Pstate, Ls, p0, p1 = well_mixed_model.infer_hmm(A,D,num_iterations=10)
        # decode ML states
        #Ls = numpy.zeros_like(A)
        #for i in xrange(0,A.shape[0]):
        #    Ls[i] = (Pstate[i,:,:]).argmax(axis=1)

        # make thing that says if filtered:

        filtered_mutations = numpy.ones_like(A)
        for i in range(0,A.shape[0]):
            if D[i,:].sum() < 0.5:
                filtered_mutations[i,:] = numpy.zeros_like(times)

        ns = numpy.zeros((Pstate.shape[2],Pstate.shape[1]))

        for l in range(0,ns.shape[0]):
            ns[l,:] = (Pstate[:,:,l]*filtered_mutations).sum(axis=0)

        #Ls = numpy.zeros_like(A)
        #for i in xrange(0,A.shape[0]):
        #    if filtered_mutations[i,0] < 0.5:
        #        Ls[i] = numpy.zeros_like(times)
        #    else:
        #        Ls[i] = (Pstate[i,:,:]).argmax(axis=1)

        hard_ns = numpy.zeros_like(ns)
        for l in range(0,ns.shape[0]):
            hard_ns[l,:] = (Ls==l).sum(axis=0)

        # Write stuff
        haplotype_filename = pt.get_path() + '/data/timecourse_final/' +("%s_well_mixed_state_timecourse.txt" % population)
        haplotype_file = open(haplotype_filename, "w")
        haplotype_file.write("%s\n" % ", ".join([str(t) for t in times]))
        haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[0,:]]))
        haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[1,:]]))
        haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[2,:]]))
        haplotype_file.write("%s\n" % ", ".join([str(nl) for nl in hard_ns[3,:]]))
        for i in range(0,Ls.shape[0]):
            haplotype_file.write("%s\n" % ", ".join([str(l) for l in Ls[i,:]]))

        haplotype_file.close()
