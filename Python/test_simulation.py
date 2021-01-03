from __future__ import division
import numpy
import random
import sys
import pickle

import phylo_tools as pt

numpy.random.seed(123456789)




def run_simulation():
    # weird sampling going on

    #4,292,969
    # mutation rate from Lynch paper, assume 10% sites are beneficial
    #mu = (3.28*10**-10 ) * 0.1
    #L =  4292969
    # keep order of magnitude for conveinance
    mu = (1.0*10**-10 )
    L =  1000000

    N = 10**6
    M = 10
    K = N/M
    c = 0.00001
    s_scale = 10**-3

    # average time in a dormant state = M

    n_active_to_dormant = int(c*N)
    n_dormant_to_active = int(c*K*M)

    if n_active_to_dormant != n_dormant_to_active:
        print("Unqueal number of individuals switching states!!")

    # rate of entering dormancy, per-capita = c
    # rate of exiting dormancy, per-capita = c*K
    #d = (c* K) / N
    #r = c / M

    # double mutants slow the simulation so we're assuming single mutants
    # e.g., the largest lineage size = 10**6, generated L*mu*N (~1000) mutants per-generation
    # probability that an individual gets two mutations ~= 10 ** -7


    generations_to_sample = [330*i for i in range(1, 11)]

    sampled_timepoints = {}

    generations = 3300


    n_clone_lineages = 0

    clone_size_dict = {}
    clone_size_dict[n_clone_lineages] = {}
    clone_size_dict[n_clone_lineages]['n_clone_active'] = N
    clone_size_dict[n_clone_lineages]['n_clone_dormant'] = M
    clone_size_dict[n_clone_lineages]['s'] = 1
    clone_size_dict[n_clone_lineages]['mutations'] = set([])

    # pre-assign fitness benefits to all sites
    all_sites = set(range(L))
    fitness_effects = numpy.random.exponential(scale=s_scale, size=L)


    # dict of what clones have a given mutation
    for generation in range(generations):
        # generate dormancy transition rates for all lineages
        # get keys and make sure they're in the same order
        #clones_active = [ clone_i for clone_i in clone_size_dict.keys() if ('n_clone_active' in clone_size_dict[clone_i]) and (clone_size_dict[clone_i]['n_clone_active'] > 0) ]
        #clones_active.sort()
        #clones_dormant = [ clone_i for clone_i in clone_size_dict.keys() if ('n_clone_dormant' in clone_size_dict[clone_i]) and (clone_size_dict[clone_i]['n_clone_dormant'] > 0)  ]
        #clones_dormant.sort()

        # get array of clone labels, the number of times each label is in the array is the size of the lineage
        clone_labels_active = [[int(clone_i)] * clone_size_dict[clone_i]['n_clone_active'] for clone_i in clone_size_dict.keys()]
        clone_labels_dormant = [[int(clone_i)] * clone_size_dict[clone_i]['n_clone_dormant'] for clone_i in clone_size_dict.keys() if ('n_clone_dormant' in clone_size_dict[clone_i]) and (clone_size_dict[clone_i]['n_clone_dormant'] > 0 )]

        clone_labels_active = numpy.concatenate(clone_labels_active).ravel()
        clone_labels_dormant = numpy.concatenate(clone_labels_dormant).ravel()

        clone_labels_active = clone_labels_active.astype(numpy.int)
        clone_labels_active = clone_labels_active.astype(numpy.int)


        # number of dormant individuals not constant???
        print(generation, len(clone_labels_active), len(clone_labels_dormant))
        active_to_dormant_sample = numpy.random.choice(clone_labels_active, size = n_active_to_dormant, replace=False)
        active_to_dormant_sample_bincount = numpy.bincount(active_to_dormant_sample)
        active_to_dormant_sample_bincount_nonzero = numpy.nonzero(active_to_dormant_sample_bincount)[0]

        dormant_to_active_sample = numpy.random.choice(clone_labels_dormant, size = n_dormant_to_active, replace=False)
        dormant_to_active_sample_bincount = numpy.bincount(dormant_to_active_sample)
        dormant_to_active_sample_bincount_nonzero = numpy.nonzero(dormant_to_active_sample_bincount)[0]

        for active_to_dormant_clone_i, active_to_dormant_n_clone_i in zip(active_to_dormant_sample_bincount_nonzero, active_to_dormant_sample_bincount[active_to_dormant_sample_bincount_nonzero]):

            clone_size_dict[active_to_dormant_clone_i]['n_clone_active'] -= active_to_dormant_n_clone_i

            if 'n_clone_dormant' not in clone_size_dict[active_to_dormant_clone_i]:
                clone_size_dict[active_to_dormant_clone_i]['n_clone_dormant'] = 0

            clone_size_dict[active_to_dormant_clone_i]['n_clone_dormant'] += active_to_dormant_n_clone_i


        for dormant_to_active_clone_i, dormant_to_active_n_clone_i in zip(dormant_to_active_sample_bincount_nonzero, dormant_to_active_sample_bincount[dormant_to_active_sample_bincount_nonzero]):

            clone_size_dict[dormant_to_active_clone_i]['n_clone_dormant'] -= dormant_to_active_n_clone_i

            if 'n_clone_dormant' not in clone_size_dict[dormant_to_active_clone_i]:
                clone_size_dict[dormant_to_active_clone_i]['n_clone_active'] = 0

            clone_size_dict[dormant_to_active_clone_i]['n_clone_active'] += dormant_to_active_n_clone_i

        # now move on to evolution
        for clone_i in list(clone_size_dict):

            if (clone_size_dict[clone_i]['n_clone_dormant'] == 0):

                if (clone_size_dict[clone_i]['n_clone_active'] == 0):
                    del clone_size_dict[clone_i]
                    continue

                else:
                    continue

            #print(clone_size_dict.keys())

            n_clone_i = clone_size_dict[clone_i]['n_clone_active']

            # mutation step#
            # lineage size can't be negative
            n_mutations_clone = min(numpy.random.poisson(mu*L*n_clone_i), n_clone_i)
            if n_mutations_clone == 0:
                continue
            # remove these individuals from the clone
            clone_size_dict[clone_i]['n_clone_active'] -= n_mutations_clone
            # all individuals in the clone have the same mutations
            # so just sample from nonmutated sites in the ancestral clone
            non_mutated_sites = all_sites - clone_size_dict[clone_i]['mutations']

            # sample without replacement
            #mutated_sites = random.sample(non_mutated_sites, n_mutations_clone)
            mutated_sites = numpy.random.choice(list(non_mutated_sites), size=n_mutations_clone, replace=False)
            #print(mutated_sites)
            #unique, counts = numpy.unique(mutated_sites, return_counts=True)
            for mutated_site in mutated_sites:

                n_clone_lineages += 1

                clone_size_dict[n_clone_lineages] = {}
                clone_size_dict[n_clone_lineages]['n_clone_active'] = 1
                clone_size_dict[n_clone_lineages]['n_clone_dormant'] = 0
                clone_size_dict[n_clone_lineages]['s'] = clone_size_dict[clone_i]['s'] + fitness_effects[mutated_site]
                clone_size_dict[n_clone_lineages]['mutations'] = clone_size_dict[clone_i]['mutations'].copy()
                clone_size_dict[n_clone_lineages]['mutations'].add(mutated_site)


            #if (clone_size_dict[clone_i]['n_clone_active'] == 0) and (clone_size_dict[clone_i]['n_clone_dormant'] == 0):
            #    del clone_size_dict[clone_i]


        #sampling_numerator = numpy.asarray( [ clone_size_dict[clone_i]['n_clone']*numpy.exp(clone_size_dict[clone_i]['s']) for clone_i in sorted(clone_size_dict.keys())] )
        sampling_numerator = numpy.asarray( [ clone_size_dict[clone_i]['n_clone_active']*numpy.exp(clone_size_dict[clone_i]['s']) for clone_i in clone_size_dict.keys()] )
        sampling_probability = sampling_numerator / sum(sampling_numerator)
        clone_sizes_after_selection = numpy.random.multinomial(N, sampling_probability)

        for clone_i_idx, clone_i in enumerate(list(clone_size_dict)):
            clone_i_size = clone_sizes_after_selection[clone_i_idx]

            #if clone_i_size == 0:
            #    del clone_size_dict[clone_i]
            #else:
            clone_size_dict[clone_i]['n_clone_active'] = clone_i_size


        if generation %100 == 0:

            sys.stderr.write("%d generations...\n" % generation)


        if generation in generations_to_sample:
            clone_size_dict_copy = clone_size_dict.copy()
            sampled_timepoints[generation] = clone_size_dict_copy



        N = sum( [ clone_size_dict[x]['n_clone_active'] for x in  clone_size_dict.keys() ] )
        M = sum( [ clone_size_dict[x]['n_clone_dormant'] for x in  clone_size_dict.keys() ] )

        print(generation, N, M)



    saved_data_file='%s/data/simulations/test2.dat' % (pt.get_path())

    with open(saved_data_file, 'wb') as outfile:
        pickle.dump(sampled_timepoints, outfile, protocol=pickle.HIGHEST_PROTOCOL)




def parse_simulation_output():

    saved_data_file='%s/data/simulations/test.dat' % (pt.get_path())
    sampled_timepoints = pickle.load( open(saved_data_file, "rb" ) )

    allele_freq_trajectory_dict = {}

    for key, value in sampled_timepoints.items():

        #print(value)

        N = sum( [ value[x]['n_clone_active'] for x in  value.keys() ] )
        M = sum( [ value[x]['n_clone_dormant'] for x in  value.keys() ] )

        print(N, M)

        #for clone, clone_dict in value.items():

        #    print(clone_dict['mutations'])




run_simulation()


#parse_simulation_output()
