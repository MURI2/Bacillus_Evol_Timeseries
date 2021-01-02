from __future__ import division
import numpy

numpy.random.seed(123456789)

#4,292,969
mu = 3.28*10**-10
L =  4292969

N = 10**6
M = 1000
generations = 10
total_fitness = 0

s_scale = 10**-5

# mutation dict with lineage size, mutation IDs, fitness effects
clone_size_dict = {}
clone_size_dict[0] = {}

# dict of what clones have a given mutation
clone_size_dict[0]['n_clone'] = N
clone_size_dict[0]['mutations'] = []
clone_size_dict[0]['fitness'] = 0

n_clone_lineages = 0
n_mutations = 0



def wrap_choice(list_to_sample, no_samples):
    list_size = len(list_to_sample)
    takes = no_samples // list_size
    samples = list_to_sample * (no_samples // list_size) + list(numpy.random.choice(list_to_sample, no_samples - takes * list_size))
    return samples



def sample_zero_truncated_poisson(_lambda):
    k=1
    t = numpy.exp(-_lambda) / (1 - numpy.exp(-_lambda)) * _lambda
    s = t
    u = numpy.random.uniform()
    while s < u:
        k += 1
        t *= _lambda / k
        s += t

    return k



for g in range(generations):

    print(g, len(clone_size_dict))

    # mutation with DFE step
    for clone_i in list(clone_size_dict):

        #clone_dict_i = clone_size_dict[clone_i]

        #n_clone_i = clone_dict_i['n_clone']
        n_clone_i = clone_size_dict[clone_i]['n_clone']

        n_mutations_clone = numpy.random.poisson(mu*L*n_clone_i)
        fitness_effects = numpy.random.exponential(scale=s_scale, size=n_mutations_clone)

        #s = numpy.random.poisson(len(fitness_effects)/clone_size )
        mu = n_mutations_clone/n_clone_i
        k_count = 0
        while k_count < n_mutations_clone:
            # get number of mutations in new clone
            k = sample_zero_truncated_poisson(mu)
            # get fitness effects
            fitnesses_clone = fitness_effects[k_count:k_count+k]
            # new clone
            n_clone_lineages += 1

            # make new dict
            clone_size_dict[n_clone_lineages] = {}
            clone_size_dict[n_clone_lineages]['n_clone'] = 1
            clone_size_dict[n_clone_lineages]['fitness'] = clone_size_dict[clone_i]['fitness'] + sum(fitnesses_clone)
            clone_size_dict[n_clone_lineages]['mutations'] = clone_size_dict[clone_i]['mutations'].copy()
            # add new mutations to dict
            for k_j in range(k):
                n_mutations += 1
                clone_size_dict[n_clone_lineages]['mutations'].append(n_mutations)

            # remove individual from ancestor
            clone_size_dict[clone_i]['n_clone'] -= 1
            # delete dictionry item if empty
            if clone_size_dict[clone_i]['n_clone'] == 0:
                del clone_size_dict[clone_i]

            k_count += k


    # multinomial sampling for selection
    # calculate denominator
    #sampling_denominator = sum([ clone_size_dict[clone_i]['n_clone']*numpy.exp(clone_size_dict[clone_i]['fitness']) for clone_i in clone_size_dict.keys()])

    sampling_numerator = numpy.asarray( [ clone_size_dict[clone_i]['n_clone']*numpy.exp(clone_size_dict[clone_i]['fitness']) for clone_i in sorted(clone_size_dict.keys())] )
    sampling_probability = sampling_numerator / sum(sampling_numerator)
    clone_sizes_after_selection = numpy.random.multinomial(N, sampling_probability)

    for clone_i_idx, clone_i in enumerate(list(clone_size_dict)):
        clone_i_size = clone_sizes_after_selection[clone_i_idx]

        if clone_i_size == 0:
            del clone_size_dict[clone_i]
        else:
            clone_size_dict[clone_i]['n_clone'] = clone_i_size


    #print(clone_size_dict)


#print( clone_size_dict.keys())

# prob sampled =( n_clone_i * (1+relative_s)) / (n_clone_i * (1+relative_s)  +  )


#print(clone_size_dict)


    #numpy.concatenate((numpy.repeat(fitness_effects, clone_size // len(fitness_effects), numpy.random.choice(fitness_effects, clone_size - clone_size // len(fitness_effects), replace=False ))))


    #print( wrap_choice(fitness_effects, 10000))

# split mutations into fitness classes
#xxxxx = numpy.random.choice(fitnesses, size=N, replace=False, p=None)
#print(xxxxx)





    # multinomial sampling for GOING dormant
    #dormancy_probability = [clone_size_dict[clone_i]['n_clone_active']/N for clone_i in clones_active]
    # multiniomial sampling for resuscitating
    #resuscitation_probability = [clone_size_dict[clone_i]['n_clone_dormant']/N for clone_i in clones_dormant]
    # number of individuals in each linege that ENTER dormancy
    #dormant_clone_changes = numpy.random.multinomial(c*N, dormancy_probability)
    # number of individuals in each linege that GET RESUSCITATED
    #resucitated_clone_changes = numpy.random.multinomial(c*K*M, resuscitation_probability)

    # might have to hard code something to make sure lineage size doesn't get < 0....we'll see
    #for dormant_clone_change_idx, dormant_clone_change in enumerate(dormant_clone_changes):

    #    clone_name = clones_active[dormant_clone_change_idx]

    #    # take out as many individuals as you can and then throw the extra
    #    # one to a random clone
    #    if dormant_clone_change > clone_size_dict[clone_name]['n_clone_active']:

    #        clone_size_dict[clone_name]['n_clone_active'] -= clone_size_dict[clone_name]['n_clone_active']


    #    else:

    #        clone_size_dict[clone_name]['n_clone_active'] -= dormant_clone_change

    #    if 'n_clone_dormant' not in clone_size_dict[clone_name]:
    #        clone_size_dict[clone_name]['n_clone_dormant'] = 0

    #    clone_size_dict[clone_name]['n_clone_dormant'] += dormant_clone_change


    #for resucitated_clone_change_idx, resucitated_clone_change in enumerate(resucitated_clone_changes):

    #    clone_name = clones_active[resucitated_clone_change_idx]
    #    clone_size_dict[clone_name]['n_clone_dormant'] -= resucitated_clone_change

    #    if 'n_clone_active' not in clone_size_dict[clone_name]:
    #        clone_size_dict[clone_name]['n_clone_active'] = 0

    #    clone_size_dict[clone_name]['n_clone_active'] += resucitated_clone_change





# assign mutations to individuals
#testtt = numpy.random.choice(fitnesses, size=3, replace=False, p=clonse_sizes/sum(clonse_sizes))
