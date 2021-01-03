############
#
# Calculates convergence matrices for different mutation identity classes
#
# These are of the form:
#
############

import numpy
import sys
import parse_file

import timecourse_utils
import phylo_tools as pt

###################
#
# Load gene information
#
###################

#taxa = pt.taxa
taxa = pt.taxa
treatments = pt.treatments
replicates = pt.replicates

for taxon in taxa:
    sys.stderr.write("Calculating gene sizes for %s...\t" % taxon)
    gene_size_map = parse_file.create_gene_size_map(taxon)
    sys.stderr.write("Done!\n")
    for treatment in treatments:
        if treatment+taxon == '1J':
            continue
        populations = [treatment+ taxon + replicate for replicate in replicates]

        sys.stderr.write('Calculating convergence matrix for mutations...\n')

        excluded_types = set(['synonymous'])

        identifier_size_map = gene_size_map

        def get_identifier_name(gene_name):
            return gene_name

        convergence_matrix = {}
        for identifier_name in sorted(identifier_size_map.keys()):

            length = identifier_size_map[identifier_name]

            convergence_matrix[identifier_name] = {'length': length, 'mutations': {population: [] for population in populations}}

        position_snv_dict = {}

        all_poly_list = []
        all_fixed_list = []
        overlap_poly_count = 0
        overlap_fixed_count = 0
        polymorphic_count = 0
        fixed_count = 0

        for population in populations:
            sys.stderr.write("Processing %s...\t" % population)

            # calculate mutation trajectories
            # load mutations
            mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)

            population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple

            state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

            num_processed_mutations = 0

            for mutation_idx in range(0,len(mutations)):

                position, gene_name, allele, var_type, codon, position_in_codon, AAs_count, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
                #Ls = haplotype_trajectories[mutation_idx]

                state_Ls = state_trajectories[mutation_idx]

                if gene_name=='intergenic':
                    continue

                if var_type in excluded_types:
                    continue

                identifier = get_identifier_name(gene_name)

                if identifier==None:
                    sys.stderr.write("No identifier for %s!\n" % gene_name)
                    continue


                good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

                freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

                masked_times = times[good_idxs]
                masked_freqs = freqs[good_idxs]
                masked_state_Ls = state_Ls[good_idxs]
                #masked_Ls = Ls[good_idxs]

                if (sum([masked_state_Ls==parse_file.POLYMORPHIC][0]) == 0) and (sum([masked_state_Ls==parse_file.FIXED][0]) == 0):
                    continue

                if sum([masked_state_Ls==parse_file.FIXED][0]) != 0:
                    fixed_count += 1
                    if (position, allele) in all_fixed_list:
                        overlap_fixed_count += 1
                    all_fixed_list.append((position, allele) )

                else:
                    polymorphic_count += 1
                    if (position, allele) in all_poly_list:
                        overlap_poly_count += 1
                    all_poly_list.append((position, allele) )


                num_processed_mutations += 1

                t = timecourse_utils.calculate_appearance_time(masked_times, masked_freqs, masked_state_Ls)

                convergence_matrix[identifier]['mutations'][population].append((t, masked_state_Ls[-1], masked_freqs[-1], max(masked_freqs)))

            sys.stderr.write("processed %d mutations!\n" % num_processed_mutations)

        # Print it out
        output_filename = pt.get_path() + '/data/timecourse_final/' + ("%s_convergence_matrix.txt" % (treatment+ taxon))

        convergence_matrix_file = open(output_filename,"w")

        # Header
        convergence_matrix_file.write(", ".join(["Identifier"]+["Size"]+[population for population in populations]))

        for identifier in sorted(convergence_matrix.keys()):

            length = convergence_matrix[identifier]['length']
            mutations = convergence_matrix[identifier]['mutations']

            convergence_matrix_file.write("\n")
            convergence_matrix_file.write(", ".join([identifier, "%0.1f" % length]+[";".join(["%d:%d:%g:%g" % (t,L,f,f_max) for t,L,f,f_max in mutations[population]]) for population in populations]))

        convergence_matrix_file.close()

        sys.stderr.write("Non-unique fraction of polymorphisms = %0.3f \n" % (overlap_poly_count/(polymorphic_count+1))  )

        sys.stderr.write("Non-unique fraction of fixations = %0.3f \n" % (overlap_fixed_count/(fixed_count+1))  )
