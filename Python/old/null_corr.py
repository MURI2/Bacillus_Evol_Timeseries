
    sys.stderr.write("Generating null distributions...\n")
    for treatment in ['1']:
        for taxon in taxa:
            r2s_null = []
            for replicate_idx_i in range(0,len(replicates)):

                population_i = treatment + taxon + replicates[replicate_idx_i]
                sys.stderr.write("Processing %s...\n" % population_i)

                mutations_i, depth_tuple_i = parse_file.parse_annotated_timecourse(population_i)
                population_avg_depth_times_i, population_avg_depths_i, clone_avg_depth_times_i, clone_avg_depths_i = depth_tuple_i
                state_times_i, state_trajectories_i = parse_file.parse_well_mixed_state_timecourse(population_i)

                times_i = mutations_i[0][9]
                Ms_i = np.zeros_like(times_i)*1.0
                fixed_Ms_i = np.zeros_like(times_i)*1.0

                for mutation_idx_i_k in range(0,len(mutations_i)):

                    location_i_k, gene_name_i_k, allele_i_k, var_type_i_k, test_statistic_i_k, pvalue_i_k, cutoff_idx_i_k, depth_fold_change_i_k, depth_change_pvalue_i_k, times_i_k, alts_i_k, depths_i_k, clone_times_i_k, clone_alts_i_k, clone_depths_i_k = mutations_i[mutation_idx_i_k]

                    state_Ls_i_k = state_trajectories_i[mutation_idx_i_k]
                    good_idx_i_k, filtered_alts_i_k, filtered_depths_i_k = timecourse_utils.mask_timepoints(times_i_k, alts_i_k, depths_i_k, var_type_i_k, cutoff_idx_i_k, depth_fold_change_i_k, depth_change_pvalue_i_k)
                    freqs_i_k = timecourse_utils.estimate_frequencies(filtered_alts_i_k, filtered_depths_i_k)

                    masked_times_i_k = times_i[good_idx_i_k]
                    masked_freqs_i_k = freqs_i_k[good_idx_i_k]
                    masked_state_Ls_i_k = state_Ls_i_k[good_idx_i_k]

                    P_idx_i_k = np.where(masked_state_Ls_i_k == 3)[0]
                    if len(P_idx_i_k) < min_trajectory_length:
                        continue
                    first_P_i_k = P_idx_i_k[0]
                    last_P_i_k = P_idx_i_k[-1]

                    masked_freqs_P_i_k = masked_freqs_i_k[first_P_i_k:last_P_i_k+1]
                    masked_times_P_i_k = masked_times_i_k[first_P_i_k:last_P_i_k+1]

                    delta_masked_freqs_P_i_k = masked_freqs_P_i_k[1:] - masked_freqs_P_i_k[:-1]
                    delta_masked_times_P_i_k = masked_times_P_i_k[:-1]

                    # now go through all the other mutations in the remaining populations
                    for replicate_idx_j in range(replicate_idx_i+1,len(replicates)):

                        population_j = treatment + taxon + replicates[replicate_idx_j]

                        mutations_j, depth_tuple_j = parse_file.parse_annotated_timecourse(population_j)
                        population_avg_depth_times_j, population_avg_depths_j, clone_avg_depth_times_j, clone_avg_depths_j = depth_tuple_j
                        state_times_j, state_trajectories_j = parse_file.parse_well_mixed_state_timecourse(population_j)

                        times_j = mutations_j[0][9]
                        Ms_j = np.zeros_like(times_j)*1.0
                        fixed_Ms_j = np.zeros_like(times_j)*1.0

                        for mutation_idx_j_l in range(0,len(mutations_i)):

                            location_j_l, gene_name_j_l, allele_j_l, var_type_j_l, test_statistic_j_l, pvalue_j_l, cutoff_idx_j_l, depth_fold_change_j_l, depth_change_pvalue_j_l, times_j_l, alts_j_l, depths_j_l, clone_times_j_l, clone_alts_j_l, clone_depths_j_l = mutations_i[mutation_idx_j_l]

                            state_Ls_j_l = state_trajectories_i[mutation_idx_j_l]
                            good_idx_j_l, filtered_alts_j_l, filtered_depths_j_l = timecourse_utils.mask_timepoints(times_j_l, alts_j_l, depths_j_l, var_type_j_l, cutoff_idx_j_l, depth_fold_change_j_l, depth_change_pvalue_j_l)
                            freqs_j_l = timecourse_utils.estimate_frequencies(filtered_alts_j_l, filtered_depths_j_l)

                            masked_times_j_l = times_i[good_idx_j_l]
                            masked_freqs_j_l = freqs_j_l[good_idx_j_l]
                            masked_state_Ls_j_l = state_Ls_j_l[good_idx_j_l]

                            P_idx_j_l = np.where(masked_state_Ls_j_l == 3)[0]
                            if len(P_idx_j_l) < min_trajectory_length:
                                continue
                            first_P_j_l = P_idx_j_l[0]
                            last_P_j_l = P_idx_j_l[-1]

                            masked_freqs_P_j_l = masked_freqs_j_l[first_P_j_l:last_P_j_l+1]
                            masked_times_P_j_l = masked_times_j_l[first_P_j_l:last_P_j_l+1]

                            delta_masked_freqs_P_j_l = masked_freqs_P_j_l[1:] - masked_freqs_P_j_l[:-1]
                            delta_masked_times_P_j_l = masked_times_P_j_l[:-1]

                            intersect_times = np.intersect1d(delta_masked_times_P_i_k, delta_masked_times_P_j_l)

                            if len(intersect_times)>=3:

                                intersect_idx_i_k = [np.where(delta_masked_times_P_i_k == intersect_time)[0][0] for intersect_time in intersect_times ]
                                intersect_delta_i_k = delta_masked_freqs_P_i_k[intersect_idx_i_k]

                                intersect_idx_j_l = [np.where(delta_masked_times_P_j_l == intersect_time)[0][0] for intersect_time in intersect_times ]
                                intersect_delta_j_l = delta_masked_freqs_P_j_l[intersect_idx_j_l]


                                if len(intersect_delta_i_k) != len(intersect_delta_i_k):
                                    print(len(intersect_delta_i_k), len(intersect_delta_j))

                                r2 = stats.pearsonr(intersect_delta_i_k, intersect_delta_j_l)[0] ** 2
                                r2s_null.append(r2)


            r2s_null_dict[treatment + taxon] = r2s_null
