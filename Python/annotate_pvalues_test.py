from __future__ import division
import sys, copy, bz2
import numpy as np


depth_filename = sys.argv[1]
snp_filename = sys.argv[2]
indel_filename = sys.argv[3]
likelihood_filename = sys.argv[4]

print(depth_filename, snp_filename, indel_filename, likelihood_filename)

# where to read the input junctions from
depth_file = bz2.open(depth_filename, "rt")
#snp_file = bz2.open(snp_filename, "rt")
indel_file = bz2.open(indel_filename, "rt")
#likelihood_file = bz2.open(likelihood_filename, "wt")

np.random.seed(123456789)

n_iter = 100
n_iter_max = 100

# 1/2 bc we're using two metrics
P_star = (0.05) ** (1/2)
trajectories_count = 0
#lte =less than or equal
#gte = grater than or equal

times_to_exclude = {'0C1': [900], '0D1':[600], '0D2':[900], '0D3':[900], '0D4':[900],
                    '1B5': [1000], '1J3': [400], '1J4': [500], '2S2':[700], '2S3':[700],
                    '2S4':[700], '2S5':[700]}

line_name = depth_filename.split('/')[-1].split('_')[0]
print(line_name)

for depth in depth_file:
    depth_split = [x.strip() for x in depth.split(',')]
    D_pt_median = depth_split[-1].split(' ')
    D_pt_median = np.asarray([int(float(x)) for x in D_pt_median])

    times = depth_split[-3].split(' ')
    times = [int(float(x)) for x in times]
    if line_name in times_to_exclude:
        times_to_exclude_idx = np.asarray([times.index(x) for x in times_to_exclude[line_name]])
        D_pt_median = np.delete(D_pt_median, times_to_exclude_idx)
    else:
        times_to_exclude_idx = None





def get_likelihood(t_star, t, d_pm, sigma):

    t_star_index = np.where(t==t_star)[0][0]
    d_pm_gt = d_pm[t_star_index+1:]
    n_gt = len(d_pm_gt)
    d_pm_lte = d_pm[:t_star_index+1]
    n_lte = len(d_pm_lte)

    lambda_0_t_star = (1/n_lte) * sum(d_pm_lte)
    r_t_star = (1/lambda_0_t_star) * (1/n_gt) *  sum(d_pm_gt)
    sigma_lt_t_star = np.sqrt( (1/n_lte)* sum((d_pm_lte - lambda_0_t_star) **2) )
    sigma_gt_t_star = np.sqrt( (1/n_gt)* sum((d_pm_gt - (lambda_0_t_star*r_t_star)) **2) )
    l_t_star = -1*n_lte * np.log(sigma_lt_t_star) - n_gt*np.log(sigma_gt_t_star)

    delta_l_t_star = n_lte*np.log(sigma/sigma_lt_t_star) + n_gt*np.log(sigma/sigma_gt_t_star)

    return (t_star, n_gt, sigma_gt_t_star, n_lte, sigma_lt_t_star, r_t_star, l_t_star, delta_l_t_star)



def get_test_statistics(t_pm, A_pm, D_pm):
    # autocorrelation
    f_pm = A_pm/D_pm
    f_bar = min((sum(A_pm) / sum(D_pm)), 0.5)
    autocorr_filter = np.where(abs(A_pm- f_bar*D_pm) < 1)[0]
    t_pm_auto = np.delete(t_pm, autocorr_filter)
    A_pm_auto = np.delete(A_pm, autocorr_filter)
    D_pm_auto = np.delete(D_pm, autocorr_filter)
    f_pm_auto = np.delete(f_pm, autocorr_filter)
    C_num_star = 0
    C_denom_star = 0
    for t_i in list(range(len(t_pm_auto)-1)):
        C_num_star += np.sqrt(D_pm_auto[t_i +1] * D_pm_auto[t_i]) * \
        (f_pm_auto[t_i +1] - f_bar) * (f_pm_auto[t_i] - f_bar) * \
        np.heaviside(abs(A_pm_auto[t_i +1] - f_bar*D_pm_auto[t_i +1])-1, 1/2) * \
        np.heaviside(abs(A_pm_auto[t_i] - f_bar*D_pm_auto[t_i])-1, 1/2)
        C_denom_star += np.sqrt(D_pm_auto[t_i +1] * D_pm_auto[t_i])


    if C_denom_star == 0:
        C_star = 0

    else:
        C_star = C_num_star / ( (f_bar**2) * C_denom_star)

    # derived allele sojourn weight
    n_0 = len(np.where(A_pm == 0)[0])
    n = len(t_pm)
    if n_0 > 0.3*n:
        f_star = f_bar / (1+ np.exp( (n_0 - 0.3*n)/5 ))
    else:
        f_star = f_bar / (1+ np.exp(-1* (0.3*n - n_0)/5 ))

    t_pm_f_star = np.delete(t_pm, np.where(f_pm <= f_star)[0])
    f_pm_f_star = np.delete(f_pm, np.where(f_pm <= f_star)[0])
    if len(t_pm_f_star) < 2:
        return np.nan,np.nan
    if (len(t_pm_f_star) == 2) and ((t_pm_f_star[1] - t_pm_f_star[0]) > 500):
        return np.nan,np.nan

    I_list = []
    for win_size in list(range(2, len(t_pm_f_star)+1)):
        for win_start in  range(len(t_pm_f_star) - win_size +1):
            t1 = t_pm_f_star[win_start]
            t2 = t_pm_f_star[win_start + win_size-1]
            f_pm_window = np.asarray(f_pm_f_star[win_start: win_start + win_size])
            I_list.append( sum(f_pm_window - f_star) )
    I = max(I_list)
    # average frequency relaxation time
    #mean_freq = sum(A_pm) / sum(D_pm)
    #if n_0 > 0.3*n:
    #    T = 0
    #else:
    #    T = 0
    #    #A_pm
    #    for t_dash in range(2, len(t_pm)):
    #        partial_mean_freq = sum(A_pm[:t_dash+1]) / sum(D_pm[:t_dash+1])
    #        if (partial_mean_freq <= 0.6*mean_freq) and (t_dash > T):
    #            T = t_dash
    # dont return T
    return C_star, I




def get_composite_p_value(t_pm, A_pm, D_pm):
    # moderate coverage timepoints
    #C_star_mut, I_mut, T_mut = get_test_statistics(t_pm, A_pm, D_pm)
    C_star_mut, I_mut = get_test_statistics(t_pm, A_pm, D_pm)
    if np.isnan(np.sum([C_star_mut, I_mut])):
        return None, None
    # error model for moderate-coverage timepoints
    f_mean = (sum(A_pm) / sum(D_pm))
    c = max(1, (1/2)*f_mean)
    f_pm = A_pm/D_pm
    f_pm[np.isnan(f_pm)] = 0
    C_star_mut_null = []
    I_mut_null = []
    #T_mut_null = []
    count = copy.copy(n_iter)
    count2 = copy.copy(n_iter)
    count_attempts = copy.copy(n_iter_max*2)
    #while (count_attempts > 0):

    print(count)

    while (count >= 0) and (count_attempts >= 0):
        #print(count, count_attempts)

        p_hat_pm = np.random.permutation(f_pm) * c
        alpha_pm = np.random.poisson(D_pm*p_hat_pm)
        #print(f_pm, D_pm*(1-p_hat_pm))
        beta_pm = np.random.poisson(D_pm*(1-p_hat_pm))
        #print((alpha_pm/ (alpha_pm+beta_pm)) * D_pm)
        A_hat_pm = np.round( (alpha_pm/ (alpha_pm+beta_pm)) * D_pm )
        if np.count_nonzero(A_hat_pm) == 0:
            count_attempts -= 1
            continue
        C_star_mut_i, I_mut_i = get_test_statistics(t_pm, A_hat_pm, D_pm)
        if np.isnan(np.sum([C_star_mut_i, I_mut_i])):
            count_attempts -= 1
            continue
        C_star_mut_null.append(C_star_mut_i)
        I_mut_null.append(I_mut_i)
        count -= 1


    if len(C_star_mut_null) >= n_iter:

        P_C_star = (len([x for x in C_star_mut_null if x > C_star_mut]) +1) /  ( n_iter +1)
        P_I = (len([x for x in I_mut_null if x > I_mut]) +1) / ( n_iter +1)
        P_list = [P_C_star, P_I]
        comp_stat = sum([np.heaviside(P_star-P_k, 1/2) * np.log(1/P_k) for P_k in P_list])
        comp_stat_null = []
        while (count2 != 0):

            p_hat_pm = np.random.permutation(f_pm) * c
            alpha_pm = np.random.poisson(D_pm*p_hat_pm)
            beta_pm = np.random.poisson(D_pm*(1-p_hat_pm))
            A_hat_pm = np.round( (alpha_pm/ (alpha_pm+beta_pm)) * D_pm )
            A_hat_pm[np.isnan(A_hat_pm)] = 0

            C_star_mut_i, I_mut_i = get_test_statistics(t_pm, A_hat_pm, D_pm)
            if np.isnan(np.sum([C_star_mut_i, I_mut_i])):
                count_attempts -= 1
                continue
            # add pseudocount of 1
            P_C_star_i = (len([x for x in C_star_mut_null if x > C_star_mut_i]) +1) / (n_iter +1)
            P_I_i = (len([x for x in I_mut_null if x > I_mut_i]) +1) / (n_iter +1)

            P_list_i = [P_C_star_i, P_I_i]
            comp_stat_i = sum([np.heaviside(P_star-P_k, 1/2) * np.log(1/P_k) for P_k in P_list_i])
            comp_stat_null.append(comp_stat_i)
            count2 -= 1

        P_comp_stat = len([x for x in comp_stat_null if x > comp_stat]) / len(comp_stat_null)

        return C_star_mut, P_comp_stat

    else:
        return None, None




def filter_trajectory(snp_split):

    t_pm = snp_split[3].split(' ')
    t_pm = np.asarray([int(float(x)) for x in t_pm])
    A_pm = snp_split[4].split(' ')
    A_pm = np.asarray([int(float(x)) for x in A_pm])
    D_pm = snp_split[5].split(' ')
    D_pm = np.asarray([int(float(x)) for x in D_pm])

    if times_to_exclude_idx != None:
        t_pm = np.delete(t_pm, times_to_exclude_idx)
        A_pm = np.delete(A_pm, times_to_exclude_idx)
        D_pm = np.delete(D_pm, times_to_exclude_idx)


    count_At_Atplus1 = 0
    for i in range(0,len(A_pm)-1):
        if (A_pm[i] != 0) and (A_pm[i+1] != 0):
            count_At_Atplus1 += 1
    if count_At_Atplus1 == 0:
        return None, None, None
    if len(t_pm) < 4:
        return None, None, None
    A_D_stack = np.stack((A_pm, D_pm, A_pm/D_pm), axis=-1 )
    # criterion 1
    A_D_stack_At_2 = [x for x in A_D_stack if (x[0] >= 2) ]
    if len(A_D_stack_At_2) < 2:
        return None, None, None
    if len([x for x in A_D_stack_At_2 if (x[1] >=10) and (x[2] >=0.05)]) <1:
        return None, None, None
    # criterion 2
    if (D_pm[0] < 10) and (A_pm[0]/D_pm[0] >0.1):
        return None, None, None
    # criteria 3
    if len([x for x in A_D_stack if (x[0] >= 2) and (x[1]>=5) ]) < 3:
        return None, None, None
    # criterion 4
    if len([x for x in A_D_stack if (x[0] >= 3) and (x[1]>=5) and (x[2]>=0.1) ]) < 1:
        return None, None, None
    # criterion 5
    if max(A_pm/D_pm) - min((sum(A_pm) / sum(D_pm)), 0.5) < 0.1:
        return None, None, None

    #if (np.mean(A_pm/D_pm) > 0.85) and (np.mean(D_pm) > 90):
    #    return None, None, None

    return t_pm, A_pm, D_pm




for indel in indel_file:
    trajectories_count += 1
    #if trajectories_count != 51:
    #    continue
    #if (trajectories_count % 100) == 0:
    print("Indel trajectory " + str(trajectories_count))
    indel_split = [x.strip() for x in indel.split(',')]
    t_pm, A_pm, D_pm = filter_trajectory(indel_split)
    if t_pm is None:
        continue

    d_pm = D_pm / D_pt_median

    sigma = np.sqrt((1/len(d_pm))*sum(d_pm**2) - (( (1/len(d_pm)) * sum(d_pm))**2))

    l_list = [get_likelihood(t_, t=t_pm, d_pm=d_pm, sigma=sigma) for t_ in t_pm[1:-1] ]
    # remove items with r < 0.5 to focus on deletions
    l_list_candidates = [ l_list_i for l_list_i in l_list if (l_list_i[-3] < 0.5)  and (np.isinf(l_list_i[-1])==False)  ]

    if len(l_list_candidates) != 0:
        # sort by increasing value
        l_list_candidates = sorted(l_list_candidates, key=lambda x: x[-1])
        delta_l_max = l_list_candidates[-1][-1]
        t_star = l_list_candidates[-1][0]
        delta_l_max_null_list = []
        count_low_cov = copy.copy(n_iter)
        total_iters = 0
        while count_low_cov != 0:
            # ranomize all but last two
            t_pm_null = np.copy(t_pm)
            np.random.shuffle(t_pm_null)
            l_list_null = [get_likelihood(t_null, t=t_pm_null, d_pm=d_pm, sigma=sigma) for t_null in t_pm_null[1:-1]]

            l_list_null = [ l_list_i for l_list_i in l_list_null if (l_list_i[-3] < 0.5)  and (np.isinf(l_list_i[-1])==False)  ]
            if len(l_list_null) == 0:
                continue
            else:
                l_list_null = sorted(l_list_null, key=lambda x: x[-1])
                delta_l_max_null_list.append(l_list_null[-1][-1])
                count_low_cov -= 1

        # r threshold of 0.5 too conservative
        P_delta_l = (len([x for x in delta_l_max_null_list if x > delta_l_max]) +1) / (len(delta_l_max_null_list) +1)
        if P_delta_l < 0.1:
            print("Deletion detected: excluding samples after t-star from rest of analysis")
            # remove samples after t_star

            after_t_star_idx = np.where(t_pm > t_star)[0]
            t_pm_t_star = np.delete(t_pm, after_t_star_idx)
            A_pm_t_star = np.delete(A_pm, after_t_star_idx)
            D_pm_t_star = np.delete(D_pm, after_t_star_idx)

            delta_l_null_list = delta_l_max

            D_pm[D_pm == 0] = 1

            C_star_mut, P_comp_stat = get_composite_p_value(t_pm_t_star, A_pm_t_star, D_pm_t_star)

        else:
            t_star = delta_l_null_list = P_delta_l = None
            D_pm[D_pm == 0] = 1
            C_star_mut, P_comp_stat = get_composite_p_value(t_pm, A_pm, D_pm)


    else:
        t_star = delta_l_null_list = P_delta_l = None
        D_pm[D_pm == 0] = 1
        C_star_mut, P_comp_stat = get_composite_p_value(t_pm, A_pm, D_pm)


    #if C_star_mut == None:
    #    continue

    indel_split.append( " ".join([ str(C_star_mut), str(P_comp_stat), str(t_star), str(delta_l_null_list), str(P_delta_l)]) )
    indel_split_merge = ", ".join(indel_split)

    print(indel_split_merge)

    #likelihood_file.write(indel_split_merge)
    #likelihood_file.write("\n")





depth_file.close()
#snp_file.close()
indel_file.close()
#likelihood_file.close()


#python ~/GitHub/Phylo_Evol_Timeseries/Python/annotate_pvalues_test.py ~/GitHub/Phylo_Evol_Timeseries/data/timecourse_depth/2B5_depth_timecourse.bz what ~/GitHub/Phylo_Evol_Timeseries/data/timecourse_indel/2B5_indel_timecourse.bz what
