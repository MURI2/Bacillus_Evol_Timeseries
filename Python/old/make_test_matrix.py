from __future__ import division
import os, imp, sys, bz2
from datetime import date
import phylo_tools as pt
import numpy as np
import parse_file as pf

storage_directory = 'cluster_output_files/'
final_directory = pt.get_path() + '/data/final_directory'

#treatments = ['0', '1', '2']
treatments = ['0', '1']

reps = ['1', '2', '3', '4', '5']
#reps = ['4']
# '1', '2', '3', '4',
#strains = ['B', 'C', 'D', 'F', 'J', 'P', 'S']
#strains = ['B', 'C', 'D']
strains = ['S']

def get_likelihood(t_star, t, d_pm):
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
    return (t_star, n_gt, sigma_gt_t_star, n_lte, sigma_lt_t_star, r_t_star, l_t_star)



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
        return np.nan,np.nan,np.nan
    if (len(t_pm_f_star) == 2) and ((t_pm_f_star[1] - t_pm_f_star[0]) > 500):
        return np.nan,np.nan,np.nan

    I_list = []
    for win_size in list(range(2, len(t_pm_f_star)+1)):
        for win_start in  range(len(t_pm_f_star) - win_size +1):
            t1 = t_pm_f_star[win_start]
            t2 = t_pm_f_star[win_start + win_size-1]
            f_pm_window = np.asarray(f_pm_f_star[win_start: win_start + win_size])
            I_list.append( sum(f_pm_window - f_star) )
    I = max(I_list)

    # average frequency relaxation time
    mean_freq = sum(A_pm) / sum(D_pm)
    if n_0 > 0.3*n:
        T = 0
    else:
        T = 0
        #A_pm
        for t_dash in range(2, len(t_pm)):
            partial_mean_freq = sum(A_pm[:t_dash+1]) / sum(D_pm[:t_dash+1])
            if (partial_mean_freq <= 0.6*mean_freq) and (t_dash > T):
                T = t_dash

    return C_star, I, T

def process_output():
    for strain in strains:
        parse_gene_list = pf.parse_gene_list(taxon=strain)
        for treatment in treatments:
            for rep in reps:
                print('%s%s%s' % (treatment, strain, rep) )
                sample = '%s%s%s' % (treatment, strain, rep)
                snp_timecourse_filename = '%s%s%s_snp_timecourse.bz' % (treatment, strain, rep)
                snp_timecourse_path = pt.get_path() + '/data/timecourse_snp/' + snp_timecourse_filename
                snp_file = bz2.open(snp_timecourse_path, "rt")
                depth_timecourse_filename = '%s%s%s_depth_timecourse.bz' % (treatment, strain, rep)
                depth_timecourse_path = pt.get_path() + '/data/timecourse_depth/' + depth_timecourse_filename
                depth_file = bz2.open(depth_timecourse_path, "rt")
                for depth in depth_file:
                    depth_split = [x.strip() for x in depth.split(',')]
                    D_pt_median = depth_split[-1].split(' ')
                    D_pt_median = np.asarray([float(x) for x in D_pt_median])

                file_mutations = []
                for snp in snp_file:
                    snp_split = [x.strip() for x in snp.split(',')]
                    t_pm = snp_split[3].split(' ')
                    t_pm = np.asarray([int(x) for x in t_pm])
                    A_pm = snp_split[4].split(' ')
                    A_pm = np.asarray([int(x) for x in A_pm])
                    D_pm = snp_split[5].split(' ')
                    D_pm = np.asarray([int(x) for x in D_pm])
                    if len(t_pm) == 1:
                        continue
                    # remove D_pmt < 5
                    remove_D_pmt = np.asarray([x for x,y in enumerate(D_pm) if y < 5])
                    D_pt_median_copy = np.empty_like(D_pt_median)
                    D_pt_median_copy[:] = D_pt_median
                    if len(remove_D_pmt) > 0:
                        t_pm = np.delete(t_pm, remove_D_pmt)
                        A_pm = np.delete(A_pm, remove_D_pmt)
                        D_pm = np.delete(D_pm, remove_D_pmt)
                        D_pt_median_copy = np.delete(D_pt_median_copy, remove_D_pmt)
                    # remove low coverage timepoints
                    d_pm = D_pm / D_pt_median_copy
                    # don't look at trajectories with fewer than four
                    if len(t_pm[1:]) < 4:
                        continue
                    l_list = [get_likelihood(t_, t=t_pm, d_pm=d_pm) for t_ in t_pm[1:-2] ]
                    l_list = sorted(l_list, key=lambda x: x[-1])
                    max_l = l_list[-1]
                    # r threshold of 0.5 too conservative
                    #if max_l[-2] >= 0.5:
                    #    continue
                    n = len(t_pm)
                    sigma = np.sqrt( (1/n)*sum(d_pm**2) - ( ((1/n)*sum(d_pm ))**2 ) )
                    if (max_l[4] == 0) or (max_l[2] == 0):
                        continue
                    delta_l = max_l[3] * np.log(sigma/max_l[4]) + max_l[1]*np.log(sigma/max_l[2])
                    # try permutation test for log likelihood?

                    # need upper threshold for delta_l, chose 20 for now
                    #if delta_l > 20:
                    #    continue

                    C_star_mut, I_mut, T_mut = get_test_statistics(t_pm, A_pm, D_pm)
                    if np.isnan(np.sum([C_star_mut, I_mut, T_mut])):
                        continue


                    f_pm = A_pm/D_pm
                    if (f_pm[0] > 0.9) and (f_pm[-1] > 0.9):
                        continue
                    #if f_pm[-1] > 0.05:
                    if (sample == '0C1') or (sample == '0C4'):
                        last_timepoint = -2
                    else:
                        last_timepoint = -1
                    if f_pm[last_timepoint] > 0.05:

                        allele_split = snp_split[2].split('->')
                        anc = allele_split[0]
                        der = allele_split[1]
                        site = int(snp_split[1])
                        genes_names_site = []
                        for k in list(range(len(parse_gene_list[0]))):
                            gene_name = parse_gene_list[0][k]
                            start = int(parse_gene_list[1][k])
                            stop = int(parse_gene_list[2][k])
                            if (site >= start ) and (site <= stop):
                                genes_names_site.append(gene_name)
                        genes_names_site_merged = ",".join(genes_names_site)

                        file_mutations.append([snp_split[0], genes_names_site_merged, snp_split[1], anc, der, str(f_pm[-1])])
                # out file
                header = ['Contig', 'Locus_tag', 'Site', 'Ancestral', 'Derived', 'Final_freq']
                mutation_filename = '%s%s%s_snp_final.txt' % (treatment, strain, rep)
                mutation_file = open(pt.get_path() + '/data/timecourse_prelim_poly/' + mutation_filename, 'w')
                mutation_file.write('\t'.join(header) + '\n')
                for mutation in file_mutations:

                    mutation_file.write("\t".join(mutation))
                    mutation_file.write("\n")

                mutation_file.close()

process_output()
