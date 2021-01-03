from __future__ import division
import os, sys
import numpy as np

import phylo_tools as pt
import parse_file
import timecourse_utils

import  matplotlib.pyplot as plt
import scipy.stats as stats




def plot_allele_freqs_all_treats(strain):

    fig = plt.figure(figsize = (12, 6))

    row_count = 0
    for treatment in pt.treatments:

        column_count = 0

        for replicate in pt.replicates:

            population = treatment + strain + replicate
            annotated_timecourse_path = pt.get_path() + "/data/timecourse_final/%s_annotated_timecourse.txt" % population
            if os.path.exists(annotated_timecourse_path) == False:
                continue
            annotated_timecourse_file = open(annotated_timecourse_path ,"r")
            print(population)


            ax_i = plt.subplot2grid((3, 5), (row_count, column_count), colspan=1)

            first_line = annotated_timecourse_file.readline()
            first_line = first_line.strip()
            first_line_items = first_line.split(",")
            times = np.asarray([float(x.strip().split(':')[1]) for x in first_line_items[16::2]])
            times = np.insert(times, 0, 0, axis=0)


            times_to_ignore = None


            if times_to_ignore != None:
                times_to_ignore_idx=[np.where(times == t)[0] for t in times_to_ignore][0]
                times = np.delete(times, times_to_ignore_idx)

            for i, line in enumerate(annotated_timecourse_file):
                line = line.strip()
                items = line.split(",")
                pass_or_fail = items[15].strip()
                if pass_or_fail == 'FAIL':
                    continue

                alt_cov = np.asarray([ float(x) for x in items[16::2]])
                total_cov = np.asarray([ float(x) for x in items[17::2]])
                # pseudocount to avoid divide by zero error
                alt_cov = alt_cov
                total_cov = total_cov + 1
                freqs = alt_cov / total_cov
                freqs = np.insert(freqs, 0, 0, axis=0)

                if times_to_ignore != None:
                    freqs = np.delete(freqs, times_to_ignore_idx)

                rgb = pt.mut_freq_colormap()
                rgb = pt.lighten_color(rgb, amount=0.5)

                if len(times) == len(freqs) + 1:
                    freqs = np.insert(freqs, 0, 0, axis=0)

                ax_i.plot(times, freqs, '.-', c=rgb, alpha=0.4)

            ax_i.set_xlim([0,max(times)])
            ax_i.set_ylim([0, 1])

            if column_count == 0:
                ax_i.set_ylabel( str(10**int(treatment)) + '-day transfers', fontsize =12  )

            ax_i.tick_params(axis="x", labelsize=8)
            ax_i.tick_params(axis="y", labelsize=8)

            column_count += 1
        row_count +=1


    fig.text(0.5, 0.04, 'Days, ' + r'$t$', ha='center', va='center', fontsize=24)
    fig.text(0.05, 0.5, 'Allele frequency, ' + r'$f(t)$', ha='center', va='center', rotation='vertical',  fontsize=24)

    fig.suptitle(pt.latex_dict[strain], fontsize=28, fontweight='bold')
    fig_name = pt.get_path() + '/figs/mut_trajectories_%s.jpg'
    #fig_name = pt.get_path() + '/figs/mut_trajectories_%s.pdf'
    fig.savefig(fig_name % strain, bbox_inches = "tight",  format='jpg',  pad_inches = 0.4)#, dpi = 600)
    #fig.savefig(fig_name % strain, bbox_inches = "tight",format='png', pad_inches = 0.4)#, dpi = 600)

    plt.close()



# ['B','D','F','J','P', 'S']:

for taxon in ['B', 'S']:
    #plot_within_taxon_paralleliism(taxon)
    plot_allele_freqs_all_treats(taxon)
