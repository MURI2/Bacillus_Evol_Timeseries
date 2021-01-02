from __future__ import division
import os, sys
import numpy as np

import  matplotlib.pyplot as plt
from matplotlib.colors import ColorConverter

import scipy.stats as stats
import statsmodels.api as sm

import phylo_tools as pt
import parse_file
import timecourse_utils
import mutation_spectrum_utils


significant_multiplicity_dict = {}

for taxon in pt.taxa:
    significant_multiplicity_dict[taxon] = {}
    for treatment_idx, treatment in enumerate(pt.treatments):

        significant_multiplicity_taxon_path = pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon)
        if os.path.exists(significant_multiplicity_taxon_path) == False:
            continue
        significant_multiplicity_taxon = open(significant_multiplicity_taxon_path, "r")
        for i, line in enumerate( significant_multiplicity_taxon ):
            if i == 0:
                continue
            line = line.strip()
            items = line.split(",")
            if items[0] not in significant_multiplicity_dict[taxon]:
                significant_multiplicity_dict[taxon][items[0]] = {}
            significant_multiplicity_dict[taxon][items[0]][treatment] = float(items[-2])


treatment_pairs = [['0','1'],['0','2'],['1','2']]
fig = plt.figure(figsize = (9, 6))
plt.axhline( y=0, color='k', lw=2.5, linestyle='--', alpha = 1, zorder=1)
ax_count=0
for treatment_pair in treatment_pairs:

    for taxon in pt.taxa:

        result = [(x[treatment_pair[0]],x[treatment_pair[1]]) for x in significant_multiplicity_dict[taxon].values() if (treatment_pair[0] in x) and (treatment_pair[1] in x)]
        if len(result) ==0:
            continue

        mult_x = np.log10([x[0] for x in result])
        mult_y = np.log10([x[1] for x in result])

        #slope, intercept, r_value, p_value, std_err = stats.linregress(mult_x, mult_y)
        mult_x = sm.add_constant(mult_x)
        mod = sm.OLS(mult_y, mult_x)
        res = mod.fit()
        slope = res.params[1]
        CI_025 = res.conf_int(0.05)[1][0]
        CI_975 = res.conf_int(0.05)[1][1]

        def flip_slope_and_CIs(slope, CI_025, CI_975,null=1):
            new_slope = abs(slope-null)
            delta_CI_025 = abs(CI_025-slope)
            delta_CI_975 = abs(CI_975-slope)

            if slope < null:
                new_CI_025 = new_slope - delta_CI_975
                new_CI_975 = new_slope + delta_CI_025
            else:
                new_CI_025 = new_slope - delta_CI_025
                new_CI_975 = new_slope + delta_CI_975

            return new_slope,new_CI_025,new_CI_975


        new_slope,new_CI_025,new_CI_975 = flip_slope_and_CIs(slope, CI_025, CI_975)

        # get mean colors
        ccv = ColorConverter()

        color_1 = np.array(ccv.to_rgb( pt.get_colors( treatment_pair[0] ) ))
        color_2 = np.array(ccv.to_rgb( pt.get_colors( treatment_pair[1] ) ))

        mix_color = 0.7 * (color_1 + color_2)
        mix_color = np.min([mix_color, [1.0, 1.0, 1.0]], 0)

        if (treatment_pair[0] == '0') and (treatment_pair[1] == '1'):
            #mix_color = pt.lighten_color(mix_color, amount=2.8)
            mix_color = 'gold'

        plt.errorbar(ax_count, new_slope, yerr = [ [new_slope-new_CI_025], [new_CI_975-new_slope]], \
                fmt = 'o', alpha = 1, barsabove = True, marker = pt.plot_species_marker(taxon), \
                mfc = 'white', mec = 'white', lw=3, c = 'k', zorder=2, ms=17)

        plt.scatter(ax_count, new_slope, marker=pt.plot_species_marker(taxon), s = 250, \
            linewidth=2, facecolors=mix_color, edgecolors='k', alpha=1, zorder=3)

        ax_count+=1

    plt.axvline( x=ax_count-0.5, color='k', lw=2, linestyle=':', alpha = 1, zorder=1)

    #plt.text(ax_count-2, 0.1, '%s-day vs. %s-day' %(str(10**int(treatment_pair[0])), str(10**int(treatment_pair[1]))),  fontsize=14)


#plt.xlabel("Transfer time (days)", fontsize = 20)
#plt.xticks((0,1,2), ('1', '10', '100'), fontsize=14  )
plt.xticks([], [])


plt.xticks([2,7.5,13], ['1-day vs. 10-days', '1-day vs. 100-days', '10-days vs. 100-days'], fontsize=13)
plt.tick_params(axis='x', labelsize=14, length = 0)



plt.ylabel("Difference between multiplicity slope\nand null expectation (" + r'$\left | \beta_{1} - 1 \right |$' +')' , fontsize = 18)

fig_name = pt.get_path() + '/figs/divergence_slopes.pdf'
fig.savefig(fig_name, bbox_inches = "tight", format='pdf',pad_inches = 0.4, dpi = 600)
plt.close()
