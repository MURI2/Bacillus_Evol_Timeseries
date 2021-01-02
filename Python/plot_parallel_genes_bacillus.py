import os, copy, sys
import matplotlib.pyplot as plt
#import matplotlib as mpl
from matplotlib import colors
from matplotlib.patches import Patch

import phylo_tools as pt
import parse_file

import scipy.stats as stats
import numpy as np


taxa = [ 'B', 'S']
#treatments=pt.treatments
treatments = ['0', '1', '2']

gene_dict = {}

gene_data = parse_file.parse_gene_list('B')
gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

locus_tag_to_gene_dict = {}
for gene_name_idx, gene_name in enumerate(gene_names):
    gene = genes[gene_name_idx]
    if gene == '':
        continue
    locus_tag_to_gene_dict[gene_name] = genes[gene_name_idx]


for taxon in taxa:
    for treatment in treatments:

        genes_significant_file_path = pt.get_path() +'/data/timecourse_final/' +  ("parallel_%ss_%s.txt" % ('gene', treatment+taxon))
        output_notsignificant_file_path = pt.get_path() +'/data/timecourse_final/' +  ("parallel_not_significant_%ss_%s.txt" % ('gene', treatment+taxon))

        if os.path.exists(genes_significant_file_path) == False:
            continue

        genes_significant_file = open(genes_significant_file_path, 'r')
        genes_notsignificant_file = open(output_notsignificant_file_path, 'r')
        first_line_significant = genes_significant_file.readline()
        first_line_notsignificant = genes_notsignificant_file.readline()

        for line in genes_significant_file:
            line_split = line.strip().split(', ')
            gene_name = line_split[0]

            if gene_name not in gene_dict:
                gene_dict[gene_name] = {}

            gene_dict[gene_name][treatment+taxon] = 2

        for line in genes_notsignificant_file:
            line_split = line.strip().split(', ')
            gene_name = line_split[0]

            if gene_name not in gene_dict:
                gene_dict[gene_name] = {}

            gene_dict[gene_name][treatment+taxon] = 1


# add zero for genes that you couldn't test
for gene, gene_dict_i in gene_dict.items():
    for taxon in taxa:
        for treatment in treatments:
            if treatment + taxon not in gene_dict_i:
                gene_dict_i[treatment + taxon] = 0


# Go back through and perform test before you remove genes you don't want
for treatment in treatments:

    count_significant_B = 0
    count_significant_S = 0
    count_significant_B_S = 0

    for gene, gene_dict_i in gene_dict.items():

        if gene_dict_i[treatment+'B'] == 2:
            count_significant_B += 1

        if gene_dict_i[treatment+'S'] == 2:
            count_significant_S += 1

        if (gene_dict_i[treatment+'B'] == 2) and (gene_dict_i[treatment+'S'] == 2):
            count_significant_B_S += 1

    print(len(gene_names))

    print(count_significant_B_S, count_significant_S, count_significant_B, len(gene_names) - count_significant_S - count_significant_B - count_significant_B_S)


    oddsratio, pvalue = stats.fisher_exact([[count_significant_B_S, count_significant_S], [count_significant_B, len(gene_names) - count_significant_S - count_significant_B - count_significant_B_S]], alternative='less')

    sys.stderr.write("%s-day Bacillus WT vs delta spo0A divergence test...\n" % treatment )

    sys.stderr.write("Fisher's exact test odds-ratio = %f, p = %f\n" %  (round(oddsratio, 3),  round(pvalue, 3) ))






gene_dict_copy = copy.deepcopy(gene_dict)
# remove genes if less than two genes are significant
for gene in list(gene_dict_copy):
    if list(gene_dict_copy[gene].values()).count(2) < 2:
        del gene_dict_copy[gene]



gene_values = []
gene_names = []
for gene, gene_dict_i in gene_dict_copy.items():

    gene_values.append(list(gene_dict_i.values()))
    if gene in locus_tag_to_gene_dict:
        gene_names.append(locus_tag_to_gene_dict[gene])
    else:
        gene_names.append(gene)






# 0 = no test
# 1 = tested, P > P*
# 2 = tested, P < P*

# define the bins and normalize

#fig = plt.figure(figsize=(6,3))
fig, ax = plt.subplots(figsize=(2,7))

data = np.asarray(gene_values)


ax.xaxis.tick_top()
ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

ax.set_yticklabels(gene_names, minor=False, fontsize=5)
ax.set_xticklabels([ r'$\mathrm{WT}$', r'$\Delta spo0A$']*3, minor=False, fontsize=5.5)



#import matplotlib.transforms as mtrans
#from matplotlib.text import TextPath
#from matplotlib.patches import PathPatch

#def curly(x,y, scale, ax=None):
#    if not ax: ax=plt.gca()
#    tp = TextPath((0, 0), "}", size=1)
#    trans = mtrans.Affine2D().scale(1, scale) + \
#        mtrans.Affine2D().translate(x,y) + ax.transData + mpl.transforms.Affine2D().rotate_deg(90)
#    pp = PathPatch(tp, lw=0, fc="k", transform=trans)
#    ax.add_artist(pp)

ax.text(0.07 , 1.05, "1-day", fontsize=6.5, transform=ax.transAxes, weight="bold")
ax.text(0.36, 1.05, "10-days", fontsize=6.5, transform=ax.transAxes, weight="bold")
ax.text(0.67, 1.05, "100-days", fontsize=6.5, transform=ax.transAxes, weight="bold")

#curly(0.5 ,40,3, ax=ax)


legend_elements = [Patch(color='lightgrey', label=r'$n_{mut} < 3$'),
                    Patch(color='orangered', label=r'$P\nless P^{*}$'),
                    Patch(color='deepskyblue', label=r'$P < P^{*}$')]


# Create the figure
#plt.legend(handles=legend_elements, bbox_to_anchor=(0., 1.02, 1., .102),mode="expand", ncol=3, loc="upper left", fontsize=8)
plt.legend(handles=legend_elements, bbox_to_anchor=(-0.65,1.12), loc="upper left", fontsize=8)

#bbox_to_anchor=(1,1.04)


cols = {0:'lightgrey',1:'orangered',2:'deepskyblue'}
cvr = colors.ColorConverter()
tmp = sorted(cols.keys())
cols_rgb = [cvr.to_rgb(cols[k]) for k in tmp]
intervals = np.asarray(tmp + [tmp[-1]+1]) - 0.5
cmap, norm = colors.from_levels_and_colors(intervals,cols_rgb)

plt.pcolor(data,cmap = cmap, norm = norm, edgecolors='k', linewidths=1.5)

fig.savefig(pt.get_path() + '/figs/genes_table.pdf', format= 'pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
