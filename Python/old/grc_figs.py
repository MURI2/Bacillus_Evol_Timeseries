from __future__ import division
import os
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy import stats
import phylo_tools as pt
import collections
from sklearn.decomposition import PCA
from Bio import SeqIO
import scipy.stats as stats
from collections import Counter


def fixed_fig():
    strains = ['B', 'S', 'C', 'D']
    reps = ['1', '2', '3', '4', '5']
    treatments = ['0', '1', '2']
    count=0
    fig = plt.figure()
    fig.subplots_adjust(bottom= 0.15)

    for i, strain in enumerate(strains):
        x_list = []
        y_list = []
        for treatment in treatments:
            for rep in reps:
                muts_filename = '%s%s%s_snp_final.txt' % (treatment, strain, rep)
                muts = pd.read_csv(pt.get_path() + '/data/timecourse_prelim_fixed/' + muts_filename, sep = '\t', header = 'infer')
                x_list.append(10 ** int(treatment))
                y_list.append(int(muts.shape[0]) +1)

        ax = fig.add_subplot(2, 2, count+1)
        #plt.scatter(x_list, y_list)

        x_y_zip = list(zip(x_list, y_list))
        x_y_0 = [x for x in x_y_zip if x[0] == 1]
        x_y_1 = [x for x in x_y_zip if x[0] == 10]
        x_y_2 = [x for x in x_y_zip if x[0] == 100]
        print(x_y_2)


        x_0 = [x[0] for x in x_y_0]
        y_0 = [x[1] for x in x_y_0]
        x_1 = [x[0] for x in x_y_1]
        y_1 = [x[1] for x in x_y_1]
        x_2 = [x[0] for x in x_y_2]
        y_2 = [x[1] for x in x_y_2]


        if i == 0:
            ax.set_title(r'$B. \, subtilis$')
            plt.scatter(x_0, y_0, marker = "o", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.9, s = 40)
            plt.scatter(x_1, y_1, marker = "o", edgecolors='#FFA500', c = '#FFA500', alpha = 0.9, s = 40)
            plt.scatter(x_2, y_2, marker = "o", edgecolors='#FF6347', c = '#FF6347', alpha = 0.9, s = 40)
        elif i == 1:
            ax.set_title(r'$B. \, subtilis\; \Delta spo0A$')
            plt.scatter(x_0, y_0, marker = "o", edgecolors='#87CEEB', alpha = 0.9, s = 40, linewidth=3, facecolors='none')
            plt.scatter(x_1, y_1, marker = "o", edgecolors='#FFA500', alpha = 0.9, s = 40, linewidth=3, facecolors='none')
            plt.scatter(x_2, y_2, marker = "o", edgecolors='#FF6347', alpha = 0.9, s = 40, linewidth=3, facecolors='none')

        elif i == 2:
            ax.set_title(r'$C. \, crescentus$')
            plt.scatter(x_0, y_0, marker = "P", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.9, s = 40)
            plt.scatter(x_1, y_1, marker = "P", edgecolors='#FFA500', c = '#FFA500', alpha = 0.9, s = 40)
            plt.scatter(x_2, y_2, marker = "P", edgecolors='#FF6347', c = '#FF6347', alpha = 0.9, s = 40)
        elif i == 3:
            ax.set_title(r'$D. \, radiodurans$')
            plt.scatter(x_0, y_0, marker = "^", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.9, s = 40)
            plt.scatter(x_1, y_1, marker = "^", edgecolors='#FFA500', c = '#FFA500', alpha = 0.9, s = 40)
            plt.scatter(x_2, y_2, marker = "^", edgecolors='#FF6347', c = '#FF6347', alpha = 0.9, s = 40)


        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)
        ax.set_xlim([0.6,150])
        ax.set_ylim([0.6,450])
        plt.subplots_adjust(wspace=0.2, hspace=0.32)

        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_list), np.log10(y_list))
        x_log10_fit_range =  np.linspace(np.log10(0.6), np.log10(150), 10000)
        y_fit_range = 10 ** (slope*x_log10_fit_range + intercept)
        ax.plot(10**x_log10_fit_range, y_fit_range, c='k')

        y_null_range = 10 ** (-1*x_log10_fit_range + intercept)
        ax.plot(10**x_log10_fit_range, y_null_range, c='darkgrey', linestyle='--')

        plt.text(30, 270, r'$\beta_{1}=$' + str(round(slope, 2)), fontsize=8)
        plt.text(30, 140, r'$r^{2}=$' + str(round(r_value**2, 2)), fontsize=8)
        # hypothetical slope of -1
        ratio = (slope - (-1)) / std_err
        pval = stats.t.sf(np.abs(ratio), len(y_list)-2)*2
        # two sided or one sided?
        print(strain, ratio, pval)

        if pval < 0.05:
            plt.text(30, 80, r'$\mathrm{p} < 0.05$', fontsize=8)
        else:
            plt.text(30, 80, r'$\mathrm{p} \nless 0.05$', fontsize=8)

        count += 1

    fig.text(0.50, 0.07, 'Transfer time (days)', ha='center', va='center', fontsize=14)
    fig.text(0.04, 0.5, 'Fixed mutations at day 1,000', ha='center', va='center', rotation='vertical', fontsize=14)
    fig_name = str(pt.get_path() + '/figs/fixed_fig.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()


def PCA_B():
    strains = ['B', 'S']
    reps = ['1', '2', '3', '4', '5']
    treatments = ['0', '1']
    sites = []
    count_dict = {}
    for i, strain in enumerate(strains):
        for treatment in treatments:
            for rep in reps:
                muts_filename = '%s%s%s_snp_final.txt' % (treatment, strain, rep)
                sample = '%s%s%s' % (treatment, strain, rep)
                count_dict[sample] = {}
                muts = pd.read_csv(pt.get_path() + '/data/timecourse_prelim_fixed/' + muts_filename, sep = '\t', header = 'infer')
                for index, row in muts.iterrows():
                    gene_name = row['Locus_tag']
                    if isinstance(gene_name, float):
                        continue
                    #if ',' in gene_name:
                    gene_names = gene_name.split(',')
                    for x in gene_names:
                        if x in count_dict[sample]:
                            count_dict[sample][x] += 1
                        else:
                            count_dict[sample][x] = 1
    count_df = pd.DataFrame.from_dict(count_dict, orient='index')
    count_df.fillna(0, inplace=True)
    count_df = (count_df.T / count_df.T.sum()).T

    pca = PCA()
    pca_out = pca.fit_transform(count_df)
    pca_df = pd.DataFrame(pca_out, index=count_df.index.values)
    pca_0B = pca_df[pca_df.index.str.contains('0B')]
    pca_1B = pca_df[pca_df.index.str.contains('1B')]
    pca_0S = pca_df[pca_df.index.str.contains('0S')]
    pca_1S = pca_df[pca_df.index.str.contains('1S')]

    eigenvalues = pca.explained_variance_

    fig = plt.figure()

    plt.scatter(pca_0B.values[:,0], pca_0B.values[:,1], marker = "o", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.8, s = 120, zorder=4)
    plt.scatter(pca_1B.values[:,0], pca_1B.values[:,1], marker = "o", edgecolors='#FFA500', c = '#FFA500', alpha = 0.8, s = 120, zorder=4)
    plt.scatter(pca_0S.values[:,0], pca_0S.values[:,1], marker = "o", edgecolors='#87CEEB', linewidth=3, facecolors='none', alpha = 0.8, s = 120, zorder=4)
    plt.scatter(pca_1S.values[:,0], pca_1S.values[:,1], marker = "o", edgecolors='#FFA500', linewidth=3, facecolors='none', alpha = 0.8, s = 120, zorder=4)

    plt.xlabel('PC 1 (' + str(round(pca.explained_variance_ratio_[0],3)*100) + '%)' , fontsize = 14)
    plt.ylabel('PC 2 (' + str(round(pca.explained_variance_ratio_[1],3)*100) + '%)' , fontsize = 14)



    plt.xlim(-0.1,0.075)
    plt.ylim(-0.08,0.01)
    fig_name = str(pt.get_path() + '/figs/PCA_B.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()




def module_pca():
    strains = ['B', 'D', 'S', 'C']
    reps = ['1', '2', '3', '4', '5']
    treatments = ['0', '1']
    sites = []
    module_count_dict = {}
    kegg_count_dict = {}

    modules_list = []
    # create kegg to module dict for all taxa before moving on
    for strain in strains:
        if strain == 'S':
            kegg_path = pt.get_path() + '/data/reference_assemblies_task2/MAPLE/MAPLE_modules/' + 'B' + '_KO_to_M.txt'
        else:
            kegg_path = pt.get_path() + '/data/reference_assemblies_task2/MAPLE/MAPLE_modules/' + strain + '_KO_to_M.txt'
        kegg_df = pd.read_csv(kegg_path, sep = '\t', header = 'infer')
        modules_list.extend(list(set(kegg_df.Pathway_ID.to_list())))

    modules_list_all = [k for k, v in Counter(modules_list).items() if v == 4]

    for strain in strains:
        kegg_module_dict = {}
        if strain == 'S':
            kegg_module_path = pt.get_path() + '/data/reference_assemblies_task2/MAPLE/MAPLE_modules/' + 'B' + '_KO_to_M.txt'
        else:
            kegg_module_path = pt.get_path() + '/data/reference_assemblies_task2/MAPLE/MAPLE_modules/' + strain + '_KO_to_M.txt'
        kegg_module_df = pd.read_csv(kegg_module_path, sep = '\t', header = 'infer')
        for index, row in kegg_module_df.iterrows():
            kegg_module_dict[row.KEGG_Orthology] = row.Pathway_ID

        protein_id_kegg_dict = {}
        locustag_protein_id_dict = {}
        if strain == 'S':
            query_path = pt.get_path() + '/data/reference_assemblies_task2/MAPLE/' + 'B' + '_MAPLE_result/query.fst.ko'
        else:
            query_path = pt.get_path() + '/data/reference_assemblies_task2/MAPLE/' + strain + '_MAPLE_result/query.fst.ko'
        query_df = pd.read_csv(query_path, sep = '\t', header = None)
        query_df.columns = ['Locus_tag', 'KEGG_id', 'strain', 'species', 'some_number']
        for index, row in query_df.iterrows():
            if row['KEGG_id'] == 'K_NA':
                continue
            kegg_ids = row['KEGG_id'].split(',')
            for kegg_id in kegg_ids:
                protein_id_kegg_dict[row['Locus_tag']] = kegg_id

        filename= pt.get_path() + '/' + pt.get_ref_gbff_dict()[strain]
        recs = [rec for rec in SeqIO.parse(filename, "genbank")]
        for rec in recs:
            for feat in rec.features:
                if 'pseudo' in list((feat.qualifiers.keys())):
                    continue
                if (feat.type == "source") or (feat.type == "gene"):
                    continue
                if 'protein_id' in feat.qualifiers:
                    locustag_protein_id_dict[feat.qualifiers['locus_tag'][0]] = feat.qualifiers['protein_id'][0]


        for treatment in treatments:
            for rep in reps:
                muts_filename = '%s%s%s_snp_final.txt' % (treatment, strain, rep)
                sample = '%s%s%s' % (treatment, strain, rep)
                print(sample)
                kegg_count_dict[sample] = {}
                module_count_dict[sample] = {}
                muts = pd.read_csv(pt.get_path() + '/data/timecourse_prelim_poly/' + muts_filename, sep = '\t', header = 'infer')

                for index, row in muts.iterrows():
                    gene_name = row['Locus_tag']
                    if isinstance(gene_name, float):
                        continue
                    gene_names = gene_name.split(',')
                    for gene_name in gene_names:
                        # get protein id
                        if gene_name in locustag_protein_id_dict:
                             protein_id = locustag_protein_id_dict[gene_name]

                             if protein_id in protein_id_kegg_dict:
                                kegg_id = protein_id_kegg_dict[protein_id]
                                if kegg_id in kegg_count_dict[sample]:
                                    kegg_count_dict[sample][kegg_id] += 1
                                else:
                                    kegg_count_dict[sample][kegg_id] = 1

                                if kegg_id in kegg_module_dict:
                                    module_id = kegg_module_dict[kegg_id]
                                    if module_id in module_count_dict:
                                        module_count_dict[sample][module_id] += 1

                                    else:
                                        module_count_dict[sample][module_id] = 1


    count_df = pd.DataFrame.from_dict(module_count_dict, orient='index')
    count_df.fillna(0, inplace=True)
    count_df = count_df[(count_df.T != 0).any()]

    count_df = (count_df.T / count_df.T.sum()).T
    count_df = count_df.drop('1B5')
    modules_list_all_intersect = list(set(modules_list_all).intersection(set(count_df.columns.to_list())))
    count_df = count_df[modules_list_all_intersect]
    #modules_list_all


    pca = PCA()
    pca_out = pca.fit_transform(count_df)
    pca_df = pd.DataFrame(pca_out, index=count_df.index.values)

    pca_0B = pca_df[pca_df.index.str.contains('0B')]
    pca_1B = pca_df[pca_df.index.str.contains('1B')]

    pca_0S = pca_df[pca_df.index.str.contains('0S')]
    pca_1S = pca_df[pca_df.index.str.contains('1S')]

    pca_0C = pca_df[pca_df.index.str.contains('0C')]
    pca_1C = pca_df[pca_df.index.str.contains('1C')]

    pca_0D = pca_df[pca_df.index.str.contains('0D')]
    pca_1D = pca_df[pca_df.index.str.contains('1D')]

    eigenvalues = pca.explained_variance_

    fig = plt.figure()

    plt.scatter(pca_0B.values[:,0], pca_0B.values[:,1], marker = "o", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.8, s = 120, zorder=4, label=r'$B. \, subtilis \; \mathrm{1-Day}$')
    plt.scatter(pca_1B.values[:,0], pca_1B.values[:,1], marker = "o", edgecolors='#FFA500', c = '#FFA500', alpha = 0.8, s = 120, zorder=4, label=r'$B. \, subtilis \; \mathrm{10-Day}$')
    plt.scatter(pca_0D.values[:,0], pca_0D.values[:,1], marker = "^", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.8, s = 120, zorder=4, label=r'$D. \, radiodurans \; \mathrm{1-Day}$')
    plt.scatter(pca_1D.values[:,0], pca_1D.values[:,1], marker = "^", edgecolors='#FFA500', c = '#FFA500', alpha = 0.8, s = 120, zorder=4, label=r'$D. \, radiodurans \; \mathrm{10-Day}$')
    plt.scatter(pca_0S.values[:,0], pca_0S.values[:,1], marker = "o", edgecolors='#87CEEB', linewidth=3, facecolors='none', alpha = 0.8, s = 120, zorder=4, label=r'$B. \, subtilis \; \mathrm{\Delta spo0A \; 1-Day}$')
    plt.scatter(pca_1S.values[:,0], pca_1S.values[:,1], marker = "o", edgecolors='#FFA500', linewidth=3, facecolors='none', alpha = 0.8, s = 120, zorder=4, label=r'$B. \, subtilis \; \mathrm{\Delta spo0A \; 10-Day}$')
    plt.scatter(pca_0C.values[:,0], pca_0C.values[:,1], marker = "P", edgecolors='#87CEEB', c = '#87CEEB', alpha = 0.8, s = 120, zorder=4, label=r'$C. \, crescentus \; \mathrm{1-Day}$')
    plt.scatter(pca_1C.values[:,0], pca_1C.values[:,1], marker = "P", edgecolors='#FFA500', c = '#FFA500', alpha = 0.8, s = 120, zorder=4, label=r'$C. \, crescentus \; \mathrm{10-Day}$')

    plt.xlabel('PC 1 (' + str(round(pca.explained_variance_ratio_[0]*100,2)) + '%)' , fontsize = 14)
    plt.ylabel('PC 2 (' + str(round(pca.explained_variance_ratio_[1]*100,2)) + '%)' , fontsize = 14)

    lgnd=plt.legend(loc='upper right', prop={'size': 8})
    #change the marker size manually for both lines
    lgnd.legendHandles[0]._sizes = [30]
    lgnd.legendHandles[1]._sizes = [30]
    lgnd.legendHandles[2]._sizes = [30]
    lgnd.legendHandles[3]._sizes = [30]
    lgnd.legendHandles[4]._sizes = [30]
    lgnd.legendHandles[5]._sizes = [30]
    lgnd.legendHandles[6]._sizes = [30]
    lgnd.legendHandles[7]._sizes = [30]


    #plt.xlim(-0.08,0.08)
    #plt.ylim(-0.08,0.08)
    fig_name = str(pt.get_path() + '/figs/PCA_kegg.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()





module_pca()
#fixed_fig()
#PCA_B()
