from __future__ import division
import os, sys, itertools
import parse_file
import phylo_tools as pt
import timecourse_utils
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

np.random.seed(123456789)

iter=10000

iter_corr= 10000

spectra_latex_dict = {'AT_GC':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{G}:\mathrm{C}$',
                    'AT_CG':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{C}:\mathrm{G}$',
                    'AT_TA':r'$\mathrm{A}:\mathrm{T} \rightarrow  \mathrm{T}:\mathrm{A}$',
                    'GC_AT':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{A}:\mathrm{T}$',
                    'GC_TA':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{T}:\mathrm{A}$',
                    'GC_CG':r'$\mathrm{G}:\mathrm{C} \rightarrow  \mathrm{C}:\mathrm{G}$'}


sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']
all_subplot_counts = 0

treatments=pt.treatments
replicates = pt.replicates
taxa = ['B','C','D','F','J','P']
iter
spectra_dict = {}
snps_indel_dict = {}
for taxon in taxa:

    for treatment in treatments:

        for replicate in replicates:

            population = treatment + taxon + replicate
            if population in pt.populations_to_ignore:
                continue

            mutation_counts = {'GC_AT': 1,
                            'GC_TA': 1,
                            'GC_CG': 1,
                            'AT_GC': 1,
                            'AT_CG': 1,
                            'AT_TA': 1}

            mutations, depth_tuple = parse_file.parse_annotated_timecourse(population)
            population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
            state_times, state_trajectories = parse_file.parse_well_mixed_state_timecourse(population)

            for mutation_idx in range(0,len(mutations)):

                location, gene_name, allele, var_type, codon, position_in_codon, AAs_count,  test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
                state_Ls = state_trajectories[mutation_idx]

                if AAs_count == 1:

                    ref, alt = allele.split('->')

                    if (ref=='A' and alt=='T') or (ref=='T' and alt=='A'):
                        mutation_counts['AT_TA']+=1

                    elif (ref=='A' and alt=='G') or (ref=='T' and alt=='C'):
                        mutation_counts['AT_GC']+=1

                    elif (ref=='A' and alt=='C') or (ref=='T' and alt=='G'):
                        mutation_counts['AT_CG']+=1

                    elif (ref=='G' and alt=='A') or (ref=='C' and alt=='T'):
                        mutation_counts['GC_AT']+=1

                    elif (ref=='G' and alt=='T') or (ref=='C' and alt=='A'):
                        mutation_counts['GC_TA']+=1

                    elif (ref=='G' and alt=='C') or (ref=='C' and alt=='G'):
                        mutation_counts['GC_CG']+=1

                    else:

                        print(ref, alt)

            mutation_counts_rel = {k: v / sum(mutation_counts.values()) for k, v in mutation_counts.items()}

            spectra_dict[population] = mutation_counts_rel


# get indel data now






df = pd.DataFrame.from_dict(spectra_dict)
X = df.values.transpose()

X_std=StandardScaler().fit_transform(X) # normalizing the features, subtract mean divide by SD

pca_ = PCA(n_components=2)
#principalComponents_all_ = pca_.fit_transform(X_std)
#principalComponents_ = principalComponents_all_[:,:2]
principalComponents_ = pca_.fit_transform(X_std)

principalComponents_df = pd.DataFrame(principalComponents_, index=df.columns, columns = ['PC1', 'PC2'])

pc_df_reinxed = []
treatment_counts = []
for treatment in pt.treatments:
    pc_df_treatment = [x for x in principalComponents_df.index.to_list() if x[0]==treatment]
    pc_df_reinxed.extend(pc_df_treatment)

    treatment_counts.append(len(pc_df_treatment))

principalComponents_df = principalComponents_df.reindex(pc_df_reinxed)
principalComponents_np = principalComponents_df.values

# get dictionary of numpy arrays
taxon_pc_dict = {}
taxon_pc_labels = []
for taxon in taxa:
    taxon_index_names = [x for x in principalComponents_df.index.to_list() if x[1]==taxon]
    taxon_pc_df = principalComponents_df.loc[taxon_index_names]
    taxon_pc_dict[taxon] = taxon_pc_df.values
    taxon_pc_labels.extend(taxon_index_names)


treatment_counts = np.asarray(treatment_counts)

F = pt.get_F_2(principalComponents_np, treatment_counts)

# permute within each species
F_nested_list = []
for i in range(iter):
    #principalComponents_np_copy = np.copy(principalComponents_np)
    taxon_pc_dict_copy = taxon_pc_dict.copy()
    taxon_pc_permute_list = []
    for taxon in taxa:
        taxon_pc_dict_copy_i = taxon_pc_dict_copy[taxon]
        np.random.shuffle(taxon_pc_dict_copy_i)
        taxon_pc_permute_list.append(taxon_pc_dict_copy_i)

    principalComponents_np_permute = np.concatenate(taxon_pc_permute_list, axis=0)
    principalComponents_df_permute = pd.DataFrame(principalComponents_np_permute, index=taxon_pc_labels, columns = ['PC1', 'PC2'])
    # now resort by treatment
    pc_df_reinxed_permute = []
    for treatment in pt.treatments:
        pc_df_permute_treatment = [x for x in principalComponents_df_permute.index.to_list() if x[0]==treatment]
        pc_df_reinxed_permute.extend(pc_df_permute_treatment)

    principalComponents_df_permute = principalComponents_df_permute.reindex(pc_df_reinxed_permute)

    F_i = pt.get_F_2(principalComponents_df_permute.values, treatment_counts)
    F_nested_list.append(F_i)

p_value = (len([f for f in F_nested_list if f > F]) +1)/(iter+1)

sys.stderr.write("PERMANOVA F = %f, P = %f\n" % (round(F, 3), round(p_value, 6)))



# examine factor loadings
# variables were standardized, so we don't need to divive loadings by S.D. of variables
# if the variables weren't standardized, we'd end up w/ something like V * sqrt(E)/ std(X)
loadings = pca_.components_.T * np.sqrt(pca_.explained_variance_)

loading_squared_matrix = pd.DataFrame(loadings**2, columns=['PC1', 'PC2'], index=df.index)

loading_1_name=loading_squared_matrix.loc[loading_squared_matrix['PC1'].idxmax()].name
loading_1_rho_2=loading_squared_matrix.loc[loading_squared_matrix['PC1'].idxmax()][0]
loading_1_rho_2_idx = df.index.to_list().index(loading_1_name)

loading_2_name=loading_squared_matrix.loc[loading_squared_matrix['PC2'].idxmax()].name
loading_2_rho_2=loading_squared_matrix.loc[loading_squared_matrix['PC2'].idxmax()][1]
loading_2_rho_2_idx = df.index.to_list().index(loading_2_name)

transition_names = df.index.to_list()
# run permutation tests
permutation_dict = {}
for transition in transition_names:
    permutation_dict[transition] = {}
    permutation_dict[transition]['PC1'] = []
    permutation_dict[transition]['PC2'] = []

for i in range(iter_corr):
    X_copy = X.copy()
    #np.random.shuffle(pca_components_T_copy)
    for i in range(X_copy.shape[1]):
        np.random.shuffle(X_copy[:,i])

    X_copy_std=StandardScaler().fit_transform(X_copy) # normalizing the features, subtract mean divide by SD

    pca_copy_ = PCA(n_components=2)
    principalComponents_copy_ = pca_copy_.fit_transform(X_copy_std)
    loadings_copy_squared = (pca_copy_.components_.T * np.sqrt(pca_copy_.explained_variance_)) **2
    for loadings_row_idx, loadings_row in enumerate(loadings_copy_squared):

        permutation_dict[transition_names[loadings_row_idx]]['PC1'].append(loadings_row[0])
        permutation_dict[transition_names[loadings_row_idx]]['PC2'].append(loadings_row[1])


# get significant loadings
for index, row in loading_squared_matrix.iterrows():
    P_PC1 = (len([k for k in permutation_dict[transition]['PC1'] if k > float(row['PC1']) ])+1)/ (iter_corr+1)
    P_PC2 = (len([k for k in permutation_dict[transition]['PC2'] if k > float(row['PC2']) ])+1)/ (iter_corr+1)

    if P_PC1 < 0.05:
        sys.stderr.write("PC1: %s, rho^2 = %f, P=%f\n" % (index, round(row['PC1'], 3), round(P_PC1, 3)))

    if P_PC2 < 0.05:
        sys.stderr.write("PC2: %s, rho^2 = %f, P=%f\n" % (index, round(row['PC2'], 3), round(P_PC2, 3)))



# perform ANOVA and spectra correlated with top 2 PCs




fig = plt.figure(figsize = (12, 9))
gs = gridspec.GridSpec(nrows=6, ncols=4)


mutation_spectra_list = [['AT_GC','AT_CG','AT_TA'],['GC_AT','GC_TA','GC_CG']]


for taxon_list_idx, taxon_list in enumerate([['B','C','D'],['F','J','P']]):

    for taxon_idx, taxon in enumerate(taxon_list):

        if taxon == 'J':
            treatments = ['0','2']
        else:
            treatments = pt.treatments

        #set_time = set_time_dict[taxon]

        ax = fig.add_subplot(gs[taxon_idx*2:(taxon_idx*2)+2, taxon_list_idx])
        ax.set_title(pt.latex_genus_bold_dict[taxon], fontsize=12, fontweight='bold' )

        for treatment in treatments:

            PCs_ = principalComponents_df[principalComponents_df.index.str.contains(treatment+taxon)]

            ax.axhline(y=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
            ax.axvline(x=0, color='k', linestyle=':', alpha = 0.8, zorder=1)
            ax.scatter(0, 0, marker = "o", edgecolors='none', c = 'darkgray', s = 120, zorder=2)

            ax.scatter(PCs_.PC1.values, PCs_.PC2.values, \
                    c=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), s = 70, \
                    edgecolors=pt.get_colors(treatment), linewidth = 0.6, alpha = 0.8, zorder=4)#, edgecolors='none'

            pt.confidence_ellipse(PCs_.PC1.values, PCs_.PC2.values, ax,
                n_std=2, edgecolor=pt.get_colors(treatment), linestyle='--', lw=4, zorder=3)


        ax.text(-0.1, 1.07, sub_plot_labels[all_subplot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax.transAxes)

        all_subplot_counts += 1



count = 0
for snp_i in range(2):
    ax_snp = fig.add_subplot(gs[snp_i*3:(snp_i+1)*3, 2:])
    ax_snp.axhline(y=0, color='k', linestyle=':', lw=2, alpha = 0.8, zorder=1)
    ax_snp.set_xticklabels([])
    ax_snp.set_xticks([])
    for spectrum_idx, spectrum in enumerate(mutation_spectra_list[snp_i]):
        for taxon in ['B','C','D','J']:
            for treatment in pt.treatments:
                spectrum_reps = []
                for replicate in replicates:
                    population = treatment+taxon+replicate
                    if population in spectra_dict:
                        spectrum_reps.append(spectra_dict[population][spectrum])

                if len(spectrum_reps) < 2:
                    continue

                spectrum_reps = np.asarray(spectrum_reps) - pt.get_ma_mutation_spectrum(taxon)[spectrum]

                ax_snp.scatter(list(itertools.repeat(count, len(spectrum_reps))), spectrum_reps,\
                        c=pt.get_colors(treatment), marker=pt.plot_species_marker(taxon), s = 40, \
                        edgecolors=pt.get_colors(treatment), linewidth = 0.6, alpha = 0.8, zorder=2)#, edgecolors='none'

                count+=1

        ax_snp.text((spectrum_idx+0.5)*(1/3), -0.05, spectra_latex_dict[spectrum], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_snp.transAxes)
        if spectrum_idx < 2:

            ax_snp.axvline(x=count-0.5, color='k', linestyle='--', alpha = 0.8, zorder=1)


    ax_snp.text(0.03, 1.07, sub_plot_labels[all_subplot_counts], fontsize=12, fontweight='bold', ha='center', va='center', transform=ax_snp.transAxes)

    all_subplot_counts += 1




fig.text(0.25, -0.01, 'PC 1 (' + str(round(pca_.explained_variance_ratio_[0]*100,2)) + '%)', ha='center', va='center', fontsize=18)
fig.text(-0.01, 0.5, 'PC 2 (' + str(round(pca_.explained_variance_ratio_[1]*100,2)) + '%)', ha='center', va='center', rotation='vertical', fontsize=18)

fig.text(0.51, 0.5, 'Observed spectra - MA spectra', ha='center', va='center', rotation='vertical', fontsize=18)


fig.subplots_adjust(hspace=0.4, wspace=0.45) #hspace=0.3, wspace=0.5
fig.tight_layout()
fig.savefig(pt.get_path() + '/figs/mutation_spectra_pca.pdf', format='pdf', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()



#    #ma_mutation_spectrum = pt.get_ma_mutation_spectrum(taxon)
