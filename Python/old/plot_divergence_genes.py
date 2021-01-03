from __future__ import division
import os, sys, pickle, random
import numpy as np

import  matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import ColorConverter

import scipy.stats as stats
import statsmodels.api as sm

import parse_file
import timecourse_utils
import mutation_spectrum_utils
import phylo_tools as pt

#import get_random_matrix

import phik

np.random.seed(123456789)

permutations_divergence = 10000




significant_multiplicity_dict = {}
significant_n_mut_dict = {}
gene_size_dict = {}
gene_mean_size_dict = {}
for taxon in pt.taxa:
    significant_multiplicity_dict[taxon] = {}
    significant_n_mut_dict[taxon] = {}
    gene_size_dict[taxon] = {}

    gene_data = parse_file.parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    convergence_matrix = parse_file.parse_convergence_matrix(pt.get_path() + '/data/timecourse_final/' +("%s_convergence_matrix.txt" % ('0'+taxon)))
    Ltot = 0
    for gene_name in sorted(convergence_matrix.keys()):
        Lmin=0
        L = max([convergence_matrix[gene_name]['length'],Lmin])
        Ltot += L
    Lavg = Ltot*1.0/len(convergence_matrix.keys())

    gene_mean_size_dict[taxon] = Lavg

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
            gene_size_dict[taxon][items[0]] = float(items[-5])
            if items[0] not in significant_multiplicity_dict[taxon]:
                significant_multiplicity_dict[taxon][items[0]] = {}

            if items[0] not in significant_n_mut_dict[taxon]:
                significant_n_mut_dict[taxon][items[0]] = {}

            significant_multiplicity_dict[taxon][items[0]][treatment] = float(items[-2])
            significant_n_mut_dict[taxon][items[0]][treatment] = float(items[-4])
