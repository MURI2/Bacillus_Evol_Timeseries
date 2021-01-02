from __future__ import division
import os
import pandas as pd
import  matplotlib.pyplot as plt
from Bio import SeqIO

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/")


def get_spore_gene_names():
    spore_gene_list = []
    IN_path = mydir + 'data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff'
    genome = SeqIO.parse(IN_path, "genbank")
    count = 0
    for record in genome:
        for f in record.features:
            if 'product' in f.qualifiers:
                if ('spore' in  f.qualifiers['product'][0]) or ('sporulation' in  f.qualifiers['product'][0]):
                    #print f.qualifiers['product']
                    count += 1
                    spore_gene_list.append(f.qualifiers['locus_tag'][0])
                #print f.qualifiers['product']
    B_path = mydir + 'data/gene_by_sample/B/sample_by_gene_Gscore.txt'
    S_path = mydir + 'data/gene_by_sample/S/sample_by_gene_Gscore.txt'
    B = pd.read_csv(B_path, sep = '\t', header = 'infer', index_col = 0)
    S = pd.read_csv(S_path, sep = '\t', header = 'infer', index_col = 0)
    B_S_merged = B.append(S)
    B_S_merged = B_S_merged.fillna(0.0)
    # only get day 100
    B_S_merged = B_S_merged[['_D100' in s for s in B_S_merged.index]]
    #df[df['A'].str.contains("hello")]
    #OUTname = mydir + 'gene_by_sample/B_S/sample_by_gene.txt'
    #B_S_merged.to_csv(OUTname, sep = '\t', index = True)
    #print list(B_S_merged.columns.values)
    spore_inter = list(set(list(B_S_merged.columns.values)).intersection(set(spore_gene_list)))
    not_spore_inter = list(set(list(B_S_merged.columns.values)).difference(set(spore_gene_list)))
    #spore_inter = list(set(spore_gene_list).intersection(set(list(B_S_merged.columns.values))))
    #not_spore_inter = list(set(spore_gene_list).difference(set(list(B_S_merged.columns.values))))
    # get just the columns
    B_S_merged_spore = B_S_merged[spore_inter]
    B_S_merged_not_spore = B_S_merged[not_spore_inter]
    #print not_spore_inter
    #print spore_inter
    spore_mean = B_S_merged_spore.mean(axis = 1).values
    not_spore_mean = B_S_merged_not_spore.mean(axis = 1).values
    diff_spore = spore_mean - not_spore_mean
    #print B_S_merged_spore.index.values




    # now make plot.


get_spore_gene_names()
