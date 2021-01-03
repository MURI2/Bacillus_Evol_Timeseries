from __future__ import division
import os
from collections import Counter

def get_path():
    return os.path.expanduser("~/GitHub/Bacillus_Evol_Timeseries")

def get_to_keep():
    return ['SNP', 'INS', 'DEL']

def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def get_bacillus_mut_bias():
    return  1.2739

def get_bacillus_mut_rate():
    return 3.35 * (10**-10)

def get_bacillus_indel_rate():
    return 1.20 * (10**-10)

def get_colors():
    return {'0':'#87CEEB', '1': '#FFA500', '2':'#FF6347'}


def get_ts_tv_dict():

    ts_tv_dict = {
    ("A", "C"):"V", ("A", "G"):"S", ("A", "T"):"V",
    ("C", "A"):"V", ("C", "G"):"V", ("C", "T"):"S",
    ("G", "A"):"S", ("G", "C"):"V", ("G", "T"):"V",
    ("T", "A"):"V", ("T", "C"):"S", ("T", "G"):"V"}

    return ts_tv_dict

def get_codon_dict():
    # translation table 11
    codon_dict = {
        "TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
        "TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
        "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
        "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",

        "CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
        "CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
        "CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
        "CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",

        "ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
        "ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
        "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
        "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",

        "GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
        "GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
        "GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
        "GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
        }

    return codon_dict

def mutations_to_exclude():
    directory = os.fsencode(get_path() + '/data/pool_pop_seq/rebreseq_annotated')
    list_muts = []
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('-100.gd'):
            in_df = open(os.path.join(str(directory, 'utf-8'), filename), 'r')
            for line in in_df:
                line_split = line.strip().split()
                if line_split[0] not in get_to_keep():
                    continue
                freq = [s for s in line_split if 'frequency=' in s][0].split('=')[1]
                if float(freq) == float(1):
                    list_muts.append(line_split[3] + '_' + line_split[4])
    dict_muts = Counter(list_muts)
    return [i for i in dict_muts if dict_muts[i] >= 10]
