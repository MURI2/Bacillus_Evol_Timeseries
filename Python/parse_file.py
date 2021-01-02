import numpy
import sys
from math import fabs
import glob, os, sys, re
from bz2 import BZ2File
from Bio import SeqIO

import phylo_tools as pt
import timecourse_utils
#data_directory='data_files/'
#figure_directory='manuscript/figures/'

default_min_depth=5

#Y = pYrimidines
#R = puRines
#S = strong ineractions, C or G
#W = weak interactions, A or T
#K = Ketones, G or T
#M = aMino groups, C or A
bases_to_skip = ['K', 'S', 'R', 'N', 'Y', 'M', 'W']
base_table = {'A':'T','T':'A','G':'C','C':'G',
            'Y':'R', 'R':'Y', 'S':'W', 'W':'S', 'K':'M', 'M':'K', 'N':'N'}

codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA':'R',
'CGG':'R', 'AGA':'R', 'AGG':'R', 'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D',
'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H',
'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L',
'CTG':'L', 'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F', 'CCT':'P', 'CCC':'P', 'CCA':'P',
'CCG':'P', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ACT':'T', 'ACC':'T',
'ACA':'T', 'ACG':'T', 'TGG':'W', 'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
'TAA':'!', 'TGA':'!', 'TAG':'!'}#, 'KTC':'F', 'KAC':'Y', 'KCC':'A', 'KGC':'D'}

# calculate number of synonymous opportunities for each codon
codon_synonymous_opportunity_table = {}
for codon in codon_table.keys():
    codon_synonymous_opportunity_table[codon] = {}
    for i in range(0,3):
        codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
        codon_list = list(codon)
        for base in ['A','C','T','G']:
            codon_list[i]=base
            new_codon = "".join(codon_list)
            if 'K' in new_codon:
                continue
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_opportunity_table[codon][i]+=1

bases = set(['A','C','T','G'])
substitutions = []
for b1 in bases:
    for b2 in bases:
        if b2==b1:
            continue

        substitutions.append( '%s->%s' % (b1,b2) )

codon_synonymous_substitution_table = {}
codon_nonsynonymous_substitution_table = {}
for codon in codon_table.keys():
    codon_synonymous_substitution_table[codon] = [[],[],[]]
    codon_nonsynonymous_substitution_table[codon] = [[],[],[]]

    for i in range(0,3):
        reference_base = codon[i]

        codon_list = list(codon)
        for derived_base in ['A','C','T','G']:
            if derived_base==reference_base:
                continue
            substitution = '%s->%s' % (reference_base, derived_base)
            codon_list[i]=derived_base
            new_codon = "".join(codon_list)
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_substitution_table[codon][i].append(substitution)
            else:
                codon_nonsynonymous_substitution_table[codon][i].append(substitution)



def get_mutation_fixation_trajectories(population):

    mutations, depth_tuple = parse_annotated_timecourse(population)
    population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths = depth_tuple
    state_times, state_trajectories = parse_well_mixed_state_timecourse(population)
    times = mutations[0][9]
    Ms = numpy.zeros_like(times)*1.0
    fixed_Ms = numpy.zeros_like(times)*1.0

    #transit_times[population] = []

    for mutation_idx in range(0,len(mutations)):

        #location, gene_name, allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]
        location, gene_name, allele, var_type, codon, position_in_codon, AAs_count,  test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths, clone_times, clone_alts, clone_depths = mutations[mutation_idx]

        state_Ls = state_trajectories[mutation_idx]

        good_idxs, filtered_alts, filtered_depths = timecourse_utils.mask_timepoints(times, alts, depths, var_type, cutoff_idx, depth_fold_change, depth_change_pvalue)

        freqs = timecourse_utils.estimate_frequencies(filtered_alts, filtered_depths)

        masked_times = times[good_idxs]
        masked_freqs = freqs[good_idxs]
        masked_state_Ls = state_Ls[good_idxs]

        t0,tf,transit_time = timecourse_utils.calculate_appearance_fixation_time_from_hmm(masked_times, masked_freqs, masked_state_Ls)
        #print(t0,tf,transit_time)
        if t0==tf==transit_time==None:
            continue

        #print(masked_times, masked_freqs)

        interpolating_function = timecourse_utils.create_interpolation_function(masked_times, masked_freqs, tmax=100000)

        fs = interpolating_function(times)
        fs[fs<0]=0

        # Record
        Ms += fs
        if masked_state_Ls[-1] in well_mixed_fixed_states:
            fixed_Ms += (times>=tf)


    return times, Ms, fixed_Ms




def create_annotation_map(taxon, gene_data=None):

    if gene_data==None:
        gene_data = parse_gene_list(taxon)

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data
    position_gene_map = {}
    gene_position_map = {}
    # new
    gene_feature_map = {}

    # then greedily annotate genes at remaining sites
    for gene_name, feature, start, end in zip(gene_names, features, gene_start_positions, gene_end_positions):
        gene_feature_map[gene_name] = feature
        for position in range(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position] = gene_name
                if gene_name not in gene_position_map:
                    gene_position_map[gene_name]=[]
                gene_position_map[gene_name].append(position)

    # remove 'partial' genes that have < 10bp unmasked sites
    for gene_name in list(sorted(gene_position_map.keys())):
        if len(gene_position_map[gene_name]) < 10:
            for position in gene_position_map[gene_name]:
                position_gene_map[position] = 'repeat'
            del gene_position_map[gene_name]

    # count up number of synonymous opportunities
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}

    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}

    for gene_name, start, end, gene_sequence, strand in zip(gene_names, gene_start_positions, gene_end_positions, gene_sequences, strands):

        if gene_name not in gene_position_map:
            continue

        if strand=='forward':
            oriented_gene_sequence = gene_sequence
        else:
            oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)

        for position in gene_position_map[gene_name]:

            if gene_name not in effective_gene_synonymous_sites:
                effective_gene_synonymous_sites[gene_name]=0
                effective_gene_nonsynonymous_sites[gene_name]=0

            if 'CDS' not in gene_feature_map[gene_name]:
                continue

            else:
                # calculate position in gene
                if strand=='forward':
                    position_in_gene = position-start
                else:
                    position_in_gene = end-position

                # calculate codon start
                codon_start = int(position_in_gene/3)*3
                if codon_start+3 > len(gene_sequence):
                    continue

                #codon = gene_sequence[codon_start:codon_start+3]
                codon = oriented_gene_sequence[codon_start:codon_start+3]
                if any(codon_i in codon for codon_i in bases_to_skip):
                    continue
                position_in_codon = position_in_gene%3

                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0

                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_synonymous_sites[substitution] += 1

                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in substitution_specific_synonymous_sites.keys()}
    # then annotate promoter regions at remaining sites
    for gene_name,start,end in zip(gene_names,promoter_start_positions,promoter_end_positions):
        for position in range(start,end+1):
            if position not in position_gene_map:
                # position hasn't been annotated yet

                if gene_name not in gene_position_map:
                    # the gene itself has not been annotated
                    # so don't annotate the promoter
                    continue
                else:
                    position_gene_map[position] = gene_name
                    gene_position_map[gene_name].append(position)

    # calculate effective gene lengths
    effective_gene_lengths = {gene_name: len(gene_position_map[gene_name])-effective_gene_synonymous_sites[gene_name] for gene_name in gene_position_map.keys()}
    effective_gene_lengths['synonymous'] = sum([effective_gene_synonymous_sites[gene_name] for gene_name in gene_position_map.keys()])
    effective_gene_lengths['nonsynonymous'] = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in gene_position_map.keys()])
    effective_gene_lengths['noncoding'] = pt.get_genome_size(taxon=taxon)-effective_gene_lengths['synonymous']-effective_gene_lengths['nonsynonymous']


    return position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction


def calculate_synonymous_nonsynonymous_target_sizes(taxon):
    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction  = create_annotation_map(taxon=taxon)
    return effective_gene_lengths['synonymous'], effective_gene_lengths['nonsynonymous'], substitution_specific_synonymous_fraction


def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])


def calculate_codon_sequence(dna_sequence):
    return "".join(codon_table[dna_sequence[3*i:3*i+3]] for i in range(0,len(dna_sequence)/3))


def get_closest_gene_name(location, gene_data):
    gene_names, start_positions, end_positions, gene_sequences, strands = gene_data

    closest_start_idx = numpy.fabs(start_positions-location).argmin()
    closest_end_idx = numpy.fabs(end_positions-location).argmin()
    if fabs(start_positions[closest_start_idx]-location) < fabs(end_positions[closest_end_idx]-location):
        return gene_names[closest_start_idx]
    else:
        return gene_names[closest_end_idx]


var_types = ['synonymous','missense','nonsense','noncoding','indel','sv']


def annotate_gene(position, position_gene_map):

    if position in position_gene_map:
        gene_name = position_gene_map[position]
    else:
        gene_name = 'intergenic'

    return gene_name

def annotate_variant(position, allele, gene_data, position_gene_map):

    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

    # get gene
    gene_name = annotate_gene(position, position_gene_map)

    if allele.startswith('Depth'):
        var_type = 'unknown'
        codon=None
        position_in_codon=None
        fold_count=None
    elif allele.startswith('MOB') or allele.startswith('junction'):
        var_type = 'sv'
        codon=None
        position_in_codon=None
        fold_count=None
    elif allele.startswith('indel'):
        var_type = 'indel'
        codon=None
        position_in_codon=None
        fold_count=None
    elif allele[1:3]=='->':
        # a SNP, so annotate it
        if gene_name=='intergenic':
            var_type = 'noncoding'
            codon=None
            position_in_codon=None
            fold_count=None
        elif gene_name=='repeat':
            var_type = 'repeat'
            codon=None
            position_in_codon=None
            fold_count=None
        else:
            # must be in a real gene
            # so get it
            i = gene_names.index(gene_name)

            gene_start_position = gene_start_positions[i]
            gene_end_position = gene_end_positions[i]
            promoter_start_position = promoter_start_positions[i]
            promoter_end_position = promoter_end_positions[i]
            gene_sequence = gene_sequences[i]
            strand = strands[i]

            if position<gene_start_position or position>gene_end_position:
                var_type='noncoding' # (promoter)
                codon=None
                position_in_codon=None
                fold_count=None
            else:

                if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                    var_type='noncoding'
                    codon=None
                    position_in_codon=None
                    fold_count=None
                else:

                    # calculate position in gene
                    if strand=='forward':
                        position_in_gene = position-gene_start_position
                        oriented_gene_sequence = gene_sequence
                        new_base = allele[3]
                    else:
                        position_in_gene = gene_end_position-position
                        oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
                        new_base = base_table[allele[3]]

                    # calculate codon start
                    codon_start = int(position_in_gene/3)*3
                    codon = oriented_gene_sequence[codon_start:codon_start+3]
                    codon_list = list(codon)
                    position_in_codon = position_in_gene%3
                    if (len(codon_list) == 0) or (len(set(codon_list) - set('ACGT')) != 0) :
                        var_type='unknown'
                        fold_count='unknown'
                        codon='unknown'
                        position_in_codon='unknown'
                        fold_count='unknown'

                    else:
                        codon_list[position_in_codon]=new_base
                        new_codon="".join(codon_list)
                        # one wrong bp in caulobacter reference genome
                        if codon_table[codon]==codon_table[new_codon]:
                            var_type='synonymous'
                        else:

                            if codon_table[new_codon]=='!':
                                var_type='nonsense'
                            else:
                                var_type='missense'

                        # count fold
                        if position_in_codon >= 3:
                            fold_count = 'unknown'

                        else:
                            amino_acids = []
                            for base in ['A','C','T','G']:

                                if position_in_codon == 0:
                                    new_fold_codon = list(base) + codon_list[1:]
                                elif position_in_codon == 1:
                                    new_fold_codon = list(codon_list[0]) + list(base) + list(codon_list[2])
                                else:
                                    new_fold_codon = codon_list[:2] + list(base)

                                amino_acids.append(codon_table["".join(new_fold_codon)])

                            fold_count=len(set(amino_acids))

    else:
        sys.stderr.write("Unknown: %s\n" % allele)
        var_type='unknown'
        codon='unknown'
        position_in_codon='unknown'
        fold_count='unknown'

    return gene_name, var_type, codon, position_in_codon, fold_count



def parse_reference_genome(taxon):
    filename= pt.get_path() + '/' + pt.get_ref_gbff_dict(taxon)

    reference_sequences = []

    # GBK file
    if filename[-3:] == 'gbk':
        file = open(filename,"r")
        origin_reached = False
        for line in file:
            if line.startswith("ORIGIN"):
                origin_reached=True
            if origin_reached:
                items = line.split()
                if items[0].isdigit():
                    reference_sequences.extend(items[1:])
        file.close()

    # FASTA file
    else:
        file = open(filename,"r")
        file.readline() # header
        for line in file:
            reference_sequences.append(line.strip())
        file.close()

    reference_sequence = "".join(reference_sequences).upper()
    return reference_sequence

def calculate_genome_length(taxon=None):
    reference_sequence = pt.classFASTA(pt.get_path() +'/'+ pt.get_ref_fna_dict()[taxon]).readFASTA()
    return sum([len(contig[1]) for contig in reference_sequence])

def print_reference_fasta(reference_sequence):
    print(">chr1")
    for i in range(0,len(reference_sequence),70):
        print(reference_sequence[i:min([len(reference_sequence),i+70])])






def create_gene_size_map(taxon, effective_gene_lengths=None):

    if effective_gene_lengths==None:
        reference_sequence = parse_reference_genome(taxon)
        gene_data = parse_gene_list(taxon)

        gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = gene_data

        position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(taxon=taxon,gene_data=gene_data)


    excluded_genes=set(['synonymous','nonsynonymous','noncoding','masked'])

    gene_size_map = {}
    for gene_name in effective_gene_lengths.keys():

        #if gene_name.startswith('tRNA'):
        #    print gene_name

        if gene_name in excluded_genes:
            continue

        gene_size_map[gene_name] = effective_gene_lengths[gene_name]

    return gene_size_map





#####################################################################
#
# Reads through the Genbank file for the reference and
# compiles a list of genes, tRNAs, etc.
#
#####################################################################
def parse_gene_list(taxon, reference_sequence=None):
    gene_names = []
    start_positions = []
    end_positions = []
    promoter_start_positions = []
    promoter_end_positions = []
    gene_sequences = []
    strands = []
    genes = []
    features = []
    protein_ids = []

    filename= pt.get_path() + '/' + pt.get_ref_gbff_dict(taxon)
    gene_features = ['CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA']
    recs = [rec for rec in SeqIO.parse(filename, "genbank")]
    count_riboswitch = 0
    for rec in recs:
        reference_sequence = rec.seq
        contig = rec.annotations['accessions'][0]
        for feat in rec.features:
            if 'pseudo' in list((feat.qualifiers.keys())):
                continue
            if (feat.type == "source") or (feat.type == "gene"):
                continue

            locations = re.findall(r"[\w']+", str(feat.location))
            if feat.type in gene_features:
                locus_tag = feat.qualifiers['locus_tag'][0]
            elif (feat.type=="regulatory"):
                locus_tag = feat.qualifiers["regulatory_class"][0] + '_' + str(count_riboswitch)
                count_riboswitch += 1
            else:
                continue
            # for frameshifts, split each CDS seperately and merge later
            # Fix this for Deinococcus, it has a frameshift in three pieces
            split_list = []
            if 'join' in locations:
                location_str = str(feat.location)
                minus_position = []
                if '-' in location_str:
                    minus_position = [r.start() for r in re.finditer('-', location_str)]
                pos_position = []

                if '+' in location_str:
                    if taxon == 'D':
                        pos_position = [pos for pos, char in enumerate(location_str) if char == '+']
                    elif taxon == 'J':
                        pos_position = [pos for pos, char in enumerate(location_str) if char == '+']
                    else:
                        pos_position = [r.start() for r in re.finditer('+', location_str)]


                if len(minus_position) + len(pos_position) == 2:
                    if len(minus_position) == 2:
                        strand_symbol_one = '-'
                        strand_symbol_two = '-'
                    elif len(pos_position) == 2:
                        strand_symbol_one = '+'
                        strand_symbol_two = '+'
                    else:
                        # I don't think this is possible, but might as well code it up
                        if minus_position[0] < pos_position[0]:
                            strand_symbol_one = '-'
                            strand_symbol_two = '+'
                        else:
                            strand_symbol_one = '+'
                            strand_symbol_two = '-'

                    start_one = int(locations[1])
                    stop_one = int(locations[2])
                    start_two = int(locations[3])
                    stop_two = int(locations[4])

                    locus_tag1 = locus_tag + '_1'
                    locus_tag2 = locus_tag + '_2'

                    split_list.append([locus_tag1, start_one, stop_one, strand_symbol_one])
                    split_list.append([locus_tag2, start_two, stop_two, strand_symbol_two])

                else:
                    if len(pos_position) == 3:
                        strand_symbol_one = '+'
                        strand_symbol_two = '+'
                        strand_symbol_three = '+'
                    start_one = int(locations[1])
                    stop_one = int(locations[2])
                    start_two = int(locations[3])
                    stop_two = int(locations[4])
                    start_three = int(locations[5])
                    stop_three = int(locations[6])

                    locus_tag1 = locus_tag + '_1'
                    locus_tag2 = locus_tag + '_2'
                    locus_tag3 = locus_tag + '_3'

                    split_list.append([locus_tag1, start_one, stop_one, strand_symbol_one])
                    split_list.append([locus_tag2, start_two, stop_two, strand_symbol_two])
                    split_list.append([locus_tag3, start_three, stop_three, strand_symbol_three])


            else:
                strand_symbol = str(feat.location)[-2]
                start = int(locations[0])
                stop = int(locations[1])
                split_list.append([locus_tag, start, stop, strand_symbol])

            for split_item in split_list:
                locus_tag = split_item[0]
                start = split_item[1]
                stop = split_item[2]
                strand_symbol = split_item[3]


                if feat.type == 'CDS':
                    #  why was a -1 there originally?
                    #gene_sequence = reference_sequence[start-1:stop]
                    gene_sequence = str(reference_sequence[start:stop])
                else:
                    gene_sequence = ""


                if 'gene' in list((feat.qualifiers.keys())):
                    gene = feat.qualifiers['gene'][0]
                else:
                    gene = ""

                if 'protein_id' in list((feat.qualifiers.keys())):
                    protein_id = feat.qualifiers['protein_id'][0]
                else:
                    protein_id = ""


                if strand_symbol == '+':
                    promoter_start = start - 100 # by arbitrary definition, we treat the 100bp upstream as promoters
                    promoter_end = start - 1
                    strand = 'forward'
                else:
                    promoter_start = stop+1
                    promoter_end = stop+100
                    strand = 'reverse'


                if gene_sequence!="" and (not len(gene_sequence)%3==0):
                    print(locus_tag, start, "Not a multiple of 3")
                    continue

                # dont need to check if gene names are unique because we're using
                # locus tags

                start_positions.append(start)
                end_positions.append(stop)
                promoter_start_positions.append(promoter_start)
                promoter_end_positions.append(promoter_end)
                gene_names.append(locus_tag)
                gene_sequences.append(gene_sequence)
                strands.append(strand)
                genes.append(gene)
                features.append(feat.type)
                protein_ids.append(protein_id)

    gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids = (list(x) for x in zip(*sorted(zip(gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes, features, protein_ids), key=lambda pair: pair[1])))

    return gene_names, numpy.array(start_positions), numpy.array(end_positions), numpy.array(promoter_start_positions), numpy.array(promoter_end_positions), gene_sequences, strands, genes, features, protein_ids


def parse_timecourse(filename):

    mutations = []

    file = open(filename,"r")
    header_line = file.readline()
    items = header_line.strip().split(",")

    times = []
    for i in range(10,len(items),2):
        times.append(int(items[i].split(":")[1]))
    times = numpy.array(times)

    for line in file:
        items = line.strip().split(",")
        location = int(items[0])
        gene_name = items[1].strip()
        ref_allele = items[2].strip()
        alt_allele = items[3].strip()
        var_type = items[4].strip()
        test_statistic = float(items[5])
        pvalue = float(items[6])
        cutoff_idx = int(items[7])
        depth_fold_change = float(items[8])
        depth_change_pvalue = float(items[9])

        alts = []
        depths = []
        for i in range(10,len(items),2):
            alts.append(int(float(items[i])))
            depths.append(int(float(items[i+1])))

        alts = numpy.array(alts)
        depths = numpy.array(depths)

        pop_times = times[times<1000000]
        pop_alts = alts[times<1000000]
        pop_depths = depths[times<1000000]

        clone_times = times[times>1000000]-1000000
        clone_alts = alts[times>1000000]
        clone_depths = depths[times>1000000]

        mutations.append((location, gene_name, ref_allele, alt_allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, pop_times, pop_alts, pop_depths, clone_times, clone_alts, clone_depths))

    file.close()

    # sort by position
    keys = [mutation[0] for mutation in mutations]
    keys, mutations = (list(t) for t in zip(*sorted(zip(keys, mutations))))
    return mutations


def parse_annotated_timecourse(population, only_passed=True, min_coverage=5):

    mutations = []

    timecourse_filename =  pt.get_path() + '/data/timecourse_final/' +("%s_annotated_timecourse.txt" % population)

    file = open(timecourse_filename, "r")

    header_line = file.readline()
    items = header_line.strip().split(",")

    times = []
    # 13
    for i in range(16,len(items),2):
        times.append(int(items[i].split(":")[1]))
    times = numpy.array(times)

    # depth line
    depth_line = file.readline()
    items = depth_line.strip().split(",")
    avg_depths = []
    for i in range(16,len(items),2):
        avg_depths.append(float(items[i+1]))
    avg_depths = numpy.array(avg_depths)
    population_avg_depth_times = times[times<1000000]
    population_avg_depths = avg_depths[times<1000000]
    clone_avg_depth_times = times[times>1000000]-1000000
    clone_avg_depths = avg_depths[times>1000000]

    for line in file:
        items = line.strip().split(",")
        location = int(items[0])
        gene_name = items[1].strip()
        allele = items[2].strip()
        var_type = items[3].strip()

        codon = items[4].strip()
        position_in_codon = items[5].strip()
        if (position_in_codon != 'None') and (position_in_codon != 'unknown'):
            position_in_codon = int(position_in_codon)
        fold_count = items[6].strip()
        if (fold_count != 'None') and (fold_count != 'unknown'):
            fold_count = int(fold_count)

        test_statistic = float(items[7])
        pvalue = float(items[8])
        cutoff_idx = int(items[9])
        depth_fold_change = float(items[10])
        depth_change_pvalue = float(items[11])

        duplication_idx = int(items[12])
        fold_increase = float(items[13])
        duplication_pvalue = float(items[14])

        passed_str = items[15]
        if passed_str.strip()=='PASS':
            passed = True
        else:
            passed = False

        alts = []
        depths = []

        for i in range(16,len(items),2):
            alts.append(int(float(items[i])))
            depths.append(int(float(items[i+1])))

        alts = numpy.array(alts)
        depths = numpy.array(depths)

        # zero out timepoints with individual coverage lower than some threshold
        alts *= (depths>=min_coverage)*(avg_depths>=min_coverage)
        depths *= (depths>=min_coverage)*(avg_depths>=min_coverage)

        pop_times = times[(times<1000000)]
        pop_alts = alts[(times<1000000)]
        pop_depths = depths[(times<1000000)]


        clone_times = times[(times>1000000)]-1000000
        clone_alts = alts[(times>1000000)]
        clone_depths = depths[(times>1000000)]

        if passed or (not only_passed):
            mutations.append((location, gene_name, allele, var_type, codon, position_in_codon, fold_count, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, pop_times, pop_alts, pop_depths, clone_times, clone_alts, clone_depths))

    file.close()
    #print(mutations[0])
    # sort by position
    keys = [mutation[0] for mutation in mutations]
    keys, mutations = (list(t) for t in zip(*sorted(zip(keys, mutations))))
    return mutations, (population_avg_depth_times, population_avg_depths, clone_avg_depth_times, clone_avg_depths)

clade_hmm_states = {'A':0,'E':1,'FB':2,'FM':3, 'Fm':4,'PB':5,'PM':6,'Pm':7,'PB*':8}
well_mixed_hmm_states = {'A':0,'E':1,'F':2,'P':3}

UNBORN = well_mixed_hmm_states['A']
EXTINCT= well_mixed_hmm_states['E']
POLYMORPHIC= well_mixed_hmm_states['P']

#ANCESTRAL_FIXED = clade_hmm_states['FB']
#MINOR_FIXED=clade_hmm_states['Fm']
#MAJOR_FIXED=clade_hmm_states['FM']
#ANCESTRAL_POLYMORPHIC=clade_hmm_states['PB']
#MINOR_POLYMORPHIC=clade_hmm_states['Pm']
#MAJOR_POLYMORPHIC=clade_hmm_states['PM']

clade_extinct_states = set([clade_hmm_states['A'],clade_hmm_states['E']])
#clade_fixed_states = set([clade_hmm_states['FB'], clade_hmm_states['FM'], clade_hmm_states['Fm'], clade_hmm_states['PB*']])

#clade_majority_states = set([clade_hmm_states['FB'], clade_hmm_states['FM'], clade_hmm_states['PB'],clade_hmm_states['PM'],clade_hmm_states['PB*']])

#clade_polymorphic_states = set([clade_hmm_states['PB'], clade_hmm_states['PM'], clade_hmm_states['Pm']])

well_mixed_extinct_states = set([well_mixed_hmm_states['A'], well_mixed_hmm_states['E']])
well_mixed_polymorphic_states = set([well_mixed_hmm_states['P']])
well_mixed_fixed_states = set([well_mixed_hmm_states['F']])

FIXED = well_mixed_hmm_states['F']
POLYMORPHIC = well_mixed_hmm_states['P']


def parse_haplotype_timecourse(population):

    haplotype_filename = data_directory+('%s_haplotype_timecourse.txt' % population)

    file = open(haplotype_filename,"r")

    times = numpy.array([float(item) for item in file.readline().split(",")])
    fmajors = numpy.array([float(item) for item in file.readline().split(",")])
    fminors = numpy.array([float(item) for item in file.readline().split(",")])
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    haplotypes = []
    for line in file:
        Ls = numpy.array([float(item) for item in line.split(",")])
        haplotypes.append(Ls)
    file.close()
    return times,fmajors,fminors,haplotypes

def parse_well_mixed_state_timecourse(population):

    haplotype_filename = pt.get_path() + '/data/timecourse_final/' +('%s_well_mixed_state_timecourse.txt' % population)

    file = open(haplotype_filename,"r")

    times = numpy.array([float(item) for item in file.readline().split(",")])
    num_unborn = numpy.array([float(item) for item in file.readline().split(",")])
    num_extinct = numpy.array([float(item) for item in file.readline().split(",")])
    num_fixed = numpy.array([float(item) for item in file.readline().split(",")])
    num_polymorphic = numpy.array([float(item) for item in file.readline().split(",")])

    states = []
    for line in file:
        Ls = numpy.array([float(item) for item in line.split(",")])
        states.append(Ls)
    file.close()
    return times, states


def print_timecourse(mutations):

    print_strs = ['Position', 'Gene', 'Ref Allele', 'Alt Allele', 'Annotation', 'Test statistic', 'P-value', 'Cutoff index', 'Fold change', 'P-value']
    times = mutations[0][10]
    for t in zip(times):
            print_strs.append('AC:%d' % t)
            print_strs.append('DP:%d' % t)

    print(", ".join(print_strs))

    for location, gene_name, ref_allele, alt_allele, var_type, test_statistic, pvalue, cutoff_idx, depth_fold_change, depth_change_pvalue, times, alts, depths in mutations:

        print_strs = [str(location),gene_name, ref_allele, alt_allele, var_type, str(test_statistic), str(pvalue), str(cutoff_idx), str(depth_fold_change), str(depth_change_pvalue)]
        for alt,depth in zip(alts,depths):
            print_strs.append(str(alt))
            print_strs.append(str(depth))

        print(", ".join(print_strs))

def parse_coverage(coverage_filename):
    coverage_file = open(coverage_filename)
    ts = numpy.array([int(item) for item in coverage_file.readline().split()])
    ds = numpy.array([float(item) for item in coverage_file.readline().split()])
    return ts,ds





all_line_colors = ['#5DA5DA', '#FAA43A', '#60BD68', '#B276B2', '#F15854', '#4D4D4D']*2


#nonmutator_line_colors = ['#d0d1e6', '#a6bddb', '#67a9cf', '#3690c0', '#02818a', '#016450']
# or reversed (blue to light)
#nonmutator_line_colors = ['#016450', '#02818a', '#3690c0', '#67a9cf', '#a6bddb', '#d0d1e6']
#mutator_line_colors = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f', '#990000']

nonmutator_line_colors = ['#4A1486','#807DBA','#084594','#4292C6','#005A32','#41AB5D']

mutator_line_colors = ['#8C2D04','#CC4C02','#B10026','#E31A16','#FC4E2A','#FD8D3C']

line_color_map = {'m5': '#4A1486', 'p2': '#807DBA', 'p4': '#084594','p1': '#4292C6', 'p5': '#005A32', 'm6': '#41AB5D', 'm1': '#8C2D04', 'm2': '#CC4C02', 'm3': '#B10026', 'm4': '#E31A16', 'p3': '#FC4E2A', 'p6': '#FD8D3C'}

def get_line_color(population):
    return line_color_map[population]

#4D4D4D (gray)
#5DA5DA (blue)
#FAA43A (orange)
#60BD68 (green)
#F17CB0 (pink)
#B2912F (brown)
#B276B2 (purple)
#DECF3F (yellow)
#F15854 (red)

def parse_simulation_trajectory(line):
    items = line.split(":")
    #sys.stderr.write(line)
    params = [float(subitem) for subitem in items[0].split(",")]
    xs = numpy.array([float(subitem) for subitem in items[1].split(",")])
    ps = numpy.array([float(subitem) for subitem in items[2].split(",")])

    if items[3].count(",") > 0:
        dxs = numpy.array([float(subitem) for subitem in items[3].split(",")])
    else:
        dxs = [0]
    #sys.stderr.write("%s \n" % items[4])
    if items[4].count(";") > 0:
        s_trajectories = []
        for subitem in items[4].split(";")[:-1]:
            s_trajectory = []
            if subitem.count(",") > 0:
                s_trajectory = [float(subsubitem) for subsubitem in subitem.split(",")]
            else:
                s_trajectory = []
            if len(s_trajectory) < 5:
                s_trajectory.extend([0]*(5-len(s_trajectory)))
            s_trajectories.append(numpy.array(s_trajectory))
    else:
        s_trajectories = [numpy.array([0,0,0,0,0])]

    if items[5].count(";") > 0:
        x_trajectories = numpy.array([[float(subsubitem) for subsubitem in subitem.split(",")] for subitem in items[5].split(";")[:-1]])
    else:
        print("No x trajectories!")
        x_trajectories = numpy.array([xs])

    if items[6].count(";") > 0:
        p_trajectories = numpy.array([[float(subsubitem) for subsubitem in subitem.split(",")] for subitem in items[6].split(";")[:-1]])
    else:
        print("No p trajectories!")
        p_trajectories = [ps]

    trajectory_data = [params,xs,ps,dxs,s_trajectories,x_trajectories,p_trajectories]
    for i in range(7,len(items)):
        trajectory_data.append([float(subitem) for subitem in items[i].split(",")])
    return trajectory_data

def parse_simulation_trajectories_bzip(filename):
    trajectories = []
    file = BZ2File(filename, "r")
    for line in file:
        #print line
        trajectories.append(parse_simulation_trajectory(line))
    return trajectories

def parse_simulation_trajectories(filename):
    trajectories = []
    file = open(filename, "r")
    for line in file:
        try:
            items = line.split(":")
            params = [float(subitem) for subitem in items[0].split(",")]
            xs = numpy.array([float(subitem) for subitem in items[1].split(",")])
            ps = numpy.array([float(subitem) for subitem in items[2].split(",")])

            if items[3].count(",") > 0:
                dxs = numpy.array([float(subitem) for subitem in items[3].split(",")])
            else:
                dxs = [0]
            if items[4].count(";") > 0:
                s_trajectories = [numpy.array([float(subsubitem) for subsubitem in subitem.split(",")]) for subitem in items[4].split(";")[:-1]]
            else:
                s_trajectories = [numpy.array([0,0,0,0,0])]
            if items[5].count(";") > 0:
                x_trajectories = numpy.array([[float(subsubitem) for subsubitem in subitem.split(",")] for subitem in items[5].split(";")[:-1]])
            else:
                print("No x trajectories!")
                x_trajectories = numpy.array([xs])
            trajectories.append([params,xs,ps,dxs,s_trajectories,x_trajectories,times])
        except ValueError:
            pass
    file.close()
    return trajectories

def parse_simulation_timecourse_trajectories(filename):
    trajectories = []
    file = open(filename, "r")
    params = [float(subitem) for subitem in file.readline().split()]
    s = float(params[0])
    Ub = float(params[1])
    times = numpy.array([float(item) for item in file.readline().split()])
    avg_xs = numpy.zeros_like(times)
    avg_ps = numpy.zeros_like(times)
    x_trajectories = []
    p_trajectories = []
    s_trajectories = []
    for line in file:
        items = line.split(":")
        xs = numpy.array([float(subitem) for subitem in items[0].split()])
        ps = numpy.array([float(subitem) for subitem in items[1].split()])
        ss = numpy.array([float(subitem.split(",")[0]) for subitem in items[2].split()])
        avg_xs += xs
        avg_ps += ps
        x_trajectories.append(xs)
        p_trajectories.append(ps)
        s_trajectories.append(ss)

    file.close()
    return [[params,xs,ps,['0'],s_trajectories,numpy.array(x_trajectories),times]]

def get_time_indices(smaller_times, bigger_times):
    indices = []
    for t in smaller_times:
        i = numpy.fabs(bigger_times-t).argmin()
        indices.append(i)
    return numpy.array(indices)


def count_differences(mutation_list_1, mutation_list_2):
    unique_mutations = set()
    unique_mutations.update(mutation_list_1)
    unique_mutations.update(mutation_list_2)
    #print unique_mutations
    #print len(unique_mutations),len(mutation_list_1),len(mutation_list_2)
    return 2*len(unique_mutations)-len(mutation_list_1)-len(mutation_list_2)





def parse_parallel_genes(filename):

    parallel_genes_file = open(filename,"r")
    line = parallel_genes_file.readline() # header
    parallel_genes = []
    for line in parallel_genes_file:
        items = line.split(",")
        parallel_genes.append( items[0].strip() )

    return parallel_genes

def parse_convergence_matrix(filename):

    convergence_matrix = {}

    convergence_matrix_file = open(filename,"r")
    # Header line
    line = convergence_matrix_file.readline()
    populations = [item.strip() for item in line.split(",")[2:]]

    for line in convergence_matrix_file:

        items = line.split(",")
        gene_name = items[0].strip()
        length = float(items[1])

        convergence_matrix[gene_name] = {'length':length, 'mutations': {population: [] for population in populations}}

        for population, item in zip(populations,items[2:]):
            if item.strip()=="":
                continue

            subitems = item.split(";")
            for subitem in subitems:
                subsubitems = subitem.split(":")
                mutation = (float(subsubitems[0]), int(subsubitems[1]), float(subsubitems[2]), float(subsubitems[3]) )
                convergence_matrix[gene_name]['mutations'][population].append(mutation)


    return convergence_matrix



#if __name__=='__main__':
#reference_sequence = parse_reference_genome()
#gene_data = parse_gene_list()
#gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data

#position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(taxon=taxon, gene_data=gene_data)

#print("Total:", len(reference_sequence))
#print("Masked:", effective_gene_lengths['masked'])
#print("Synonymous sites:", effective_gene_lengths['synonymous'])
#print("Nonsynonymous sites:", effective_gene_lengths['nonsynonymous'])
#print("Noncoding sites:", effective_gene_lengths['noncoding'])

#print("Nonsynonymous:synonymous ratio:", effective_gene_lengths['nonsynonymous']/effective_gene_lengths['synonymous'])
#print("Noncoding:synonymous ratio:", effective_gene_lengths['noncoding']/effective_gene_lengths['synonymous'])
#print(len(gene_names), "genes")
