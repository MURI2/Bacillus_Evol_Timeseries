from Bio import SeqIO
import phylo_tools as pt
import re
import numpy



def parse_gene_list(taxon='B', reference_sequence=None):
    gene_names = []
    start_positions = []
    end_positions = []
    promoter_start_positions = []
    promoter_end_positions = []
    gene_sequences = []
    strands = []
    genes = []

    filename= pt.get_path() + '/' + pt.get_ref_gbff_dict()[taxon]
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
            strand_symbol = str(feat.location)[-2]
            start = int(locations[0])
            stop = int(locations[1])

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



            if feat.type in gene_features:
                locus_tag = feat.qualifiers['locus_tag'][0]

            elif (feat.type=="regulatory"):
                locus_tag = feat.qualifiers["regulatory_class"][0] + '_' + str(count_riboswitch)
                count_riboswitch += 1

            else:
                continue


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

    gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes = (list(x) for x in zip(*sorted(zip(gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands, genes), key=lambda pair: pair[1])))

    return gene_names, numpy.array(start_positions), numpy.array(end_positions), numpy.array(promoter_start_positions), numpy.array(promoter_end_positions), gene_sequences, strands, genes




parse_gene_list(taxon='B')
