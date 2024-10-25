import pandas as pd
import argparse
from collections import defaultdict
import pyranges as pr


parser = argparse.ArgumentParser()
parser.add_argument('--gtf', help = 'path to gtf file', required = True)
parser.add_argument('--gtf_format', help = 'format of the annotation file: "gtf", "gff3"', required = True)
parser.add_argument('--species', help = '', required = True)


args = parser.parse_args()
if args.gtf_format == 'gff3':
    genes = pr.read_gff3(args.gtf).df
else:
    genes = pr.read_gtf(args.gtf).df
genes = genes[genes['Feature'] == 'gene']

# p vulgaris genes: 
# /blue/foster/l.kondratova/pvulgaris/ref/annotationPvulgaris_442_v2.1.gene.gff3
# genes = read_gtf('/blue/foster/l.kondratova/pvulgaris/ref/annotation/Pvulgaris_442_v2.1.gene_exons.gff3')
#genes = read_gtf('/blue/foster/l.kondratova/pvulgaris/conservation_test/references/Glycine_max.Glycine_max_v2.1.52.gtf')

def extract_promoters(chromosome, strand, start, end, gene, promoter_start, promoter_end, bed):
    #add all promotes positions to the list
    if strand == '+':
        promoters_bed = [chromosome, start - promoter_start, start + promoter_end, gene, 1, '+']
    if strand == '-':
        promoters_bed = [chromosome, end - promoter_end, end + promoter_start, gene, 1, '-']
    bed.append(promoters_bed)
    return(bed)

def set_to_zero(start):
    if start < 0:
        start = 0
    return(start)

bed=[]
if args.gtf_format == 'gff3':
    genes.apply(lambda gene: extract_promoters(gene['Chromosome'], gene['Strand'], gene['Start'], gene['End'],
        gene['Name'], 2000, 200, bed), axis = 1)
else:
    genes.apply(lambda gene: extract_promoters(gene['Chromosome'], gene['Strand'], gene['Start'], gene['End'],
        gene['gene_id'], 2000, 200, bed), axis = 1)
bed = pd.DataFrame(bed, columns = ['seqname', 'start', 'end', 'gene_id', 'score', 'strand'])

#set start to 0 if it is < 0 
bed['start']=bed['start'].map(lambda x: set_to_zero(x))

bed.to_csv('{}_promoters_2000_200.bed'.format(args.species), header = False, index = False, sep = '\t')
                                                                                                             
