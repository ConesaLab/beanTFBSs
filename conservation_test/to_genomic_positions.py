import argparse
import pandas as pd
from collections import defaultdict
from gtfparse import read_gtf

parser = argparse.ArgumentParser()
parser.add_argument('--conserved_motifs', help = "Concatenated output of conservation_test.py", required = True)
parser.add_argument('--q_gtf', help = "Query species gtf", required = True)
parser.add_argument('--species1', help = 'First subject species name', required = True)
parser.add_argument('--gtf1', help = "", required = True)
parser.add_argument('--species2', help = 'Second subject species name', required = True)
parser.add_argument('--gtf2', help = "", required = True)
parser.add_argument('--species3', help = 'Third subject species name', required = True)
parser.add_argument('--gtf3', help = "", required = True)
parser.add_argument('--tf_info', help = "Table with 'Gene_id' and 'Family' columns, where Gene_id is transcription factor and Family is its family", required = True)


class gene_info:
    def __init__(self, gene, strand, start, end, chrom):
        self.gene = gene
        self.strand = strand
        self.start = start
        self.end = end
        self.chrom = chrom
        
def add_gene_info(gene, strand, start, end, chrom):
    globals()[gene] = gene_info(gene, strand, start, end, chrom)
    
def create_dict(key, value, d):
    d[key] = value
        
def q_gene_to_genomic_position(gene, gene_info, tf, family, tfbs_seq, relative_start, final_dict):
    identifier = gene+tfbs_seq+str(relative_start)
    if identifier not in final_dict:
        if gene_info.strand == '-':  #relative position is negative if upstream TSS and positive if downstream
            genomic_start = int(gene_info.start) + int(relative_start)
            genomic_end = genomic_start + len(tfbs_seq)
        else:
            genomic_end = int(gene_info.end) - int(relative_start) + 2 
            genomic_start = genomic_end - len(tfbs_seq) 
        final_dict[identifier] = [gene, tf, family, tfbs_seq, gene_info.chrom, genomic_start, genomic_end,
                                 None, None, None, None, None,
                                 None, None, None, None, None,
                                 None, None, None, None, None, None, None]
    else:
        pass
    
def s_gene_genomic_position(q_gene, tfbs_seq, q_relative_start, q_relative_end,
                            s_gene, s_gene_info, s_relative_start, final_dict,species, pos1, pos2, pos3, pos4, pos5):
    identifier = q_gene+tfbs_seq+str(q_relative_start)
    if s_gene_info.strand == '-':  #relative position is negative if upstream TSS and positive if downstream
        genomic_start = int(s_gene_info.start) + int(s_relative_start)
        genomic_end = genomic_start + len(tfbs_seq)
    else:
        genomic_end = int(s_gene_info.end) - int(s_relative_start) + 2
        genomic_start = genomic_end - len(tfbs_seq)  
    final_dict[identifier][pos1] = species
    final_dict[identifier][pos2] = s_gene
    final_dict[identifier][pos3] = s_gene_info.chrom
    final_dict[identifier][pos4] = int(genomic_start)
    final_dict[identifier][pos5] = int(genomic_end)
    final_dict[identifier][22] = q_relative_start
    final_dict[identifier][23] = q_relative_end
        
args = parser.parse_args()
gtf = read_gtf(args.q_gtf)
s1 = args.species1
gtf1 = read_gtf(args.gtf1)
s2 = args.species2
gtf2 = read_gtf(args.gtf2)
s3 = args.species3
gtf3 = read_gtf(args.gtf3)
tf_info = pd.read_csv(args.tf_info, sep = '\t')

motifs = pd.read_csv(args.conserved_motifs, sep = '\t', low_memory = False) #gm_tfbs
gtf[gtf['feature'] == 'gene'].apply(lambda x: add_gene_info(x['gene_id'], x['strand'], x['start'], x['end'], x['seqname']), axis = 1) 
gtf1[gtf1['feature'] == 'gene'].apply(lambda x: add_gene_info(x['gene_id'], x['strand'], x['start'], x['end'], x['seqname']), axis = 1) 
gtf2[gtf2['feature'] == 'gene'].apply(lambda x: add_gene_info(x['gene_id'], x['strand'], x['start'], x['end'], x['seqname']), axis = 1) 
gtf3[gtf3['feature'] == 'gene'].apply(lambda x: add_gene_info(x['gene_id'], x['strand'], x['start'], x['end'], x['seqname']), axis = 1) 

tf_fam_d = {}
tf_info.apply(lambda x: create_dict(x['Gene_id'], x['Family'], tf_fam_d), axis = 1)

#identifier is a combination of gene+motif+motif_start
#transform coordinated for subject species and add columns to final dictionary values
     
final_dict = {}    
motifs.apply(lambda x: q_gene_to_genomic_position(x['q_gene'], globals()[x['q_gene']], 
                                                  x['tf'], tf_fam_d[x['tf']],
                                                  x['tfbs'], x['tfbs_start'], 
                                                  final_dict), axis = 1)

motifs[motifs['species'] == s1].apply(lambda x: s_gene_genomic_position(x['q_gene'], x['tfbs'], 
                                                                              x['tfbs_start'], x['tfbs_end'], x['hom_gene'], 
                                                                              globals()[x['hom_gene']], 
                                                                              x['hom_tfbs_start'], 
                                                                              final_dict,s1, 
                                                                              7, 8, 9, 10, 11), axis = 1)
motifs[motifs['species'] == s2].apply(lambda x: s_gene_genomic_position(x['q_gene'], x['tfbs'], 
                                                                              x['tfbs_start'], x['tfbs_end'], x['hom_gene'], 
                                                                              globals()[x['hom_gene']], 
                                                                              x['hom_tfbs_start'], 
                                                                              final_dict, s2, 
                                                                              12, 13, 14, 15, 16), axis = 1)
motifs[motifs['species'] == s3].apply(lambda x: s_gene_genomic_position(x['q_gene'], x['tfbs'], 
                                                                              x['tfbs_start'], x['tfbs_end'], 
                                                                              x['hom_gene'], 
                                                                              globals()[x['hom_gene']], 
                                                                              x['hom_tfbs_start'], 
                                                                              final_dict,
                                                                              s3, 
                                                                              17, 18, 19, 20, 21), axis = 1)
                                            
df = pd.DataFrame(final_dict.values())
df.columns = ['q_gene', 'tf', 'tf_family', 'conserved_motif', 'q_chrom', 'q_start', 'q_end',
              s1, '{}_gene'.format(s1), '{}_chrom'.format(s1),'{}_start'.format(s1), '{}_end'.format(s1),
              s2, '{}_gene'.format(s2), '{}_chrom'.format(s2),'{}_start'.format(s2), '{}_end'.format(s2),
              s3, '{}_gene'.format(s3), '{}_chrom'.format(s3),'{}_start'.format(s3), '{}_end'.format(s3),
             'motif_relative_start', 'motif_relative_end']

def add_conservation(s1, s2, s3, record1, record2, record3):
    counter = 0
    for i in [record1, record2, record3]:
        if i in [s1, s2, s3]:
            counter += 1
    return(counter)

df['conservation_level'] = df.apply(lambda x: add_conservation(s1, s2, s3, x[s1], x[s2], x[s3]), axis = 1)

#overlapping_genes.bed contains studied promoter positions overlapping with the body of upstream genes.
bed = pd.read_csv('overlapping_genes.bed', sep = '\t', low_memory = False)
bed.columns = ['source', 'start', 'end', 'score', 'strand', 'gene']

def fill_overlap_dict(gene, overlap_start, overlap_end, d):
    if not gene in d:
        d[gene] = []
        d[gene].append((int(overlap_start)+1, int(overlap_end)+1))
    else: 
        d[gene].append((int(overlap_start)+1, int(overlap_end)+1))
    
#relative positions of motifs are 1-based and relative to tss
def add_overlap(rel_motif_start, rel_motif_end, overlap_d, gene):
    start = rel_motif_start+2000 #make positions relative to the promoter start/end
    end = rel_motif_end+2000
    half_len = int((rel_motif_end - rel_motif_start)/2) #calculate half length of a motif. 
    within_up_gene = False
    if gene in overlap_d:
        for i in overlap_d[gene]:
            overlap_range = range(i[0], i[-1]+1)
            if start + half_len in overlap_range and end - half_len in overlap_range:
                within_up_gene = True
    return(within_up_gene)
    
overlap_d = {}
bed.apply(lambda x: fill_overlap_dict(x['gene'], x['start'], x['end'], overlap_d), axis = 1)
df['within_upstream_gene'] = f.apply(lambda x: add_overlap(x['motif_relative_start'], x['motif_relative_end'], overlap_d, x['q_gene']), axis = 1)

motifs_bed = pd.DataFrame()
motifs_bed['cl'] = f['conservation_level']
motifs_bed['start'] = f['motif_relative_start'] + 2000
motifs_bed['end'] = f['motif_relative_end'] + 2000
motifs_bed['gene'] = f['q_gene']
motifs_bed['n'] = f['within_upstream_gene']

f.to_csv("conserved_sites_genomic_position.txt", sep = '\t', header=True, index=False)  
motifs_bed.to_csv("conserved_motifs.bed", sep = '\t', header=False, index=False)

