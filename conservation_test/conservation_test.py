import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import re


"""
homologous promoters are obtained by aligning orthologous genes promoters to target genes promoters with minimap2 in .paf format:
minimap2 -x splice --secondary=yes -C5 $ref $reads
alignments in .sam format can be converted to .paf with minimap paftools
.paf file must contain columns: 
q_gene (pv_gene in old version)
q_start
q_end
orth_gene (hom_gene in old version)
orth_start
orth_end
"""

parser = argparse.ArgumentParser()
parser.add_argument('--hom_promoters', help = 'Alignment of orthologous promoters to target genes promoters (file2)', required = True)
parser.add_argument('--orth_seq', help = '.fasta sequences of orthologous genes promoters that were used for homologous promoters search', required = True)
parser.add_argument('--tfbs', help = '.tsv FIMO output file with predicted TFBSs in query sequences.', required = True)
parser.add_argument('--species', help = 'Subject species name which will be used in output files and column names', required = False, default = 'out')


args = parser.parse_args()
hom_promoters = pd.read_csv(args.hom_promoters, sep = '\t') #==g_max
orth_seq = SeqIO.parse(open(args.orth_seq),'fasta') #g_max_seq
tfbs = pd.read_csv(args.tfbs, names = ['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand', 'score', 'p-value', 'q-value', 'matched_sequence'], sep = '\t', low_memory=False)
species = args.species

prom_dict = {} #==g_max_prom


def promoters_info (q_gene, q_start, q_end, orth_gene, hom_start, hom_end, dictionary): #(pv_gene, pv_start, pv_end, hom_gene, hom_start, hom_end, dictionary):
    if q_gene not in dictionary:
        dictionary[q_gene] = [[q_start, q_end, orth_gene, hom_start, hom_end]]
    else:
        dictionary[q_gene].append([q_start, q_end, orth_gene, hom_start, hom_end])
        
    
hom_promoters.apply(lambda row: promoters_info(row['q_gene'], row['q_start'], row['q_end'], row['ortholog'], row['s_start'], row['s_end'], prom_dict), axis = 1)

 
orth_dict = {} #g_max_fasta = {}
for fasta in orth_seq:
    orth_dict[fasta.id] = str(fasta.seq)    
#define alignment region from start – end in .paf files
#Check if predicted TFBS is within the alignment region, if
#select region on homologous region -100 : len(matched_sequence) : + 100 


def bs_region(q_start, q_end, motif_start, motif_end, strand, tfbs):
    bs = True if q_start <= motif_start and q_end >= motif_end else False
    tfbs1 = Seq(tfbs)
    tfbs = str(tfbs1.reverse_complement()) if bs and strand == '-' else tfbs
    return(bs, tfbs)

final = []

def find_tfbs(prom_info, orth_dict, tfbs, tfbs_start, tfbs_end, tf, species, q_gene, df, strand):
    try:
        for i in prom_info:
            q_start = i[0]
            q_end = i[1]
            orth_gene = i[2]
            hom_start = i[3]
            hom_end = i[4]
            hom_seq = orth_dict[i[2]]
            bs, corrected_tfbs = bs_region(q_start, q_end, tfbs_start, tfbs_end, strand, tfbs)
            if bs:
                #define window in homologous region of an orthologous gene for a conserved tfbs search
                #+/- 100 nucleotides around the motif in the query
                hom_tf_start = tfbs_start - q_start + hom_start - 100
                hom_tf_end = hom_tf_start + 200 + len(tfbs)
                hom_tf_start = hom_tf_start if hom_tf_start > hom_start else hom_start
                hom_tf_end = hom_tf_end if hom_tf_end < hom_end else hom_end
                region_of_interest = hom_seq[hom_tf_start:hom_tf_end]
                if corrected_tfbs in region_of_interest:
                    hom_cons_start = re.search(corrected_tfbs, region_of_interest).start()
                    hom_tfbs_start = hom_tf_start + hom_cons_start+1
                    #FIMO was ran on promoters of 5000 bp upstream the TSS, therefore to retreave the relative position we subtract 5000 from a tfbs's start and end
                    df.append([q_gene, tf, tfbs_start-5000, tfbs_end-5000, tfbs, species, orth_gene, hom_tfbs_start-5000])
    except KeyError:
        pass  

        
genes_to_keep = list(prom_dict.keys())  
tfbs_orth = tfbs[tfbs['sequence_name'].isin(genes_to_keep)]

tfbs_orth.apply(lambda row: find_tfbs(prom_dict[row['sequence_name']], orth_dict, row['matched_sequence'], int(row['start']), int(row['stop']), row['motif_id'], species, row['sequence_name'], final, row['strand']), axis = 1)  

df = pd.DataFrame(final, columns = ['q_gene', 'tf', 'tfbs_start', 'tfbs_end', 'tfbs', 'species', 'hom_gene', 'hom_tfbs_start'])
df.to_csv(species + "_conserved_sites.txt", sep = '\t', header=True, index=False) 
