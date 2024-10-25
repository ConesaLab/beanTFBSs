import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import pandas as pd
from gtfparse import read_gtf

parser = argparse.ArgumentParser()
parser.add_argument('--gtf', help = 'Reference gtf file', required = True)
parser.add_argument('--fasta', help = 'Reference assembly', required = True)
parser.add_argument('--positions', help = '.tsv file with motifs genomic positions', required = True)
parser.add_argument('--motif', help = 'Name of a column with motifs to search in reference assembly', 
                    required = True)
parser.add_argument('--record_id', help = 'Region identifier column (such as gene name), not required', required = False)
parser.add_argument('--start', help = 'Column name in positions file with start of genomic region to test (default "start").', 
                    default = 'start')
parser.add_argument('--end', help = 'Column name in positions file with end of genomic region to test (default "end")',
                    default = 'end')
parser.add_argument('--chromosome', help = 'Name of column with reagions chromosomes (default "chromosome")', 
                   default = 'chromosome')
parser.add_argument('--same_strand', help = 'Look for matches on the same strand only if True (default False)',
                   default = False)

args = parser.parse_args()
gtf = read_gtf(args.gtf)
ref = SeqIO.parse(open(args.fasta),'fasta')
positions = pd.read_csv(args.positions, header=0, low_memory = False, sep = '\t')

def default_d():
    return(0)

def find_match(chrom, start, end, ref_dict, gene, region, out):
    if start.isdigit():
        region2 = ref_dict[chrom][int(start)-1:int(end)-1]
        match = False
        if args.same_strand and region == region2:
            match = True
        if not args.same_strand:
            if region == region2 or region == Seq(region2).reverse_complement():
                match = True
        else:
            pass
        if args.record_id:
            out.append([gene, region, start, end, region2, match])
        else:
            out.append([region, start, end, region2, match])
        
ref_dict = defaultdict(default_d)
for fasta in ref:
    ref_dict[fasta.id] = str(fasta.seq) 
        
if args.record_id:
    out = [['region_id', 'region', 'start', 'end', 'genomic_region', 'match']]
else:
    out = [['region', 'start', 'end', 'genomic_region', 'match']]

positions.apply(lambda x: find_match(x[args.chromosome], x[args.start], x[args.end], 
                                     ref_dict, x[args.record_id], x[args.motif], out), axis = 1)

out = pd.DataFrame(out)

out.to_csv("test_genomic_positions.txt", sep = '\t', header=False, index=False)  
