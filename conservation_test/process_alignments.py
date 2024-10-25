import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--alignment', help = 'Minimap alignment in .paf format', required = True)
parser.add_argument('--orthologs', help = 'Ensembl orthologs file', required = True)
parser.add_argument('--pairs', help = 'orthologs.txt obtained with combine_ortholog_files.py', required = True)
parser.add_argument('--query_species', help = 'Query species name as it is in the pairs file (pv_gene)', required = True)
parser.add_argument('--s_species', help = 'Subject species name as it is in the pairs file (v_angularis, v_radiata, or g_max)', required = True)

class ortholog:
    def __init__(self, s_gene):
        self.s_gene = s_gene #a list of one or more genes with the highest gene homology
        self.length = 1
        self.percent_identity = None
        self.homology = None
        self.hom_type = None
        self.matches_pos = {}
        self.s_pos = {}
        
    def add_homology_info(self, homology, homology_type):
        self.homology = homology
        self.hom_type = homology_type
        
    def add_alignment_info (self, length, su_gene, matches_start, matches_end, matches, s_start, s_end):
        if self.s_gene == su_gene:
            self.matches_pos[matches] = (matches_start, matches_end)
            self.length = length
            self.s_pos[matches] = (s_start, s_end) #start and end of homology in orthologous gene

def add_ortholog(q_gene, group, gene_list, homology_d, pairs_d):
    gene_list.append(q_gene)
    if len(set(group['s_gene'])) == 1:
        globals()[q_gene] = ortholog(s_gene = list(group['s_gene'])[0])
        #gene_list.append(q_gene)
    if len(set(group['s_gene'])) > 1:
        group['homology'] = group['s_gene'].apply(lambda x: homology_d[x] if x in homology_d else 0)
        highest_h = group.loc[group['homology'] == group['homology'].max()]
        globals()[q_gene] = ortholog(s_gene = list(highest_h['s_gene'])[0]) 
        #gene_list.append(q_gene)
    if len(group['s_gene']) == 0:
        non_aligned = pd.DataFrame()
        non_aligned['genes'] = pairs_d[q_gene].split(',')
        non_aligned['hom'] = non_aligned['genes'].apply(lambda x: homology_d[x])
        highest_h = non_aligned.loc[non_aligned['hom'] == non_aligned['hom'].max()]
        globals()[q_gene] = ortholog(s_gene = list(highest_h['genes'])[0])

    
def subset_orthologs(q_gene, s_gene, homology_perc, homology_type, orth_subset, gene_list):
    if q_gene in gene_list and globals()[q_gene].s_gene == s_gene:
        orth_subset.append([q_gene, s_gene, homology_perc, homology_type])
        
        
def create_dict(key, value, d):
    #!!!value to list
    if isinstance(value, str):
        value = value.split(',')
        d[key] = list(value) 

def create_dict2(key, value, d):
    d[key] = value

args = parser.parse_args()
alignment = pd.read_csv(args.alignment, sep = '\t')
alignment.columns = ['q_gene', 'q_len', 'q_start', 'q_end', 'same_strand',
                     's_gene', 's_len', 's_start', 's_end', 'matches', 
                     'matches_with_gaps', 'q', 't1', 't2', 't3', 't4', 't5', 't6', 
                     't7', 't8', 't9', 't10', 't11']
#filter off all alignments that are on the opposite strands
alignment = alignment.loc[alignment['same_strand'] == '+']

orthologs = pd.read_csv(args.orthologs, sep = '\t')

#in Ensembl orthologs file:
#column[0] is P. vulgaris gene name
#column[3] orthology type (one-to-one, one-to-many)
#column[4] protein homology percentage between the two orthologs
#column[10] is matching orhologous gene
orthologs = orthologs.rename(columns={orthologs.columns[0]:'q_gene', orthologs.columns[10]:'orth_gene', 
                                      orthologs.columns[3]:'homology_type', orthologs.columns[4]:'homology_perc'})
orthologs = orthologs[['q_gene', 'orth_gene', 'homology_perc', 'homology_type']]
pairs = pd.read_csv(args.pairs, sep='\t')
#pvGene  v_angularis     v_raidata       g_max

q_species = args.query_species
s_species = args.s_species

#creatre dictionary of all orthologs and their homology level
homology_d = {}
orthologs.apply(lambda x: create_dict2(x['orth_gene'], x['homology_perc'], homology_d), axis = 1) 
#in the alignemnts files, if one promoter aligns to multiple orthologs, select one with the highest gene homology
#group by q_gene
alignment_groupped = alignment.groupby('q_gene')
#subset by orthologs
#if shape[1] > 1
#apply another column with homology
#select rows with max homology
#copy the rows to subsetted alignment df
pairs_d = {}
all_orthologs_column = 'all_{s}'.format(s=s_species)
pairs.apply(lambda x: create_dict2(x[q_species], x[all_orthologs_column], pairs_d), axis = 1)

#initiate object for each gene
gene_list = [] #list of genes with known orthologs 

for key in pairs_d:
#for name, group in alignment_groupped:
    if pairs_d[key] == pairs_d[key]:
        name = key
        try:
            group = alignment_groupped.get_group(name)
            name_orth = pairs_d[name].split(',')
            group = group[group['s_gene'].isin(name_orth)]
        except KeyError:
            group = pd.DataFrame()
            group['s_gene'] = []
        #group = group[group['s_gene'] == group['s_gene']]
        add_ortholog(name, group, gene_list, homology_d, pairs_d)
        

#create dataframe with information on othology type for each pv gene and its orhtolog
orth_subset = []
orthologs.apply(lambda x: subset_orthologs(x['q_gene'], x['orth_gene'], x['homology_perc'], x['homology_type'], orth_subset, gene_list), 
                  axis = 1)
orth_subset = pd.DataFrame(orth_subset)
orth_subset.columns = ['q_gene', 'orth_gene', 'homology_perc', 'homology_type']

#add homology and homology types to orhtologs objects
orth_subset.apply(lambda x: globals()[x['q_gene']].add_homology_info(x['homology_perc'], x['homology_type']), axis = 1)

#subset alignments for pv genes that have orhtologs
alignment = alignment[alignment['q_gene'].isin(gene_list)]

#add alignment info orthologs objects
alignment.apply(lambda x: globals()[x['q_gene']].add_alignment_info(x['q_len'], x['s_gene'], x['q_start'], x['q_end'], 
                                                                    x['matches'], x['s_start'], x['s_end']), axis = 1)
file1 = [['q_gene', 'q_len', 'ortholog', 'homology_perc', 'matches_combined', 'pi_combined', 
          'homologous_prom', 'orthology_type', 'number_of_alignments', 'species']]
file2 = [['q_gene', 'ortholog', 'homology_perc', 'matches', 'q_start', 'q_end', 's_start', 's_end',
          'orthology_type', 'number_of_alignments', 'species']]
def create_files(gene_list, file1, file2, s_species):
    for i in gene_list:
        align_num = len(globals()[i].matches_pos.keys())
        ortholog_gene = globals()[i].s_gene
        q_len = globals()[i].length
        homology_perc = globals()[i].homology
        matches_combined = sum(globals()[i].matches_pos.keys()) if align_num > 0 else 0
        pi_combined = int(matches_combined)/int(q_len)
        homologous_prom = True if align_num > 0 else False
        orthology_type = globals()[i].hom_type
        file1.append([i, q_len, ortholog_gene, homology_perc, matches_combined, pi_combined,
                     homologous_prom, orthology_type, align_num, s_species])
        for key in globals()[i].matches_pos:
            matches = key
            start = globals()[i].matches_pos[key][0]
            end = globals()[i].matches_pos[key][1]
            s_start = globals()[i].s_pos[key][0]
            s_end = globals()[i].s_pos[key][1]
            file2.append([i, ortholog_gene, homology_perc, matches, start, end, s_start, s_end,
                          orthology_type, align_num, s_species])
            
create_files(gene_list, file1, file2, s_species)
                 
df1 = pd.DataFrame(file1)
df2 = pd.DataFrame(file2)


df1.to_csv('{s}_orthologs_combined_alignment.txt'.format(s = s_species), header = False, index = False, sep = '\t')
df2.to_csv('{s}_orthologs_alignment_info.txt'.format(s = s_species), header = False, index = False, sep = '\t')
