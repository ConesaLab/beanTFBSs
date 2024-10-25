import pandas as pd
#files1,2,3 are orhtologs file downloaded from Ensembl Plants
file1='p_vulgaris_vs_vigna_angularis.txt'
file2='p_vulgaris_vs_vigna_radiata.txt'
file3='pv_vs_gm.txt'
f1 = pd.read_csv(file1, sep = '\t')
f2 = pd.read_csv(file2, sep = '\t')
f3 = pd.read_csv(file3, sep = '\t')
f1 = f1.rename(columns={'Gene stable ID':'pvGene', 'Vigna angularis gene stable ID':'homolog_id','Vigna angularis homology type':'homology_type', '%id. target Vigna angularis gene identical to query gene':'homology_perc'})
f2 = f2.rename(columns={'Gene stable ID':'pvGene', 'Vigna radiata gene stable ID':'homolog_id', 'Vigna radiata homology type':'homology_type', '%id. target Vigna radiata gene identical to query gene':'homology_perc'})
f3 = f3.rename(columns={'Gene stable ID':'pvGene', 'Glycine max gene stable ID':'homolog_id', 'Glycine max homology type':'homology_type', '%id. target Glycine max gene identical to query gene':'homology_perc'})
#for each homolog file group by pvGene
f1 = f1.groupby('pvGene')
f2 = f2.groupby('pvGene')
f3 = f3.groupby('pvGene')

Create column with all genes from p.vulgaris (unique in 'Gene stable ID')
df = pd.DataFrame()
df['pvGene'] = f1.pvGene.unique()
def select_ortholog(gene, orthologs_file, max_orth_l, all_orth_l):
    print(gene)
    h_gene = None
    all_genes = orthologs_file.get_group(gene[0])
    all_orth = None
    #all_genes = all_genes.drop(all_genes[all_genes.homology_perc < 50].index)
    if all_genes.shape[0] > 0:
        one2one = all_genes.loc[all_genes['homology_type'] == 'ortholog_one2one']
        if one2one.shape[0] == 1:
            h_gene = list(one2one['homolog_id'])
            all_orth = list(one2one['homolog_id'])
        else:
            one2many = all_genes.loc[all_genes['homology_type'] == 'ortholog_one2many']
            if one2many.shape[0] > 0:
                max_perc = one2many.loc[one2many['homology_perc'] == one2many['homology_perc'].max()]
                h_gene = list(max_perc['homolog_id'])
                all_orth = list(one2many['homolog_id'])
    max_orth_l.append(h_gene)
    all_orth_l.append(all_orth)
    print(h_gene, all_orth)
    #return(h_gene, all_orth)

def replace(string):
    if isinstance(string, list):
        #string = string.replace("'", "").replace("[", "").replace("]", "")
        string = (',').join(string)
        return(string)
      
max_orth = []
all_orth = []
df['pvGene'].map(lambda x: select_ortholog(x, f1, max_orth, all_orth))
df['v_angularis'] = max_orth
df['v_angularis'] = df['v_angularis'].apply(lambda x: replace(x))
df['all_v_angularis'] = all_orth
df['all_v_angularis'] = df['all_v_angularis'].apply(lambda x: replace(x))

max_orth = []
all_orth = []
df['pvGene'].map(lambda x: select_ortholog(x, f2, max_orth, all_orth))
df['v_radiata'] = max_orth
df['all_v_radiata'] = all_orth
df['v_radiata'] = df['v_radiata'].apply(lambda x: replace(x))
df['all_v_radiata'] = df['all_v_radiata'].apply(lambda x: replace(x))

max_orth = []
all_orth = []
df['pvGene'].map(lambda x: select_ortholog(x, f3, max_orth, all_orth))
df['g_max'] = max_orth
df['all_g_max'] = all_orth
df['g_max'] = df['g_max'].apply(lambda x: replace(x))
df['all_g_max'] = df['all_g_max'].apply(lambda x: replace(x))

df.to_csv("orthologs.txt", sep = '\t', header=True)
