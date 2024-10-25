This repository is designed for parsing conserved transcription factor binding sites (TFBSs) in common bean (Phaseolus vulgaris). The provided script can process a .txt file, where a single column houses a user-specified list of genes. For proper formatting, refer to the example file at examples/gene_list.txt. The output is a comprehensive report, motifs_summary.pdf, containing essential information about the requested set of genes, facilitating a detailed analysis of TFBSs.

The arguments required to run the tfbs_report.R:
--genes_path /path/to/the/list/of/genes/of/interest
--overlap_path /path/to/overlapping_genes.bed
--sites_path /path/to/S_Table_1_V2_Conserved_Motifs.txt (make sure the file is unzipped)

The components of the report are the next:
1.	This heatmap visually represents the presence or absence of conserved motifs on gene promoters. In this visualization, each row corresponds to a gene, while each column represents different motifs grouped by their respective transcription factor families. Red boxes signify a relatively high number of motifs present on a promoter.
2.	This visual representation illustrates the relationships between transcription factor families, highlighting their co-occurrence on gene promoters. Nodes in the plot represent individual transcription factor families, while edges denote the co-occurrence of putative sites on genes' promoters.
For each gene individually:
3.	This table presents the count of conserved motifs associated with each transcription factor family on a gene's promoter. The "Binomial Statistics" column quantifies the probability of observing the given number of motifs. This probability is calculated based on the expected number derived from the overall frequency of each family at the genome level and the total number of motifs observed on the gene's promoter.
4.	This plot illustrates the locations of observed motifs, with a shaded area indicating the annotated upstream gene.

