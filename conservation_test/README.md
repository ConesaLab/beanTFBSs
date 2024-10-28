This folder contains scripts used in the paper "Profiling conserved transcription factor binding motifs in *Phaseolus vulgaris* through comparative genomics" to identify conserved transcription factor binding motifs (TFBSs) in promoter regions among *P. vulgaris*, *V. angularis*, *V. radiata*, and *G. max*. Below, we describe the steps performed and provide supporting intermediate files where relevant. 
1. Download orthologs from the Ensembl Plants BioMart (https://plants.ensembl.org/):
   - `p_vulgaris_vs_g_max.txt`
   - `p_vulgaris_vs_vigna_angularis.txt`
   - `p_vulgaris_vs_vigna_radiata.txt`
     
2. For convenience, extract orthologous gene names into one file using `combine_orholog_files.py`. The output of the script is:
   - `orthologs.txt`, which lists all common bean genes, the ortholog with the highest protein percent identity for each corresponding species, and all orthologs for that species.
     
3. Extract all analyzed species promoters using `extract_up2000_down200_promoter_positions.py`:
   
    `python extract_up2000_down200_promoter_positions.py --gtf species_gtf_file --gtf_format gtf/gff3 --species p_vulgaris/v_angularis/v_radiata/g_max`
  
   The script generates a bed file of all promoters, spanning 2000 upstream to 200 downstream of TSS. Further, bed files were converted to fasta files using bedtools.
   
4. Use extracted *P. vulgaris* promoters and known transcription factor binding sites in MEME format to run FIMO:

    `fimo -text -bfile p_vulgaris_background_markov_zero.txt Pvu_TF_binding_motifs.meme p_vulgaris_promoters_2000_200.fasta > pv_tfbs_fimo.tsv`

6. The common bean promoters were separately aligned with *V. angularis*, *V. radiata*, and *G. max* promoters using minimap2:

   `minimap2 -c -secondary=yes compared_species+promoters.fasta p_vulgaris_promoters.fasta > output.paf`
   
8. The alignment files were processed using `process_alignments.py` to add information about the orthology type, protein homology between the two orthologs, and other useful information for downstream analysis. The output files of this step are as follows:
   - `*species*_orthologs_combined_alignment.txt` contains information for combined alignments. If minimap2 identified more than one alignment region between the two promoters, the information is combined. The file contains the following columns: q_gene (query gene), q_len (the length of the query alignment), ortholog (name of orthologous gene), homology_perc (homology percentage between the orthologous genes), matches_combined (number of matching nucleotides across all alignment regions between the two promoters), pi_combined (percent identity between the two promoters), homologous_prom (True if alignment exists between the two promoters; False if no alignment), orthology_type, number_of_alignments (number of alignment regions spanning the promoter), species (*V. angularis, V. radiata,* or *G. max*).
   - `*species*_orthologs_alignment_info.txt` with the following columns: q_gene (query gene), ortholog (name of orthologous gene), homology_perc (homology percentage between orthologous genes), matches (number of matching nucleotides), q_start (start of the alignment span on the *P. vulgaris* promoter), q_end (end of the alignment span on the *P. vulgaris* promoter), s_start (start of the alignment span on the orthologous gene's promoter), s_end (end of the alignment span on the orthologous gene's promoter), orthology_type, number_of_alignments (indicates whether this is the only alignment span between the two promoters or if there are more), species (*V. angularis*, *V.radiata*, or *G. max*).

9. Use the '*species*_orthologs_alignment_info.txt', fasta sequences of the orthologous genes' promoters, and FIMO output to run the `conservation_test.py` script:

     `python conservation_test.py --hom_promoters *species*_orthologs_alignment_info.txt --orth_seq orthologous_genes_promoters.fasta --tfbs pv_tfbs_fimo.tsv --species name_of_analyzed_species`
   
The conservation test output contains the next columns: q_gene (*P. vulgaris* gene), tf (transcription factor gene), tfbs_start (relative start of a TFBS), tfbs_end (relative end of a TFBS), tfbs (sequence of a TFBS), species (*V. angularis*, *V. radiata*, *G. max*), orth_gene (orthologous gene), hom_tfbs_start (start of the conserved TFBS on the orthologous gene's promoter).

11. Convert the relative positions of TFBSs to genomic positions using `to_genomic_positions.py`.

    `python to_genomic_positions.py --conserved_motifs conservation_test_output --q_gtf p_vulgaris_gtf_file --species1 v_angularis --gtf1 v_angularis_gtf_file --species2 v_radiata --gtf2 v_radiata_gtf_file --species3 g_max --gtf3 g_max_gtf_file --tf_info tf_information`

The TF information file contains two columns: 'Gene_id' (transcription factor gene) and 'Family' (transcription factor family).
   The output file is located in this repository:
   `beanTFBSs/parsing_script/S_Table_2_V1_Conserved_Motifs.txt`.
   
13. The genomic positions were tested with `test_genomic_positions.py`.

*All output files are in a tsv format.*
