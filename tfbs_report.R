library(ggplot2)
library(ggplotify)
library(pheatmap)
library(dplyr)
library(RColorConesa)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggraph)
library(igraph)
library(gridBase)

args <- commandArgs(trailingOnly = TRUE)
genes_path <- args[1]
overlap_path <- args[2]
sites_path <- args[3]


mytheme <- theme_linedraw(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) +
  theme(axis.title.x = element_text(size=15),
        axis.text.x  = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y  = element_text(vjust=0.5, size=15) ) +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=11), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=16)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) 

write_to_file = function(){
  pdf(file='motifs_summary.pdf', paper="a4")
  if (length(plots.list) > 0){
    # Overall title for the figure
    overall_title <- textGrob("Presence of Conserved Motifs\non Gene Promoters", 
                              gp=gpar(fontface="italic", fontsize=25, col="blue4", vjust = -21, hjust = 0))
    
    grid.draw(overall_title)
    grid.newpage()
    
    heatmap_title <- textGrob("Heatmap of motif presence on gene promoters", gp=gpar(fontsize=18), vjust=-18)
    heatmap_description <- textGrob("The heatmap shows the presence or absence of conserved motifs on gene promoters.\nEach row represents a gene, and each column represents motifs.\nMotifs are grouped by transcription factor family. \nRed boxes indicate relatively high number of motifs on a promoter.", 
                                    gp=gpar(fontsize=12), vjust=-1)
    
    grid.draw(heatmap_title)
    grid.draw(heatmap_description)
    plot.new()
    print(heatmap)
    
    # Network plot
    network_title <- textGrob("Network plot of transcription factor families", gp=gpar(fontsize=18), vjust=-18)
    network_description <- textGrob("The network plot shows the relationships between transcription factor\nfamilies based on their co-occurrence on gene promoters.\nNodes represent transcription factor families,\nand edges represent co-occurrence.", 
                                    gp=gpar(fontsize=12), vjust=-1)
    grid.arrange(gList(network_title, network_description), ncol=1)
    print(network_plot)
    
    # Plots and tables for individual genes
    for (i in 1:length(plots.list)){
      print(plots.list[[i]])

      table.gene = tableGrob(as.data.frame(tables.list[i]), rows=NULL)
      #table_description <- textGrob("The table shows the number of conserved motifs belonging\nto each transcription factor family on the promoter of a gene.\nThe Binomial statistics column indicates the probability of \nobserving the observed number of motifs, given the expected number \nbased on the overall frequency of each family on the genome level.", 
      #                             gp=gpar(fontsize=12), vjust=-2)
      #gt.table = gTree(children = gList(table.gene, table_description))
      grid.arrange(table.gene, ncol=1)
    }
  }
  dev.off()
  print("Summary file generated.")
}


# Load data
genes = read.table(genes_path, header = F, as.is = T)
colnames(genes) = c('gene')

overlap = read.table(overlap_path, header = F, as.is = T, sep = '\t')

colnames(overlap) = c('region', 'start', 'end', 'score', 'gene.strand', 'q_gene')

sites = read.table(sites_path, header = T, as.is = T, sep = '\t')
sites$uniq_id = paste(sites$tf, sites$q_start, sites$q_end)
uniq.sites.n = length(unique(sites$uniq_id))

uniq_sites = sites[!duplicated(sites$uniq_id), ]
fam.counts = as.data.frame(table(uniq_sites$tf_family), row.names='Var1') 

#subset motifs that are on requested genes
motifs.subset = subset(sites, sites$q_gene %in% genes$gene)

#from the gene list subset only genes that have conserved motifs
genes = subset(genes, genes$gene %in% motifs.subset$q_gene)

#generate network of tf families
network = unique(motifs.subset[ , c("q_gene", "tf_family")]) 
network = crossprod(table(network$q_gene, network$tf_family))
diag(network) <- 0
#network <- as.data.frame(network)

p1 <- graph_from_adjacency_matrix(network, 
                                  mode='plus', 
                                  diag=F, weighted = T)

V(p1)$degree <- degree(p1)

network_plot = ggraph(p1, 'igraph', algorithm = 'kk') + 
  geom_edge_link2(aes(width = weight), edge_alpha = 0.1) + 
  geom_node_point(aes(size = 6, colour = name), show.legend = F) +
  scale_color_conesa(palette = "main") +
  geom_node_text(aes(label = name), color = 'black', 
                 size = 3) + 
  theme_void() +
  theme(legend.position = "none")

#generate heatmap
sub_counts_table = as.data.frame.matrix(table(motifs.subset[, c("q_gene", "tf_family")]))

sub_counts_table = sub_counts_table %>% select_if(colSums(.) != 0)
sub_counts_scaled = scale(sub_counts_table)
heatmap = pheatmap(sub_counts_scaled, legend = F)

#gene = 'PHAVU_L007800g' --> could be an example
#gene = 'PHAVU_001G219900g'

# Define function for binomial test
function_PV <- function(x, n, p)
{
  binom.test(x = x,
             n = n,
             p = p,
             alternative = "greater",
             conf.level = 0.95)$p.value
}

titles = genes$gene

# Initialize lists to store results
plots.list = list()
tables.list = list()

# Loop through genes
for (i in 1:length(titles)){
  gene = titles[i]
  
  # Subset sites for current gene
  subs = subset(sites, sites$q_gene == gene)
  
  # Check if any sites were found for current gene
  if (!dim(subs)[1] == 0){
    
    # Subset count table for current gene
    gene.table = sub_counts_table[rownames(sub_counts_table) == gene, ] %>% select_if(colSums(.) != 0)
  
    # Transpose table and rename columns
    gene.table = as.data.frame(t(gene.table))
    colnames(gene.table) = c('motifs_count')
    
    # Count number of promoter conserved sites
    prom.cons.sites = dim(subs)[1]
    gene.table$n = prom.cons.sites
    
    # Merge with family counts
    gene.table <- merge(gene.table, fam.counts,
                        by = 'row.names', all = FALSE)
    
    # Calculate frequency of TF family on promoters
    gene.table['p'] = gene.table$Freq / uniq.sites.n 
    
    # Calculate binomial test for each TF family
    gene.table$binom_statistics <- mapply(FUN = function_PV, x = gene.table$motifs_count, n = gene.table$n,
                                          p = gene.table$p)
    
    # Subset columns for output table
    gene.table = gene.table[c('Row.names', 'motifs_count', 'binom_statistics')]
    colnames(gene.table) = c('TF family', 'Motifs count', 'Binomial statistics')
    
    # Add gene name to table
    gene.table['Gene'] = gene
    
    # Generate plot
    g1 = ggplot(subs, aes(x=motif_relative_start, y=tf_family, color = tf_family)) + 
      geom_point(size=3, stroke = 1.3, shape=21) +
      xlab("Motif relative start") +
      geom_vline(xintercept=0, linetype="dashed") +
      scale_color_conesa(palette = "complete") +
      xlim(-2000, 200) +
      ylab("TF family") +
      mytheme +
      theme(legend.position = "none") +
      labs(title = paste(gene,"\n\n\n"))
    ov_sub = subset(overlap, overlap$q_gene == gene)
    
    # Shadow upstream genes regions
    for(j in 1:nrow(ov_sub)){
      rowDF<-ov_sub[j,]
      dfstart <-rowDF$start - 2000
      dfend<-rowDF$end - 2000
      g1<-g1+annotate("rect", xmin = dfstart, xmax = dfend, ymin = 0, 
                      ymax = Inf, alpha = .2)
    }
    tables.list[[i]] = gene.table
    plots.list[[i]] = g1
    g1  
  }
}

invisible(write_to_file())


