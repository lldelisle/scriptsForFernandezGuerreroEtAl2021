# From Sofía Zdral
# Written 17/11/2020
# Unfortunately due to update in databases the results cannot be reproduced
# The results (selected terms) are in Figure 1 (PCA global), Figure 2 (PCA MES), Figure 5 (PCA ECT), Table 1 and 2. The top 20 GO enriched terms for each comparasion are in additional files. 

### Prepare the libraries

library("goseq")
library("geneLenDataBase")
library(org.Mm.eg.db)

n.genes.in.hm <- 100 # How many genes to retrieve for each PCA
tableWithNormalizedExpression <- "output/inputs_pc/FPKM_onlypc.txt"
expressionDF<-read.delim(tableWithNormalizedExpression)

# First for the PCA:
for (tissue in c("all", "MES", "ECT")){
  pcaDF <- read.delim(paste0("output/PCA_", tissue, "/geneContributionToPCA.txt"))
  for (pc in 1:2){
    # For the first 2 PC
    # Extract the top 100 genes
    top.id <- pcaDF$gene_id[order(pcaDF[, paste0("PC", pc)], decreasing = T)[1:n.genes.in.hm]]
    gene.vector = as.integer(expressionDF$gene_id %in% top.id)
    names(gene.vector) = expressionDF$gene_id
    pwf=nullp(gene.vector,"mm10","ensGene")
    GO.BP = goseq(pwf, "mm10", "ensGene", test.cats=c("GO:BP"))[1:20, c("term", "over_represented_pvalue")]
    GO.MF = goseq(pwf,"mm10","ensGene",test.cats=c("GO:MF"))[1:20, c("term", "over_represented_pvalue")]
  }
}
# Then for the selected terms:
xmax <- 25
# For figure 2B:
# xmax <- 30
ggplot(GO.MF, aes(x = -log10(over_represented_pvalue), y = term)) +
  geom_col(fill = "#ff7f50ff", color = "#e2623dff") +
  ggtitle("GO Term Molecular Function") +
  ylab("") +
  xlab("-log10(P-value)") +
  coord_cartesian(xlim = c(0, xmax))
ggplot(GO.BP, aes(x = -log10(over_represented_pvalue), y = term)) +
  geom_col(fill = "#6395ecff", color = "#5d69e1ff") +
  ggtitle("GO Term Biological Process") +
  ylab("") +
  xlab("-log10(P-value)") +
  coord_cartesian(xlim = c(0, xmax))
# Then for each of the cluster:
for (tissue in c("MES", "ECT")){
  df <- read.delim(paste0("output/DESeq2/", tissue, "_Clusters_DESeq2.txt"))
  for (cluster in unique(df$cluster)){
    selected.id <- df$gene_id[df$cluster == cluster]
    gene.vector = as.integer(expressionDF$gene_id %in% selected.id)
    names(gene.vector) = expressionDF$gene_id
    pwf=nullp(gene.vector,"mm10","ensGene")
    GO.BP = goseq(pwf, "mm10", "ensGene", test.cats=c("GO:BP"))[1:20, c("term", "over_represented_pvalue")]
    GO.MF = goseq(pwf,"mm10","ensGene",test.cats=c("GO:MF"))[1:20, c("term", "over_represented_pvalue")]
  }
}
