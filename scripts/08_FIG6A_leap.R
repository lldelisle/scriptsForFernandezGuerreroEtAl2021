################ stage pairwise comparisons ###############
######### Author: Alejandro Castilla-Ibeas
options(stringsAsFactors = F)
rm(list = ls())
# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("eulerr")

#### LOADING RESULTS FROM DESEQ2 ####

# comparisons from script "4_DESeqHeatmapAndClusteringPerTissue.R" Wald test, all stages pairwise comparisons for each tissue:

input.folder <- "output/DESeq2/"

ECT9.5vs10.5 <- read.table(file.path(input.folder, "ECT_E9.5VSE10.5_DESeq2significant.txt"), header = T)
ECT10.5vs11.5 <- read.table(file.path(input.folder, "ECT_E10.5VSE11.5_DESeq2significant.txt"), header = T)
ECT11.5vs12.5 <- read.table(file.path(input.folder, "ECT_E11.5VSE12.5_DESeq2significant.txt"), header = T)

MES9.5vs10.5 <- read.table(file.path(input.folder, "MES_E9.5VSE10.5_DESeq2significant.txt"), header = T)
MES10.5vs11.5 <- read.table(file.path(input.folder, "MES_E10.5VSE11.5_DESeq2significant.txt"), header = T)
MES11.5vs12.5 <- read.table(file.path(input.folder, "MES_E11.5VSE12.5_DESeq2significant.txt"), header = T)

#### EULERR ####
if (! dir.exists("figure")){
  dir.create("figure")
}
pdf("figure/Fig6Aa.pdf", width = 9)
p <- plot(eulerr::euler(list("E9.5vsE10.5"=MES9.5vs10.5$gene_short_name,
                             "E10.5vsE11.5"=MES10.5vs11.5$gene_short_name,
                             "E11.5vsE12.5"=MES11.5vs12.5$gene_short_name)), labels =F, #labels were manually added w/ Ps
          quantities = list(cex=2.5), fills = list(fill = c("#c2abfe", "#d44a5e", "#5dd1aa")),
          legend = T)
print(p)
dev.off()
pdf("figure/Fig6Ab.pdf", width = 9)
p <- plot(eulerr::euler(list("E9.5vsE10.5"=ECT9.5vs10.5$gene_short_name,
                             "E10.5vsE11.5"=ECT10.5vs11.5$gene_short_name,
                             "E11.5vsE12.5"=ECT11.5vs12.5$gene_short_name)), labels = F,#labels were manually added w/ Ps
          quantities = list(cex=2.5), fills = list(fill = c("#c2abfe", "#d44a5e", "#5dd1aa")),
          legend = T)
print(p)
dev.off()

