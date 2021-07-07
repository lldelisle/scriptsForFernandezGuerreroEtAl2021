options(stringsAsFactors = F)
rm(list = ls())
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")

# Fixed variables:
tableWithNormalizedExpression <- "output/inputs_pc/FPKM_onlypc.txt"
fixedColors<- list('tissue'=c('ECT'="#FC7C76",'MES'="#17B7FF"),
                   'stage'=c('E9.5'="#E3A208",'E10.5'="#D87FFF",'E11.5'="#1CD570",'E12.5'="#86C306"),
                   'Replicate'=c('1'="mediumturquoise",'2'="orchid")) 
samplesPlan <- "input/samplesplan.txt"

# Preparation for plots
samplesPlanDF <- read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
samplesToPlot <- samplesPlanDF$sample
factorizedSP <- samplesPlanDF
for (cn in colnames(factorizedSP)){
  uniqVal <- unique(factorizedSP[,cn])
  factorizedSP[,cn] <- factor(factorizedSP[,cn], levels=uniqVal)
}
cols.hm <- intersect(names(fixedColors), colnames(factorizedSP))
annotForHM <- factorizedSP[, cols.hm]
expressionDF<-read.delim(tableWithNormalizedExpression)
metaCols<-which(sapply(colnames(expressionDF),function(cn){class(expressionDF[,cn])!="numeric"}))
colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))

genesToPlot <- readLines("input/geneLists/Fig7_individual_genes.txt")

if (! dir.exists("figure")){
  dir.create("figure")
}

for (gene in genesToPlot){
  df.sub <- factorizedSP
  df.sub$gene <- unname(t(expressionDF[expressionDF$gene_short_name == gene, rownames(factorizedSP)]))
  ggplot(df.sub, aes(x = stage, y = gene)) +
    geom_dotplot(aes(fill = tissue), binaxis = "y", stackdir = "center") +
    scale_fill_manual(values = fixedColors[["tissue"]]) +
    ylab(paste("FPKM of", gene))
  ggsave(paste0("figure/Fig7_", gene, ".pdf"), width = 5, height= 5)
}

