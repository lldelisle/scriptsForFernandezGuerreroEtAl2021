options(stringsAsFactors = F)
rm(list = ls())
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)

safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")

# Fixed variables:
tableWithNormalizedExpression <- "output/inputs_pc/FPKM_onlypc.txt"
fixedColors<- list('tissue'=c('ECT'="#FC7C76",'MES'="#17B7FF"),
                   'stage'=c('E9.5'="#E3A208",'E10.5'="#D87FFF",'E11.5'="#1CD570",'E12.5'="#86C306"),
                   'Replicate'=c('1'="mediumturquoise",'2'="orchid")) 
samplesPlan <- "input/samplesplan.txt"
epsilon <- 0.0000001
cc <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(103)
geneListsFolder <- "input/geneLists/"

# Preparation for heatmap
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
data <- expressionDF[, -metaCols]
ldata <- log2(1 + data)

if (! dir.exists("figure")){
  dir.create("figure")
}

# First I plot heatmaps with no annotation for genes:
genelists <- grep("Fig7", list.files(geneListsFolder), value = T, invert = T)
for (figName in gsub(".txt$", "", genelists)){
  genesToPlot <- readLines(file.path(geneListsFolder, paste0(figName, ".txt")))
  df.sub <- ldata[match(genesToPlot, expressionDF$gene_short_name), ]
  rownames(df.sub) <- genesToPlot
  fig.width <- 7.4
  if (grepl("Fig[34]", figName)){
    df.sub <- df.sub[, grep("MES", colnames(df.sub))]
  }
  if (grepl("Fig6", figName)){
    fig.width <- 5
    df.sub <- df.sub[, grep("ECT", colnames(df.sub))]
  }
  pdf(paste0("figure/", figName, ".pdf"),
      title=figName,
      onefile = TRUE, width= fig.width, height= 5 + 0.22 * length(genesToPlot))
  breaksListAbs<-c(seq(0, 10, length.out = 100), max(max(df.sub), 10) + epsilon)
  pheatmap(df.sub,
           cluster_rows= FALSE,
           cluster_cols= FALSE ,
           annotation_col = annotForHM,
           annotation_colors = fixedColors,
           breaks = breaksListAbs,
           main=paste("log2(1+FPKM)", figName),
           cellheight = 16,
           color = cc)
  dev.off()
}

# Then I plot what is in figure 7
# First I get annotations for genes:
summary.df <- read.delim("output/DESeq2_MESECT_eachStage/summary_significant.txt")
ann_color_row <- list()
for (stage in rev(unique(samplesPlanDF$stage))){
  summary.df[, paste0(stage, "_signif")] <- as.character(! is.na(summary.df[, paste0(stage, "_padj")]))
  ann_color_row[[paste0(stage, "_signif")]] <- c('FALSE'='white', 'TRUE'=unname(fixedColors[["stage"]][stage]))
}
genelists <- grep("Fig7[A-Z]", list.files(geneListsFolder), value = T)
for (figName in gsub(".txt$", "", genelists)){
  genesToPlot <- readLines(file.path(geneListsFolder, paste0(figName, ".txt")))
  df.sub <- ldata[match(genesToPlot, expressionDF$gene_short_name), ]
  rownames(df.sub) <- genesToPlot
  annR <- summary.df[summary.df$gene_short_name %in% genesToPlot, grep("_sign", colnames(summary.df))]
  rownames(annR) <- summary.df$gene_short_name[summary.df$gene_short_name %in% genesToPlot]
  fig.width <- 7.4
  pdf(paste0("figure/", figName, ".pdf"),
      title=figName,
      onefile = TRUE, width= fig.width, height= 5 + 0.22 * length(genesToPlot))
  breaksListAbs<-c(seq(0, 10, length.out = 100), max(max(df.sub), 10) + epsilon)
  pheatmap(df.sub,
           cluster_rows= FALSE,
           cluster_cols= FALSE ,
           annotation_col = annotForHM,
           annotation_row = annR,
           annotation_colors = c(fixedColors, ann_color_row),
           breaks = breaksListAbs,
           main=paste("log2(1+FPKM)", figName),
           cellheight = 16,
           color = cc)
  dev.off()
}