options(stringsAsFactors = F)
rm(list = ls())
# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("DESeq2")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("writexl")
source("scripts/functions.R")

# Fixed variables:
tableWithNormalizedExpression <- "output/inputs_pc/FPKM_norm_onlypc.txt"
tableWithCounts <- "output/inputs_pc/Counts_onlypc.txt"
tableWithTF <- "output/inputs_pc/Ravasi_TF_ensID.txt"
fixedColors<- list('tissue'=c('ECT'="#FC7C76",'MES'="#17B7FF"),'stage'=c('E9.5'="#E3A208",'E10.5'="#D87FFF",'E11.5'="#1CD570",'E12.5'="#86C306"),'Replicate'=c('1'="mediumturquoise",'2'="orchid")) 
samplesPlan <- "input/samplesplan.txt"
epsilon <- 0.0000001
cols.hm <- names(fixedColors)
cc <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(103)
factorToStudy <- "stage"
pathForDESeq2 <- "output/DESeq2/"

# Prepare the inputs
samplesPlanDF<-read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
count.table <- read.delim(tableWithCounts)
rownames(count.table)<-count.table$gene_id
expressionDF<-read.delim(tableWithNormalizedExpression)
metaCols<-which(sapply(colnames(expressionDF),function(cn){class(expressionDF[,cn])!="numeric"}))
colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
data <- expressionDF[, -metaCols]
ldata <- log2(1 + data)
tf.id <- readLines(tableWithTF)

if (!dir.exists(pathForDESeq2)){
  dir.create(pathForDESeq2, recursive = T)
}
if (! dir.exists("figure")){
  dir.create("figure")
}
if (! dir.exists("additional_files")){
  dir.create("additional_files", recursive = T)
}

# Set figure names
fig.names <- list("MES" = c('7' = "Fig2C"), "ECT" = c('6' = "Fig5C"))

for (tis in names(fig.names)){
  # For each tissue
  # Select the samples
  new.samples.plan <- subset(samplesPlanDF, tissue == tis)
  samplesToPlot <- new.samples.plan$sample
  # Compute DESeq2 with all stages LRT
  if ( ! file.exists(paste0(pathForDESeq2,"/",tis,"_allStages_DESeq2significant.txt"))){
    signif <- simpleDeseqAna(count.table, factorToStudy, paste0(pathForDESeq2,"/",tis,"_allStages_"),
                             new.samples.plan,
                             LRT = T, lfcT = 0)
  } else {
    signif <- read.delim(paste0(pathForDESeq2,"/",tis,"_allStages_DESeq2significant.txt"))
  }
  # Annotate the significant genes
  signif$tf <- signif$gene_id %in% tf.id
  # Check that variance is not equal to 0
  genesToPlot <- signif$gene_id
  varOfGenes <- apply(expressionDF[match(genesToPlot, expressionDF$gene_id), samplesToPlot], 1, var)
  if (any(varOfGenes == 0)){
    print("Genes", paste(genesToPlot[which(varOfGenes == 0)], collapse = ", "), "were significant but FPKM values is not variable.")
    genesToPlot <- genesToPlot[varOfGenes > 0]
  }
  # Convert the columns of samplesPlan to factors
  factorizedSP <- samplesPlanDF[samplesToPlot, ]
  for (cn in colnames(factorizedSP)){
    uniqVal <- unique(factorizedSP[,cn])
    factorizedSP[,cn] <- factor(factorizedSP[,cn], levels=uniqVal)
  }
  # Prepare a data frame with annotations for samples
  annotForHM <- factorizedSP[, cols.hm]
  # Subset the genes:
  df.sub <- ldata[match(genesToPlot, expressionDF$gene_id), samplesToPlot]
  rownames(df.sub) <- signif$gene_short_name[match(genesToPlot, signif$gene_id)]
  # Scale the expression and perform clustering
  sub.df.scale <- t(apply(t(df.sub), 2, function(v){(v - mean(v, na.rm=T)) / sd(v, na.rm=T)}))
  d <- as.dist(1 - cor(t(sub.df.scale)))
  hc <- hclust(d)
  # Cut the clustering in nClusters part
  nClusters <- as.numeric(names(fig.names[[tis]]))
  annR <- data.frame(k=letters[1:nClusters][cutree(hc,nClusters)])
  rownames(annR)<-hc$labels
  # Give labels from top to bottom in plotting
  pretty.names <- paste0(firstup(tolower(tis)), 1:nClusters)
  names(pretty.names) <- unique(annR[hc$order, "k"])
  annR$cluster <- pretty.names[annR$k]
  annR$k <- NULL
  # Put rainbow colors
  clusterCol <- rainbow(nClusters)
  names(clusterCol) <- pretty.names
  temp.fixedColors <- c(fixedColors, list(cluster = clusterCol))
  # Plot heatmap
  pdf(paste0("figure/", fig.names[[tis]], ".pdf"),
      title=paste0("HeatmapWithSignificantGenes_", tis, "_scale_", nClusters,"clusters"),
      onefile = TRUE, width= 7, height= 15)
  breaksListAbs<-c(min(apply(df.sub, 1, scale), -1.5) - epsilon,
                   seq(-1.5, 1.5, 0.03),
                   max(apply(df.sub, 1, scale), 1.5) + epsilon)
  pheatmap(df.sub,
           # labels_row=expressionDF[match(genesToPlot, expressionDF$gene_id), "gene_short_name"],
           labels_row=rep("", nrow(df.sub)),
           cluster_rows= TRUE,
           cluster_cols= FALSE ,
           annotation_col = annotForHM,
           annotation_row = annR,
           annotation_colors = temp.fixedColors,
           breaks = breaksListAbs,
           main=paste("log2(1+FPKM)", tis, "genes significant in DESeq2 analysis\neach gene is centered/scaled"),
           clustering_distance_rows="correlation",
           #              clustering_method="ward.D2",
           scale="row",
           color = colorRampPalette(c("blue","white","red"))(103))
  dev.off()
  # Then I do the summary where I merge the replicates
  xfrow<-round(sqrt(nClusters))
  yfrow<-ceiling(nClusters/xfrow)
  pdf(paste0("figure/", fig.names[[tis]], "_summary.pdf"),
      title=paste0("Summary_", tis, "_", nClusters,"clusters"),
      onefile = TRUE, width= 3*yfrow, height= 4.5*xfrow)
  par(mfrow=c(xfrow,yfrow),mar=c(5, 4, 4, 2) + 0.1)
  for(i in 1:nClusters){
    plot(x=factorizedSP[colnames(sub.df.scale), "stage"],
         y=rep(0, ncol(df.sub)), col="white", border="white", 
         main=paste0(pretty.names[i], "\n", length(which(annR$cluster == pretty.names[i])),
                     " genes\n", sum(signif$tf & annR$cluster == pretty.names[i]), " TF"),
         ylim=c(min(sub.df.scale), max(sub.df.scale)))#,las=2)
    meanPerSample <- apply(sub.df.scale[which(annR$cluster == pretty.names[i]), ], 2, mean)
    points(x=factorizedSP[colnames(sub.df.scale), "stage"], y=meanPerSample)
    meanVal <- aggregate(meanPerSample, by=list(stage=factorizedSP[colnames(sub.df.scale), "stage"]), FUN=mean)
    lines(x=meanVal$stage, y=meanVal$x)
  }
  par(mfrow=c(1,1))
  dev.off()
  # I add annotations before exporting
  annR$gene_name <- rownames(annR)
  annR$gene_id <- expressionDF$gene_id[match(annR$gene_name, expressionDF$gene_short_name)]
  annR$isTF <- annR$gene_id %in% tf.id
  # I will export the FPKM values for each cluster:
  new.expressionDF <- expressionDF[match(genesToPlot, expressionDF$gene_id), ]
  new.expressionDF$TF <- "no"
  new.expressionDF$TF[new.expressionDF$gene_id %in% tf.id] <- "yes"
  new.expressionDF$cluster <- annR[new.expressionDF$gene_short_name, "cluster"]
  new.expressionDF <- new.expressionDF[, c(names(metaCols), "TF", "cluster", samplesToPlot)]
  colnames(new.expressionDF)[(3 + length(metaCols)):ncol(new.expressionDF)] <-
    paste0("FPKM_", colnames(new.expressionDF)[(3 + length(metaCols)):ncol(new.expressionDF)])
  excel_list <- split(new.expressionDF, f = new.expressionDF$cluster)
  write_xlsx(excel_list, paste0("additional_files/", tis, "_clusters_with_FPKM_values.xlsx"))
  # Now DESeq2 consecutive stages:
  for (i in 1:(length(fixedColors[[factorToStudy]]) - 1)){
    # Select the samples
    new.samples.plan <- samplesPlanDF[samplesPlanDF$tissue == tis & 
                                        samplesPlanDF[, factorToStudy] %in% names(fixedColors[[factorToStudy]])[i + 0:1], ]
    samplesToPlot <- new.samples.plan$sample
    # Compute DESeq2 with Wald log2FC threshold 1.5
    expe <- paste0(names(fixedColors[[factorToStudy]])[i], "VS", 
                   names(fixedColors[[factorToStudy]])[i + 1])
    deseq2.file <- paste0(pathForDESeq2,"/",tis,"_", expe, "_DESeq2significant.txt")
    if ( ! file.exists(deseq2.file)){
      signif <- simpleDeseqAna(count.table, factorToStudy, gsub("DESeq2significant.txt", "", deseq2.file),
                               new.samples.plan,
                               LRT = F, lfcT = 1.5)
    } else {
      signif <- read.delim(deseq2.file)
    }
    # Add the annotations on annR
    annR[, paste0(expe, "_l2fc")] <- signif$log2FoldChange[match(annR$gene_id, signif$gene_id)]
    annR[, paste0(expe, "_padj")] <- signif$padj[match(annR$gene_id, signif$gene_id)]
  }
  # Export
  write.table(annR, paste0(pathForDESeq2, "/", tis, "_Clusters_DESeq2.txt"), quote = F, row.names = F, sep = "\t")
}

