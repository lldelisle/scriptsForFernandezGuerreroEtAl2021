options(stringsAsFactors = F)
rm(list = ls())
# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("writexl")
source("scripts/functions.R")

# Fixed variables:
tableWithNormalizedExpression <- "output/inputs_pc/FPKM_onlypc.txt"
fixedColors<- list('tissue'=c('ECT'="#FC7C76",'MES'="#17B7FF"),'stage'=c('E9.5'="#E3A208",'E10.5'="#D87FFF",'E11.5'="#1CD570",'E12.5'="#86C306"),'Replicate'=c('1'="mediumturquoise",'2'="orchid")) 
samplesPlan <- "input/samplesplan.txt"
n.genes.in.hm <- 100 # How many genes to retrieve for each PCA
epsilon <- 0.0000001
nb.genes <- 500 # How many genes to use for PCA analysis
cols.hm <- names(fixedColors)
cc <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(103)

# Prepare the inputs
samplesPlanDF<-read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
expressionDF<-read.delim(tableWithNormalizedExpression)
metaCols<-which(sapply(colnames(expressionDF),function(cn){class(expressionDF[,cn])!="numeric"}))
colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
data <- expressionDF[, -metaCols]
ldata <- log2(1 + data)
sumperline <- apply(ldata,1,sum)
lnzdata <- ldata[sumperline != 0,]
if (! dir.exists("figure")){
  dir.create("figure")
}
if (! dir.exists("additional_files")){
  dir.create("additional_files", recursive = T)
}
# Give the name of the final figures
fig.names <- list("Global" = c("Fig1B", "Fig1C", "Fig1E"),
                  "MES" = c("Fig2A", NA, NA),
                  "ECT" = c("Fig5A", NA, NA)
)

for (tissue in names(fig.names)){
  if (tissue == "Global"){
    samplesToPlot <- intersect(colnames(expressionDF),samplesPlanDF$sample)
  } else {
    samplesToPlot <- grep(tissue, intersect(colnames(expressionDF),samplesPlanDF$sample), value = T)
  }
  # Will transform each column as factor
  factorizedSP <- samplesPlanDF[samplesToPlot, ]
  for (cn in colnames(factorizedSP)){
    uniqVal <- unique(factorizedSP[,cn])
    factorizedSP[,cn] <- factor(factorizedSP[,cn], levels=uniqVal)
  }
  # Build a dataframe with annotations for samples
  annotForHM <- factorizedSP[, cols.hm]
  outputFolder <- paste0("output/PCA_", tissue, "/")
  if ( ! dir.exists(outputFolder)){
    dir.create(outputFolder, recursive = T)
  }
  # Select samples
  rldata <- lnzdata[, samplesToPlot]
  # Select genes
  rldata <- rldata[order(apply(rldata,1,var),decreasing = T)[1:min(nrow(rldata), 500 )],]
  # Compute pca
  sample.pca <- prcomp(t(rldata),
                       center = TRUE,
                       scale. = FALSE)
  # Retrieve loadings
  new.df <- data.frame(factorizedSP,sample.pca$x[samplesToPlot,])
  # Evaluate variance explained
  var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)
  current.figure.names <- fig.names[[tissue]]
  # Plot first 2 PC:
  g <- ggplot(new.df, aes(PC1,PC2)) +
    geom_point( aes(fill=tissue, color=Replicate, shape=stage), stroke=2 ,size=3) +
    theme_grey(base_size = 20) +
    xlab(paste0("PC1: ",var[1],"% variance"))+
    ylab(paste0("PC2: ",var[2],"% variance")) +
    scale_fill_manual(values=fixedColors[["tissue"]])+
    scale_color_manual(values=fixedColors[["Replicate"]])+
    scale_shape_manual(values=c(21,22,23,24)) +
    guides(fill=guide_legend(override.aes = list(shape = 21)),alpha=guide_legend(override.aes = list(shape = 21)), color=guide_legend(override.aes = list(shape = 21)))  +
    ggtitle(paste("PC1-PC2", tissue, nb.genes, "genes most variable"))
  ggsave(paste0("figure/", current.figure.names[1], ".pdf"), g, width= 7.26, height= 7.26)
  # Extract the genes driving the pc
  pcaDF <- cbind(expressionDF[rownames(sample.pca$rotation), metaCols], sample.pca$rotation * sample.pca$rotation)
  write.table(pcaDF, file.path(outputFolder, "geneContributionToPCA.txt"), sep = "\t", row.names = F, quote=F)
  for (pc in 1:2){
    # For the first 2 PC
    # Extract the top 100 genes
    top.id <- pcaDF$gene_id[order(pcaDF[, paste0("PC", pc)], decreasing = T)[1:n.genes.in.hm]]
    if (! is.na(current.figure.names[3 * pc - 1])){
      # Plot a heatmap for those genes
      df.sub <- ldata[match(top.id, expressionDF$gene_id), samplesToPlot]
      pdf(paste0("figure/", current.figure.names[3 * pc - 1], ".pdf"),
          title=paste0("HeatmapWithTop", nb.genes, "PC", pc, "driving_genes Complete clustering"),
          onefile = TRUE, width= 7.4, height= 5 + 0.22 * n.genes.in.hm)
      breaksListAbs<-c(seq(0, 10, length.out = 100), max(max(df.sub), 10) + epsilon)
      pheatmap(df.sub,
               labels_row=expressionDF[match(top.id, expressionDF$gene_id), "gene_short_name"],
               cluster_rows= TRUE,
               cluster_cols= FALSE ,
               annotation=annotForHM,
               annotation_colors = fixedColors,
               breaks = breaksListAbs,
               main=paste("log2(1+FPKM)", tissue, nb.genes, "genes most variable driving PC", pc),
               clustering_distance_rows="correlation",
               cellheight = 16,
               color = cc)
      dev.off()
    }
    writexl::write_xlsx(expressionDF[match(top.id, expressionDF$gene_id), ], paste0("additional_files/", tissue, "_PC", pc, "_with_FPKM_values.xlsx"))
  }
}

