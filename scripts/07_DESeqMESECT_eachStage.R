options(stringsAsFactors = F)
rm(list = ls())
# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("DESeq2")
safelyLoadAPackageInCRANorBioconductor("writexl")
source("scripts/functions.R")

# Fixed variables:
tableWithCounts <- "output/inputs_pc/Counts_onlypc.txt"
samplesPlan <- "input/samplesplan.txt"
factorToStudy <- "tissue"
pathForDESeq2 <- "output/DESeq2_MESECT_eachStage/"

# Prepare inputs
samplesPlanDF<-read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
factorizedSP <- samplesPlanDF
for (cn in colnames(factorizedSP)){
  uniqVal <- unique(factorizedSP[,cn])
  factorizedSP[,cn] <- factor(factorizedSP[,cn], levels=uniqVal)
}
samplesToPlot <- samplesPlanDF$sample
count.table <- read.delim(tableWithCounts)
rownames(count.table)<-count.table$gene_id

if (!dir.exists(pathForDESeq2)){
  dir.create(pathForDESeq2, recursive = T)
}
# Prepare a big table with the results of all DESeq2
big.annot <- count.table[, c("gene_id", "gene_short_name", "locus")]
# For each stage
for (cur.stage in unique(samplesPlanDF$stage)){
  # Select the samples
  new.samples.plan <- subset(samplesPlanDF, stage == cur.stage)
  # Run or read DESeq2 results with Wald test threshold of log2FC at 1.5
  if ( ! file.exists(paste0(pathForDESeq2,"/",cur.stage,"_ECTvsMES_DESeq2significant.txt"))){
    signif <- simpleDeseqAna(count.table, factorToStudy, paste0(pathForDESeq2,"/",cur.stage,"_ECTvsMES_"),
                             new.samples.plan,
                             LRT = F, lfcT = 1.5, writeRLOG = F)
  } else {
    signif <- read.delim(paste0(pathForDESeq2,"/",cur.stage,"_ECTvsMES_DESeq2significant.txt"))
  }
  # Add results to the dataframe
  big.annot[, paste0(cur.stage, "_l2fc")] <- signif$log2FoldChange[match(big.annot$gene_id, signif$gene_id)]
  big.annot[, paste0(cur.stage, "_padj")] <- signif$padj[match(big.annot$gene_id, signif$gene_id)]
}
# Summary the results
big.annot$is.ECT <- apply(big.annot[, grep("l2fc", colnames(big.annot))], 1, function(v){!all(is.na(v)) & all(na.omit(v) < 0)})
big.annot$is.MES <- apply(big.annot[, grep("l2fc", colnames(big.annot))], 1, function(v){!all(is.na(v)) & all(na.omit(v) > 0)})
big.annot$is.never.signif <- apply(big.annot[, grep("l2fc", colnames(big.annot))], 1, function(v){all(is.na(v))})

big.annot$status <- NA
big.annot$status[big.annot$is.ECT] <- "ECT_spe"
big.annot$status[big.annot$is.MES] <- "MES_spe"
big.annot$status[big.annot$is.never.signif] <- "None"
big.annot$status[is.na(big.annot$status)] <- "Both"
print(table(big.annot$status))
# 67 genes are both
# 1899 are ECT
# 1486 are MES

write.table(big.annot, file.path(pathForDESeq2, "summary_significant.txt"), sep = "\t", quote = F, row.names = F)
write_xlsx(big.annot, "additional_files/Summary_MESECT_eachstage.xlsx")
