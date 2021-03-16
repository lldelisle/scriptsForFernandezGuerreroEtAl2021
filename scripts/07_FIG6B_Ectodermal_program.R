######### Obtaining the ectodermal program #########
####### Author: Alejandro Castilla-Ibeas
options(stringsAsFactors = F)
rm(list = ls())
# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("biomaRt")
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("eulerr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("DESeq2")
safelyLoadAPackageInCRANorBioconductor("writexl")
source("scripts/functions.R")

if (! dir.exists("figure")){
  dir.create("figure")
}
if (! dir.exists("additional_files")){
  dir.create("additional_files", recursive = T)
}

### load data ####

# comparisons from script "05_DESeqMESECT_eachStage.R" Wald test, all stages pairwise comparisons for each tissue:

input.folder <- "output/DESeq2_MESECT_eachStage/"

e9_1 <- read.table(file.path(input.folder, "E9.5_ECTvsMES_DESeq2Results.txt"), header = T)
e10_1 <- read.table(file.path(input.folder, "E10.5_ECTvsMES_DESeq2Results.txt"), header = T)
e11_1 <- read.table(file.path(input.folder, "E11.5_ECTvsMES_DESeq2Results.txt"), header = T)
e12_1 <- read.table(file.path(input.folder, "E12.5_ECTvsMES_DESeq2Results.txt"), header = T)

### 

### Subset significant genes for the ectoderm ####

E9_dif1.5_ECT<-subset(e9_1, log2FoldChange < -1.5 & padj < 0.05)

E10_dif1.5_ECT<-subset(e10_1, log2FoldChange < -1.5 & padj < 0.05)

E11_dif1.5_ECT<-subset(e11_1, log2FoldChange < -1.5 & padj < 0.05)

E12_dif1.5_ECT<-subset(e12_1, log2FoldChange < -1.5 & padj < 0.05)

### Subset significant genes for the mesoderm ####

E9_dif1.5_MES<-subset(e9_1, log2FoldChange > 1.5 & padj < 0.05)

E10_dif1.5_MES<-subset(e10_1, log2FoldChange > 1.5 & padj < 0.05)

E11_dif1.5_MES<-subset(e11_1, log2FoldChange > 1.5 & padj < 0.05)

E12_dif1.5_MES<-subset(e12_1, log2FoldChange > 1.5 & padj < 0.05)

### Make euler diagram of the 4 comparisons ####

pdf("figure/Fig6Ba.pdf", width = 9)
p <- plot(eulerr::venn(list("E9.5"=E9_dif1.5_ECT$gene_short_name,
                       "E12.5"=E12_dif1.5_ECT$gene_short_name,
                       "E10.5"=E10_dif1.5_ECT$gene_short_name,
                       "E11.5"=E11_dif1.5_ECT$gene_short_name)), labels=F, #labels are manually added with Ps
     fills = list(fill = c("#d8f1c5" ,"#d8d8d8", "#d77f7f","#faf37a")),
     quantities = list(cex=1.9), legend = T)
print(p)
dev.off()

### Merge datasets and retrieve genes significant at all stages ####
dir.create("output/panect/", showWarnings = F)
#My solution:
ECT_SPE1 <- plyr::join_all(list(E9_dif1.5_ECT, E10_dif1.5_ECT, E11_dif1.5_ECT, E12_dif1.5_ECT), by= "gene_id", type = "inner")
#save list
write.table(ECT_SPE1$gene_short_name, "output/panect/panect_1.5threshhold_617genes.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")

### Intersect gene list with count values to run DESeq2 ####
# fpkm_ECTSPE_1.5L2FC <- fpkm_onlypc[fpkm_onlypc$gene_short_name %in% ECT_SPE1genes$x,]
# write.table(fpkm_ECTSPE_1.5L2FC, "threshold_1.5L2FC_ect_fpkms.txt", row.names = F, col.names = T, quote = F)

counts_ect <- read.table("output/inputs_pc/Counts_onlypc.txt", header = T)
counts_ECTSPE_1.5L2FC <- counts_ect[counts_ect$gene_short_name %in% ECT_SPE1$gene_short_name,]
write.table(counts_ECTSPE_1.5L2FC, "output/panect/threshold_1.5L2FC_ect_counts.txt", row.names = F, col.names = T, quote = F)

### DESQE2 run ########
samplesPlan <- "input/samplesplan.txt"
samplesPlanDF<-read.delim(samplesPlan)
rownames(samplesPlanDF) <- samplesPlanDF$sample
samplesPlanDF <- subset(samplesPlanDF, tissue == "ECT")

rownames(counts_ECTSPE_1.5L2FC)<-counts_ECTSPE_1.5L2FC$gene_id
simpleDeseqAna(counts_ECTSPE_1.5L2FC, "stage", "output/panect/ECT_allStages_", samplesPlanDF,
               LRT = T)

### filters for ectoderm program #######

#load results from deseq2 script

ectsigdseq2 <- read.table("output/panect/ECT_allStages_DESeq2Results.txt", header = T)
fpkm_onlypc <- read.table("output/inputs_pc/FPKM_onlypc.txt", header = T)
ectsigdseq2b <- fpkm_onlypc[fpkm_onlypc$gene_short_name %in% ectsigdseq2$gene_short_name,]
write.table(ectsigdseq2b, "output/panect/ECT_Signature_All_1.5L2fc_FPKM.xlsx", sep = "\t", row.names = F, quote = F)  #fpkm values

#apply filters

## potentially stable genes following Lucille guideline

ectdseq2stable <- subset(ectsigdseq2,ectsigdseq2$pvalue > 0.5) 
ectdseq2stableb <- fpkm_onlypc[fpkm_onlypc$gene_short_name %in% ectdseq2stable$gene_short_name,]
writexl::write_xlsx(ectdseq2stableb, "additional_files/ECT_Signature_Stable.xlsx")  #fpkm values


## undetermined
ectsig_undetermined <- subset(ectsigdseq2, ectsigdseq2$pvalue < 0.5 & ectsigdseq2$padj >0.05) 
ectsig_undetermined <- fpkm_onlypc[fpkm_onlypc$gene_short_name %in% ectsig_undetermined$gene_short_name,]
writexl::write_xlsx(ectsig_undetermined, "additional_files/ECT_Signature_undetermined.xlsx")  #fpkm values

# statistically significant
ectsig <- subset(ectsigdseq2, ectsigdseq2$padj < 0.05) 
ectsig <- fpkm_onlypc[fpkm_onlypc$gene_short_name %in% ectsig$gene_short_name,]
writexl::write_xlsx(ectsig, "additional_files/ECT_Signature_Variable.xlsx")  #fpkm values

### pie chart with ggplot2 #####

dfpie <- data.frame(
  group = c("Pot.Stable", "Undetermined", "Sig.Variable"),
  count = c(nrow(ectdseq2stable), nrow(ectsig_undetermined), nrow(ectsig))
)
print(dfpie)
#          group count
# 1   Pot.Stable    66
# 2 Undetermined   185
# 3 Sig.Variable   366

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

sum(nrow(ectdseq2stable), nrow(ectsig_undetermined), nrow(ectsig))

bp<- ggplot2::ggplot(dfpie, aes(x="617 ect signature", y=count, fill=group))+
  geom_bar(width = 0.5, stat = "identity")
bp + coord_polar("y", start=0) + blank_theme + scale_fill_brewer("blues")

ggsave("figure/Fig6Bb.pdf", width = 7, height = 7)
