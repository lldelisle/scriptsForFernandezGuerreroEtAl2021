options(stringsAsFactors = F)
rm(list = ls())
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("readxl")
safelyLoadAPackageInCRANorBioconductor("httr")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("biomaRt")

# Get the mouse TF list:
GET("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2836267/bin/NIHMS177825-supplement-02.xls", write_disk(tf <- tempfile(fileext = ".xls")))
ravasiTF <- as.data.frame(read_excel(tf, range = "G2:H1990"))
colnames(ravasiTF) <- c("EntrezGene", "Symbol")
# Remove the uninformative lines
ravasiTF <- unique(subset(ravasiTF, EntrezGene != "---"))

# Preparation for mart requests:
# I use the same version as gtf (93):
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://jul2018.archive.ensembl.org")
mart <- useDataset("mmusculus_gene_ensembl", mart)

entrezgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene"),
                 mart = mart)

# Then using the entrezgene
entrezgene$is_TF <- entrezgene$entrezgene %in% ravasiTF$EntrezGene

ravasiTF$alreadyReported <- ravasiTF$EntrezGene %in% entrezgene$entrezgene

sum(! ravasiTF$alreadyReported)

# 23:
cat("Genes which are ignored:\n")
ravasiTF[! ravasiTF$alreadyReported, ]

# We ignore them and write found ensembl ids in a file.
cat(entrezgene[entrezgene$is_TF, c("ensembl_gene_id")], file = "output/inputs_pc/Ravasi_TF_ensID.txt", sep = "\n")
