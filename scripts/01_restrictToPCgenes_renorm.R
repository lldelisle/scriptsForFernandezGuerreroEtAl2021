options(stringsAsFactors = F)
rm(list = ls())
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
source("scripts/functions.R")

input.folder <- "input/"
# Check there is only one file per pattern
gtf.file <- list.files(input.folder, pattern="gtf", full.names = T)
stopifnot("Was expecting exactly one file containing 'gtf' in its name."=length(gtf.file) == 1)
fpkm.file <- list.files(input.folder, pattern = "cufflinks|fpkm", ignore.case = T, full.names = T)
stopifnot("Was expecting exactly one file containing 'cufflinks' or 'fpkm' in its name."=length(fpkm.file) == 1)
counts.file <- list.files(input.folder, pattern = "count", ignore.case = T, full.names = T)
stopifnot("Was expecting exactly one file containing 'count' in its name."=length(counts.file) == 1)

# Load gtf
gtf.df <- import.gff(gtf.file, format="gtf")
# Get protein.coding.genes
protein.coding.genes <- unique(gtf.df$gene_id[gtf.df$gene_biotype=="protein_coding"])

output.folder <- file.path(input.folder, "../output/inputs_pc/")
if (! dir.exists(output.folder)){
  dir.create(output.folder, recursive = T)
}
# Load fpkm
fpkm.df <- read.delim(fpkm.file)
stopifnot("Was expecting a column named gene_id in fpkm file"="gene_id" %in% colnames(fpkm.df))
sub.fpkm.df <- subset(fpkm.df, gene_id %in% protein.coding.genes)
# Write fpkm
write.table(sub.fpkm.df, file = file.path(output.folder, "FPKM_onlypc.txt"), sep = "\t",
            quote = F, row.names = F)
# Load counts
counts.df <- read.delim(counts.file)
stopifnot("Was expecting a column named gene_id in fpkm file"="gene_id" %in% colnames(counts.df))
sub.counts.df <- subset(counts.df, gene_id %in% protein.coding.genes)
# Write counts
write.table(sub.counts.df, file = file.path(output.folder, "Counts_onlypc.txt"), sep = "\t",
            quote = F, row.names = F)

# Now we will normalize the FPKM using the House Keeping rank method
# We normalize the full fpkm
resl <- normalizeHKRank(fpkm.df, 1000)
# We restrict to protein coding genes
sub.fpkm.norm.df <- subset(resl[["normData"]], gene_id %in% protein.coding.genes)
# Write fpkm norm
write.table(sub.fpkm.norm.df, file = file.path(output.folder, "FPKM_norm_onlypc.txt"), sep = "\t",
            quote = F, row.names = F)
