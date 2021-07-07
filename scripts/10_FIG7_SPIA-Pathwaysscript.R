####### SPIA analysis FIG 7 from Marc et al 2021 ######
####### author: Alejandro Castilla-Ibeas
options(stringsAsFactors = F)
rm(list = ls())
# Load dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("biomaRt")
safelyLoadAPackageInCRANorBioconductor("SPIA")
safelyLoadAPackageInCRANorBioconductor("writexl")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("ggrepel")

# we need to make a database for SPIA with all pathways:
# At the time of analysis, KEGG pathways were obtained by:
# curl -s http://rest.kegg.jp/list/pathway/mmu | awk '{split($1,a,":"); print "curl -s http://rest.kegg.jp/get/"a[2]"/kgml -o extdata/keggxml/mmu/"a[2]".xml"}' | bash
# They are now in input/KEGGPATHWAYS
makeSPIAdata(kgml.path="input/KEGGPATHWAYS/", organism="mmu", out.path=system.file("extdata/",package="SPIA"))

# first we need a vector of all pc genes with entrez ids
expressionDF <- read.table("output/inputs_pc/FPKM_onlypc.txt", header=T) #Expression matrix here
# We use the same version as gtf (93):
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://jul2018.archive.ensembl.org")
mart <- useDataset("mmusculus_gene_ensembl", mart)
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene"),
  values=expressionDF$gene_id,
  mart=mart)
allgenesvect <- as.vector(genes$entrezgene) 


## INPUT
# inputs are results obtained from deseq2 script, wald test ectvsmes at each stage, and filtered by log2FC and pvalue
files <- list.files("output/DESeq2_MESECT_eachStage/", 
                    pattern = "*_DESeq2significant.txt",
                    full.names = T)
names(files) <- gsub("_DESeq2significant.txt", "", basename(files))
# Order:
file.stage <- sapply(strsplit(names(files), "_"), "[", 1)
names(files) <- file.stage
file.stage.order <- paste0("E", sort(as.numeric(gsub("E", "", file.stage))))
files <- files[file.stage.order]
names(files) <- gsub("_DESeq2significant.txt", "", basename(files))

SPIAoutputpath <- "output/SPIA"
if (!dir.exists(SPIAoutputpath)){
  dir.create(SPIAoutputpath, recursive = T)
}
if (! dir.exists("figure")){
  dir.create("figure")
}
if (! dir.exists("additional_files")){
  dir.create("additional_files", recursive = T)
}

excel_list <- list()
for (i in 1:length(files)) {
  form= sprintf('output/SPIA/results_%s.RDS', names(files)[i])
  if (! file.exists(form)){
    ## read strings as character
    tmp <- read.table(files[i], header = T)
    tmp1 <- getBM(filters="ensembl_gene_id", ## retrieve entrez ids
                  attributes=c("ensembl_gene_id", "entrezgene"),
                  values=tmp$gene_id,
                  mart = mart)
    tmp1 <- tmp1[!duplicated(tmp1$entrezgene),]
    tmp1 <- tmp1[!is.na(tmp1$entrezgene),] #remove duplicates and missing genes
    interect <- merge(tmp,tmp1, by.x="gene_id", by.y="ensembl_gene_id") #merge entrez ids to results
    input <- as.vector(interect$log2FoldChange) #extract log2FC
    names(input) <- as.vector(interect$entrezgene) #assing entrezid to log2FC 
    res=spia(de=input,all=allgenesvect, organism="mmu",
             nB=2000,plots=TRUE,beta=NULL,combine="fisher",verbose=FALSE) #Run SPIA
    excel_list[[names(files)[i]]] <- res
    saveRDS(res, form) # save RDS
  } else {
    res = readRDS(form)
    excel_list[[names(files)[i]]] <- res
  }
}
writexl::write_xlsx(excel_list, "additional_files/SPIA_raw_results.xlsx") #save xlsx
all.res <- do.call(rbind, lapply(names(excel_list), function(na){
  df <- excel_list[[na]]
  df$analysis <- na
  return(df)
}))
all.res$analysis <- factor(all.res$analysis, levels = names(excel_list))

pathways_to_check <- c("Wnt", "MapK", "Ras", "Hippo")
names(pathways_to_check) <- c("04310", "04010", "04014", "04390")

my.labeller <- sapply(strsplit(levels(all.res$analysis), "_"), "[[", 1)
names(my.labeller) <- levels(all.res$analysis)

all.res$status <- "NS"
all.res$status[all.res$pGFWER < 0.05 & all.res$tA > 0] <- "MES"
all.res$status[all.res$pGFWER < 0.05 & all.res$tA < 0] <- "ECT"
fixedColors <- list('status'=c('ECT'="#FC7C76",'MES'="#17B7FF",'NS'="grey"))
ggplot(all.res, aes(x = tA, y = - log(pGFWER))) +
  geom_point(aes(color = status)) +
  facet_grid(. ~ analysis, labeller = labeller(analysis = my.labeller)) +
  geom_label_repel(data = subset(all.res, ID %in% names(pathways_to_check)),
                   aes(label = pathways_to_check[ID],
                       color = status),
                   min.segment.length = 0,
                   show.legend = F) +
  theme(legend.position = "bottom") +
  scale_color_manual("",
                     values = fixedColors[["status"]])
ggsave("figure/Fig7_SPIA.pdf", width = 9, height = 3)
