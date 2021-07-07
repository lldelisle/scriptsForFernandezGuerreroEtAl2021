options(stringsAsFactors = F)

# Dependencies
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("DEXSeq")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("plyr")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("writexl")

input.folder <- "input/"
input.dexseq <- file.path(input.folder, "DEXseq")

if (! dir.exists("additional_files")){
  dir.create("additional_files", recursive = T)
}
if (! dir.exists("figure")){
  dir.create("figure")
}

# Get gtfs as GRanges
gtf.file <- list.files(input.folder, pattern="gtf", full.names = T)
stopifnot("Was expecting exactly one file containing 'gtf' in its name."=length(gtf.file) == 1)
gtf <- readGFF(gtf.file)
gtf <- subset(gtf, type == "exon")
gtf.gr <- makeGRangesFromDataFrame(gtf, keep.extra.columns = T)
gtf.dexseq.file <- list.files(input.dexseq, pattern="gtf", full.names = T)
stopifnot("Was expecting exactly one file containing 'gtf' in its name."=length(gtf.dexseq.file) == 1)
gtf.dexseq <- readGFF(gtf.dexseq.file)
gtf.dexseq.gr <- makeGRangesFromDataFrame(subset(gtf.dexseq, type == "exonic_part"), keep.extra.columns = T)

# Reduce the gtf
temp <- reduce(gtf.gr)
ov.red <- as.data.frame(findOverlaps(temp, gtf.gr))

# Compute overlap
ov <- as.data.frame(findOverlaps(gtf.dexseq.gr, gtf.gr))
# Annotate ov
ov$exon.coo <- paste0(seqnames(gtf.gr[ov$subjectHits]), ":",
                      GenomicRanges::start(gtf.gr[ov$subjectHits]),
                      "-", GenomicRanges::end(gtf.gr[ov$subjectHits]))
ov$gene <- gtf.gr$gene_name[ov$subjectHits]
ov$transcript <- gtf.gr$transcript_name[ov$subjectHits]
ov$exon <- gtf.gr$exon_number[ov$subjectHits]
ov$biotype <- gtf.gr$gene_biotype[ov$subjectHits]
ov$exonic_part <- gtf.dexseq.gr$exonic_part_number[ov$queryHits]
ov$gene_id <- gtf.dexseq.gr$gene_id[ov$queryHits]
ov$id <- paste0(gtf.dexseq.gr$gene_id, ":E", gtf.dexseq.gr$exonic_part_number)[ov$queryHits]
ov$reduced.exon <- ov.red$queryHits[match(ov$subjectHits, ov.red$subjectHits)]

# Subset ov to protein_coding:
ov <- subset(ov, biotype == "protein_coding")
countFiles = list.files(input.dexseq, pattern="dexseq", full.names=TRUE)
sampleTable <- read.delim(file.path(input.folder, "samplesplan.txt"))
rownames(sampleTable) <- sampleTable$sample

dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + tissue:exon,
  flattenedfile=gtf.dexseq.file )
# Keep only protein_coding genes:
dxd <- dxd[which(rownames(dxd) %in% ov$id), ]
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="tissue")
dxr1 = DEXSeqResults( dxd )
# Convert output to dataframe
dxr1.sorted <- as.data.frame(dxr1)
dxr1.sorted <- dxr1.sorted[order(dxr1.sorted$padj), ]
dxr1.sorted$rank <- 1:nrow(dxr1.sorted)
# Annotate for the condition
dxr1.sorted$tissue <- "MES"
dxr1.sorted$tissue[dxr1.sorted$log2fold_MES_ECT < 0] <- "ECT"
dxr1.signif <- dxr1.sorted[dxr1.sorted$padj < 0.05 & !is.na(dxr1.sorted$padj), ]
dxr1.signif$transcripts <- sapply(dxr1.signif$transcripts, paste, collapse = ",")
ov.by.id <- ddply(ov, .(id), summarize,
                  transcript.names = paste(unique(sort(transcript)), collapse = ","),
                  gene.names = paste(unique(sort(gene)), collapse = ","),
                  exon.nb = paste(unique(sort(exon)), collapse = ","))
dxr1.signif <- cbind(dxr1.signif, ov.by.id[match(rownames(dxr1.signif), ov.by.id$id), ])
write_xlsx(dxr1.signif[order(dxr1.signif$groupID), ], file.path("additional_files", "DEexonicparts_pc.xlsx"))
# Propagate the annotation to ov
ov$signif <- ov$id %in% dxr1.signif$id
ov$tissue <- dxr1.signif[match(ov$id, dxr1.signif$id), "tissue"]

all.signif <- list()
my.units <- rev(c("gene", "transcript", "reduced.exon", "id"))
my.units.pretty <- rev(c("gene", "transcript", "reduced exon", "exonic part"))
for (unit in my.units){
  all.signif[[unit]] <- ddply(subset(ov, signif), unit, summarize,
                              signif.tissue = paste(sort(unique(tissue)), collapse = ","))
}
my.units.pretty <- paste0(my.units.pretty, "\n(", unlist(lapply(all.signif, nrow))[my.units], ")")
summary <- do.call(rbind, lapply(names(all.signif), function(unit){
  df <- as.data.frame(table(all.signif[[unit]]$signif.tissue))
  df$unit <- unit
  return(df)
}))
# Total number of exonic parts:
nrow(dxr1.sorted)
# [1] 285924
# percentage of significant exonic part:
nrow(all.signif[["id"]]) / nrow(dxr1.sorted) * 100
# [1] 0.4973315

# Order:
summary$unit <- factor(summary$unit, levels = my.units)
summary$Var1 <- factor(summary$Var1, levels = c("ECT", "MES", "ECT,MES"))
summary <- summary[order(summary$Var1, decreasing = T), ] 
summary <- ddply(summary, "unit",
                 transform, label_ypos=cumsum(Freq))
ggplot(summary, aes(x = unit, y = Freq)) +
  geom_bar(aes(fill = Var1), stat = "identity") +
  geom_text(aes(label = Freq, y = label_ypos), vjust=-0.3, size=3.5) +
  ylab("Number") +
  xlab("") +
  scale_fill_manual("tissue",
                    values = c('ECT'="#FC7C76", 'MES'="#17B7FF", 'ECT,MES'='purple'),
                    labels = c("ECT", "MES", "both")) + 
  scale_x_discrete(labels = my.units.pretty) +
  ylim(0, max(summary$label_ypos) * 1.05)

ggsave(file.path("figure", "Fig2A_DEXseq.pdf"), width = 5, height = 5)

pdf(file.path("figure","Fig2B_Fgfr2.pdf"))
plotDEXSeq( dxr1, unique(dxr1.signif$groupID[dxr1.signif$gene.names == "Fgfr2"]), legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , fitExpToVar = "tissue",
            splicing = T, displayTranscripts = T)
dev.off()

saveRDS(object = dxr1, "dexseq.RDS")
