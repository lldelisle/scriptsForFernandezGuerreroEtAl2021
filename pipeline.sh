# First download the data:
# gtf:
wget https://zenodo.org/record/3820860/files/mergeGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.93_ExonsOnly_UCSC.gtf.gz?download=1 -P input/
# FPKM:
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150701&format=file&file=GSE150701%5FAllCufflinks%5FSimplified%2Etxt%2Egz" -O input/GSE150701_AllCufflinks_Simplified.txt.gz
# Counts:
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150701&format=file&file=GSE150701%5FCounts%2Etxt%2Egz" -O input/GSE150701_Counts.txt.gz
# All dexseq_counts:
i=4556930
for tissue in MES ECT; do
  for stage in 09 10 11 12; do
    for rep in 1 2; do
      wget "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4556nnn/GSM${i}/suppl/GSM${i}%5FdexseqCount%5F${tissue}${stage}5%5FR${rep}%2Etxt%2Egz" -O input/DEXseq/dexseqCount_${tissue}${stage}5_R${rep}.txt.gz
      gunzip input/DEXseq/dexseqCount_${tissue}${stage}5_R${rep}.txt.gz
      i=$((i + 1))
    done
  done
done

# Prepare the inputs
Rscript scripts/01_restrictToPCgenes_renorm.R
Rscript scripts/02_curatingRavasiTF.R

# Global analysis
# Generate figures 1, 3A, 5A
Rscript scripts/03_plotPCA_Heatmap_GO.R
# Generate figure 2A and 2B
Rscript scripts/04_DEXseq.R
# Generate figures 3C, 5C
Rscript scripts/05_DESeqHeatmapAndClusteringPerTissue.R
# Generate supplementary tables and figures 2B and 5B
Rscript scripts/06_GO_analysis.R
# Run DESeq2 at each stage between MES and ECT:
Rscript scripts/07_DESeqMESECT_eachStage.R
# Plot euler diagrams of figure 6A
Rscript scripts/08_FIG6A_leap.R
# Generate Figure 6B
Rscript scripts/09_FIG6B_Ectodermal_program.R
# Generate SPIA analysis summarized in Figure 7A
Rscript scripts/10_FIG7_SPIA-Pathwaysscript.R
mv SPIAPerturbationPlots.pdf output/SPIA/
# Plot heatmaps of list of genes
Rscript scripts/11_plotHeatmap_for_list_of_gene.R
# Plot the individual genes of Figure 7
Rscript scripts/12_plotIndividualGenes_from_list_of_genes.R
