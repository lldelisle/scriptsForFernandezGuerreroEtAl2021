# First download the data:
# gtf:
wget https://zenodo.org/record/3820860/files/mergeGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.93_ExonsOnly_UCSC.gtf.gz?download=1 -P input/
# FPKM:
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150701&format=file&file=GSE150701%5FAllCufflinks%5FSimplified%2Etxt%2Egz" -O input/GSE150701_AllCufflinks_Simplified.txt.gz
# Counts:
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE150701&format=file&file=GSE150701%5FCounts%2Etxt%2Egz" -O input/GSE150701_Counts.txt.gz
# Prepare the inputs
Rscript scripts/01_restrictToPCgenes_renorm.R
Rscript scripts/02_curatingRavasiTF.R
# Global analysis
# Generate figures 1, 2A, 5A
Rscript scripts/03_plotPCA_Heatmap.R
# Generate figures 2C, 5C
Rscript scripts/04_DESeqHeatmapAndClusteringPerTissue.R
# Run DESeq2:
Rscript scripts/05_DESeqMESECT_eachStage.R
# Plot euler diagrams of figure 6A
Rscript scripts/06_FIG6A_leap.R
# Generate Figure 6B
Rscript scripts/07_FIG6B_Ectodermal_program.R
# Generate SPIA analysis summarized in Figure 7 top left
Rscript scripts/08_FIG7_SPIA-Pathwaysscript.R
mv SPIAPerturbationPlots.pdf output/SPIA/
# Plot heatmaps of list of genes
Rscript scripts/09_plotHeatmap_for_list_of_gene.R
# Plot the individual genes of Figure 7
Rscript scripts/10_plotIndividualGenes_from_list_of_genes.R
# Generate supplementary tables and figures 2B and 5B
Rscript scripts/11_GO_analysis.R
