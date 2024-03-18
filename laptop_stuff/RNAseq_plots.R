#!/usr/bin/env Rscript

# check if user requested help menu
args = commandArgs(trailingOnly=TRUE)
#if (length(args) != 1) {
#    stop("Improper number of arguments provided. See 'RNAseq_plots.R -h'", call.=FALSE)
if (args[1] == ("-h")) {
  cat(paste("v1.0.0",
      "RNAseq_plots.R graphs dataframes of interest for RNAseq data mapped by Salmon. Note that a 'mapping' directory containing each sample's Salmon output directory (including the quant.sf file) is required.",
      "",
      "RNAseq_plots requires the following packages, which will be installed upon running:",
      "- BiocManager (from bioconductor.org), which contains the following for differential expression analysis:",
         "- tximport",
         "- edgeR",
         "- limma",
      "For output data management:",
        "- ggplot2",
        "- ggrepel",
        "- tidyr",
      "",
      "Syntax: RNAseq_plots.R <libraries.tsv> <-graph_option> <baseline_name & assessed_name>",
      "",
      "RNAseq_plots.R requires 3 user arguments:",
      "     1. <libraries.tsv> must be a 2-column dataframe of library names and categories (tissue types, developmental stages, etc.) seperated by a tab. Don't include headers.",
      "     A minimum of 3 replicates per category is required. Example dataframe:",
      "         MG_SK1    skin",
      "         MG_SK2    skin",
      "         MG_SK3    skin",
      "         MG_B1    barbel",
      "         MG_B2    barbel",
      "         MG_B3    barbel",
      "",
      "     2. A RNAseq_plots.R graph option (given in <>):",
      "     All graphs <-all>",
      "     Read (CPM) Distribution <-CPMd>",
      "     Multidimensional Scaling <-MDS>",
      "     Biological Coefficient of Variation <-BCV>",
      "     MA plot <-MA>",
      "     Volcano plot of DGE <-Vol>",
      "",
      "     3. A string denoting the baseline library category and the library category for which differential gene expression will be assessed. From the barbel/skin example, 'skin barbel' would make skin the baseline to which up/down regulation of barbel expression is compared.",
      "\n", sep="\n"))
  stop(call.=FALSE)
}

# mapping directory in cwd is required
if (file.exists("mapping") == FALSE) {
  stop("Either your current working directory is incorrect or you did not make a 'mapping' directory. Ensure cwd contains RNAseq_plots.R and a 'mapping' directory; see 'RNAseq_plots.R -h'")
} else { 
    dir <- getwd()
    dir.create("RNAseqPlots_outdata")
    result_dir <- paste0(dir,"/RNAseqPlots_outdata")
}

# Load packages and install as needed
if (require("BiocManager") == FALSE) {source("https://bioconductor.org/install")}
if (require("tximport") == FALSE) {BiocManager::install("tximport"); library(tximport)}
if (require("edgeR") == FALSE) {BiocManager::install("edgeR"); library(edgeR)}
if (require("ggplot2") == FALSE) {install("ggplot2"); library(ggplot2)}
if (require("ggrepel") == FALSE) {install("ggrepel"); library(ggrepel)}
if (require("tidyr") == FALSE) {install("tidyr"); library(tidyr)}
if (require("cowplot") == FALSE) {install("cowplot"); library(tidyr)}
# Print session info so package and R versions are documented
paste("RNAseq_plots.R v1.0.0")
print(sessionInfo())

# Define each plotting function
CPMdist_plot <- function(DGE_object, categories){
  print("-------------------------------------------------")
  print("--------Plotting Library Distributions-----------")
  print("-------------------------------------------------")
  # Calculate log10 CPM for kept transcripts, rename columns with libraries, write long format dataframe
  log10CPM_dist <- data.frame(read_counts=log10(DGE_object$counts))
  colnames(log10CPM_dist) <- DGE_object$samples$samples
  long_format <- gather(log10CPM_dist, library, log10cpm, colnames(log10CPM_dist), factor_key=TRUE)
  # Add in category column HOW???????????
        #conversion <- data.frame(libraries=colnames(log10CPM_dist), groups=categories)
  # Plot library distributions & save plot data table
  write_path <- paste0(result_dir,"/CPMdist_plot.tsv")
  write.table(log10CPM_dist, file=write_path, sep="\t")
  ggplot(long_format, aes(library, log10cpm)) + geom_violin(scale="count") + labs(x="Library", y="log10(CPM)") + theme_cowplot()
  ggsave(filename="CPMdist_plot.jpeg", path=result_dir, device="jpeg", dpi=300)
}

MDS_plot <- function(DGE_object, categories){
  print("-------------------------------------------------")
  print("-------Plotting Multidimensional Scaling---------")
  print("-------------------------------------------------")
  MDS_result <- plotMDS.DGEList(DGE_object, plot=FALSE)
  print("Multidimensional scaling plotted via EdgeR plotMDS.DGEList()")
  mds_data <- data.frame(dimension1=MDS_result$x, dimension2=MDS_result$y, samples=row.names(DGE_object$samples), categories=as.factor(categories))
  # Visualize the MDS plot & save data table
  write_path <- paste0(result_dir,"/mds_data_table.tsv")
  write.table(mds_data, file=write_path, sep="\t")
  ggplot(mds_data, aes(dimension1, dimension2, color=categories)) +  geom_point() + geom_text_repel(size=3, aes(label=samples)) + labs(x="Leading LogFC dim1", y="Leading LogFC dim2", color="Category") + theme_cowplot()
  ggsave(filename="MDS_plot.jpeg", path=result_dir, device="jpeg", dpi=300)
}

BCV_plot <- function(DGE_object, DGEEX_object){
  print("----------------------------------------------------")
  print("----Plotting Biological Coefficient of Variation----")
  print("----------------------------------------------------")
  BCV_result <- data.frame(log2_aveCPM=DGE_object$AveLogCPM, tagWdis_bcv=sqrt(DGE_object$tagwise.dispersion))
  row.names(BCV_result) <- row.names(DGEEX_object$table)
  #Visualize the BCV graph & save data table
  write_path <- paste0(result_dir,"/BCV_result_table.tsv")
  write.table(BCV_result, file=write_path, sep="\t")
  ggplot(BCV_result, aes(log2_aveCPM, tagWdis_bcv)) + geom_hline(yintercept=sqrt(DGE_object$common.dispersion), col="blue") + labs(x="log2(mean CPM)", y="Biological Coefficient of Variation") + geom_point(size=0.50) + theme_cowplot()
  ggsave(filename="BCV_plot.jpeg", path=result_dir, device="jpeg", dpi=300)
}

MA_plot <- function(DGEEX_object, de_id, write_de){
  print("------------------------------------------------")
  print("-----------------Plotting MA--------------------")
  print("------------------------------------------------")
  MA_data <- data.frame(log2FC=DGEEX_object$table$logFC, log2_aveCPM=DGEEX_object$table$logCPM, Pvalues=DGEEX_object$table$PValue, deGenes=as.factor(de_id[,1]))
  row.names(MA_data) <- row.names(DGEEX_object$table)
  # Visualize the MA plot & save data table
  write_path <- paste0(result_dir,"/MA_plot_table.tsv")
  write.table(MA_data, file=write_path, sep="\t")
  ggplot(MA_data, aes(log2_aveCPM, log2FC, color=deGenes)) + labs(x="log2(MeanCPM)", y="Log2(FC)") + geom_point(size=0.50) + geom_hline(yintercept=c(-1, 1), col="black") + scale_color_manual(values=c("red", "black", "blue")) + theme_cowplot() + theme(legend.position="none")
  ggsave(filename="MA_plot.jpeg", path=result_dir, device="jpeg", dpi=300)

  # Write TSV table of differentially expressed transcripts
  tab <- summary(write_de)
  sig_DGE <- tab[1]+tab[3]
  write_path <- paste0(result_dir,"/DGE_table.tsv")
  write.table(topTags(DGEEX_object, n=sig_DGE), file=write_path, sep="\t")
}

Volcano_plot <- function(DGEEX_object, de_id){
  print("---------------------------------------------------")
  print("---------------Plotting Volcano Plot---------------")
  print("---------------------------------------------------")
  adjusted_P <- -log10(DGEEX_object$table$PValue)
  print("Exact Test P values adjusted by log10")
  Volc_data <- data.frame(log2FC=DGEEX_object$table$logFC, log10_pval=adjusted_P, deGenes=as.factor(de_id[,1]))
  print("log2FC taken from logFC column in exactTest() object.")
  print("log10(Pvalue) taken from log10 of Exact Test() object P values.")
  row.names(Volc_data) <- row.names(DGEEX_object$table)
  # Visualize the Volcano plot & save data table
  write_path <- paste0(result_dir,"/Volc_data_table.tsv")
  write.table(Volc_data, file=write_path, sep="\t")
  ggplot(Volc_data, aes(log2FC, log10_pval, color=deGenes)) + labs(x="Log2(FC)", y="-Log10(Pvalue)") + geom_point(size=0.50) + geom_vline(xintercept=c(-1, 1), col="black") + scale_color_manual(values=c("red", "black", "blue")) + theme_cowplot() + theme(legend.position="none")
  ggsave(filename="Volcano_plot.jpeg", path=result_dir, device="jpeg", dpi=300)
}

print("------------------------------")
print("---Dataframe Reading & Prep---")
print("------------------------------")

#verify library-category datafile and identify unique library categories
lib_cats<-read.table(args[1], header=F, row.names=1)
if (ncol(lib_cats) != 1) {
  stop("Improper dataframe submitted; there must be 2 columns only, with no headers. See './RNAseq_plots.R -h'", call.=FALSE)
}
print(paste("Submitted TSV file is:", args[1]))
print(lib_cats)
unique_cat <- unique(lib_cats$V2)
print(unique_cat)
cat_groupfiles <- vector("character")
vec_counter <- 0

# Identify the libraries in each unique category
for (cat in 1:length(unique_cat)) {
  cat_group <- vector("character")
  cat_group <- rownames(lib_cats)[which(lib_cats$V2==unique_cat[cat])]
  # For each library in a given category, assign its file path name to cat_groupfiles
  for (library in 1:length(cat_group)) {
    vec_counter <- vec_counter + 1
    cat_groupfiles[vec_counter] <- file.path("mapping", cat_group[library], "quant.sf")
  }
}
cat_groupfiles <- unlist(cat_groupfiles)

print("--------------------------------------------------")
print("------Assessing Differential Gene Expression------")
print("--------------------------------------------------")

##  read in files with tximport
txi.salmon <- tximport(cat_groupfiles, type="salmon", txOut=TRUE, dropInfReps=TRUE)
sal_counts <- txi.salmon$counts
paste("Dimensions of CPMs per sample table:", dim(sal_counts))
print("Sum of all CPMs per sample:")
print(colSums(sal_counts))

# Filter transcripts according to > 10.0 CPM
kept_transcripts <- rowSums(cpm(sal_counts)>10.0) >= length(lib_cats$V2)
sal_counts <- sal_counts[kept_transcripts,]
paste("Dimensions of filtered CPMs per sample table:", dim(sal_counts))
print("Sum of CPMs per sample after transcripts with < 10.0 CPM are removed:")
print(colSums(sal_counts))

# Make DGEList object, make normalized libraries
DGE_obj <- DGEList(counts=sal_counts, group=lib_cats$V2, samples=row.names(lib_cats))
row.names(DGE_obj$samples) <- row.names(lib_cats)
DGE_norm <- calcNormFactors(DGE_obj)
print("Libraries' CPM size normalized via EdgeR calcNormFactors()")
# calculate tagwise dispersion. Note prior.n = 50/(#libraries - #categories)
DGE_comDis <- estimateCommonDisp(DGE_norm)
print("Common dispersion calculated with EdgeR estimateCommonDisp()")
DGE_tagWdisp <- estimateTagwiseDisp(DGE_comDis, prior.n=(50 / (length(lib_cats$V2) - length(unique(lib_cats$V2)))))
print("Tagwise dispersion calculated with EdgeR estimateTagwiseDisp() using prior.n = 50 / (# of samples — # of categories)")
  
# Identify differentially expressed sequences
et_result <- exactTest(DGE_tagWdisp, pair=c(args[3],args[4]))
print("Exact test calculated with EdgeR exactTest() using tagwise dispersion object and user's category string")
#QUOTE on exactTest() comparisons:
#"Note that the first group listed in the pair is the baseline for the comparison—so if the pair is c("A","B") then the comparison is B - A, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change)."
de <- decideTestsDGE(et_result, p=0.05, adjust="BH")
print("Up/down regulation interpreted from exact test via EdgeR decideTestsDGE(p=0.05, adjust='BH')")
df_de <- as.data.frame(de)
paste("As per user input, positive/negative fold changes will be up/down regulated in", args[4], "relative to", args[3])

# Run plotting functions according to user input
if (args[2] == ("-all")) {
  CPMdist_plot(DGE_tagWdisp, lib_cats$V2)
  MDS_plot(DGE_tagWdisp, lib_cats$V2)
  BCV_plot(DGE_tagWdisp, et_result)
  MA_plot(et_result, df_de, de)
  Volcano_plot(et_result, df_de)
  # Code for single-graph requests
} else if (args[2] == ("-CPMd")) {
  CPMdist_plot(DGE_tagWdisp, lib_cats$V2)
} else if (args[2] == ("-MDS")) {
  MDS_plot(DGE_tagWdisp, lib_cats$V2)
} else if (args[2] == ("-BCV")) {
  BCV_plot(DGE_tagWdisp, et_result)
} else if (args[2] == ("-MA")) {
  MA_plot(et_result, df_de, de)
} else if (args[2] == ("-Vol")) {
  Volcano_plot(et_result, df_de)
} 
