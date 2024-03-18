args = commandArgs(trailingOnly=TRUE)

if (args[1] == ("-h")) {
  cat(paste("1. current directory",
            "2. species of interest",
            "3. x-axis label cutoff",
            "4. optional (T/F): no plot outputs (default F)",
            "\n", 
            sep="\n"))
  stop(call.=FALSE)
}

print("------------------------------------------------------")
print("------------------------setup-------------------------")
print("------------------------------------------------------")

if (require("BiocManager") == FALSE) {source("https://bioconductor.org/install")}

if (require("biomartr") == FALSE) {BiocManager::install("biomartr"); library(biomartr)}
if (require("cogeqc") == FALSE) {BiocManager::install("cogeqc"); library(cogeqc)}
if (require("dplyr") == FALSE) {BiocManager::install("dplyr"); library(dplyr)}
if (require("readr") == FALSE) {BiocManager::install("readr"); library(readr)}
if (require("stringr") == FALSE) {BiocManager::install("stringr"); library(stringr)}
if (require("ggplot2") == FALSE) {BiocManager::install("ggplot2"); library(ggplot2)}
if (require("ggrepel") == FALSE) {BiocManager::install("ggrepel"); library(ggrepel)}
if (require("cowplot") == FALSE) {BiocManager::install("cowplot"); library(cowplot)}

if (require("data.table") == FALSE) {BiocManager::install("data.table")}

# setting globals
JUST_STATS <- ifelse(args[4] == "" | args[4] != "T", FALSE, TRUE)
cd <- args[1]
# cd <- getwd()
SPECIES_INTEREST = args[2]
# SPECIES_INTEREST = "E_deani.TRINITY.hit.fasta.transdecoder_fixed_headers.pep"

print("------------------------------------------------------")
print("-----------------grabbing orthogroups-----------------")
print("------------------------------------------------------")

alpha_seqs <- read.csv("~/Desktop/Orthofinder/alpha_seqs.txt", header = FALSE)
# use cogeqc to grab orthogroups
#orthogroups <- read_orthogroups("~/Desktop/Orthofinder/Orthogroups_OLD_RUN.tsv")
orthogroups <- read_orthogroups("~/Desktop/Orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

# select the human orthogroups
human_orthogroups <- orthogroups |> 
  filter(Species == "Homo_sapiens.pep")

# select current species orthogroups
species_orthogroups <- orthogroups |> 
    filter(Species == SPECIES_INTEREST) |> 
    mutate(Sequence = sub("_[^_]*$", "", Gene))
remove(orthogroups)


print("------------------------------------------------------")
print("----------grabbing gene symbols from biomart----------")
print("------------------------------------------------------")

# ### Helper functions to find correct mart/ attributess ###
# # show all elements of the data.frame
# options(tibble.print_max = Inf)
# # return available attributes for "Homo sapiens"
# marts <- biomartr::organismAttributes("Homo sapiens")
# 
# marts |> 
#     filter(name == "pdb")

BiomartResults <- biomartr::biomart( genes      = as.vector(human_orthogroups$Gene), # genes were retrieved using biomartr::getGenome()
                                     mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                     dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                     attributes = c("pfam","pdb", "hgnc_symbol"), # attributes were selected with biomartr::getAttributes()
                                     filters    = "ensembl_peptide_id_version")

colnames(human_orthogroups)[3] <- "ensembl_peptide_id_version"
human_orthogroups <-  left_join(human_orthogroups, BiomartResults, by="ensembl_peptide_id_version")
remove(BiomartResults)

# orthogroups_groupd <- orthogroups |> 
#     dplyr::filter(Species == SPECIES_INTEREST) |> 
#     group_by(Orthogroup, Species) |> 
#     summarise(genes = list(Gene)) 

print("------------------------------------------------------")
print("------output of RNAseqPlots and orthogroup data-------")
print("------------------------------------------------------")

#RNASeq_plots headers are set weird, this removes them and adds them back nicely
eg_deg <- read_tsv(paste0(cd, "/RNAseqPlots_outdata/MA_plot_table.tsv"), skip=1, col_names=FALSE)
colnames(eg_deg) <- c("Sequence", "log2FC", "log2_aveCPM", "Pvalues", "deGenes")

# output de + orthogroup species of interest tsv with formatted name
species_out <- left_join(eg_deg, species_orthogroups, by = "Sequence")
remove(eg_deg)
remove(species_orthogroups)
species_out <- left_join(species_out, human_orthogroups, by = "Orthogroup")
remove(human_orthogroups)
temp_name <- sub("\\..*", "", SPECIES_INTEREST)

# output sequences will all info
print("outputing sequences with labeling data")
write_tsv(species_out, paste0(cd, "/", temp_name, "_diff_expression_orthogroup.tsv"))

# output missing genes
only_missing_out <- species_out |> 
    filter(is.na(hgnc_symbol) | hgnc_symbol == "") |>
    select(Sequence)
print("outputing sequences with missing gene symbols")
write_tsv(only_missing_out, paste0(cd, "/", temp_name, "_missing_genes.tsv"), col_names = FALSE)
remove(temp_name)

# format dataframe for correct labeling
plot_frame <- species_out |> 
    group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |> 
    summarize(gene_symbol = paste(hgnc_symbol, collapse = ", "))
plot_frame$deGenes <- as.factor(plot_frame$deGenes)


# end program here if only stats are wanted
if (JUST_STATS) {
  print("Skipping plot output")
  quit("no", status = 0)
}
print("------------------------------------------------------")
print("---------------labeled MA plot output-----------------")
print("------------------------------------------------------")

VALUES = c("-1"="red", "0"="grey", "1"="blue", "3"="black", "4"="purple", "5"="green")
X_LAB = as.integer(args[3])

# confine labels to specific range and length or else will block out screen
species_out |> 
    group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |> 
    summarize(gene_symbol = paste(hgnc_symbol, collapse = ", ")) |> 
    mutate(deGenes = ifelse((gene_symbol != "NA" && deGenes == 1), 4, deGenes)) |> # take out later
    mutate(deGenes = ifelse(Sequence %in% alpha_seqs$V1, 5, deGenes)) |> # uses MSA to label DE genes
    mutate(deGenes = as.factor(deGenes)) |> 
    mutate(gene_symbol = ifelse(log2_aveCPM >= X_LAB, gene_symbol, NA)) |> 
    mutate(gene_symbol = ifelse(deGenes == 1 || deGenes == 4, gene_symbol, NA)) |> # take out later
    #mutate(gene_symbol = ifelse(deGenes != 0, gene_symbol, NA)) |> 
    #mutate(gene_symbol = ifelse(nchar(gene_symbol) > 15, paste0(substr(gene_symbol, 1, 20), "..."), gene_symbol)) |> 
    mutate(gene_symbol = strsplit(gene_symbol, ",")[[1]][1]) |> 
    mutate(gene_symbol = ifelse(gene_symbol == "NA", NA, gene_symbol)) |> # take out later
    ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) + 
        geom_point(size=0.50, alpha=0.75) + 
        #geom_vline(xintercept=c(-1, 1), col="black") +
        geom_text_repel(aes(label=gene_symbol),
                        size=2.5,
                        max.overlaps = Inf,
                        segment.color = NA,
                        na.rm = TRUE) +
        labs(x="Mean Transcript Abundance (log CPM)", y="Differential Expression (Log Fold Change)") + 
        scale_color_manual(values=VALUES) + 
        theme_cowplot() + 
        theme(legend.position="none")
ggsave(filename=paste0(cd, "/Labeled_MA.pdf"), device="pdf", dpi=300)

print("------------------------------------------------------")
print("----------sequence labeled MA plot output-------------")
print("------------------------------------------------------")

species_out |> 
  group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |> 
  summarize(gene_symbol = paste(hgnc_symbol, collapse = ", ")) |> 
  mutate(deGenes = ifelse((gene_symbol != "NA" && deGenes == 1), 4, deGenes)) |> # take out later
  mutate(deGenes = ifelse(Sequence %in% alpha_seqs$V1, 5, deGenes)) |> # uses MSA to label DE genes
  mutate(deGenes = as.factor(deGenes)) |> 
  mutate(Sequence = ifelse(deGenes == 1 | deGenes == 4, Sequence, NA)) |> 
  mutate(Sequence = ifelse(log2_aveCPM >= X_LAB, Sequence, NA)) |> 
  #mutate(gene_symbol = ifelse(deGenes == 1 || deGenes == 4, gene_symbol, NA)) |> # take out later
  #mutate(gene_symbol = ifelse(deGenes != 0, gene_symbol, NA)) |> 
  #mutate(gene_symbol = ifelse(nchar(gene_symbol) > 15, paste0(substr(gene_symbol, 1, 20), "..."), gene_symbol)) |> 
  #mutate(gene_symbol = strsplit(gene_symbol, ",")[[1]][1]) |> 
  #mutate(gene_symbol = ifelse(gene_symbol == "NA", NA, gene_symbol)) |> # take out later
  ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) + 
  geom_point(size=0.50, alpha=0.75) + 
  #geom_vline(xintercept=c(-1, 1), col="black") +
  geom_text_repel(aes(label=Sequence),
                  size=2.5,
                  max.overlaps = Inf,
                  segment.color = NA,
                  na.rm = TRUE) +
  labs(x="Mean Transcript Abundance (log CPM)", y="Differential Expression (Log Fold Change)") + 
  scale_color_manual(values=VALUES) + 
  theme_cowplot() + 
  theme(legend.position="none")
ggsave(filename=paste0(cd, "/Sequence_Labeled_MA.pdf"), device="pdf", dpi=300)

print("------------------------------------------------------")
print("---------labeled MA plot only slime output------------")
print("------------------------------------------------------")

species_out |> 
    filter(log2FC > 0) |> 
    group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |> 
    summarize(gene_symbol = paste(hgnc_symbol, collapse = ", ")) |> 
    mutate(deGenes = as.factor(deGenes)) |> 
    mutate(gene_symbol = ifelse(log2_aveCPM >= X_LAB, gene_symbol, NA)) |> 
    mutate(gene_symbol = ifelse(deGenes != 0, gene_symbol, NA)) |> 
    mutate(gene_symbol = ifelse(nchar(gene_symbol) > 15, paste0(substr(gene_symbol, 1, 20), "..."), gene_symbol)) |> 
    ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) + 
        geom_point(size=0.50) + 
        #geom_vline(xintercept=c(-1, 1), col="black") + 
        geom_text_repel(aes(label=gene_symbol), 
                        size=1,
                        max.overlaps = Inf,
                        segment.color = NA,
                        na.rm = TRUE) +
        #labs(x="Log2(FC)", y="log2 average CPM") + 
        scale_color_manual(values=VALUES) + 
        theme_cowplot() + 
        theme(legend.position="none")
ggsave(filename=paste0(cd, "/Labeled_MA_slime.pdf"), device="pdf", dpi=300)

print("------------------------------------------------------")
print("-----------labeled MA plot int fil output-------------")
print("------------------------------------------------------")

# HGNC intermediate filament list
int_fils <- read_delim("~/Desktop/Orthofinder/group-607.csv", delim=",")
int_fils <- int_fils |> 
    select(`Approved symbol`)

# filaments
species_out |> 
    #filter(log2FC > 0) |> 
    mutate(deGenes = ifelse(hgnc_symbol %in% int_fils$`Approved symbol` & deGenes == 0, 2, deGenes)) |>
    mutate(deGenes = ifelse(hgnc_symbol %in% int_fils$`Approved symbol` & deGenes != 0, 3, deGenes)) |>
    group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |> 
    summarize(gene_symbol = paste(hgnc_symbol, collapse = ", ")) |> 
    mutate(deGenes = as.factor(deGenes)) |> 
    #mutate(gene_symbol = ifelse(log2_aveCPM >= X_LAB, gene_symbol, NA)) |> 
    mutate(gene_symbol = ifelse(deGenes == 2 | deGenes == 3, gene_symbol, NA)) |> 
    mutate(gene_symbol = ifelse(nchar(gene_symbol) > 15, paste0(substr(gene_symbol, 1, 20), "..."), gene_symbol)) |> 
    ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) + 
        geom_point(size=0.50) + 
        #geom_vline(xintercept=c(-1, 1), col="black") + 
        geom_text_repel(aes(label=gene_symbol), 
                        size=1,
                        max.overlaps = Inf,
                        segment.color = NA,
                        na.rm = TRUE) +
        #labs(x="Log2(FC)", y="log2 average CPM") + 
        scale_color_manual(values=VALUES) + 
        theme_cowplot() 
        #theme(legend.position="none")
ggsave(filename=paste0(cd, "/Labeled_MA_inter_fils.pdf"), device="pdf", dpi=300)

print("------------------------------------------------------")
print("--------labeled MA plot glycoprotein output-----------")
print("------------------------------------------------------")

# glycoproteins https://beta.glycosmos.org/glycoproteins/index?page=2
glyco <- read_delim("~/Desktop/Orthofinder/glycoproteinDB.csv", delim=",")
glyco <- glyco |> 
    select(`Protein Name`, `Gene Symbol`, Species, `No. of Glycosylation Sites`) |> 
    filter(grepl("Homo sapiens", Species)) |> 
    rename(hgnc_symbol = `Gene Symbol`) |> 
    rename(Prot_Name = `Protein Name`) |> 
    rename(glycosylation_sites = `No. of Glycosylation Sites`) |> 
    group_by(hgnc_symbol, Species, glycosylation_sites) |> 
    summarize(Prot_Name = paste(Prot_Name, collapse = ", "))

# removes " [ ] characters
glyco$hgnc_symbol <- gsub("\\[|\\]|\"", "", glyco$hgnc_symbol)

test <- left_join(species_out, glyco, by = "hgnc_symbol")

test |> 
    mutate(deGenes = as.factor(deGenes)) |> 
    group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |>
    summarize(glycosylation_sites = paste(glycosylation_sites, collapse = ", ")) |> 
    mutate(glycosylation_sites = ifelse(log2_aveCPM >= X_LAB, glycosylation_sites, NA)) |> 
    mutate(glycosylation_sites = ifelse(deGenes == -1 | deGenes == 1, glycosylation_sites, NA)) |> 
    #mutate(deGenes = ifelse(is.na(Prot_Name), deGenes, 2)) |> 
    mutate(glycosylation_sites = ifelse(nchar(glycosylation_sites) > 10, 
                                        paste0(substr(glycosylation_sites, 1, 10), "..."), 
                                        glycosylation_sites)) |> 
    ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) + 
        geom_point(size=0.50) + 
        #geom_vline(xintercept=c(-1, 1), col="black") + 
        geom_text_repel(aes(label=glycosylation_sites), 
                        size=1,
                        max.overlaps = Inf,
                        segment.color = NA,
                        na.rm = TRUE) +
        #labs(x="Log2(FC)", y="log2 average CPM") + 
        scale_color_manual(values=VALUES) + 
        theme_cowplot() 
        #theme(legend.position="none")
ggsave(filename=paste0(cd, "/Labeled_MA_glycosylation_sites.pdf"), device="pdf", dpi=300)

mucins <- read_delim("~/Desktop/Orthofinder/group-648.csv", delim=",")
mucins <- mucins |>
    select(`Approved symbol`)

print("------------------------------------------------------")
print("------------labeled MA plot mucins output-------------")
print("------------------------------------------------------")

# mucins
species_out |>
    #filter(log2FC > 0) |>
    mutate(deGenes = ifelse(hgnc_symbol %in% mucins$`Approved symbol` & deGenes == 0, 2, deGenes)) |>
    mutate(deGenes = ifelse(hgnc_symbol %in% mucins$`Approved symbol` & deGenes != 0, 3, deGenes)) |>
    group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |>
    summarize(gene_symbol = paste(hgnc_symbol, collapse = ", ")) |>
    mutate(deGenes = as.factor(deGenes)) |>
    #mutate(gene_symbol = ifelse(log2_aveCPM >= X_LAB, gene_symbol, NA)) |>
    mutate(gene_symbol = ifelse(deGenes == 2 | deGenes == 3, gene_symbol, NA)) |>
    mutate(gene_symbol = ifelse(nchar(gene_symbol) > 15, paste0(substr(gene_symbol, 1, 20), "..."), gene_symbol)) |>
    ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) +
        geom_point(size=0.50) +
        #geom_vline(xintercept=c(-1, 1), col="black") +
        geom_text_repel(aes(label=gene_symbol),
                        size=1,
                        max.overlaps = Inf,
                        segment.color = NA,
                        na.rm = TRUE) +
        #labs(x="Log2(FC)", y="log2 average CPM") +
        scale_color_manual(values=VALUES) +
        theme_cowplot()
        #theme(legend.position="none")
ggsave(filename=paste0(cd, "/Labeled_MA_mucins.pdf"), device="pdf", dpi=300)

# # pfam
# species_out |> 
#     filter(log2FC > 0) |> 
#     mutate(deGenes = ifelse(pfam == "PF00038", 2, deGenes)) |> 
#     group_by(Sequence, log2FC, log2_aveCPM, deGenes, Orthogroup) |> 
#     summarize(pfam = paste(pfam, collapse = ", ")) |> 
#     mutate(deGenes = as.factor(deGenes)) |> 
#     mutate(pfam = ifelse(log2_aveCPM >= 0, pfam, NA)) |> 
#     mutate(pfam = ifelse(deGenes != 0, pfam, NA)) |> 
#     mutate(pfam = ifelse(nchar(pfam) > 40, paste0(substr(pfam, 1, 20), "..."), pfam)) |> 
#     ggplot(aes(y=log2FC, x=log2_aveCPM, color=deGenes)) + 
#         geom_point(size=0.50) + 
#         #geom_vline(xintercept=c(-1, 1), col="black") + 
#         geom_text_repel(aes(label=pfam), 
#                         size=1,
#                         max.overlaps = Inf,
#                         segment.color = NA) +
#         #labs(x="Log2(FC)", y="log2 average CPM") + 
#         scale_color_manual(values=c("black", "blue", "purple")) + 
#         theme_cowplot() + 
#         theme(legend.position="none")
# #ggsave(filename="Labeled_Volcano_slime.pdf", device="pdf", dpi=300)
