#!/usr/bin/env Rscript

### ARGUMENT PARSING ###
if (suppressMessages(require("argparser")) == FALSE) {
    install.packages("argparser")
    library(argparser, quietly = TRUE)
}

parser <- arg_parser(
    "RNAseq Differential Expression Analysis and Plotting",
    name = "RNAdeap",
    hide.opts = TRUE
)

parser <- add_argument(
    parser,
    "-i",
    help = "input directory that contains: 1. libraries.tsv sorted by condition 2. mapping directory with each quant.sf file stored in its own directory",
    type = "character"
)
parser <- add_argument(
    parser,
    "-f",
    help = "filter libraries.tsv by this imput and determine differential expression. Example: 'skin,slime' would make skin the baseline to which up/down regulation of slime expression is compared",
    type = "character"
)
parser <- add_argument(parser, "-p", help = "plot output", flag = TRUE)
parser <- add_argument(
    parser,
    "-e",
    help = "EnTAP input from \"entap_results.tsv\"",
    type = "character"
)
parser <- add_argument(
    parser,
    "-l",
    help = "label plots based on Entap Input",
    flag = TRUE
)
parser <- add_argument(
    parser,
    "-x",
    help = "x axis label cutoff if -l selected",
    type = "integer",
    default = 0
)

args <- parse_args(parser)
rm(parser)

### REQUIRED LIBRARIES ###
cat("\nLoading Libraries\n")
# Bioconductor #
if (suppressMessages(require("BiocManager")) == FALSE) {
    source("https://bioconductor.org/install")
}
# Differential Analysis #
if (suppressMessages(require("tximport")) == FALSE) {
    BiocManager::install("tximport")
    library(tximport)
}
if (suppressMessages(require("edgeR")) == FALSE) {
    BiocManager::install("edgeR")
    library(edgeR)
}
# Plotting/Labeling #
if (suppressMessages(require("glue")) == FALSE) {
    BiocManager::install("glue")
    library(glue)
}
if (suppressMessages(require("dplyr")) == FALSE) {
    BiocManager::install("dplyr")
    library(dplyr)
}
if (suppressMessages(require("readr")) == FALSE) {
    BiocManager::install("readr")
    library(readr)
}
if (suppressMessages(require("tidyr")) == FALSE) {
    BiocManager::install("tidyr")
    library(tidyr)
}
if (suppressMessages(require("purrr")) == FALSE) {
    BiocManager::install("purrr")
    library(purrr)
}
if (suppressMessages(require("stringr")) == FALSE) {
    BiocManager::install("stringr")
    library(stringr)
}
if (suppressMessages(require("ggplot2")) == FALSE) {
    BiocManager::install("ggplot2")
    library(ggplot2)
}
if (suppressMessages(require("cowplot")) == FALSE) {
    BiocManager::install("cowplot")
    library(cowplot)
}
if (suppressMessages(require("ggrepel")) == FALSE) {
    BiocManager::install("ggrepel")
    library(ggrepel)
}
if (suppressMessages(require("knitr")) == FALSE) {
    BiocManager::install("knitr")
    library(knitr)
}
if (suppressMessages(require("reshape2")) == FALSE) {
    BiocManager::install("reshape2")
    library(reshape2)
}
if (suppressMessages(require("extrafont")) == FALSE) {
    BiocManager::install("extrafont")
    library(extrafont)
}

# ###TODO: Testing remove
# args <- data.frame("1"=FALSE,
#                    "help"=FALSE,
#                    "p"=TRUE,
#                    "l"=TRUE,
#                    "i"="/Users/kale/Desktop/salmon_expression_filtered/M_phantasma",
#                    "e"="/Users/kale/Desktop/EnTAP/M_phantasma/final_results/entap_results.tsv",
#                    "x"=9,
#                    f="skin,slime")

custom_theme <- function(
    base_size = 12,
    base_family = "Arial",
    line_size = 0.5
) {
    ret_theme <- theme_minimal(
        base_size = base_size,
        base_family = base_family
    ) +
        theme(
            text = element_text(family = base_family, size = base_size),
            axis.title = element_text(size = base_size * 1.2),
            axis.text = element_text(size = base_size * 0.8),
            axis.line = element_line(size = line_size, colour = "black"),
            axis.ticks = element_line(linewidth = line_size),
            panel.grid.major = element_blank(), #element_line(linewidth = line_size, colour = "grey80"),
            panel.grid.minor = element_blank(), #element_line(linewidth = line_size / 2, colour = "grey90"),
            legend.title = element_text(size = base_size),
            legend.text = element_text(size = base_size * 0.8)
        ) +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title.position = "left",
            legend.title.align = 0.5,
            legend.text.position = "top",
            panel.border = element_blank(),

            #plot.title=element_text(size=20, vjust=2, family="Arial"),
            axis.title.x = element_text(size = 18, family = "Arial"),
            axis.title.y = element_text(size = 18, family = "Arial"),
            axis.text.x = element_text(size = 12, family = "Arial"),
            axis.text.y = element_text(size = 12, family = "Arial"),
            axis.line = element_line(
                colour = "black",
                size = 1,
                linetype = "solid",
                lineend = "round"
            )
            #axis.text = element_text(size = 12)
        )
    return(ret_theme)
}

### OUTPUT DIRECTORY CREATION ###
cat("Creating Output Directories\n")
arg_filter <- str_split(args$f, ",")[[1]]
output_dir <- paste0(
    args$i,
    "/RNAdeap_output_",
    arg_filter[1],
    "-",
    arg_filter[2],
    "/"
)
output_dir_fig <- paste0(output_dir, "figures/")
output_dir_data <- paste0(output_dir, "data/")
output_dir_list <- paste0(output_dir, "list/")
dir.create(output_dir, showWarnings = FALSE)
dir.create(output_dir_fig, showWarnings = FALSE)
dir.create(output_dir_data, showWarnings = FALSE)
dir.create(output_dir_list, showWarnings = FALSE)
file_start <- paste0(basename(args$i), "_")
rm(output_dir, arg_filter)

### DATA IMPORT ###
cat("\nImporting Data\n")
# libraries.tsv input
cat("Reading libraries.tsv\n")
lib <- list.files(path = args$i, pattern = "libraries.tsv", full.names = TRUE)
lib <- read_tsv(lib, col_names = FALSE, show_col_types = FALSE)
colnames(lib) <- c("sample", "condition")
if (ncol(lib) != 2) {
    stop(
        "impropper libraries.tsv file, see help message for propper formatting",
        call. = FALSE
    )
}
lib$ind <- sub("_.*", "", lib$sample)
filter <- str_split(args$f, ",")[[1]]
lib <- lib |>
    filter(condition %in% filter)
rm(filter)

# table of input data
cat("Libraries Information:\n")
lib_tab <- lib |>
    group_by(condition) |>
    summarise(n = n(), conditions = paste(sample, collapse = ","))
names(lib_tab) <- c("Condition", "Count", "Reads")
kable(lib_tab)
cat("\n")
rm(lib_tab)

# quant.sf input
cat("Finding quant.sf files\n")
lib["path"] <- NA
for (i in 1:nrow(lib)) {
    directory <- dir(
        path = args$i,
        pattern = lib$sample[i],
        full.names = TRUE,
        recursive = TRUE,
        include.dirs = TRUE
    )
    file <- list.files(
        path = directory,
        pattern = "quant.sf",
        full.names = TRUE
    )
    if (is.na(file)) {
        stop(
            glue("Could not find quant.sf file for {lib$sample[i]}"),
            call. = FALSE
        )
    }
    lib$path[i] <- file
}
rm(i, directory, file)

### DIFFERENTIAL EXPRESSION ANALYSIS ###
cat("\nImporting Expression Data from quant.sf files using tximport\n")
# tximport salmon in
files <- lib$path
names(files) <- lib$sample
dge_import <- tximport(files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
rm(files)

# DGE #
# melt(dge_import$counts) |>
#     ggplot(aes(x=value)) +
#     geom_histogram(bins=200) +
#     facet_wrap(~Var2, scales="free_x")

# edgeR input
cat("Differential Expression Analysis with edgeR\n")
DGE <- DGEList(
    counts = dge_import$counts,
    group = lib$condition,
    samples = lib$sample
)
rm(dge_import)
DGE_full <- DGE

# filter by expression
cat("Filtering by Expression\n")
#filter <- rowSums(cpm(DGE$counts)>10.0) >= length(lib$condition) # extreme filter
filter <- filterByExpr(DGE, group = DGE$samples$group)
DGE <- DGE[filter, , keep.lib.sizes = FALSE]
rm(filter)

# print stats
stats_file <- paste0(output_dir_data, "stats.txt")
cat("Outputing Stats to: ", stats_file, "\n", sep = "")
cat(
    "Full Dimensions: ",
    nrow(DGE_full$counts),
    "\n",
    sep = "",
    file = stats_file,
    append = FALSE
)
cat(
    "Filtered Dimensions:",
    nrow(DGE$counts),
    "\n",
    sep = "",
    file = stats_file,
    append = TRUE
)
cat("Full Library sizes:\n", file = stats_file, append = TRUE)
cat(
    paste0(kable(colSums(DGE_full$counts)), collapse = "\n"),
    "\n\n",
    file = stats_file,
    append = TRUE
)
cat("Filtered Library sizes:\n", file = stats_file, append = TRUE)
cat(
    paste0(kable(colSums(DGE$counts)), collapse = "\n"),
    "\n\n",
    file = stats_file,
    append = TRUE
)
rm(DGE_full)

# normalize
cat("Normalizing\n")
DGE <- calcNormFactors(DGE)
cat("Final Samples:\n", file = stats_file, append = TRUE)
cat(
    paste0(kable(DGE$samples), collapse = "\n"),
    "\n\n",
    file = stats_file,
    append = TRUE
)


# experimental design
cat("Experimental Design Setup\n")
cat("Experimental Design:\n", file = stats_file, append = TRUE)
cat(
    paste0(
        lib |> dplyr::select(sample, ind, condition) |> kable(),
        collapse = "\n"
    ),
    "\n\n",
    file = stats_file,
    append = TRUE
)
design <- model.matrix(~ lib$ind + lib$condition)
rownames(design) <- colnames(DGE)
write.csv(design, paste0(output_dir_data, file_start, "design.csv"))

# dispersion estimation
cat("Dispersion Estimation\n")
DGE <- estimateDisp(DGE, design) # common and tagwise in one
et <- exactTest(DGE, pair = str_split(args$f, ",")[[1]])
de <- decideTests.DGEExact(et, p.value = 0.005, lfc = 1)
de <- data.frame(de)
top_degs = topTags(object = et, n = "Inf")
cat("Diffential Gene Expression:\n", file = stats_file, append = TRUE)
cat(
    paste0(
        summary(decideTests(object = et, p.value = 0.005, lfc = 1)) |>
            kable(),
        collapse = "\n"
    ),
    "\n\n",
    file = stats_file,
    append = TRUE
)

### OUTPUT LISTS ###
if (!is.na(args$e)) {
    # EnTAP Labeling
    cat("Labeling with EnTAP\n")
    entap <- read_tsv(args$e, col_names = TRUE, show_col_types = FALSE)
    cat(
        "EnTAP total raw sequences: ",
        nrow(entap),
        "\n",
        sep = "",
        file = stats_file,
        append = TRUE
    )
    cat(
        "EnTAP In Frame sequences: ",
        nrow(entap |> filter(!is.na(Frame))),
        "\n",
        sep = "",
        file = stats_file,
        append = TRUE
    )
    cat(
        "EggNOG predicted genes: ",
        nrow(entap |> filter(!is.na(`EggNOG Predicted Gene`))),
        "\n",
        sep = "",
        file = stats_file,
        append = TRUE
    )

    entap_filtered <- entap |>
        filter(!is.na(Frame)) |>
        filter(!is.na(`EggNOG Predicted Gene`)) |>
        dplyr::select(`Query Sequence`, `EggNOG Predicted Gene`)

    de_genes <- data.frame(
        names = rownames(et$table),
        log2_aveCPM = et$table$logCPM,
        log2FC = et$table$logFC,
        Pvalues = et$table$PValue,
        deGenes = de[, 1]
    )
    de_genes <- left_join(de_genes, entap, by = c("names" = "Query Sequence"))

    # output sequence lists
    cat("Outputing Differential Lists\n")
    splits <- str_split(args$f, ",")[[1]]
    neg_one <- splits[1]
    one <- splits[2]
    # up DE sequences
    de_genes |>
        filter(deGenes == 1) |>
        dplyr::select(names) |>
        filter(!is.na(names)) |>
        write.table(
            file = paste0(
                output_dir_list,
                file_start,
                one,
                "_de_sequences_list.txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
    # down DE sequences
    de_genes |>
        filter(deGenes == -1) |>
        dplyr::select(names) |>
        write.table(
            file = paste0(
                output_dir_list,
                file_start,
                neg_one,
                "_de_sequences_list.txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
    # up DE genes
    de_genes |>
        filter(deGenes == 1) |>
        filter(!is.na(`EggNOG Predicted Gene`)) |>
        dplyr::select(names) |>
        write.table(
            file = paste0(
                output_dir_list,
                file_start,
                one,
                "_de_gene_list.txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
    # down DE genes
    de_genes |>
        filter(deGenes == -1) |>
        filter(!is.na(`EggNOG Predicted Gene`)) |>
        dplyr::select(names) |>
        write.table(
            file = paste0(
                output_dir_list,
                file_start,
                neg_one,
                "_de_gene_list.txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
    # up DE seqs/genes
    de_genes |>
        filter(deGenes == 1) |>
        dplyr::select(names, `EggNOG Predicted Gene`) |>
        write.csv(
            file = paste0(output_dir_list, file_start, one, "_de_list.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            na = ""
        )
    # down DE seqs/genes
    de_genes |>
        filter(deGenes == -1) |>
        dplyr::select(names, `EggNOG Predicted Gene`) |>
        write.csv(
            file = paste0(output_dir_list, file_start, neg_one, "_de_list.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            na = ""
        )
    # all genes
    entap |>
        dplyr::select(`Query Sequence`, `EggNOG Predicted Gene`) |>
        filter(!is.na(`EggNOG Predicted Gene`)) |>
        write.table(
            file = paste0(output_dir_list, file_start, "all_genes_list.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE
        )
    rm(splits, de_genes)
}
### PLOTTING ###
# logCPM distribution
cat("Plotting CPM distribution\n")
logCPM_plot <- data.frame(read_counts = log10(DGE$counts + 1))
colnames(logCPM_plot) <- DGE$samples$samples
logCPM_plot <- melt(logCPM_plot)
colnames(logCPM_plot) <- c("sample", "log10cpm")
logCPM_plot <- merge(logCPM_plot, lib[, 1:3], by = "sample")
logCPM_plot |>
    ggplot(aes(x = ind, y = log10cpm, color = condition, fill = condition)) +
    geom_violin(scale = "count") +
    labs(x = "Sample", y = "log10(CPM)") +
    theme_cowplot() +
    custom_theme()
ggsave(
    filename = paste0(file_start, "CPM_violin.png"),
    path = output_dir_fig,
    device = "png",
    dpi = "retina"
)
cat("Saved at: ", paste0(output_dir_fig, file_start, "CPM_violin.png"), "\n")
write.csv(logCPM_plot, paste0(output_dir_data, file_start, "CPM_violin.csv"))
rm(logCPM_plot)

# MDS plot
cat("Plotting MDS\n")
mds <- plotMDS.DGEList(DGE, plot = FALSE)
mds_plot <- data.frame(
    x_val = mds$x,
    y_val = mds$y,
    sample = row.names(DGE$samples),
    Condition = lib$condition,
    Individual = lib$ind
)
varience1 <- mds$var.explained[1] * 100
varience2 <- mds$var.explained[2] * 100
labels <- str_split(args$f, ",")[[1]]
mds_plot |>
    ggplot(aes(x = x_val, y = y_val, color = Condition, shape = Individual)) +
    geom_point() +
    geom_text_repel(size = 4, aes(label = sample)) +
    labs(
        x = paste0("Leading LogFC dim1 (", round(varience1, 2), "%)"),
        y = paste0("Leading LogFC dim2 (", round(varience2, 2), "%)"),
        color = "Condition"
    ) +
    scale_color_manual(
        name = "Sample Type", #testing
        labels = c(labels[1], labels[2]),
        values = c("red", "blue")
    ) +
    #scale_color_manual(values = c("skin" = "red", "slime" = "blue")) +

    theme_cowplot() +
    custom_theme()
ggsave(
    filename = paste0(file_start, "LogFC_MDS.png"),
    path = output_dir_fig,
    device = "png",
    dpi = "retina",
    width = 10,
    height = 10
)
cat("Saved at: ", paste0(output_dir_fig, file_start, "LogFC_MDS.png"), "\n")
write.csv(mds_plot, paste0(output_dir_data, file_start, "LogFC_MDS.csv"))
rm(mds, mds_plot, varience1, varience2)

# BCV plot
cat("Plotting BCV\n")
bcv <- data.frame(
    log2_aveCPM = DGE$AveLogCPM,
    tagWdis_bcv = sqrt(DGE$tagwise.dispersion)
)
row.names(bcv) <- row.names(et$table)
common_bcv <- sqrt(DGE$common.dispersion)
trend_bcv <- sqrt(DGE$trended.dispersion)
bcv |>
    ggplot(aes(log2_aveCPM, tagWdis_bcv)) +
    geom_point(aes(color = "Tagwise"), size = 0.50, alpha = 0.5) +
    geom_line(aes(y = common_bcv, color = "Common")) +
    geom_line(aes(y = trend_bcv, color = "Trend")) +
    labs(x = "log2(mean CPM)", y = "Biological Coefficient of Variation") +
    scale_color_manual(
        name = "BCV Type",
        values = c("Tagwise" = "black", "Common" = "red", "Trend" = "blue")
    ) +
    theme_cowplot() +
    custom_theme()
ggsave(
    filename = paste0(file_start, "BCV.png"),
    path = output_dir_fig,
    device = "png",
    dpi = "retina"
)
cat("Saved at: ", paste0(output_dir_fig, file_start, "BCV.png"), "\n")
write.csv(bcv, paste0(output_dir_data, file_start, "BCV.csv"))
rm(bcv, common_bcv, trend_bcv)

# MA plot
cat("Plotting MA\n")
labels <- str_split(args$f, ",")[[1]]
ma <- data.frame(
    log2_aveCPM = et$table$logCPM,
    log2FC = et$table$logFC,
    Pvalues = et$table$PValue,
    deGenes = as.factor(de[, 1])
)
row.names(ma) <- row.names(et$table)
ma |>
    ggplot(aes(log2_aveCPM, log2FC, color = deGenes)) +
    labs(x = "log2(MeanCPM)", y = "Log2(FC)") +
    geom_point(size = 0.50, alpha = 0.5) +
    geom_hline(yintercept = c(-1, 1), col = "black") +

    scale_color_manual(
        name = "DE Genes",
        labels = c(labels[1], "none", labels[2]),
        values = c("red", "black", "blue")
    ) +
    theme_cowplot() +
    custom_theme()
ggsave(
    filename = paste0(file_start, "MA.png"),
    path = output_dir_fig,
    device = "png",
    dpi = "retina"
)
cat("Saved at: ", paste0(output_dir_fig, file_start, "MA.png"), "\n")
write.csv(ma, paste0(output_dir_data, file_start, "MA.csv"))
rm(ma)

# Volcano plot
cat("Plotting Volcano\n")
adjusted_P <- -log10(et$table$PValue)
volc <- data.frame(
    log2FC = et$table$logFC,
    log10_pval = adjusted_P,
    deGenes = as.factor(de[, 1])
)
row.names(volc) <- row.names(et$table)
volc |>
    ggplot(aes(log2FC, log10_pval, color = deGenes)) +
    labs(x = "Log2(FC)", y = "-Log10(Pvalue)") +
    geom_point(size = 0.50, alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), col = "black") +
    scale_color_manual(
        name = "DE Genes",
        labels = c(labels[1], "none", labels[2]),
        values = c("red", "black", "blue")
    ) +
    theme_cowplot() +
    custom_theme()
ggsave(
    filename = paste0(file_start, "Volcano.png"),
    path = output_dir_fig,
    device = "png",
    dpi = "retina"
)
cat("Saved at: ", paste0(output_dir_fig, file_start, "Volcano.png"), "\n")
write.csv(volc, paste0(output_dir_data, file_start, "Volcano.csv"))
rm(volc, adjusted_P)

if (!is.na(args$e)) {
    # Labeled Volcano plot
    cat("Plotting Labeled Volcano\n")
    adjusted_P <- -log10(et$table$PValue)
    volc <- data.frame(
        log2FC = et$table$logFC,
        log10_pval = adjusted_P,
        deGenes = de[, 1],
        labels = rownames(et$table)
    )
    volc <- left_join(volc, entap_filtered, by = c("labels" = "Query Sequence"))
    row.names(volc) <- row.names(et$table)

    VALUES = c(
        "-1" = "red",
        "0" = "grey",
        "1" = "blue",
        "3" = "black",
        "4" = "purple",
        "5" = "green"
    )

    for (i in 1:nrow(volc)) {
        if (!is.na(volc$`EggNOG Predicted Gene`[i]) && volc$deGenes[i] == 1) {
            volc$deGenes[i] <- 4
        }
    }
    rm(i)

    cat("Labeled DE genes above", args$x, "log2FC:\n")
    volc |>
        mutate(deGenes = ifelse((log2FC > args$x & deGenes == 4), 4, 0)) |>
        group_by(deGenes) |>
        summarise(n = n()) |>
        kable()
    cat("\n")

    volc |>
        #mutate(`EggNOG Predicted Gene` = ifelse(is.na(`EggNOG Predicted Gene`), "NA", `EggNOG Predicted Gene`)) |>
        mutate(
            `EggNOG Predicted Gene` = ifelse(
                (deGenes != 4),
                NA,
                `EggNOG Predicted Gene`
            )
        ) |>
        mutate(
            `EggNOG Predicted Gene` = ifelse(
                log2FC >= args$x,
                `EggNOG Predicted Gene`,
                NA
            )
        ) |>
        mutate(deGenes = as.factor(deGenes)) |>
        ggplot(aes(x = log2FC, y = log10_pval, color = deGenes)) +
        labs(x = "Log2(FC)", y = "-Log10(Pvalue)") +
        geom_point(size = 0.50, alpha = 0.5) +
        geom_vline(xintercept = c(-1, 1), col = "black") +
        geom_text_repel(
            aes(label = `EggNOG Predicted Gene`),
            size = 1,
            max.overlaps = Inf,
            segment.color = NA,
            na.rm = TRUE
        ) +
        scale_color_manual(
            name = "DE Genes",
            labels = c("Skin", "not DE", "Slime", "Top Labeled "),
            values = VALUES
        ) +
        theme_cowplot() +
        custom_theme()

    ggsave(
        filename = paste0(file_start, "Volcano_labeled.png"),
        path = output_dir_fig,
        device = "png",
        dpi = "retina"
    )
    cat(
        "Saved at: ",
        paste0(output_dir_fig, file_start, "Volcano_labeled.png"),
        "\n"
    )
    write.csv(volc, paste0(output_dir_data, file_start, "Volcano_labeled.csv"))
    rm(volc, adjusted_P, labels)

    ### von Willebrand MA ###
    cat("von Willebrand MA\n")
    search_strings <- c("VWA", "VWD", "VWC")
    entap_VWF <- entap |>
        dplyr::select(
            `Query Sequence`,
            `EggNOG Predicted Gene`,
            `EggNOG Description`,
            `EggNOG Protein Domains`
        ) |>
        #filter(!is.na(`EggNOG Predicted Gene`)) |>
        mutate(
            PFAM = ifelse(
                grepl("PFAM", `EggNOG Protein Domains`),
                str_split(
                    gsub(
                        ".*PFAM \\((.*?)\\).*",
                        "\\1",
                        `EggNOG Protein Domains`
                    ),
                    ", "
                ),
                NA
            ),
            SMART = ifelse(
                grepl("SMART ", `EggNOG Protein Domains`),
                str_split(
                    gsub(
                        ".*SMART \\((.*?)\\).*",
                        "\\1",
                        `EggNOG Protein Domains`
                    ),
                    ", "
                ),
                NA
            )
        ) |>
        filter(
            map_lgl(PFAM, ~ any(.x %in% search_strings)) |
                map_lgl(SMART, ~ any(.x %in% search_strings))
        )

    labels <- str_split(args$f, ",")[[1]]
    VALUES = c(
        "-1" = "red",
        "0" = "grey",
        "1" = "blue",
        "3" = "black",
        "4" = "purple",
        "5" = "green"
    )
    #VALUES = c("-1"="#ef476f", "0"="#373F51", "1"="#4cc9f0", "3"="grey", "4"="#7209b7", "5"="#06d6a0")
    ma <- data.frame(
        names = rownames(et$table),
        log2_aveCPM = et$table$logCPM,
        log2FC = et$table$logFC,
        Pvalues = et$table$PValue,
        deGenes = de[, 1],
        adjusted_P = -log10(et$table$PValue)
    )
    ma <- left_join(ma, entap_VWF, by = c("names" = "Query Sequence"))

    write.csv(
        ma |>
            filter(deGenes >= 1 & !is.na(`EggNOG Predicted Gene`)) |>
            select(`names`, `EggNOG Predicted Gene`),
        paste0(output_dir_list, file_start, "vonWillebrand_slime.txt"),
        na = "",
        quote = FALSE,
        row.names = FALSE
    )

    ylim_slime <- c(1, NA)
    ylim_skin <- c(NA, -1)
    ma_plot <- ma |>
        mutate(deGenes = ifelse(log2FC >= 1 | log2FC <= -1, deGenes, 0)) |>
        mutate(
            `EggNOG Predicted Gene` = ifelse(
                deGenes == 0,
                NA,
                `EggNOG Predicted Gene`
            )
        ) |>
        #mutate(deGenes = ifelse(!is.na(`EggNOG Predicted Gene`), 5, deGenes)) |>
        mutate(
            pred_gene_slime = ifelse(deGenes == 1, `EggNOG Predicted Gene`, NA)
        ) |>
        mutate(
            pred_gene_skin = ifelse(deGenes == -1, `EggNOG Predicted Gene`, NA)
        ) |>
        #mutate(`pred_gene_slime` = ifelse(grepl("^MUC", `pred_gene_slime`), `pred_gene_slime`, NA)) |>
        #mutate(`pred_gene_skin` = ifelse(grepl("^MUC", `pred_gene_skin`), `pred_gene_skin`, NA)) |>
        mutate(
            deGenes = ifelse(
                !is.na(`EggNOG Protein Domains`) & deGenes != 0,
                5,
                deGenes
            )
        ) |>
        mutate(
            deGenes = ifelse(
                is.na(`EggNOG Predicted Gene`) & deGenes == 5,
                4,
                deGenes
            )
        ) |>

        mutate(
            deGenes = ifelse(
                (!is.na(pred_gene_skin) | !is.na(pred_gene_slime)),
                5,
                deGenes
            )
        ) |>
        mutate(
            deGenes = ifelse(
                is.na(`EggNOG Predicted Gene`) & deGenes == 5,
                4,
                deGenes
            )
        ) |>
        # mutate(`EggNOG Predicted Gene` = ifelse(is.na(`EggNOG Predicted Gene`) & deGenes == 5, "NA", `EggNOG Predicted Gene`)) |>
        mutate(deGenes = factor(deGenes)) |>
        mutate(alpha = ifelse(deGenes == 5 | deGenes == 4, 1, 0.5))

    ma_plot |>
        ggplot(aes(log2_aveCPM, log2FC, color = deGenes)) +
        labs(x = "log2(MeanCPM)", y = "Log2(FC)") +
        #geom_point(data = subset(ma, deGenes == 4 | deGenes == 5), size=0.50, aes(alpha=alpha)) +
        geom_point(
            data = function(x)
                subset(x, deGenes == -1 | deGenes == 0 | deGenes == 1),
            size = 0.50,
            aes(log2_aveCPM, log2FC, color = deGenes, alpha = alpha)
        ) +
        geom_point(
            data = function(x) subset(x, deGenes == 4 | deGenes == 5),
            size = 0.50,
            aes(log2_aveCPM, log2FC, color = deGenes, alpha = alpha)
        ) +
        geom_hline(yintercept = c(-1, 1), col = "black") +
        geom_text_repel(
            aes(label = pred_gene_slime),
            ylim = ylim_slime,
            color = "black",
            segment.color = "gray",
            na.rm = TRUE,
            size = 2,
            force = 20,
            box.padding = 0.2,
            segment.size = 0.2,
            min.segment.length = 0.1,
            max.overlaps = Inf,
            point.padding = 0.01,
            max.time = 1000
        ) +
        geom_text_repel(
            aes(label = pred_gene_skin),
            ylim = ylim_skin,
            color = "black",
            segment.color = "grey",
            na.rm = TRUE,
            size = 2,
            force = 20,
            box.padding = 0.2,
            segment.size = 0.2,
            min.segment.length = 0.1,
            max.overlaps = Inf,
            point.padding = 0.01,
            max.time = 1000
        ) +
        theme_cowplot() +
        scale_alpha_continuous(guide = "none") +
        scale_color_manual(
            name = "DE Genes",
            values = VALUES,
            labels = c(
                neg_one,
                "not DE",
                one,
                "VW no Prediction",
                "Predicted Protein VW domains"
            )
        ) +
        #custom_theme(base_family="fira sans",line_size=0.4) +
        custom_theme()

    ggsave(
        filename = paste0(file_start, "vonWillebrand_MA.png"),
        path = output_dir_fig,
        device = "png",
        dpi = "retina"
    )
    cat(
        "Saved at: ",
        paste0(output_dir_fig, file_start, "vonWillebrand_MA.png"),
        "\n"
    )
    #write.csv(ma_plot, file=paste0(output_dir_data, file_start, "vonWillebrand_MA.csv"))

    xlim_slime <- c(1, NA)
    xlim_skin <- c(NA, -1)
    ma_plot |>
        ggplot(aes(x = log2FC, y = adjusted_P, color = deGenes)) +
        labs(x = "log2FC", y = "-log10(Pvalue)") +
        geom_point(
            data = function(x)
                subset(x, deGenes == -1 | deGenes == 0 | deGenes == 1),
            size = 0.50,
            aes(log2FC, adjusted_P, color = deGenes, alpha = alpha)
        ) +
        geom_point(
            data = function(x) subset(x, deGenes == 4 | deGenes == 5),
            size = 0.50,
            aes(log2FC, adjusted_P, color = deGenes, alpha = alpha)
        ) +
        geom_vline(xintercept = c(-1, 1), col = "black") +
        geom_text_repel(
            aes(label = pred_gene_slime),
            xlim = xlim_slime,
            color = "black",
            segment.color = "gray",
            na.rm = TRUE,
            size = 2,
            force = 20,
            box.padding = 0.2,
            segment.size = 0.2,
            min.segment.length = 0.1,
            max.overlaps = Inf,
            point.padding = 0.01,
            max.time = 1000
        ) +
        geom_text_repel(
            aes(label = pred_gene_skin),
            xlim = xlim_skin,
            color = "black",
            segment.color = "gray",
            na.rm = TRUE,
            size = 2,
            force = 20,
            box.padding = 0.2,
            segment.size = 0.2,
            min.segment.length = 0.1,
            max.overlaps = Inf,
            point.padding = 0.01,
            max.time = 1000
        ) +
        theme_cowplot() +
        scale_alpha_continuous(guide = "none") +
        scale_color_manual(
            name = "DE Genes",
            values = VALUES,
            labels = c(
                neg_one,
                "not DE",
                one,
                "VW no Prediction",
                "Predicted Protein VW domains"
            )
        ) +
        #custom_theme(base_family="fira sans",line_size=0.4) +
        custom_theme()

    ggsave(
        filename = paste0(file_start, "vonWillebrand_volc.png"),
        path = output_dir_fig,
        device = "png",
        dpi = "retina"
    )
    cat(
        "Saved at: ",
        paste0(output_dir_fig, file_start, "vonWillebrand_volc.png"),
        "\n"
    )
    rm(ma, ma_plot, entap_VWF, search_strings, ylim_slime, ylim_skin)
}
#
# entap |>
#     select(`Query Sequence`, `EggNOG Predicted Gene`) |>
#
# # GO analysis
# library(clusterProfiler)
# library(org.Hs.eg.db)
# ego2 <- enrichGO(gene         = ,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENSEMBL',
#                  ont           = "CC",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)
#
# de$sample <- rownames(de)
# thing2 <- merge(entap, de, by.x="Query Sequence", by.y="sample")
# things <- thing2 |>
#     filter(slime.skin == 1) |>
#     filter(`EggNOG Predicted Gene` != "NA") |>
#     dplyr::select(`EggNOG Predicted Gene`)
# write.csv(things, file="~/Desktop/M_phantasma_slime.csv", row.names = FALSE, quote = FALSE)
# merge(entap)
# thing <- entap |>
#     dplyr::select(`EggNOG Protein Domains`, `EggNOG Predicted Gene`, ) |>
#     filter(!is.na(`EggNOG Protein Domains`))
