#! /usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) != 1){
    stop("Usage: pathway.r yaml_file")
}
if (! file.exists(args[1])){
    stop("can NOT open yaml_file")
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(SummarizedExperiment)
    library(clusterProfiler)
})

# load data
msigdb <- "~/../Desktop/reference/msigdb"
hallmark <- read.gmt(file.path(msigdb, "h.all.v7.2.entrez.gmt"))
gobp <- read.gmt(file.path(msigdb, "c5.go.bp.v7.2.entrez.gmt"))
gomf <- read.gmt(file.path(msigdb, "c5.go.mf.v7.2.entrez.gmt"))

# function
get_readable_result <- function (res){
    if (!is.null(res)){
        DOSE::setReadable(res, "org.Hs.eg.db", "ENTREZID")@result
    }
}

gsea2se <- function (df, metadata = NULL){
    SummarizedExperiment(
        assays = SimpleList(
            GSEA = dplyr::select(df, !c(ID, Description)) %>%
                as_tibble() %>%
                as("DataFrame")
        ),
        rowData = dplyr::select(df, ID, Description) %>%
            as_tibble() %>%
            as("DataFrame"),
        metadata = metadata
    )
}

gsea_dump <- function (se, prefix = "gsea", beta = "logHR", cutoff = 1){
    gls <- tibble(
            entrez = rowData(se)$entrez,
            beta = assay(se)[[beta]]
        ) %>%
        dplyr::filter(!is.na(entrez)) %>%
        arrange(desc(beta)) %>%
        distinct(entrez, .keep_all = TRUE) %>%
        with(setNames(beta, entrez))
    gseKEGG(gls, pvalueCutoff = cutoff, keyType = "ncbi-geneid", seed = 1) %>%
        get_readable_result() %>%
        gsea2se(metadata = SimpleList(db = "KEGG pathway")) %>%
        saveRDS(paste0(prefix, "-kegg.rda"))
    GSEA(gls, pvalueCutoff = cutoff, TERM2GENE = hallmark, seed = 1) %>%
        get_readable_result() %>%
        gsea2se(metadata = SimpleList(db = "Hallmark gene set")) %>%
        saveRDS(paste0(prefix, "-hallmark.rda"))
    GSEA(gls, pvalueCutoff = cutoff, TERM2GENE = gobp, seed = 1) %>%
        get_readable_result() %>%
        gsea2se(metadata = SimpleList(db = "GO biological process")) %>%
        saveRDS(paste0(prefix, "-gobp.rda"))
    GSEA(gls, pvalueCutoff = cutoff, TERM2GENE = gomf, seed = 1) %>%
        get_readable_result() %>%
        gsea2se(metadata = SimpleList(db = "GO molecular function")) %>%
        saveRDS(paste0(prefix, "-gomf.rda"))
}

# main
config <- yaml::read_yaml(args[1])

names(config) %>%
    str_subset("^-", negate = TRUE) %>%
    lapply(function (prefix){
        se <- readRDS(config[[prefix]]$file)
        modifyList(config[[prefix]], list(se = se, prefix = prefix)) %>%
            `[`(intersect(names(.), formalArgs(gsea_dump))) %>%
            do.call(gsea_dump, .)
    })
