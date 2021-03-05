#! /usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) != 1){
    stop("Usage: format.r yaml_file")
}
if (! file.exists(args[1])){
    stop("can NOT open yaml_file")
}

library(tidyverse)
library(SummarizedExperiment)

# load data
gmap <- list(
    human = "~/../Desktop/reference/Homo_sapiens/gencode/gene_map.rda",
    mouse = "~/../Desktop/reference/Mus_musculus/gencode/gene_map2human.rda"
) %>%
    lapply(function (f){
        readRDS(f) %>%
            rowData() %>%
            as_tibble()
    })

# function

# main
config <- yaml::read_yaml(args[1])

names(config) %>%
    str_subset("^-", negate = TRUE) %>%
    lapply(function (prefix){
        cfg <- config[[prefix]]
        df <- data.table::fread(cfg$file) %>%
            as_tibble() %>%
            dplyr::select(all_of(cfg$column %>% unlist()))
        kv <- intersect(gmap[[cfg$organism]] %>% colnames(), colnames(df))
        df <- distinct(df, across(all_of(kv)), .keep_all = TRUE)
        rowdata <- dplyr::select(df, all_of(kv)) %>%
            left_join(gmap[[cfg$organism]]) %>%
            distinct(across(all_of(kv)), .keep_all = TRUE)
        if (nrow(rowdata) != nrow(df)){
            str_glue("rowData is not compatible to assay in {prefix}") %>%
                stop()
        }
        SummarizedExperiment(
            as(select(df, !all_of(kv)), "DataFrame"),
            rowData = as(rowdata, "DataFrame"),
            metadata = cfg$metadata
        ) %>%
            saveRDS(str_glue("{prefix}.rda"))
    })