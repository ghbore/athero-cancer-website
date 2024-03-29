#! /usr/bin/env Rscript

suppressPackageStartupMessages({
    library(shiny)
    library(shinydashboard)
    library(shinydashboardPlus)
    library(tidyverse)
    library(SummarizedExperiment)
    library(plotly)
    library(DT)
})

source("helper.r")

# load data
GENESET <- yaml::read_yaml("data/expression/geneset.yml")
GENESET_DEFAULT <- "hallmark"

DATASET <- yaml::read_yaml("data/expression/dataset.yml")

## load athero dataset and pathway
ATHERO_DATASET <- DATASET$athero %>%
    unlist() %>%
    setNames(., ifelse("" == names(.), ., names(.)))
ATHERO_DATASET_DEFAULT <- "STARNET_case_control"
ATHERO <- load_or_save("data/athero.rda", {
    lapply(ATHERO_DATASET, function (s)
        str_glue("data/expression/{s}.rda") %>%
            readRDS()
    ) %>%
    setNames(ATHERO_DATASET)
})
ATHERO_PATHWAY <- load_or_save("data/athero-pathway.rda", {
    lapply(GENESET, function (db){
        lapply(ATHERO_DATASET, function (ds){
            str_glue("data/expression/{ds}-{db}.rda") %>%
                readRDS()
        }) %>%
            setNames(ATHERO_DATASET)
    }) %>%
        setNames(GENESET)
})

## load cancer dataset and pathway
CANCER_DATASET <- DATASET$cancer %>%
    unlist() %>%
    setNames(., ifelse("" == names(.), ., names(.)))
CANCER_DATASET_DEFAULT <- "STAD"
CANCER <- load_or_save("data/cancer.rda", {
    lapply(CANCER_DATASET, function (s)
        str_glue("data/expression/{s}.rda") %>%
            readRDS()
    ) %>%
    setNames(CANCER_DATASET)
})
CANCER_PATHWAY <- load_or_save("data/cancer-pathway.rda", {
    lapply(GENESET, function (db){
        lapply(CANCER_DATASET, function (ds){
            str_glue("data/expression/{ds}-{db}.rda") %>%
                readRDS()
        }) %>%
            setNames(CANCER_DATASET)
    }) %>%
        setNames(GENESET)
})

logging("finish loading data")

# inner name reverse to display name
DATASET_NAME_MAP <- list(
    inner = c(ATHERO_DATASET, CANCER_DATASET),
    display = setNames(
        c(names(ATHERO_DATASET), names(CANCER_DATASET)),
        c(ATHERO_DATASET, CANCER_DATASET)
    )
)

# similarity matrix and dimension reduction
GENE_MATRIX <- load_or_save("data/gene-matrix.rda", {
    combine_gene_all(c(ATHERO, CANCER))
})
GENE_DIM <- load_or_save("data/gene-dim.rda", {
    scale_cluster_redim(GENE_MATRIX)
})

PATHWAY_MATRIX <- load_or_save("data/pathway-matrix.rda", {
    lapply(setNames(GENESET, GENESET), function (db){
        combine_pathway_all(
            c(ATHERO_PATHWAY[[db]], CANCER_PATHWAY[[db]])
        )
    })
})
PATHWAY_DIM <- load_or_save("data/pathway-dim.rda", {
    lapply(setNames(GENESET, GENESET), function (db){
        scale_cluster_redim(PATHWAY_MATRIX[[db]])
    })
})

logging("finish calculating similarity matrix")