#! /usr/bin/env Rscript

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(tidyverse)
library(SummarizedExperiment)
library(plotly)
library(DT)

# function
pathway_comparison <- function (athero, cancer){
    inner_join(
        cbind(rowData(athero), assay(athero)[, c("NES", "pvalue")]) %>%
            as_tibble() %>%
            dplyr::rename("NES:athero" = NES, "p:athero" = pvalue),
        cbind(rowData(cancer), assay(cancer)[, c("NES", "pvalue")]) %>%
            as_tibble() %>%
            dplyr::rename("NES:cancer" = NES, "p:cancer" = pvalue)
    )
}

# load data
GENESET <- yaml::read_yaml("data/expression/geneset.yml")
GENESET_DEFAULT <- "hallmark"

DATASET <- yaml::read_yaml("data/expression/dataset.yml")

ATHERO_DATASET <- DATASET$athero %>% 
    unlist() %>% 
    setNames(., ifelse(""==names(.), ., names(.)))
ATHERO_DATASET_DEFAULT <- "BiKE"
ATHERO <- lapply(ATHERO_DATASET, function (s)
        str_glue("data/expression/{s}.rda") %>%
            readRDS()
    ) %>%
    setNames(ATHERO_DATASET)
ATHERO_PATHWAY <- lapply(GENESET, function (db){
    lapply(ATHERO_DATASET, function (ds){
        str_glue("data/expression/{ds}-{db}.rda") %>%
            readRDS()
    }) %>%
        setNames(ATHERO_DATASET)
}) %>%
    setNames(GENESET)

CANCER_DATASET <- DATASET$cancer %>% 
    unlist() %>% 
    setNames(., ifelse("" == names(.), ., names(.)))
CANCER_DATASET_DEFAULT <- "STAD"
CANCER <- lapply(CANCER_DATASET, function (s)
        str_glue("data/expression/{s}.rda") %>%
            readRDS()
    ) %>%
    setNames(CANCER_DATASET)
CANCER_PATHWAY <- lapply(GENESET, function (db){
    lapply(CANCER_DATASET, function (ds){
        str_glue("data/expression/{ds}-{db}.rda") %>%
            readRDS()
    }) %>%
        setNames(CANCER_DATASET)
}) %>%
    setNames(GENESET)

cat(file = stderr(), "finish loading data\n")

### TODO: gene level similarity matrix
DATASET_SIMILARITY <- list(
    pathway = lapply(GENESET, function (db, nms){
        ds <- c(ATHERO_PATHWAY[[db]], CANCER_PATHWAY[[db]])
        name <- lapply(ds, function (d) rowData(d)$Description) %>%
            Reduce(union, .)
        mtx <- lapply(ds, function (d)
                rep(0, length(name)) %>%
                    `[<-`(match(rowData(d)$Description, name),
                        value = assay(d)$NES
                    )
            ) %>%
            Reduce(cbind, .) %>%
            `colnames<-`(nms) %>%
            `rownames<-`(str_remove(name, "^GO_|^HALLMARK_") %>%
                str_replace_all("_", " ") %>%
                str_trunc(60)
            ) %>%
            as.matrix()
        mtx[is.na(mtx)] <- 0
        rh <- t(mtx) %>% scale() %>%
            `[<-`(is.na(.), value = 0) %>%
            t() %>% dist() %>% hclust()
        cd <- scale(mtx) %>%
            `[<-`(is.na(.), value = 0) %>%
            t() %>% dist()
        ch <- hclust(cd)
        list(
            dist = as.matrix(cd)[ch$order, ch$order],
            matrix = mtx[rh$order, ch$order] %>% as.matrix()
        )
    }, nms = c(names(ATHERO_DATASET), names(CANCER_DATASET))) %>%
        setNames(GENESET)
)

cat(file = stderr(), "finish calculating similarity matrix\n")