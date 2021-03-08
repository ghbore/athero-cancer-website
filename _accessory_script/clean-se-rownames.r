#! /usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0){
    stop("Usage: clean-se-rownames.r se.rda ...")
}
if (! file.exists(args[1])){
    stop("can NOT open yaml_file")
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(SummarizedExperiment)
})

for (i in seq(args)){
    se <- readRDS(args[i])

    rownames(se) <- NULL

    saveRDS(se, args[i])
    cat(file = stderr(), paste0("finish cleaning ", args[i], "\n"))
}