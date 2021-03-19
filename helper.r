#! /usr/bin/env Rscript

#' load RData or eval and save to RData
#'
#' @param file Filepath to RData
#' @param expr Expression to eval once RData does not exist.
#'        The return value of `expr` will be saved to a RData
#' @return Value from RData or evaluated `expr`
load_or_save <- function (file, expr){
    if (file.exists(file)){
        # cat(file = stderr(), "load from file\n")
        obj <- readRDS(file)
    }else {
        # cat(file = stderr(), "eval and save to file\n")
        env <- rlang::env(rlang::caller_env())
        obj <- rlang::eval_tidy(rlang::enexpr(expr))
        saveRDS(obj, file)
    }
    return(obj)
}

select_gene_statistics <- function (se, tag = "athero"){
    bind_cols(
        rowData(se) %>%
            as_tibble() %>%
            dplyr::select(ensembl, symbol),
        assay(se) %>%
            as_tibble() %>%
            dplyr::select(any_of(
                c("logHR", "logFC", "logOR", "cor", "pval")
            )) %>%
            dplyr::rename(p = pval) %>%
            dplyr::rename_with(function (nm){
                str_glue("{nm}:{tag}")
            })
    )
}
combine_gene <- function (athero, cancer){
    inner_join(
        select_gene_statistics(athero, "athero"),
        select_gene_statistics(cancer, "cancer")
    )
}
combine_gene_all <- function (lists, new_names = names(lists)){
    kv <- lapply(lists, function (o){
        rowData(o)[, c("ensembl", "symbol")] %>%
            as_tibble()
    }) %>%
        bind_rows() %>%
        distinct() %>%
        dplyr::mutate(id = paste(ensembl, symbol, sep="-"))
    mtx <- bind_cols(
        dplyr::select(kv, !id),
        matrix(NA_real_, ncol = length(lists), nrow = nrow(kv)) %>%
            `colnames<-`(new_names) %>%
            as_tibble()
    )
    for (i in seq_along(lists)){
        k <- rowData(lists[[i]])[, c("ensembl", "symbol")] %>%
            apply(1, paste, collapse = "-") %>%
            match(kv$id)
        stopifnot(!any(is.na(k)))
        v <- assay(lists[[i]])[[
            intersect(
                c("logHR", "logFC", "logOR", "cor"),
                colnames(lists[[i]])
            )[1]
        ]]
        si <- !duplicated(k)
        mtx[k[si], 2+i] <- v[si]
    }
    return(mtx)
}


select_pathway_statistics <- function (se, tag = "athero"){
    bind_cols(
        rowData(se) %>% as_tibble(),
        assay(se)[, c("NES", "pvalue")] %>%
            as_tibble() %>%
            dplyr::rename(p = pvalue) %>%
            dplyr::rename_with(function (nm){
                str_glue("{nm}:{tag}")
            })
    )
}
combine_pathway <- function (athero, cancer){
    inner_join(
        select_pathway_statistics(athero, "athero"),
        select_pathway_statistics(cancer, "cancer")
    )
}
combine_pathway_all <- function (lists, new_names = names(lists)){
    kv <- lapply(lists, function (o){
        rowData(o) %>% as_tibble()
    }) %>%
        bind_rows() %>%
        distinct()
    kv$id <- apply(kv, 1, paste, collapse = "-")
    mtx <- bind_cols(
        dplyr::select(kv, !id),
        matrix(NA_real_, ncol = length(lists), nrow = nrow(kv)) %>%
            `colnames<-`(new_names) %>%
            as_tibble()
    )
    for (i in seq_along(lists)){
        k <- rowData(lists[[i]]) %>%
            apply(1, paste, collapse = "-") %>%
            match(kv$id)
        stopifnot(!any(is.na(k)))
        v <- assay(lists[[i]])[["NES"]]
        si <- !duplicated(k)
        mtx[k[si], 2+i] <- v[si]
    }
    return(mtx)
}

combine_se <- function (lists,
    row_keys = c("ensembl", "symbol", "entrez", "ID", "Description"),
    col_keys = list(
        stat = c("logFC", "logHR", "logOR", "cor", "NES"),
        pval = c("pval", "pvalue")
    ),
    new_names = names(lists)
){
    if (is.null(new_names)){
        new_names <- as.character(seq(lists))
    }
    kv <- lapply(lists, function (o) rowData(o) %>%
        as_tibble() %>%
        dplyr::select(any_of(row_keys))
    ) %>%
        bind_rows() %>%
        distinct() %>%
        as.data.frame() %>%
        `rownames<-`(apply(., 1, paste, collapse = "-"))
    mtx <- lapply(col_keys, function (k){
        lapply(lists, function (o){
            v <- rep(NA, nrow(kv))
            k0 <- rowData(o) %>%
                as_tibble() %>%
                dplyr::select(any_of(row_keys)) %>%
                apply(1, paste, collapse = "-")
            v0 <- assay(o) %>%
                as_tibble() %>%
                dplyr::select(any_of(k)) %>%
                `[`(i =, j = 1, drop = TRUE)
            v[match(k0, rownames(kv))] <- v0
            return(v)
        }) %>%
            as.data.frame() %>%
            `colnames<-`(new_names)
    })

    rownames(kv) <- NULL
    SummarizedExperiment(
        lapply(mtx, as, Class = "DataFrame"),
        rowData = as(kv, "DataFrame"),
        colData = tibble(name = new_names) %>%
            as("DataFrame")
    )
}

scale_cluster_redim <- function (df, .scale = TRUE){
    x <- dplyr::select_if(df, is.numeric) %>%
        `[<-`(is.na(.), value = 0) %>%
        `[`(i =, j = apply(., 2, function (x) any(x != x[1])), drop = FALSE) %>%
        scale(center = FALSE, scale = .scale)
    d <- dist(t(x))
    cluster <- hclust(d)
    pca <- prcomp(x, retx = FALSE, center = FALSE) %>%
        `$`("rotation") %>%
        `[`(i =, j = 1:2) %>%
        `colnames<-`(c("PC_1", "PC_2")) %>%
        as.data.frame() %>%
        rownames_to_column("id") %>%
        as_tibble()
    set.seed(1)
    umap <- uwot::umap(t(x), n_neighbors = max(2, floor(ncol(x)/2))) %>%
        `colnames<-`(c("UMAP_1", "UMAP_2")) %>%
        as_tibble() %>%
        mutate(id = colnames(x), .before = 1)
    list(scaled = x, dist = d, cluster = cluster, pca = pca, umap = umap)
}

logging <- function (...){
    cat(file = stderr(),
        "[",
        format(Sys.time()),
        "] ",
        ...,
        sep = "",
        fill = TRUE
    )
}
