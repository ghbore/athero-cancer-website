#! /usr/bin/env Rscript

function (input, output, session){
    output$athero_dataset_volcano <- renderPlotly({
        se <- ATHERO[[input$athero_dataset_sel]]
        beta_name <- intersect(c("logHR", "logFC", "logOR", "cor"), colnames(se))
        tibble(
            ensembl = rowData(se)$ensembl,
            symbol = rowData(se)$symbol,
            beta = assay(se)[[beta_name]],
            pval = assay(se)[["pval"]]
        ) %>%
            plot_ly(x = ~beta, y = ~-log10(pval),
                type = "scatter", mode = "markers",
                hoverinfo = "text",
                text = ~str_glue("{symbol} [{ensembl}]"),
                source = "athero_dataset_volcano"
            ) %>%
            layout(
                title = DATASET_NAME_MAP$display[input$athero_dataset_sel],
                xaxis = list(title = beta_name)
            )
    })

    output$athero_dataset <- renderDataTable({
        se <- ATHERO[[input$athero_dataset_sel]]
        cbind(
            rowData(se)[, c("ensembl", "symbol")],
            assay(se)
        ) %>%
            as_tibble()
    })

    output$athero_pathway_volcano <- renderPlotly({
        se <- ATHERO_PATHWAY[[input$pathway_db_sel]][[input$athero_dataset_sel]]
        tibble(
            Description = rowData(se)$Description %>%
                str_remove("^GO_|^HALLMARK_"),
            NES = assay(se)$NES,
            pval = assay(se)$pvalue,
            genes = assay(se)$core_enrichment
        ) %>%
            plot_ly(x = ~NES, y = ~-log10(pval),
                type = "scatter", mode = "markers",
                hoverinfo = "text",
                text = ~str_glue("{Description}<br />{genes}"),
                source = "athero_pathway_volcano"
            ) %>%
            layout(
                title = str_glue(
                    DATASET_NAME_MAP$display[input$athero_dataset_sel],
                    " at ",
                    names(GENESET)[GENESET == input$pathway_db_sel]
                )
            )
    })

    output$cancer_dataset_volcano <- renderPlotly({
        se <- CANCER[[input$cancer_dataset_sel]]
        beta_name <- intersect(c("logHR", "logFC", "logOR", "cor"), colnames(se))
        tibble(
            ensembl = rowData(se)$ensembl,
            symbol = rowData(se)$symbol,
            beta = assay(se)[[beta_name]],
            pval = assay(se)[["pval"]]
        ) %>%
            plot_ly(x = ~beta, y = ~-log10(pval),
                hoverinfo = "text", type = "scatter", mode = "markers",
                text = ~str_glue("{symbol} [{ensembl}]"),
                source = "cancer_dataset_volcano"
            ) %>%
            layout(
                title = DATASET_NAME_MAP$display[input$cancer_dataset_sel],
                xaxis = list(title = beta_name)
            )
    })

    output$cancer_dataset <- renderDataTable({
        se <- CANCER[[input$cancer_dataset_sel]]
        cbind(
            rowData(se)[, c("ensembl", "symbol")],
            assay(se)
        ) %>%
            as_tibble()
    })

    output$cancer_pathway_volcano <- renderPlotly({
        se <- CANCER_PATHWAY[[input$pathway_db_sel]][[input$cancer_dataset_sel]]
        tibble(
            Description = rowData(se)$Description %>%
                str_remove("^GO_|^HALLMARK_"),
            NES = assay(se)$NES,
            pval = assay(se)$pvalue,
            genes = assay(se)$core_enrichment
        ) %>%
            plot_ly(x = ~NES, y = ~-log10(pval),
                type = "scatter", mode = "markers",
                hoverinfo = "text",
                text = ~str_glue("{Description}<br />{genes}"),
                source = "cancer_pathway_volcano"
            ) %>%
            layout(
                title = str_glue(
                    DATASET_NAME_MAP$display[input$cancer_dataset_sel],
                    " at ",
                    names(GENESET)[GENESET == input$pathway_db_sel]
                )
            )
    })

    output$gene_matrix_heatmap <- renderPlotly({
        mtx <- GENE_DIM$scaled %>%
            `rownames<-`(with(GENE_MATRIX, paste(ensembl, symbol, sep="-"))) %>%
            `colnames<-`(DATASET_NAME_MAP$display[colnames(.)])
        plot_ly(
            x = colnames(mtx),
            y = rownames(mtx),
            z = mtx,
            colors = colorRamp(c("blue", "yellow", "red")),
            type = "heatmap"
        ) %>%
            layout(
                xaxis = list(title = ""),
                yaxis = list(title = "")
            )
    })

    output$gene_dist_heatmap <- renderPlotly({
        mtx <- GENE_DIM$dist %>%
            as.matrix() %>%
            `colnames<-`(DATASET_NAME_MAP$display[colnames(.)]) %>%
            `rownames<-`(DATASET_NAME_MAP$display[rownames(.)]) %>%
            `[`(GENE_DIM$cluster$order, GENE_DIM$cluster$order)
        odr <- as_tibble(mtx) %>%
            dplyr::mutate(across(everything(), rank)) %>%
            dplyr::mutate(x = colnames(.)) %>%
            pivot_longer(!x,
                names_to = "y",
                values_to = "rank"
            ) %>%
            dplyr::mutate(rank = rank - 1) %>%
            dplyr::filter(rank > 0, rank <= 3)
        plot_ly(
            x = colnames(mtx),
            y = rownames(mtx),
            z = mtx,
            colors = colorRamp(c("blue", "yellow", "red")),
            type = "heatmap"
        ) %>%
            add_trace(x = ~x, y = ~y, text = ~rank,
                data = odr,
                hoverinfo = "none",
                type = "scatter", mode = "text"
            ) %>%
            layout(
                xaxis = list(title = ""),
                yaxis = list(title = "")
            )
    })
    
    output$gene_pca <- renderPlotly({
        dplyr::mutate(GENE_DIM$pca,
            id = DATASET_NAME_MAP$display[id]
        ) %>%
            plot_ly(x = ~PC_1, y = ~PC_2, text = ~id,
                hoverinfo = "skip", type = "scatter", mode = "markers"
            ) %>%
            add_annotations(text = ~unname(id)) %>%
            config(editable = TRUE)
    })

    output$gene_umap <- renderPlotly({
        dplyr::mutate(GENE_DIM$umap,
            id = DATASET_NAME_MAP$display[id]
        ) %>%
            plot_ly(x = ~UMAP_1, y = ~UMAP_2, text = ~id,
                hoverinfo = "skip", type = "scatter", mode = "markers"
            ) %>%
            add_annotations(text = ~unname(id)) %>%
            config(editable = TRUE)
    })

    output$gene_quad <- renderPlotly({
        dplyr::select(GENE_MATRIX,
            ensembl, symbol,
            x = all_of(input$athero_dataset_sel),
            y = all_of(input$cancer_dataset_sel)
        ) %>%
            plot_ly(x = ~x, y = ~y,
                hoverinfo = "text", type = "scatter", mode = "markers",
                text = ~str_glue("{symbol} [{ensembl}]"),
                source = "gene_quad"
            ) %>%
            layout(
                xaxis = list(title = DATASET_NAME_MAP$display[input$athero_dataset_sel]),
                yaxis = list(
                    title = DATASET_NAME_MAP$display[input$cancer_dataset_sel],
                    scaleanchor = "x"
                )
            )
    })

    output$pathway_matrix_heatmap <- renderPlotly({
        mtx <- PATHWAY_DIM[[input$pathway_db_sel]]$scaled %>%
            `rownames<-`(
                PATHWAY_MATRIX[[input$pathway_db_sel]]$Description %>%
                    str_remove("^GO_|^HALLMARK_")
            ) %>%
            `colnames<-`(DATASET_NAME_MAP$display[colnames(.)])
        plot_ly(
            x = colnames(mtx),
            y = rownames(mtx),
            z = mtx,
            colors = colorRamp(c("blue", "yellow", "red")),
            type = "heatmap"
        ) %>%
            layout(
                xaxis = list(title = ""),
                yaxis = list(title = "")
            )
    })
    
    output$pathway_dist_heatmap <- renderPlotly({
        mtx <- PATHWAY_DIM[[input$pathway_db_sel]]$dist %>%
            as.matrix() %>%
            `colnames<-`(DATASET_NAME_MAP$display[colnames(.)]) %>%
            `rownames<-`(DATASET_NAME_MAP$display[rownames(.)]) %>%
            `[`(
                PATHWAY_DIM[[input$pathway_db_sel]]$cluster$order,
                PATHWAY_DIM[[input$pathway_db_sel]]$cluster$order
            )
        odr <- as_tibble(mtx) %>%
            dplyr::mutate(across(everything(), rank)) %>%
            dplyr::mutate(x = colnames(.)) %>%
            pivot_longer(!x,
                names_to = "y",
                values_to = "rank"
            ) %>%
            dplyr::mutate(rank = rank - 1) %>%
            dplyr::filter(rank > 0, rank <= 3)
        plot_ly(
            x = colnames(mtx),
            y = rownames(mtx),
            z = mtx,
            colors = colorRamp(c("blue", "yellow", "red")),
            type = "heatmap"
        ) %>%
            add_trace(x = ~x, y = ~y, text = ~rank,
                data = odr,
                hoverinfo = "none",
                type = "scatter", mode = "text"
            ) %>%
            layout(
                xaxis = list(title = ""),
                yaxis = list(title = "")
            )
    })
    
    output$pathway_pca <- renderPlotly({
        dplyr::mutate(PATHWAY_DIM[[input$pathway_db_sel]]$pca,
            id = DATASET_NAME_MAP$display[id]
        ) %>%
            plot_ly(x = ~PC_1, y = ~PC_2, text = ~id,
                hoverinfo = "skip", type = "scatter", mode = "markers"
            ) %>%
            add_annotations(text = ~unname(id)) %>%
            config(editable = TRUE)
    })
    
    output$pathway_umap <- renderPlotly({
        dplyr::mutate(PATHWAY_DIM[[input$pathway_db_sel]]$umap,
            id = DATASET_NAME_MAP$display[id]
        ) %>%
            plot_ly(x = ~UMAP_1, y = ~UMAP_2, text = ~id,
                hoverinfo = "skip", type = "scatter", mode = "markers"
            ) %>%
            add_annotations(text = ~unname(id)) %>%
            config(editable = TRUE)
    })
    
    pathway_sel_all <- reactiveVal()
    pathway_sel_cur <- reactiveVal()
    observeEvent({
        input$athero_dataset_sel
        input$cancer_dataset_sel
        input$pathway_db_sel
        event_data("plotly_doubleclick", "pathway_quad")
    }, {
        pathway_sel_cur(NULL)
        pathway_sel_all(NULL)
    })
    observeEvent(event_data("plotly_click", "pathway_quad"), {
        rst <- event_data("plotly_click", "pathway_quad")$customdata
        if (!is.null(rst)){
            pathway_sel_cur(rst)
        }
    })
    observeEvent(pathway_sel_cur(), {
        if (length(pathway_sel_cur()) == 0){
            pathway_sel_all(NULL)
            logging("clean ... # nevel run here")
        }else if (pathway_sel_cur() %in% pathway_sel_all()){
            pathway_sel_all(setdiff(
                pathway_sel_all(),
                pathway_sel_cur()
            ))
        }else {
            pathway_sel_all(c(
                pathway_sel_all(),
                pathway_sel_cur()
            ))
        }
        logging("pathway_sel_all = ", paste(pathway_sel_all(), collapse = ","))
    })
    output$pathway_quad <- renderPlotly({
        df <- dplyr::select(PATHWAY_MATRIX[[input$pathway_db_sel]],
            ID,
            Description,
            x = all_of(input$athero_dataset_sel),
            y = all_of(input$cancer_dataset_sel)
        ) %>%
            dplyr::filter(!is.na(x), !is.na(y)) %>%
            dplyr::mutate(Description =
                str_remove(Description, "^GO_|^HALLMARK_")
            )
        p <- plot_ly(df, x = ~x, y = ~y, text = ~Description,
            hoverinfo = "text", type = "scatter", mode = "markers",
            marker = list(
                color = ifelse(df$ID %in% pathway_sel_all(), "red", "black")
            ),
            source = "pathway_quad", customdata = df$ID
        ) %>%
            layout(
                xaxis = list(title = DATASET_NAME_MAP$display[input$athero_dataset_sel]),
                yaxis = list(
                    title = DATASET_NAME_MAP$display[input$cancer_dataset_sel],
                    scaleanchor = "x"
                )
            ) %>%
            config(editable = TRUE)
        if (length(pathway_sel_all()) != 0){
            p <- add_annotations(p,
                x = ~x, y = ~y,
                text = ~Description,
                data = dplyr::filter(df, ID %in% pathway_sel_all()),
                captureevents = TRUE
            )
        }
        p
    })
    output$pathway_summary1 <- renderPlotly({
        if (length(pathway_sel_cur()) == 0)
            return(NULL)
        df <- PATHWAY_MATRIX[[input$pathway_db_sel]] %>%
            dplyr::filter(ID == pathway_sel_cur()) %>%
            pivot_longer(!c(ID, Description),
                names_to = "dataset",
                values_to = "NES"
            ) %>%
            dplyr::filter(!is.na(NES)) %>%
            dplyr::mutate(
                label = DATASET_NAME_MAP$display[dataset],
                class = case_when(
                    dataset %in% ATHERO_DATASET ~ "athero",
                    dataset %in% CANCER_DATASET ~ "cancer",
                    TRUE ~ "NA"
                ),
                class = factor(class)
            )
        plot_ly(df, x = ~as.numeric(class), y = ~NES, split = ~class,
            type = "box", hoverinfo = "skip",
            showlegend = FALSE
        ) %>%
            add_markers(x = jitter(as.numeric(df$class)),
                text = ~label,
                hoverinfo = "text",
                showlegend = FALSE
            ) %>%
            layout(
                title = pathway_sel_cur(),
                xaxis = list(
                    title = "",
                    ticktext = list("athero", "cancer"),
                    tickvals = list(1, 2)
                )
            )
    })
}