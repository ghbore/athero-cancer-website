#! /usr/bin/env Rscript

function (input, output, session){
    output$athero_dataset_title <- renderText({
        input$athero_dataset_sel
    })
    output$athero_dataset_info <- renderText({
        input$athero_dataset_sel
    })

    output$cancer_dataset_title <- renderText({
        input$cancer_dataset_sel
    })
    output$cancer_dataset_info <- renderText({
        input$cancer_dataset_sel
    })


    output$overview_dist <- renderPlotly({
        cat(file=stderr(), 'rendering overview dist plot for', input$pathway_db_sel, '\n')
        mtx <- DATASET_SIMILARITY[['pathway']][[input$pathway_db_sel]][['dist']]
        hi <- tibble(y = rownames(mtx),
            x = names(ATHERO_DATASET)[ATHERO_DATASET == input$athero_dataset_sel],
            val = mtx[, x[1]]
        ) %>%
            dplyr::filter(val > 0) %>%
            arrange(val) %>%
            mutate(idx = 1:n()) %>%
            dplyr::filter(idx <= 3)
        plot_ly(x = colnames(mtx), y = rownames(mtx), z = mtx,
            colors = colorRamp(c('blue', 'yellow', 'red')),
            type = 'heatmap'
        ) %>%
            add_trace(x = ~x, y= ~y, text =~ idx, data = hi, type = 'scatter', mode = 'text')
    })
    output$overview_matrix <- renderPlotly({
        cat(file=stderr(), 'rendering overview matrix plot for', input$pathway_db_sel, '\n')
        DATASET_SIMILARITY[['pathway']][[input$pathway_db_sel]][['matrix']] %>%
            plot_ly(x = colnames(.), y = rownames(.), z = .,
                colors = colorRamp(c('blue', 'yellow', 'red')),
                type = 'heatmap'
            )
    })

    pathway_comparison_tbl <- reactive({
        cat(file=stderr(), 'calculating pathway comparison table for', input$athero_dataset_sel, input$cancer_dataset_sel, input$pathway_db_sel, '\n')
        pathway_comparison(
            ATHERO_PATHWAY[[input$pathway_db_sel]][[input$athero_dataset_sel]],
            CANCER_PATHWAY[[input$pathway_db_sel]][[input$cancer_dataset_sel]]
        )
    })
    output$quad_pathway_tbl <- renderDataTable({ pathway_comparison_tbl() })
    output$quad_pathway <- renderPlotly({
        pathway_comparison_tbl() %>%
            transmute(
                x = `NES:athero`, 
                y = `NES:cancer`,
                label = Description %>% str_remove('^GO_|^HALLMARK_') %>% str_replace_all('_', ' ')
            ) %>%
            plot_ly(x = ~x, y = ~y, text = ~paste('term:', label), type='scatter', mode='markers') %>%
            layout(
                xaxis = list(title = names(ATHERO_DATASET[ATHERO_DATASET == input$athero_dataset_sel])), 
                yaxis = list(title = names(CANCER_DATASET[CANCER_DATASET == input$cancer_dataset_sel]), 
                    scaleanchor = "x")
            )
    })
}