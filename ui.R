#! /usr/bin/env Rscript

dashboardPage(
    header = dashboardHeader(
        title = "Athero ~ Cancer"
    ),
    sidebar = dashboardSidebar(
        sidebarMenu(
            menuItem("Overview", tabName = "overview", icon = icon("globe")),
            menuItem("Datasets",
                menuSubItem("Athero", tabName = "atheroDataSet"),
                menuSubItem("Cancer", tabName = "cancerDataSet"),
                icon = icon("database"),
                startExpanded = TRUE
            ),
            menuItem("Gene level", tabName = "geneLevel", icon = icon("dna")),
            menuItem("Pathway level", tabName = "pathwayLevel",
                icon = icon("network-wired")
            ),
            NULL
        )
    ),
    body = dashboardBody(
        tags$style(type = "text/css",
            "
            #athero_dataset_volcano,#athero_dataset,#athero_pathway_volcano,
            #cancer_dataset_volcano,#cancer_dataset,#cancer_pathway_volcano,
            #gene_matrix_heatmap,#gene_dist_heatmap,
            #gene_pca,#gene_umap,#gene_quad,
            #pathway_matrix_heatmap,#pathway_dist_heatmap,
            #pathway_pca,#pathway_umap,#pathway_quad
                {height: calc(100vh - 120px) !important;}
            "
        ),

        tabItems(
            tabItem("overview",
                span("here goes to be three-layer network")
            ),

            tabItem("atheroDataSet",
                tabsetPanel(
                    tabPanel("summary", span("here goes to be a summary")),
                    tabPanel("volcano", plotlyOutput("athero_dataset_volcano")),
                    tabPanel("table", DTOutput("athero_dataset")),
                    tabPanel("pathway", plotlyOutput("athero_pathway_volcano")),
                    type = "tabs"
                )
            ),

            tabItem("cancerDataSet",
                tabsetPanel(
                    tabPanel("summary", span("here goes to be a summary")),
                    tabPanel("volcano", plotlyOutput("cancer_dataset_volcano")),
                    tabPanel("table", DTOutput("cancer_dataset")),
                    tabPanel("pathway", plotlyOutput("cancer_pathway_volcano")),
                    type = "tabs"
                )
            ),

            tabItem("geneLevel",
                tabsetPanel(
                    tabPanel("dist", plotlyOutput("gene_dist_heatmap")),
                    tabPanel("PCA", plotlyOutput("gene_pca")),
                    tabPanel("UMAP", plotlyOutput("gene_umap")),
                    tabPanel("quadrant", plotlyOutput("gene_quad")),
                    type = "tabs"
                )
            ),

            tabItem("pathwayLevel",
                tabsetPanel(
                    tabPanel("matrix", plotlyOutput("pathway_matrix_heatmap")),
                    tabPanel("dist", plotlyOutput("pathway_dist_heatmap")),
                    tabPanel("PCA", plotlyOutput("pathway_pca")),
                    tabPanel("UMAP", plotlyOutput("pathway_umap")),
                    tabPanel("quadrant", 
                        column(width = 8, plotlyOutput("pathway_quad")),
                        column(width = 4, plotlyOutput("pathway_summary1"))
                    ),
                    type = "tabs"
                )
            )
        )
    ),
    controlbar = dashboardControlbar(
        selectInput("athero_dataset_sel", "athero dataset",
            ATHERO_DATASET, ATHERO_DATASET_DEFAULT),
        selectInput("cancer_dataset_sel", "cancer dataset",
            CANCER_DATASET, CANCER_DATASET_DEFAULT),
        selectInput("pathway_db_sel", "pathway database",
            GENESET, GENESET_DEFAULT),
        NULL
    ),
    title = "Athero ~ Cancer",
    skin = "black",
    NULL
)
