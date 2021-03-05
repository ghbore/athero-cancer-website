#! /usr/bin/env Rscript

dashboardPagePlus(
    header = dashboardHeaderPlus(
        title = 'Athero ~ Cancer',
        enable_rightsidebar = TRUE
    ),
    sidebar = dashboardSidebar(
        sidebarMenu(
            menuItem('Overview', tabName = 'overview', icon = icon('globe')),
            menuItem('Datasets',
                menuSubItem('Athero', tabName = 'atheroDataSet'),
                menuSubItem('Cancer', tabName = 'cancerDataSet'),
                icon = icon('database'),
                startExpanded = TRUE
            ),
            menuItem('Gene level', tabName = 'geneLevel', icon = icon('dna')),
            menuItem('Pathway level', tabName = 'pathwayLevel', icon = icon('network-wired')),
            NULL
        )
    ),
    body = dashboardBody(
        tags$style(type = "text/css", 
            "#overview_dist,#overview_matrix,#quad_pathway {height: calc(100vh - 120px) !important;}"
        ),

        tabItems(
            tabItem('overview',
                tabsetPanel(
                    tabPanel('distant', plotlyOutput('overview_dist')),
                    tabPanel('term-matrix', plotlyOutput('overview_matrix')),
                    type = 'tabs'
                )
            ),

            tabItem('atheroDataSet',
                h3(textOutput('athero_dataset_title')),
                textOutput('athero_dataset_info'),
                # plotlyOutput('athero_dataset_volcano'),
                NULL
            ),

            tabItem('cancerDataSet',
                h3(textOutput('cancer_dataset_title')),
                textOutput('cancer_dataset_info'),
                # plotlyOutput('cancer_dataset_volcano'),
                NULL
            ),

            tabItem('geneLevel',
                # plotlyOutput('quad_gene'),
                NULL
            ),

            tabItem('pathwayLevel',
                tabsetPanel(
                    tabPanel('plot', plotlyOutput('quad_pathway')),
                    tabPanel('table', dataTableOutput('quad_pathway_tbl')),
                    type = 'tabs'
                )
            )
        )
    ),
    rightsidebar = rightSidebar(
        .items = list(
            selectInput('athero_dataset_sel', 'athero dataset', ATHERO_DATASET, ATHERO_DATASET_DEFAULT),
            selectInput('cancer_dataset_sel', 'cancer dataset', CANCER_DATASET, CANCER_DATASET_DEFAULT),
            selectInput('pathway_db_sel', 'pathway database', GENESET, GENESET_DEFAULT),
            NULL
        )
    ),
    title = 'Athero ~ Cancer',
    skin = 'black',
    NULL
)
