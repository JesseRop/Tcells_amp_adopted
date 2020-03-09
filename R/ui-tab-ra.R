tabPanel(
  title = "Effect of Beta 2 integrin in T Cell differentiation",
  value = "data",
  
  titlePanel("Single-cell RNA-seq in Beta 2 integring Knockouts compared to wildtype cells"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("UMAP clusters", "Top cluster markers conserved between conditions", "Label populations", "Dotplot for gene expression between conditions and across cell types", "Differential gene expression between conditions per cell type", "Gene expression visualization between conditions across cell types")),
      
      # Display depending on which graph type is selected
      
      conditionalPanel(condition = "input.graph_type == 'UMAP clusters'",
                       # Select PC to plot
                       sliderInput(inputId = "umap_dim", label = strong("Number of PCs from 1"), value = dim1, min = dim1, max = dim2, step = diff_dim),
                       sliderInput(inputId = "neighbours_dim", label = strong("Number of PCs from 1 (KNN clustering in FindNeighbours function)"), value = dim1, min = dim1, max = dim2, step = diff_dim),
                       sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res),
                       actionButton(inputId = "go", label = "Run")
                       
      ),
      conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between conditions'",
                       
                       sliderInput(inputId = "marker_genes_no", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
                       #selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = , multiple = F),
                       uiOutput("dyn_clusters"),
                       actionButton(inputId = "go_marker", label = "Search"),
                       selectInput(inputId = "select_markers_umap", label = strong("Select conserved cluster markers:"), choices = all_genes_final, multiple = T)
                       
      ),
      conditionalPanel(condition = "input.graph_type == 'Label populations'",
                       #Select PC to plot
                       uiOutput("cluster_annot"),
                       actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
      ),
      conditionalPanel(condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'",
                       #Select PC to plot
                       selectInput(inputId = "select_markers_dotplot", label = strong("Select markers for dotplot:"),choices = all_genes_final, multiple = T)
      ),
      conditionalPanel(condition = "input.graph_type == 'Differential gene expression between conditions per cell type'",
                       #Select PC to plot
                       #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                       uiOutput("cluster_ids"),
                       radioButtons(inputId = "condition1", label = "1st Condition", choiceNames = c("Wildtype 1", "Wildtype 2", "Knock-out 1"), choiceValues = (c("WT1","WT2", "KO1"))),
                       radioButtons(inputId = "condition2", label = "2nd Condition", choiceNames = c("Wildtype 1", "Wildtype 2", "Knock-out 1"), choiceValues = (c("WT1","WT2", "KO1"))),
                       uiOutput("ct_de_plot_dyn"),   
                       sliderInput(inputId = "no_of_de_genes", label = strong("Number of top DE genes:"), value = 10, min = 1, max = 100, step = 1),
                       uiOutput("ct_de_table_dyn")
                       
      ),
      conditionalPanel(condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'",
                       #Select PC to plot
                       selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_final, multiple = F),
                       #actionButton(inputId = "go_heatmap", label = "Run"),
                       #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
      )
    ),
    
    # Output: Description, lineplot, and reference
    mainPanel(  
      conditionalPanel(
        condition = "input.graph_type == 'UMAP clusters'", 
        plotOutput("all_groups"), plotOutput("stim_vs_ctrl")),
      conditionalPanel(
        condition = "input.graph_type == 'Top cluster markers conserved between conditions'", 
        tableOutput("top_conserved_genes"), plotOutput("conserved_markers_umap")), 
      conditionalPanel(
        condition = "input.graph_type == 'Label populations'", 
        plotOutput("labelled_umap")),
      conditionalPanel(
        condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'", 
        plotOutput("marker_dotplot")),
      conditionalPanel(
        condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
        plotOutput("cell_type_plot", click = clickOpts(id ="plot_click")), verbatimTextOutput("click_info")),
      conditionalPanel(
        condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
        tableOutput("top_de_genes")),
      conditionalPanel(
        condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'", uiOutput("tb") )
      
      
    )
  )
)
  # #tabsetPanel(selected = "RNA-seq",
  #   
  #   tabPanel(
  #     "RNA-seq",
    #     
    #     # Show a plot of the generated distribution
    #     fluidPage(
    #       fluidRow(
    # 
    #         column(
    #           width = 12,
    #           div(
    #             radioButtons(
    #               inputId  = "rnaseq_cell_type",
    #               inline   = TRUE,
    #               label    = "Cell type:",
    #               choices  = possible_cell_types_rna,
    #               selected = "all"
    #             ),
    #             htmlOutput("tnse_marker_plot") #, height = "400px"),
    #             #style = "height: 400px;"
    #           )
    #         )
    #         
    #       )
    #         
    #     ),
    #     hr(),
    #     fluidRow(
    #       column(
    #         width = 5,
    #         htmlOutput("box_marker_plot_all", height = "500px")
    #       ),
    #       column(
    #         width = 7,
    #         DT::dataTableOutput("dg_table")#, height = "350px")
    #       )
    #     ),
    #     hr(),
    #     fluidRow(
    #       column(
    #         width = 10,
    #         h4("Bulk RNA-seq"),
    #         htmlOutput("bulk_dots", height = "500px")
    #       )
    #     ),
    #   hr(),
    #   fluidRow(column(
    #     width = 12,
    #     h2("Cluster annotations"),
    #     HTML("
    #          <p>Fibroblasts (CD45<sup>-</sup> Podoplanin<sup>+</sup>)</p>
    #          <ul>
    #          <li>SC-F1: CD34+ sublining fibroblasts</li>
    #          <li>SC-F2: HLA+ sublining fibroblasts</li>
    #          <li>SC-F3: DKK3+ sublining fibroblasts</li>
    #          <li>SC-F4: CD55+ lining fibroblasts</li>
    #          </ul>
    #          
    #          <p>Monocytes (CD45<sup>+</sup> CD14<sup>+</sup>)</p>
    #          <ul>
    #          <li>SC-M1: IL1B+ pro-inflammatory monocytes</li>
    #          <li>SC-M2: NUPR1+ monocytes</li>
    #          <li>SC-M3: C1QA+ monocytes</li>
    #          <li>SC-M4: IFN-activated monocytes</li>
    #          </ul>
    #          
    #          <p>T cells (CD45<sup>+</sup> CD3<sup>+</sup>)</p>
    #          <ul>
    #          <li>SC-T1: CCR7+ CD4+ T cells</li>
    #          <li>SC-T2: FOXP3+ Tregs</li>
    #          <li>SC-T3: PD-1+ Tph/Tfh</li>
    #          <li>SC-T4: GZMK+ CD8+ T cells</li>
    #          <li>SC-T5: GNLY+ GZMB+ CTLs</li>
    #          <li>SC-T6: GZMK+/GZMB+ T cells</li>
    #          </ul>
    #          
    #          <p>B cells (CD45<sup>+</sup> CD3<sup>-</sup> CD19<sup>+</sup>)</p>
    #          <ul>
    #          <li>SC-B1: IGHD+ CD270 naive B cells</li>
    #          <li>SC-B2: IGHG3+ CD27- memory B cells</li>
    #          <li>SC-B3: autoimmune-associated cells (ABC)</li>
    #          <li>SC-B4: Plasmablasts</li>
    #          </ul>
    #          ")
    #     ))
    # 
    # ), # tabPanel
# 
#     tabPanel(
#       "Mass cytometry",
#       
#       # Show a plot of the generated distribution
#       fluidPage(
#         fluidRow(
#           column(
#             width = 12,
#             div(
#               radioButtons(
#                 inputId  = "cytof_cell_type",
#                 inline   = TRUE,
#                 label    = "Cell type:",
#                 choices  = possible_cell_types_cytof,
#                 selected = "Fibroblast"
#               ),
#               selectInput(
#                 inputId = "cytof_marker",
#                 label = "Protein",
#                 choices = sort(c("CD45", "CD19", "RANKL", "CD64", "CD16",
#                                  "CD8a", "FAP", "CD20", "CD45RO", "CD38",
#                                  "PD.1", "CD14", "CD69", "CXCR5", "CD4",
#                                  "Podoplanin", "CD3", "CD11c", "FcRL4",
#                                  "CD138", "CD90", "CCR2", "Cadherin.11",
#                                  "FoxP3", "CD34", "CD146", "IgA", "TCRgd",
#                                  "ICOS", "CD66b", "IgM", "VE.Cadherin",
#                                  "HLA.DR", "IgD", "VCAM.1")),
#                 selected = "CD90",
#                 multiple = FALSE,
#                 selectize = TRUE,
#                 width = NULL,
#                 size = NULL
#               ),
#               htmlOutput("tnse_cytof")
#             )
#           )
#         ), # fluidRow
#         hr(),
#         fluidRow(
#           column(
#             width = 12,
#             DT::dataTableOutput("cytof_table", height = "350px")
#           )
#         )
#       ) # fluidPage
#     ) # tabPanel
#     
#   ) # tabsetPanel
#   
# ) # tabPanel
