
source("tcell_tdo_changes_libs_raw_objects.R")
##Creating the ui object for the user interface
##Creating the ui object for the user interface

ui = fluidPage(theme = shinytheme("united"),
               navbarPage("Role of CD18 in γδ T cells",
                          source(file.path("introduction_page.R"), local = TRUE)$value,
                          
                          ##The first tab which is for the differential expression visualization
                          tabPanel("Differential expression",
                                   sidebarLayout(
                                     sidebarPanel(
                                       
                                       ##The options (subtitles) for the first tab addressing differential expression visualization)
                                       radioButtons(inputId = "de_panel", label = strong("single cell RNAseq analysis of wild type (WT1/WT2) and CD18 knockout (KO) lung gamma delta T cells"), choices = c("Gene expression visualization across cell types and clusters", "Dotplot of gene expression across cell types and clusters", "Differential gene expression between cell types in individual clusters" )),
                                       
                                       
                                       #Option/subtitle for comparison of gene expression between WT and KO across cell types using umap and violin plots
                                       conditionalPanel(condition = "input.de_panel == 'Gene expression visualization across cell types and clusters'",
                                                        #Select PC to plot
                                                        selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_common_in_all_groups, multiple = F),
                                                        #actionButton(inputId = "go_heatmap", label = "Run")
                                                        #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                                                        fluidPage(
                                                          tags$hr(),
                                                          tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                                 #titlePanel("includeText"),
                                                                 #fluidRow(
                                                                 #column(12,
                                                                 "Comparison of gene expression between WT and KO across cell types using umap and violin plots"
                                                                 #))
                                                          ),
                                                          tags$hr()
                                                        )
                                                        
                                       ),
                                       
                                       #Option/subtitle for comparison of gene expression between WT and KO across cell types using dotplots
                                       conditionalPanel(condition = "input.de_panel == 'Dotplot of gene expression across cell types and clusters'",
                                                        #Select PC to plot
                                                        selectInput(inputId = "select_markers_dotplot", label = strong("Select markers for dotplot:"),choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes),
                                                        fluidPage(
                                                          tags$hr(),
                                                          tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                                 "Comparison of gene expression between WT and KO across cell types using dotplots. Red if for Wildtype and Blue for CD18KO with the increasing in intensity correlating with the level of gene expression in the population"
                                                          ),
                                                          tags$hr()
                                                        )
                                       ),
                                       
                                       #Option/subtitle for comparison of gene expression between WT and KO across cell types using scatter plots and tables looking at all the genes in a cell type
                                       conditionalPanel(condition = "input.de_panel == 'Differential gene expression between cell types in individual clusters'",
                                                        #Select PC to plot
                                                        #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                                        uiOutput("cluster_ids"),
                                                        #radioButtons(inputId = "condition2", label = "WT group to compare to KO", choiceNames = c("Wildtype 1", "Wildtype 2"), choiceValues = (c("WT1","WT2"))),
                                                        conditionalPanel(condition = "input.tabslctd == 'gg'",
                                                                         uiOutput("ct_de_plot_dyn"),
                                                                         fluidPage(
                                                                           tags$hr(),
                                                                           tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                                           "Comparison of gene expression between WT and CD18 KO across cell types using scatter plots looking at all the genes in a cell type"
                                                                             ),
                                                                           tags$hr()
                                                                         )
                                                        ),
                                                        conditionalPanel(condition = "input.tabslctd == 'tbl'",
                                                                         sliderInput(inputId = "no_of_de_genes", label = strong("Number of top DE genes:"), value = 10, min = 1, max = 100, step = 1),
                                                                         downloadButton("downloadData", "Download table of significantly DE genes"),
                                                                         fluidPage(
                                                                           tags$hr(),
                                                                           tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                                           "Listing differential expressed genes (adjusted P value <0.05) WT and CD18 KO across cell populations"
                                                                             ),
                                                                           tags$hr()
                                                                         )
                                                        )
                                       )
                                       
                                     ),
                                     
                                     ##Output allocations to contain graphs and outputs generated by server for tab 1
                                     mainPanel(  
                                       conditionalPanel(
                                         condition = "input.de_panel == 'Gene expression visualization across cell types and clusters'", 
                                         uiOutput("tb") ),
                                       conditionalPanel(
                                         condition = "input.de_panel == 'Dotplot of gene expression across cell types and clusters'", 
                                         plotOutput("marker_dotplot")),
                                       conditionalPanel(
                                         condition = "input.de_panel == 'Differential gene expression between cell types in individual clusters'", 
                                         uiOutput("de_outputs"))
                                        
                                     )
                                   )
                          ),
                          
                          ##The second tab which is for cluster adjustment
                          tabPanel("Cluster adjustment",
                                   sidebarLayout(
                                     sidebarPanel(
                                       
                                       ##The options (subtitles) for the second tab clustering
                                       radioButtons(inputId = "graph_type", label = strong("Cluster adjustment and visualization"), choices = c("UMAP clusters", "Top cluster markers conserved between WT and KO", "Label clusters")),
                                       
                                       #2.1 Option/subtitle for adjusting resolution to give more clusters
                                       conditionalPanel(condition = "input.graph_type == 'UMAP clusters'",
                                                        # Select PC to plot
                                                        # sliderInput(inputId = "umap_dim", label = strong("Number of PCs from 1"), value = dim1, min = dim1, max = dim2, step = diff_dim),
                                                        # sliderInput(inputId = "neighbours_dim", label = strong("Number of PCs from 1 (KNN clustering in FindNeighbours function)"), value = dim1, min = dim1, max = dim2, step = diff_dim),
                                                        sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F),
                                                        #actionButton(inputId = "go", label = "Run")
                                                        fluidPage(
                                                          tags$hr(),
                                                          tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                                 "Adjusting resolution (0.15, 0.25, 0.35, 0.45 or 0.55) to set the clustering granularity. This options allow the subdivision of clusters further into sub-populations to obtain rare cell types and assess differential expression in these"
                                                          ),
                                                          tags$hr()
                                                        )
                                       ),
                                       
                                       #2.2 Option/subtitle for listing top cluster markers for which can then be used in the labelling of the generated populations
                                       conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between WT and KO'",
                                                        uiOutput("dyn_clusters"),
                                                        sliderInput(inputId = "marker_genes_no", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
                                                        
                                                        #selectInput(inputId = "select_markers_umap", label = strong("select marker to visualize in clusters:"), choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes),
                                                        uiOutput("top_markers_umap"),
                                                        
                                                        fluidPage(
                                                          tags$hr(),
                                                          tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                                 "Listing top cluster markers which can then be used in the labelling of the generated populations"
                                                          ),
                                                          tags$hr()
                                                        )
                                       #)
                                                        # conditionalPanel(condition = "input.cm_tabslctd == 'marker_tbl'",
                                                        #                  uiOutput("top_markers_umap"),
                                                        #                  fluidPage(
                                                        #                    tags$hr(),
                                                        #                    tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                        #                    "Listing top cluster markers which can then be used in the labelling of the generated populations"
                                                        #                    ),
                                                        #                    tags$hr()
                                                        #                  )
                                                        # ),
                                                        # conditionalPanel(condition = "input.cm_tabslctd == 'marker_umap'",
                                                        #                  selectInput(inputId = "select_markers_umap", label = strong("select marker to visualize in clusters:"), choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes),
                                                        #                  fluidPage(
                                                        #                    tags$hr(),
                                                        #                    tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                        #                    "Visualize the top cluster markers to confirm that they are expressed in the appropriate populations"
                                                        #                    ),
                                                        #                    tags$hr()
                                                        #                  )
                                                        # )

                                       ),
                                       
                                       #2.3 Option/subtitle for allowing the labelling of the clusters
                                       conditionalPanel(condition = "input.graph_type == 'Label clusters'",
                                                        #Select PC to plot
                                                        uiOutput("cluster_annot"),
                                                        #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                                        fluidPage(
                                                          tags$hr(),
                                                          tags$p(style = "font-family:Comic Sans MS;color:navy",
                                                          "Labelling populations using the top cluster markers"
                                                          ),
                                                          tags$hr()
                                                        )
                                       )
                                     ),
                                     
                                     ##Output allocations to contain graphs and outputs generated by server for tab 2
                                     mainPanel(
                                       
                                       conditionalPanel(
                                         condition = "input.graph_type == 'UMAP clusters'", 
                                         withSpinner(plotOutput("all_groups")), withSpinner(plotOutput("stim_vs_ctrl"))),
                                       conditionalPanel(
                                         condition = "input.graph_type == 'Top cluster markers conserved between WT and KO'", 
                                         uiOutput("cluster_markers_outputs")),
                                       conditionalPanel(
                                         condition = "input.graph_type == 'Label clusters'", 
                                         withSpinner(plotOutput("labelled_umap")))
                                       
                                       
                                     )
                                   )
                          )
                          
               )
)

##server function to compute the outputs
server = function(input, output) {
  
  
  ######################
  ### Selection of seurat object and plotting UMAP plots under subtitle 2.1
  ######################
  ##Selecting the Seurat object
  umap_clusters = reactive({
    
    tcells_combined_umap_list_res[[(input$clusters_res * 10)-0.5]]
  })
  
  ##Plotting UMAP plots for clustering
  output$all_groups = renderPlot({
    
    p1 = DimPlot(umap_clusters(), reduction = "umap", group.by = "sample", label.size = 6)
    p2 = DimPlot(umap_clusters(), reduction = "umap", label = T, label.size = 6)
    plot_grid(p1, p2)
    
  })
  
  ##Plotting UMAP plots showing the stimulated and control umaps side by side
  output$stim_vs_ctrl = renderPlot({
    DimPlot(umap_clusters(), reduction = "umap", split.by = "sample")
  })
  ######################
  ### END Selection of seurat object and plotting UMAP plots under subtitle 2.1
  ######################
  
  
  ######################
  ### Cluster markers plots and tables under subtitle 2.2
  ######################
  ##Dynamic input field for selecting cluster to plot table of markers
  output$dyn_clusters <- renderUI({
    umap_clusters = umap_clusters()
    selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_clusters$seurat_clusters), multiple = F)
  })
  
  ##Displaying table of cluster markers for annotating cell types
  cluster_markers = reactive({
    head(tcells_combined_clusters_tables_res[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$marker_genes_cluster) + 1)]] %>% rownames_to_column(var = 'genes'), n = input$marker_genes_no)

  })
  
  output$top_conserved_genes = renderDataTable({
    datatable(cluster_markers() %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)))
  }, options = list(scrollX = TRUE))
  
  ##Preparing and plotting UMAP cluster markers for annotating cell types
  umap_cluster_modified_rna = reactive({
    umap_cluster_modified_ul = umap_clusters()
    DefaultAssay(umap_cluster_modified_ul) = "RNA"
    umap_cluster_modified_ul
  })
  
  # output$conserved_markers_umap = renderPlot({
  #   FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
  # })
  
  ##Conditional tabsets for the DE table and ggplot UI change for input triggered by selecting tab 
  # output$cluster_markers_outputs <- renderUI({
  #   tabsetPanel(tabPanel("Top cluster markers table", value = "marker_tbl", 
  #                        withSpinner(dataTableOutput("top_conserved_genes"))),
  #               tabPanel("Cluster marker UMAPs", value = "marker_umap",
  #                        withSpinner(plotOutput("conserved_markers_umap"))), id = "cm_tabslctd")
  # })
  
  output$top_markers_umap <- renderUI({
    a = head(cluster_markers()[,1])
    a = factor(a, levels = c(a))
    print(a)
    selectInput(inputId = "select_markers_umap", label = strong("select marker to visualize in clusters:"), choices = all_genes_common_in_all_groups, multiple = T, selected = a)
    
    })
  output$conserved_markers_umap = renderPlot({
    FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
  })
  
  #Conditional tabsets for the DE table and ggplot UI change for input triggered by selecting tab
  output$cluster_markers_outputs <- renderUI({
    fluidPage(
      fluidRow(
        h4("Cluster markers table"),
        withSpinner(dataTableOutput("top_conserved_genes"))
      ),
      fluidRow(
        h4("UMAP of 6 top markers"),
        withSpinner(plotOutput("conserved_markers_umap"))

      )
    )

  })
  
  ######################
  ### END Cluster markers plots and tables under subtitle 2.2
  ######################
  
  
  ######################
  ### Labelling clusters under subtitle 2.3
  ######################
  ##Generating dynamic fields for labelling UMAP clusters and initializing the fields with placeholders
  output$cluster_annot <- renderUI({
    umap_cluster_modified1 = umap_cluster_modified_rna()
    
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(cluster_names)){
      umap_clusters = umap_clusters()
      lapply(0:(length(unique(umap_clusters$seurat_clusters))-1), function(x) {
        textInput(inputId = paste("labeller", x), label = strong(paste("Label the cluster", x, "based on marker genes")), value = cluster_names[x+1])
      })
      
    } else {
      umap_clusters = umap_clusters()
      lapply(0:(length(unique(umap_clusters$seurat_clusters))-1), function(x) {
        textInput(inputId = paste("labeller", x), label = strong(paste("Label the cluster", x, "based on marker genes")), value = x)
      })
    }
  })
  
  ##Storing names decided on by the researchers based on optimal clustering to start off the differential expression visualization
  annotation = reactiveValues(annot = cluster_names)
  
  ##Observer to allow updating of cluster names dynamically as typed in
  observe({
    umap_cluster_modified = umap_cluster_modified_rna()
    
    req(unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    })))
    annotation$annot = unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    }))
  })
  
  ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
  output$cluster_ids <- renderUI({
    umap_names = annotation$annot
    #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
    umap_cluster_modified1 = umap_cluster_modified_rna()
    
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
      umap_names = annotation$annot
      selectInput(inputId = "select_cell_type", label = strong("Select cell type to compare gene expression across conditions:"),choices = umap_names, multiple = F)
      
    } else {
      selectInput(inputId = "select_cell_type", label = strong("Select cell type to compare gene expression across conditions:"),choices = unique(umap_cluster_modified1$seurat_clusters), multiple = F)
    }
    
  })
  
  ##Renaming clusters
  umap_cluster_modified_ren_reo = reactive({
    umap_cluster_modified1 = umap_cluster_modified_rna()
    umap_names = annotation$annot
    
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
      names(umap_names) <- levels(umap_cluster_modified1)
      umap_cluster_modified1 <- RenameIdents(umap_cluster_modified1, umap_names)
      
    } else {
      umap_cluster_modified1
    }
    
  })
  
  umap_cluster_modified_umap = reactive({
    DimPlot(umap_cluster_modified_ren_reo(), label = TRUE, label.size = 6)
    
  })
  output$labelled_umap = renderPlot({
    umap_cluster_modified_umap()
  })
  ######################
  ###END Labelling clusters under subtitle 2.3
  ######################
  
  
  ######################
  ### Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
  ######################
  stim_markers = reactive({
    
    umap_cluster_modified = umap_cluster_modified_ren_reo()
    umap_cluster_modified$celltype.group <- paste(Idents(umap_cluster_modified), umap_cluster_modified$group, sep = "_")
    umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
    Idents(umap_cluster_modified) <- "celltype.group"
    umap_cluster_modified
  })
  
  
  #Functions to update differentially expressed genes
  output$de_stim_vs_ctrl_um = renderPlot({
    
    FeaturePlot(stim_markers(), features = input$de_genes, split.by = "sample", max.cutoff = 3,cols = c("grey", "red"))
  })
  
  output$de_stim_vs_ctrl_vp = renderPlot({
    
    VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE)  
    
  })
  
  ## Define two panels for UMAP and violin plots
  output$tb <- renderUI({
    tabsetPanel(tabPanel("UMAP plot", 
                         withSpinner(plotOutput("de_stim_vs_ctrl_um"))),
                tabPanel("Violin plot", 
                         withSpinner(plotOutput("de_stim_vs_ctrl_vp")))) 
  })
  ######################
  ### END Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
  ######################
  
  
  ######################
  ### Differential expression using dotplot under subtitle 1.2
  ######################
  ##Dotplot for DE comparison between KO and WT across cell types
  output$marker_dotplot = renderPlot({
    DotPlot(umap_cluster_modified_ren_reo(), features = input$select_markers_dotplot, cols = c("blue", "red"), dot.scale = 6, split.by = "group") + RotatedAxis()
  })
  ######################
  ### Differential expression using dotplot under subtitle 1.2
  ######################
  
  
  ######################
  ### Differential expression using ggplot and tables under subtitle 1.3
  ######################
  
  ##Retrieving table for DE expression from precomputed list
  genes_in_de_order = reactive({
    umap_names = annotation$annot
    umap_cluster_modified1 = umap_cluster_modified_rna()
    
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
      tcells_combined_de_tables[[(input$clusters_res * 10)-0.5]][[as.numeric(match(input$select_cell_type,umap_names))]]
      
    } else {
      tcells_combined_de_tables[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$select_cell_type) + 1)]]
    }
    
    
  })
  
  output$ct_de_plot_dyn <- renderUI({
    actionButton(inputId = "ct_de", label = paste("Plot differential expression (DE) of all genes in ", input$select_cell_type , sep = ""))
  })
  
  ##Retrieving table for DE expression from precomputed list
  top_de_g = reactive({
    head(genes_in_de_order()  %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05), n = input$no_of_de_genes)  %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
  })
  
  ##plotting DE table
  output$top_de_genes = renderDataTable({
    datatable(top_de_g())
  })
  
  ##Allowing for download of DE table
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("differentially_expressed_genes_in",input$select_cell_type, "KO","VS","WT", ".csv", sep = "_")
      
    },
    content = function(file) {
      write.csv(genes_in_de_order() %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05)  %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)), file)
    }
  )
  
  ##Generating a table for each of the celltype for average DE expression of all genes for ggplot
  cell_type_de = eventReactive(input$ct_de,{
    
    cells_type <- subset(umap_cluster_modified_ren_reo(), idents = input$select_cell_type)
    #Idents(cells_type) <- "sample"
    Idents(cells_type) <- "group"
    avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
    avg.cells$gene <- rownames(avg.cells)
    avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene)) %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = "Gene names  (primary )")) %>% dplyr::distinct(., gene, .keep_all = T)%>% select(gene, KO, WT , `Protein names`, uniprot)
    avg.cells
    
  })
  
  ##preparing ggplot for average DE expression for genes above
  cell_type_de_plot = eventReactive(input$ct_de,{
    theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=9, fontface="italic")))
    ggplot(data=cell_type_de(), aes_string("KO", "WT")) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob)
    
    
  })
  
  ##plotting ggplot for average DE expression for genes above
  output$cell_type_plot = renderPlot({
    theme_set(theme_cowplot())
    cell_type_de_plot()
  })
  
  ##Conditional tabsets for the DE table and ggplot UI change for input triggered by selecting tab 
  output$de_outputs <- renderUI({
    tabsetPanel(tabPanel("DE ggplot", value = "gg", 
                         withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))), dataTableOutput("click_info")),
                tabPanel("DE table", value = "tbl",
                         withSpinner(dataTableOutput("top_de_genes"))), id = "tabslctd")
  })
  
  ##Displaying further details upon clicking points
  displayed_text <- reactive({
    req(input$plot_click)
    nearPoints(cell_type_de(), input$plot_click)
    
  })
  
  ##Displaying table with gene details upon click of point in DE ggplot
  output$click_info <- renderDataTable({
    req(displayed_text())
    
    #cat("Name\n")
    displayed_text = displayed_text()
    displayed_text 
  }, escape = F)
  ######################
  ### END Differential expression using ggplot and tables under subtitle 1.3
  ######################
  
  
}


shinyApp(ui = ui, server = server)
