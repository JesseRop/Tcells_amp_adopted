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
                                       radioButtons(inputId = "de_panel", label = strong("scRNA-Seq analysis of wild type (WT1/WT2) VS CD18 knockout (KO) lung γδ T cells"), choices = c("Gene expression visualization across cell types and clusters", "Dotplot of gene expression across cell types and clusters", "Differential gene expression between cell types in individual clusters" )),
                                       
                                       
                                       #1.1 Option/subtitle for comparison of gene expression between WT and CD18 KO across cell types using umap and violin plots
                                       conditionalPanel(condition = "input.de_panel == 'Gene expression visualization across cell types and clusters'",
                                                        wellPanel(
                                                          #Select PC to plot
                                                          selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_common_in_all_groups, multiple = F, selected = "Cd163l1")
                                                          #actionButton(inputId = "go_heatmap", label = "Run")
                                                          #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                                                        ),
                                                        uiOutput("box_1_1")
                                       ),
                                       
                                       #1.2 Option/subtitle for comparison of gene expression between WT and KO across cell types using dotplots
                                       conditionalPanel(condition = "input.de_panel == 'Dotplot of gene expression across cell types and clusters'",
                                                        #Select PC to plot
                                                        wellPanel(
                                                          selectInput(inputId = "select_markers_dotplot", label = strong("Select markers for dotplot:"),choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes)
                                                        ),
                                                        uiOutput("box_1_2")
                                       ),
                                       
                                       #1.3 Option/subtitle for comparison of gene expression between WT and KO across cell types using scatter plots and tables looking at all the genes in a cell type
                                       conditionalPanel(condition = "input.de_panel == 'Differential gene expression between cell types in individual clusters'",
                                                        #Select PC to plot   
                                                        uiOutput("cluster_ids"),
                                                        conditionalPanel(condition = "input.tabslctd == 'gg'",
                                                                         uiOutput("box_1_3a")
                                                        ),
                                                        conditionalPanel(condition = "input.tabslctd == 'tbl'",
                                                                         wellPanel(
                                                                           #uiOutput("cluster_ids"),
                                                                           sliderInput(inputId = "no_of_de_genes", label = strong("Number of top DE genes:"), value = 10, min = 1, max = 100, step = 1),
                                                                           downloadButton("downloadData", "Download table of significantly DE genes")
                                                                         ),
                                                                         uiOutput("box_1_3b")
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
                                                        wellPanel(
                                                          sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F)
                                                        ),
                                                        wellPanel(style = "background:Teal",
                                                                  tags$hr(),
                                                                  tags$p(style = "font-family:Comic Sans MS;color:white",
                                                                         "Adjusting resolution (0.15, 0.25, 0.35, 0.45 or 0.55) to set the clustering granularity. This option allows the subdivision of clusters further into sub-populations and the subsequent interrogation of differential expression."
                                                                  ),
                                                                  tags$hr()
                                                        )
                                       ),
                                       
                                       #2.2 Option/subtitle for listing top cluster markers for which can then be used in the labelling of the generated populations
                                       conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between WT and KO'",
                                                        wellPanel(
                                                          uiOutput("dyn_clusters"),
                                                          sliderInput(inputId = "marker_genes_no", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
                                                          uiOutput("top_markers_umap")
                                                        ),
                                                        uiOutput("box_2_2") 
                                                        

                                       ),
                                       
                                       #2.3 Option/subtitle for allowing the labelling of the clusters
                                       conditionalPanel(condition = "input.graph_type == 'Label clusters'",
                                                        wellPanel(
                                                          #Select PC to plot
                                                          uiOutput("cluster_annot")
                                                          #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                                        ),
                                                        wellPanel(style = "background:Teal",
                                                                  tags$hr(),
                                                                  tags$p(style = "font-family:Comic Sans MS;color:white",
                                                                         "Labelling populations using the top cluster markers."
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
                                         withSpinner(plotOutput("labelled_umap")),
                                         tags$hr(),
                                         HTML(
                                           "<div class='text-center'><button class='btn-success btn-lg' onclick = 'fakeClick(\"Differential expression\")'> Explore differential expression in generated clusters</button></div>"
                                         ))
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
    selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_clusters()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
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
  

  
  output$top_markers_umap <- renderUI({
    selectInput(inputId = "select_markers_umap", label = strong("select marker to visualize in clusters:"), choices = cluster_markers()[,1], multiple = T, selected = head(cluster_markers()[,1], n = 4))
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
  
  ##information box
  output$box_2_2 =renderUI({
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     paste("Listing top cluster marker genes which can subsequently be used in labelling the cluster.")
              ),
              tags$hr()
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

    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster_names)){
      lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
        textInput(inputId = paste("labeller", x), label = strong(paste("Label cluster", x, "based on marker genes")), value = cluster_names[x+1])
      })
      
    } else {
      lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
        textInput(inputId = paste("labeller", x), label = strong(paste("Label cluster", x, "based on marker genes")), value = x)
      })
    }
  })
  
  ##Storing names decided on by the researchers based on optimal clustering to start off the differential expression visualization
  annotation = reactiveValues(annot = cluster_names)
  
  ##Observer to allow updating of cluster names dynamically as typed in
  observe({
    req(unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    })))
    annotation$annot = unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    }))
  })
  
  ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
  output$cluster_ids <- renderUI({
    umap_names = annotation$annot
    #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)

    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      umap_names = annotation$annot
      selectInput(inputId = "select_cell_type", label = strong("Select cell type to compare gene expression across conditions:"),choices = umap_names, multiple = F)
      
    } else {
      selectInput(inputId = "select_cell_type", label = strong("Select cell population to compare gene expression across conditions:"),choices = unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
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
  
  ##Information box
  output$box_1_1 <- renderUI({
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     #titlePanel("includeText"),
                     #fluidRow(
                     #column(12,
                     paste("Comparison of", input$de_genes, "expression between", conditions[1], "and", conditions[2], "using umap and violin plots.")
                     #))
              ),
              tags$hr()
    )  })
  
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
  
  ##Information box
  output$box_1_2 <- renderUI({
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     
                     paste("Comparison of gene expression between", conditions[1], "and", conditions[2], "across clusters using a dotplot. The genes are on the x-axis and the clusters on the y-axis. Red is for", conditions[1], "and Blue for", conditions[2], "cells with the increase in intensity of the respective colour (from grey to blue/red) correlating with the increase in the average level of gene expression across all cells in the cluster. The size of the dot corresponds to the percentage of cells in the cluster expressing the gene.")
                     
              ),
              tags$hr()
    )})
  ######################
  ###END Differential expression using dotplot under subtitle 1.2
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
  
  ##Retrieving table for DE scatterplotfrom precomputed list
  cell_type_de = reactive({
    umap_names = annotation$annot
    umap_cluster_modified1 = umap_cluster_modified_rna()
    
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
      tcells_combined_de_ggplots_table[[(input$clusters_res * 10)-0.5]][[as.numeric(match(input$select_cell_type,umap_names))]]
      
    } else {
      tcells_combined_de_ggplots_table[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$select_cell_type) + 1)]]
    }
    
  })
  
  ##preparing ggplot for average DE expression for genes above
  cell_type_de_plot = reactive({
    theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=9, fontface="italic")))
    ggplot(data=cell_type_de(), aes_string("KO", "WT")) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob)
    
    
  })
  
  ##plotting ggplot for average DE expression for genes above
  output$cell_type_plot = renderPlot({
    #theme_set(theme_cowplot())
    cell_type_de_plot()
  })
  
  ##Conditional tabsets for the DE table and ggplot UI change for input triggered by selecting tab 
  output$de_outputs <- renderUI({
    tabsetPanel(tabPanel("DE scatterplot", value = "gg", 
                         withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))), dataTableOutput("click_info")),
                tabPanel("DE table", value = "tbl",
                         withSpinner(dataTableOutput("top_de_genes"))), id = "tabslctd")
  })
  
  ##Displaying further details upon clicking points
  displayed_text <- reactive({
    req(input$plot_click)
    nearPoints(cell_type_de(), input$plot_click)
    
  })
  
  ##Displaying table with gene details upon click of point in DE scatterplot
  output$click_info <- renderDataTable({
    req(displayed_text())
    
    #cat("Name\n")
    displayed_text = displayed_text()
    displayed_text 
  }, escape = F)
  
  ##Information box
  output$box_1_3a <- renderUI({
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     paste("Comparison of gene expression between", conditions[1], "and", conditions[2], "in", input$select_cell_type, "using a scatter plot.")
              ),
              tags$hr()
    )
  })
  
  output$box_1_3b <- renderUI({
    
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     paste("Listing differentially expressed genes (adjusted P value <0.05) in" , conditions[1], "as compared to", conditions[2], "among", input$select_cell_type,".")
              ),
              tags$hr()
    )
  })
  ######################
  ### END Differential expression using ggplot and tables under subtitle 1.3
  ######################
  
  
}


shinyApp(ui = ui, server = server)
