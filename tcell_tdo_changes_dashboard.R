
setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/Tcells_amp_adopted")
source("tcell_tdo_changes_libs_raw_objects.R")
##Creating the ui object for the user interface

library(shinydashboard)


which_numeric_cols = function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}

ui <- dashboardPage( skin = "black",
  dashboardHeader(title = "Role of CD18 in γδ T cells"),
  dashboardSidebar(
    sidebarMenu(
      id = 'sidebar',
      menuItem("Differential expression", tabName = "vis", icon = icon("dashboard"),
               menuSubItem("Gene view", tabName = "ge_vis_gene", icon = icon("dashboard")),
               menuSubItem("Cell view", tabName = "ge_vis_cell", icon = icon("dashboard"))
      ),
      
      menuItem("Cluster adjustment", tabName = "cluster_res", icon = icon("dashboard"))
     
    )
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    tabItems(
      # First tab content
      tabItem(tabName = "ge_vis_gene",
              #fluidRow(
              wellPanel(
                fluidRow(
                  #h3("Gene expression"),
                  column(6,
                         wellPanel(
                           fluidRow(
                             column(
                               12,
                               h4("Single gene visualization"),
                               selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_common_in_all_groups, multiple = F, selected = "Cd163l1"),
                             )
                           ),
                           fluidRow(
                             tabBox(
                               
                               #title = "Single gene differential expression",
                               # The id lets us use input$tabset1 on the server to find the current tab
                               id = "vis_de", 
                               #height = "250px",
                               width = 12,
                               tabPanel("Violin plot", withSpinner(plotOutput("de_stim_vs_ctrl_vp"))),
                               tabPanel("UMAP feature plot", withSpinner(plotOutput("de_stim_vs_ctrl_um")))
                             ),
                             
                             
                           ),
                           tags$hr(),
                           uiOutput("box_1_1")
                           
                         )
                  ),
                  
                  column(
                    6,
                    wellPanel(
                      h4("Dotplot for multiple gene visualization"),
                      
                      selectInput(inputId = "select_markers_dotplot", label = strong("Choose genes:"),choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes),
                      withSpinner(plotOutput("marker_dotplot")),
                      tags$hr(),
                      uiOutput("box_1_2")
                      
                      # )
                    )
                  )
                )
              )
      ),
      tabItem(tabName = "ge_vis_cell",
              
              wellPanel(
                fluidRow(
                  column(
                    4,
                    uiOutput("cluster_ids"),
                  )
                ),
                fluidRow(
                  
                  column(6,
                         wellPanel(
                           h4("DE scatterplot"),
                           fluidRow(
                             withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))), dataTableOutput("click_info"),
                             
                             
                           ),
                           tags$hr(),
                           uiOutput("box_1_3a")
                           
                         )
                  ),
                  
                  column(
                    6,
                    wellPanel(
                      h4("DE table"),
                      
                      # sliderInput(inputId = "no_of_de_genes", label = strong("Number of top DE genes:"), value = 10, min = 1, max = 100, step = 1),
                      withSpinner(dataTableOutput("top_de_genes")),
                      tags$hr(),
                      uiOutput("box_1_3b")
                      
                      # )
                    )
                  )
                )
              )
              
              
              
              
              
              #)
      ),
      
      tabItem(tabName = "cluster_res",
              
              wellPanel(fluidRow(
                
                column(
                  7,
                  tabBox(
                    title = "UMAP Clustering",
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "vis_de", 
                    #height = "250px",
                    width = 12,
                    tabPanel("Labelling", 
                             wellPanel(
                               
                               sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F),
                               withSpinner(plotOutput("labelled_umap"))
                               
                               # )
                               ,
                               
                               uiOutput("cluster_annot"),
                               # )
                               tags$hr(),
                               wellPanel(style = "background:Teal",
                                         tags$hr(),
                                         tags$p(style = "font-family:Comic Sans MS;color:white",
                                                paste("Adjusting resolution (0.15, 0.25, 0.35, 0.45 or 0.55) to set the clustering granularity. This option allows the sub-division of clusters further into sub-populations and the subsequent interrogation of differential expression. Labelling of cell populations can then be done based on top cluster markers in the 'Cluster marker' table on the left.")
                                                
                                         ),
                                         tags$hr()
                               )
                             )),
                    tabPanel("KO VS WT1 VS WT2", withSpinner(plotOutput("all_groups")))

                  ),
                  
                ),
                
                column(
                  5, 
                  tabBox(
                    title = "Cluster markers",
                    # The id lets us use input$tabset1 on the server to find the current tab
                    id = "vis_de", 
                    #height = "250px",
                    width = 12,
                    tabPanel("Marker table", 
                             wellPanel(
                               # box(
                               h4("Cluster markers"),
                               # width = NULL,
                               # solidHeader = TRUE,
                               uiOutput("dyn_clusters"),
                               withSpinner(DTOutput("top_conserved_genes")),
                               tags$hr(),
                               uiOutput("box_2_2")
                               )
                               
                               # )
                             ),
                    tabPanel("Marker feature plots", 
                             uiOutput("top_markers_umap"),
                             withSpinner(plotOutput("conserved_markers_umap")))
                  )
                  
                  
                )
              )
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
    
    DimPlot(umap_cluster_modified_ren_reo(), reduction = "umap", split.by = "sample", label.size = 6)

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
    tcells_combined_clusters_tables_res[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$marker_genes_cluster) + 1)]] %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
    
  })
  
  
  output$top_conserved_genes = DT::renderDataTable({
    numeric_cols =
      colnames(data.frame(cluster_markers()))[which_numeric_cols(data.frame(cluster_markers()))]
    
    #   # Javascript-enabled table.
    #   datatable(
    DT::datatable(
      cluster_markers(),
      # colnames = c(
      #   "Gene",
      #   "scRNA-seq cluster",
      #   "Avg_logFC",
      #   "% in interest cluster",
      #   "% in other clusters",
      #   "Adjusted P"
      # ),
      selection = "single",
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller"),
      options = list(
        deferRender = TRUE,
        scrollY = 350,
        #scroller = TRUE,
        #lengthMenu = FALSE,
        autoWidth = FALSE,
        dom = "Blfrtip",
        buttons = 
          list(list(
            extend = "collection",
            buttons = c("csv", "pdf"),
            text = "Download"
          )),  # end of buttons customization
        
        # customize the length menu
        lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
        pageLength = 10
      ), fillContainer = TRUE
    ) %>%
      DT::formatSignif(columns = numeric_cols, digits = 3)
  }, server = TRUE)
  
  
  ##Preparing and plotting UMAP cluster markers for annotating cell types
  umap_cluster_modified_rna = reactive({
    umap_cluster_modified_ul = umap_clusters()
    DefaultAssay(umap_cluster_modified_ul) = "RNA"
    umap_cluster_modified_ul
  })
  
  
  output$top_markers_umap <- renderUI({
    a = cluster_markers()[,1]
    # a = factor(a, levels = c(a))
    selectInput(inputId = "select_markers_umap", label = strong("Select marker to visualize in clusters:"), choices = a, multiple = T, selected = head(a, n=4))
    
  })
  output$conserved_markers_umap = renderPlot({
    FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
  })
  
 
  
  ##information box
  output$box_2_2 =renderUI({
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     paste("Listing top cluster marker genes irrespective of cell state which can subsequently be used in labelling the cluster. Here, markers are genes highly expressed in a cluster as compared to all other clusters in both", conditions[1], "and", conditions[2],".")
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
      umap_clusters = umap_clusters()
      do.call(flowLayout, 
              lapply(0:(length(unique(umap_clusters$seurat_clusters))-1), function(x) {
                textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = cluster_names[x+1])
              })
      )
      
    } else {
      umap_clusters = umap_clusters()
      do.call(flowLayout,
              lapply(0:(length(unique(umap_clusters$seurat_clusters))-1), function(x) {
                textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = x)
              })
      )
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
      selectInput(inputId = "select_cell_type", label = strong(paste("Select cell population to compare gene expression between",conditions[1],"and", conditions[2],":")),choices = umap_names, multiple = F)
      
    } else {
      selectInput(inputId = "select_cell_type", label = strong("Select cell population to compare gene expression across conditions:"),choices = unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
    }
    
  })
  
  ##Renaming clusters
  umap_cluster_modified_ren_reo = reactive({
    umap_names = annotation$annot
    umap_cluster_modified1 = umap_cluster_modified_rna()
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      names(umap_names) <- levels(umap_cluster_modified1)
      umap_cluster_modified1 <- RenameIdents(umap_cluster_modified1, umap_names)
      
    } else {
      umap_cluster_modified1
    }
    
  })
  
  # umap_cluster_modified_umap = reactive({
  #   DimPlot(umap_cluster_modified_ren_reo(), label = TRUE, label.size = 6)
  #   
  # })
  output$labelled_umap = renderPlot({
    DimPlot(umap_cluster_modified_ren_reo(), label = TRUE, label.size = 6)
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
  
  
  ##Information box
  output$box_1_1 <- renderUI({
    wellPanel(style = "background:Teal",
              tags$hr(),
              tags$p(style = "font-family:Comic Sans MS;color:white",
                     #titlePanel("includeText"),
                     #fluidRow(
                     #column(12,
                     paste("Comparison of", input$de_genes, "expression between", conditions[1], "and", conditions[2], "using violin plots and umap.")
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
                     #paste("Comparison of", input$select_markers_dotplot[1], input$select_markers_dotplot[2], "expression between", conditions[1], "and", conditions[2], "across clusters using a dotplot. Red for Wildtype and Blue for CD18 KO cells with the increasing in intensity of the respective colour correlating with the level of gene expression in the cluster"),
                     # paste(unlist(for (i in 1:length(input$select_markers_dotplot)) {
                     #   print(paste(input$select_markers_dotplot[i]))
                     # }
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

    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      tcells_combined_de_tables[[(input$clusters_res * 10)-0.5]][[as.numeric(match(input$select_cell_type,umap_names))]]
      
    } else {
      tcells_combined_de_tables[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$select_cell_type) + 1)]]
    }
    
    
  })
  
  ##Retrieving table for DE expression from precomputed list
  top_de_g = reactive({
    genes_in_de_order()  %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
  })
  
  
  output$top_de_genes = DT::renderDataTable({
    numeric_cols =
      colnames(data.frame(top_de_g()))[which_numeric_cols(data.frame(top_de_g()))]
    
    #   # Javascript-enabled table.
    #   datatable(
    DT::datatable(
      top_de_g(),
      # colnames = c(
      #   "Gene",
      #   "scRNA-seq cluster",
      #   "Avg_logFC",
      #   "% in interest cluster",
      #   "% in other clusters",
      #   "Adjusted P"
      # ),
      selection = "single",
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller"),
      options = list(
        deferRender = TRUE,
        scrollY = 350,
        #scroller = TRUE,
        #lengthMenu = FALSE,
        autoWidth = FALSE,
        dom = "Blfrtip",
        buttons = 
          list(list(
            extend = "collection",
            buttons = c("csv", "pdf"),
            text = "Download"
          )),  # end of buttons customization
        
        # customize the length menu
        lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
        pageLength = 10
      )
    ) %>%
      DT::formatSignif(columns = numeric_cols, digits = 3)
  }, server = TRUE)
  
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

    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      tcells_combined_de_ggplots_table[[(input$clusters_res * 10)-0.5]][[as.numeric(match(input$select_cell_type,umap_names))]]
      
    } else {
      tcells_combined_de_ggplots_table[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$select_cell_type) + 1)]]
    }
    
  })
  
  ##preparing ggplot for average DE expression for genes above
  cell_type_de_plot = reactive({
    #theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=9, fontface="italic")))
    ggplot(data=cell_type_de(), aes_string("KO", "WT")) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob) + theme_bw()
    
    
  })
  
  ##plotting ggplot for average DE expression for genes above
  output$cell_type_plot = renderPlot({
    #theme_set(theme_cowplot())
    cell_type_de_plot()
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
