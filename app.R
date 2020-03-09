#library(SeuratData)
library(shiny)
library(shinythemes)
library(plotly)
library(Seurat)
#install.packages("magrittr")
library(magrittr)
library(dplyr)
library(cowplot)
#install.packages("ggplot2")
library("ggplot2")
library(rlang)
library(tidyr)
library(tibble)
library(tidyverse)
library(ggrepel)
library(grid)


# adapt to your path
setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/Tcells_amp_adopted/")


#USE_CACHE <- TRUE

# Dependencies ----------------------------------------------------------------
#source("R/install-packages.R")
#source("R/processing_loading_data.R")



##visualization


ui = fluidPage(theme = shinytheme("lumen"),
               tabsetPanel(
                 source(file.path("R", "ui-tab-about.R"), local = TRUE)$value,
                 source(file.path("R", "ui-tab-ra.R"), local = TRUE)$value,
               ))
#                titlePanel("Effect of Beta 2 integrin in T Cell differentiation"),
#                sidebarLayout(
#                  sidebarPanel(
#                    radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("UMAP clusters", "Top cluster markers conserved between conditions", "Label populations", "Dotplot for gene expression between conditions and across cell types", "Differential gene expression between conditions per cell type", "Gene expression visualization between conditions across cell types")),
#                    
#                    # Display depending on which graph type is selected
#                    
#                    conditionalPanel(condition = "input.graph_type == 'UMAP clusters'",
#                                     # Select PC to plot
#                                     sliderInput(inputId = "umap_dim", label = strong("Number of PCs from 1"), value = dim1, min = dim1, max = dim2, step = diff_dim),
#                                     sliderInput(inputId = "neighbours_dim", label = strong("Number of PCs from 1 (KNN clustering in FindNeighbours function)"), value = dim1, min = dim1, max = dim2, step = diff_dim),
#                                     sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res),
#                                     actionButton(inputId = "go", label = "Run")
#                                     
#                    ),
#                    conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between conditions'",
#                                     
#                                     sliderInput(inputId = "marker_genes_no", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
#                                     #selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = , multiple = F),
#                                     uiOutput("dyn_clusters"),
#                                     actionButton(inputId = "go_marker", label = "Search"),
#                                     selectInput(inputId = "select_markers_umap", label = strong("Select conserved cluster markers:"), choices = all_genes_final, multiple = T)
#                                     
#                    ),
#                    conditionalPanel(condition = "input.graph_type == 'Label populations'",
#                                     #Select PC to plot
#                                     uiOutput("cluster_annot"),
#                                     actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
#                    ),
#                    conditionalPanel(condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'",
#                                     #Select PC to plot
#                                     selectInput(inputId = "select_markers_dotplot", label = strong("Select markers for dotplot:"),choices = all_genes_final, multiple = T)
#                    ),
#                    conditionalPanel(condition = "input.graph_type == 'Differential gene expression between conditions per cell type'",
#                                     #Select PC to plot
#                                     #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
#                                     uiOutput("cluster_ids"),
#                                     radioButtons(inputId = "condition1", label = "1st Condition", choiceNames = c("Wildtype 1", "Wildtype 2", "Knock-out 1"), choiceValues = (c("WT1","WT2", "KO1"))),
#                                     radioButtons(inputId = "condition2", label = "2nd Condition", choiceNames = c("Wildtype 1", "Wildtype 2", "Knock-out 1"), choiceValues = (c("WT1","WT2", "KO1"))),
#                                     uiOutput("ct_de_plot_dyn"),   
#                                     sliderInput(inputId = "no_of_de_genes", label = strong("Number of top DE genes:"), value = 10, min = 1, max = 100, step = 1),
#                                     uiOutput("ct_de_table_dyn")
#                                     
#                    ),
#                    conditionalPanel(condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'",
#                                     #Select PC to plot
#                                     selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_final, multiple = F),
#                                     #actionButton(inputId = "go_heatmap", label = "Run"),
#                                     #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
#                    )
#                  ),
#                  
#                  # Output: Description, lineplot, and reference
#                  mainPanel(  
#                    conditionalPanel(
#                      condition = "input.graph_type == 'UMAP clusters'", 
#                      plotOutput("all_groups"), plotOutput("stim_vs_ctrl")),
#                    conditionalPanel(
#                      condition = "input.graph_type == 'Top cluster markers conserved between conditions'", 
#                      tableOutput("top_conserved_genes"), plotOutput("conserved_markers_umap")), 
#                    conditionalPanel(
#                      condition = "input.graph_type == 'Label populations'", 
#                      plotOutput("labelled_umap")),
#                    conditionalPanel(
#                      condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'", 
#                      plotOutput("marker_dotplot")),
#                    conditionalPanel(
#                      condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
#                      plotOutput("cell_type_plot", click = clickOpts(id ="plot_click")), verbatimTextOutput("click_info")),
#                    conditionalPanel(
#                      condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
#                      tableOutput("top_de_genes")),
#                    conditionalPanel(
#                      condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'", uiOutput("tb") )
#                    
#                    
#                  )
#                )
# )


server = function(input, output) {
  
  umap_clusters = eventReactive(input$go, {
    tcells.combined.umap.list[[((input$umap_dim/5) - 1)]][[((input$umap_dim/5) - 1)]][[round(input$clusters_res/0.25)]]
  })
  
  output$all_groups = renderPlot({
    
    
    p1 = DimPlot(umap_clusters(), reduction = "umap", group.by = "sample")
    p2 = DimPlot(umap_clusters(), reduction = "umap", label = T)
    plot_grid(p1, p2)
    
  })
  
  ##Showing the stimulated and control umaps side by side
  output$stim_vs_ctrl = renderPlot({
    DimPlot(umap_clusters(), reduction = "umap", split.by = "sample")
  })
  
  
  output$dyn_clusters <- renderUI({
    umap_clusters = umap_clusters()
    selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_clusters$seurat_clusters), multiple = F)
  })
  
  cluster_markers = eventReactive(input$go_marker, {
    head(tcells.combined.clusters.tables[[((input$umap_dim/5) - 1)]][[((input$umap_dim/5) - 1)]][[round(input$clusters_res/0.25)]][[(as.numeric(input$marker_genes_cluster) + 1)]], n = input$marker_genes_no)
  })
  
  ##Finding conserved genes in clusters in both conditions to annotate cell types
  output$top_conserved_genes = renderTable({
    #isolate({
    cluster_markers() %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
  })
  #})
  #})
  
  umap_cluster_modified_rna = reactive({
    umap_cluster_modified_ul = umap_clusters()
    DefaultAssay(umap_cluster_modified_ul) = "RNA"
    umap_cluster_modified_ul
  })
  
  output$conserved_markers_umap = renderPlot({
    FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
  })
  
  output$cluster_annot <- renderUI({
    umap_clusters = umap_clusters()
    annotation = c("Il171", "Ccr7", "Lyc2", "Gzma", "Cdk6")
    lapply(0:(length(unique(umap_clusters$seurat_clusters))-1), function(x) {
      textInput(inputId = paste("labeller", x), label = strong(paste("Label the cluster", x, "based on marker genes")), value = annotation[x+1])
      
    })
  })
  
  # umap_cluster_modified_renamed = reactive({
  #   umap_cluster_modified = umap_cluster_modified_rna()
  #   #old_cluster_name = 0:(length(unique(umap_cluster_modified$seurat_clusters))-1)
  #   umap_cluster_modified_r = lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
  #     new_cluster_name = input[[paste("labeller",x)]]
  #     expr(`x`)
  #     print(new_cluster_name)
  #     RenameIdents(umap_cluster_modified, expr(`x`) = new_cluster_name)
  # 
  #   })
  # })
  
  umap_names = reactive({
    umap_cluster_modified = umap_cluster_modified_rna()
    new_cluster_name = unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    }))
    #print(new_cluster_name)
    
  })
  
  output$cluster_ids <- renderUI({
    umap_names = umap_names()
    selectInput(inputId = "select_cell_type", label = strong("Select cell type to compare gene expression across conditions:"),choices = umap_names, multiple = F)
    
  })
  
  umap_cluster_modified_renamed = reactive({
    umap_cluster_modified = umap_cluster_modified_rna()
    umap_names = umap_names()
    print(umap_names)
    # old_cluster_name = unlist(as.character(0:(length(unique(umap_cluster_modified$seurat_clusters))-1)))
    # new_cluster_name = unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
    #   new_cluster_name = input[[paste("labeller",x)]]
    #   }))
    # print(old_cluster_name)
    # print(new_cluster_name)
    
    #umap_cluster_modified <- RenameIdents(umap_cluster_modified, `0` = "Il171",`1` = "Ccr7", `2` = "Lyc2", `3` = "Gzma", `4`  = "Cdk6")
    names(umap_names) <- levels(umap_cluster_modified)
    umap_cluster_modified <- RenameIdents(umap_cluster_modified, umap_names)
    
    #umap_cluster_modified@active.ident <- plyr::mapvalues(umap_cluster_modified@active.ident, from = old_cluster_name, to = new_cluster_name)
    #umap_cluster_modified = RenameIdents(umap_cluster_modified, old.ident.na = x, new.ident.name = new_cluster_name)
    
    #})
  })
  
  umap_cluster_modified_umap = eventReactive(input$go_labelled_umap, {
    DimPlot(umap_cluster_modified_renamed(), label = TRUE)
  })
  output$labelled_umap = renderPlot({
    umap_cluster_modified_umap()
  })
  
  umap_cluster_modified_ren_reo = reactive({
    umap_cluster_modified_reo = umap_cluster_modified_renamed()
    umap_names = umap_names()
    Idents(umap_cluster_modified_reo) <- factor(Idents(umap_cluster_modified_reo), levels = umap_names)
    umap_cluster_modified_reo
    
  })
  
  output$marker_dotplot = renderPlot({
    DotPlot(umap_cluster_modified_ren_reo(), features = input$select_markers_dotplot, cols = c("blue", "red", "yellow"), dot.scale = 6, split.by = "sample") + RotatedAxis()
    
  })
  
  stim_markers = reactive({
    umap_cluster_modified = umap_cluster_modified_ren_reo()
    umap_cluster_modified$celltype.sample <- paste(Idents(umap_cluster_modified), umap_cluster_modified$sample, sep = "_")
    umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
    Idents(umap_cluster_modified) <- "celltype.sample"
    umap_cluster_modified
  })
  
  genes_in_de_order = reactive({
    umap_names = umap_names()
    FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, input$condition1, sep = "_"), ident.2 = paste(input$select_cell_type, input$condition2, sep = "_"), verbose = FALSE)
    
  })
  
  output$ct_de_plot_dyn <- renderUI({
    actionButton(inputId = "ct_de", label = paste("Plot differential expression (DE) of all genes in ", input$select_cell_type , sep = ""))
  })
  
  output$ct_de_table_dyn <- renderUI({
    actionButton(inputId = "ct_de_table", label = paste("Display table of top DE genes in ", input$select_cell_type , sep = ""))
  })
  
  top_de_g = eventReactive(input$ct_de_table,{
    
    ##Finding conserved genes in clusters in both conditions to annotate cell types
    head(genes_in_de_order(), n = input$no_of_de_genes) %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
  })
  
  output$top_de_genes = renderTable({
    top_de_g()
  })
  
  cell_type_de = eventReactive(input$ct_de,{
    
    cells_type <- subset(umap_cluster_modified_ren_reo(), idents = input$select_cell_type)
    Idents(cells_type) <- "sample"
    avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
    avg.cells$gene <- rownames(avg.cells)
    avg.cells
    
  })
  
  cell_type_de_plot = eventReactive(input$ct_de,{
    theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=9, fontface="italic")))
    ggplot(data=cell_type_de(), aes_string(input$condition1, input$condition2)) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob)
    
    
  })
  
  output$cell_type_plot = renderPlot({
    theme_set(theme_cowplot())
    cell_type_de_plot()
    #ggplotly(p)
  })
  
  
  # output$cell_type_plot = renderPlot({
  #   # p1 = ggplot(data=cell_type_de(), aes(KO1, WT1)) + geom_point() + ggtitle(input$select_cell_type)
  #   # p2 = ggplot(data=cell_type_de(), aes(KO1, WT2)) + geom_point() + ggtitle(input$select_cell_type)
  #   # p3 = ggplot(data=cell_type_de(), aes(WT1 , WT2)) + geom_point() + ggtitle(input$select_cell_type)
  #   # plot_grid(p1,p2,p3)
  #   
  #   ggplot(data=cell_type_de(), aes_string(input$condition1, input$condition2)) + geom_point() + ggtitle(input$select_cell_type)
  #   
  #   #ggplotly(p)
  # })
  
  
  displayed_text <- reactive({
    req(input$plot_click)
    nearPoints(cell_type_de(), input$plot_click)
    
  })
  
  output$click_info <- renderPrint({
    req(displayed_text())
    
    cat("Name\n")
    displayed_text()
  })
  
  output$de_stim_vs_ctrl_um = renderPlot({
    
    FeaturePlot(stim_markers(), features = input$de_genes, split.by = "sample", max.cutoff = 3,cols = c("grey", "red"))
  })
  
  output$de_stim_vs_ctrl_vp = renderPlot({
    
    VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE)  
    
  })
  
  output$tb <- renderUI({
    tabsetPanel(tabPanel("UMAP plot", 
                         plotOutput("de_stim_vs_ctrl_um")),
                tabPanel("Violin plot", 
                         plotOutput("de_stim_vs_ctrl_vp"))) 
  })
  
  
}


shinyApp(ui = ui, server = server)
