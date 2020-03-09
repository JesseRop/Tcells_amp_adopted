library(Seurat)
library(cowplot)
#library(SeuratData)
library(ggplot2)
library(shiny)
library(shinythemes)
library(plotly)
library(tidyr)
library(tibble)
library(tidyverse)
library(ggrepel)
library(grid)


#setwd("D:/GCRF_UoG/Vicky_JCR_Shiny")
all_genes_final_inf = readRDS("D:/GCRF_UoG/tuitorials/git_seurat/seurat_practice_files/all_immune_genes.rds")
immune.combined.umap.list_ifn = readRDS("D:/GCRF_UoG/tuitorials/git_seurat/seurat_practice_files/immune.combined.umap.list.rds")
immune.combined.clusters.tables_ifn = readRDS("D:/GCRF_UoG/tuitorials/git_seurat/seurat_practice_files/immune.combined.clusters.tables.rds")

dim1_ifn = 19
dim2_ifn = 20
diff_dim_ifn = 1


res1_ifn = 0.4
res2_ifn = 0.5 
diff_res_ifn = 0.1

if (!(all(exists("all_genes_final_inf"), exists("immune.combined.umap.list_ifn"), exists("immune.combined.clusters.tables_ifn")))) {
  data("ifnb") 
  
  #data("ifnb.SeuratData::ifnb")
  ifnb.list = SplitObject(ifnb.SeuratData::ifnb, split.by = "stim")
  
  ifnb.list = lapply(ifnb.list, function(x) { 
    x = NormalizeData(x)
    x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  } )
  
  immune.anchors = FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
  
  immune.combined = IntegrateData(anchorset = immune.anchors, dims = 1:20)
  
  #Running a single integrated analysis on all cells from the different conditions
  
  DefaultAssay(immune.combined) = "integrated"
  
  ##Running the standard workflow
  immune.combined = ScaleData(immune.combined, verbose = F)
  immune.combined = RunPCA(immune.combined, npcs = 30, verbose = F)
  
  all_genes_final_inf = row.names(immune.combined)
  saveRDS(all_genes_final_inf, "/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/all_immune_genes.rds")
  
  
  ##tSNE and clustering
  immune.combined.umap.list_ifn = lapply(dim1_ifn:dim2_ifn, function(x) {
    a = RunUMAP(immune.combined, reduction = "pca", dims = 1:x)
    lapply(dim1_ifn:dim2_ifn, function(x) {
      b = FindNeighbors(a, reduction = "pca", dims = 1:x)
      lapply(seq(res1_ifn, res2_ifn, by = 0.1), function(x) FindClusters(b, resolution = x))
    })
  })
  
  saveRDS(immune.combined.umap.list_ifn, "/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/immune.combined.umap.list.rds")
  
  immune.combined.clusters.tables_ifn = lapply(immune.combined.umap.list_ifn, function(x) { 
    lapply(x, function(x) {
      lapply(x, function(x){
        DefaultAssay(x) = "RNA"
        #cluster.markers = lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
        cluster.markers = lapply(6:7, function(y) {
          FindConservedMarkers(x, ident.1 = y, grouping.var = "stim")
          
        })
      }) 
    }) 
  })
  
  saveRDS(immune.combined.clusters.tables_ifn, "/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/immune.combined.clusters.tables.rds")
  
}

#cell_types = c( "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T")

##visualization



ui = fluidPage(theme = shinytheme("lumen"),
               titlePanel("Interferron beta stimulation of PBMCs"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("UMAP clusters", "Top cluster markers conserved between conditions", "Label populations", "Dotplot for gene expression between conditions and across cell types", "Differential gene expression between conditions per cell type", "Gene expression visualization between conditions across cell types")),
                   # Display depending on which graph type is selected
                   
                   conditionalPanel(condition = "input.graph_type == 'UMAP clusters'",
                                    # Select PC to plot
                                    sliderInput(inputId = "umap_dim_ifn", label = strong("Number of PCs from 1"), value = dim1_ifn, min = dim1_ifn, max = dim2_ifn, step = diff_dim_ifn),
                                    sliderInput(inputId = "neighbours_dim_ifn", label = strong("Number of PCs from 1 (KNN clustering in FindNeighbours function)"), value = dim1_ifn, min = dim1_ifn, max = dim2_ifn, step = diff_dim_ifn),
                                    sliderInput(inputId = "clusters_res_ifn", label = strong("Louvain algorithm resolution"), value = res1_ifn, min = res1_ifn, max = res2_ifn, step = diff_res_ifn),
                                    actionButton(inputId = "go_ifn", label = "Run")
                                    
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between conditions'",
                                    
                                    sliderInput(inputId = "marker_genes_no_ifn", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
                                    #selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = , multiple = F),
                                    uiOutput("dyn_clusters_ifn"),
                                    actionButton(inputId = "go_ifn_marker", label = "Search"),
                                    selectInput(inputId = "select_markers_umap_ifn", label = strong("Select conserved cluster markers:"), choices = all_genes_final_inf, multiple = T)
                                    
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Label populations'",
                                    #Select PC to plot
                                    uiOutput("cluster_annot_ifn"),
                                    actionButton(inputId = "go_ifn_labelled_umap", label = "View labelled clusters"),
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'",
                                    #Select PC to plot
                                    selectInput(inputId = "select_markers_dotplot_ifn", label = strong("Select markers for dotplot:"),choices = all_genes_final_inf, multiple = T)
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Differential gene expression between conditions per cell type'",
                                    #Select PC to plot
                                    #actionButton(inputId = "go_ifn_labelled_umap", label = "View labelled clusters"),
                                    uiOutput("cluster_ids_ifn"),
                                    uiOutput("ct_de_plot_dyn_ifn"),   
                                    sliderInput(inputId = "no_of_de_genes_ifn", label = strong("Number of top differentially expressed (DE) genes:"), value = 10, min = 1, max = 100, step = 1),
                                    uiOutput("ct_de_table_dyn_ifn"),   
                                    
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'",
                                    #Select PC to plot
                                    selectInput(inputId = "de_genes_ifn", label = strong("Choose gene:"),choices = all_genes_final_inf, multiple = F),
                                    #actionButton(inputId = "go_ifn_heatmap", label = "Run"),
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
                     plotOutput("cell_type_plot", hover = hoverOpts(id ="plot_hover")), verbatimTextOutput("hover_info")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
                     tableOutput("top_de_genes")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'", uiOutput("tb") )
                   
                   
                 )
               )
)

server = function(input, output) {
  
  umap_clusters = eventReactive(input$go_ifn, {
    immune.combined.umap.list_ifn[[(input$umap_dim_ifn - 18)]][[(input$neighbours_dim_ifn - 18)]][[((input$clusters_res_ifn * 10)-3)]]
  })
  
  output$all_groups = renderPlot({
    
    p1 = DimPlot(umap_clusters(), reduction = "umap", group.by = "stim")
    p2 = DimPlot(umap_clusters(), reduction = "umap", label = T)
    plot_grid(p1, p2)
    
  })
  
  ##Showing the stimulated and control umaps side by side
  output$stim_vs_ctrl = renderPlot({
    DimPlot(umap_clusters(), reduction = "umap", split.by = "stim")
  })
  
  output$dyn_clusters_ifn <- renderUI({
    
    umap_clusters = umap_clusters()
    selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for (cluster 6 & 7 only)"), choices = unique(umap_clusters$seurat_clusters), multiple = F)
  })
  
  cluster_markers = eventReactive(input$go_ifn_marker, {
    head(immune.combined.clusters.tables_ifn[[(input$umap_dim_ifn - 18)]][[(input$neighbours_dim_ifn - 18)]][[((input$clusters_res_ifn * 10)-3)]][[ifelse(as.numeric(input$marker_genes_cluster) - 5 < 0, 0, as.numeric(input$marker_genes_cluster) - 5)]], n = input$marker_genes_no_ifn)
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
    FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap_ifn, min.cutoff = "q9")
  })
  
  
  output$cluster_annot_ifn <- renderUI({
    umap_clusters = umap_clusters()
    annotation = c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", "CD16 Mono",  "B", "CD8 T", "T activated", "NK", "DC", "B Activated",  "Mk",  "pDC",  "Eryth", "Mono/Mk Doublets")
    lapply(0:(length(unique(umap_clusters$seurat_clusters))-1), function(x) {
      textInput(inputId = paste("labeller", x), label = strong(paste("Label cluster", x, "based on marker genes")), value = annotation[x+1])
      
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
  
  output$cluster_ids_ifn <- renderUI({
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
  
  umap_cluster_modified_umap = eventReactive(input$go_ifn_labelled_umap, {
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
    DotPlot(umap_cluster_modified_ren_reo(), features = input$select_markers_dotplot_ifn, cols = c("blue", "red"), dot.scale = 6, split.by = "stim") + RotatedAxis()
    
  })
  
  stim_markers = reactive({
    umap_cluster_modified = umap_cluster_modified_ren_reo()
    umap_cluster_modified$celltype.stim <- paste(Idents(umap_cluster_modified), umap_cluster_modified$stim, sep = "_")
    umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
    Idents(umap_cluster_modified) <- "celltype.stim"
    umap_cluster_modified
  })
  
  genes_in_de_order = reactive({
    umap_names = umap_names()
    FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "_STIM", sep = ""), ident.2 = paste(input$select_cell_type, "_CTRL", sep = ""), verbose = FALSE)
    
  })
  
  output$ct_de_plot_dyn_ifn <- renderUI({
    actionButton(inputId = "ct_de", label = paste("Plot differential expression (DE) of all genes in ", input$select_cell_type , sep = ""))
  })
  
  output$ct_de_table_dyn_ifn <- renderUI({
    actionButton(inputId = "ct_de_table", label = paste("Display table of top DE genes in ", input$select_cell_type , sep = ""))
  })
  
  top_de_g = eventReactive(input$ct_de_table,{
    
    ##Finding conserved genes in clusters in both conditions to annotate cell types
    head(genes_in_de_order(), n = input$no_of_de_genes_ifn) %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
  })
  
  output$top_de_genes = renderTable({
    top_de_g()
  })
  
  
  cell_type_de = eventReactive(input$ct_de,{
    cells_type <- subset(umap_cluster_modified_ren_reo(), idents = input$select_cell_type)
    Idents(cells_type) <- "stim"
    avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
    avg.cells$gene <- rownames(avg.cells)
    avg.cells
    
    
  })
  
  
  cell_type_de_plot = eventReactive(input$ct_de,{
    theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Place pointer on dots to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=9, fontface="italic")))
    ggplot(data=cell_type_de(), aes(CTRL , STIM)) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob)
    
    
  })
  output$cell_type_plot = renderPlot({
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
    req(input$plot_hover)
    nearPoints(cell_type_de(), input$plot_hover)
    
  })
  
  output$hover_info <- renderPrint({
    req(displayed_text())
    
    cat("Name\n")
    displayed_text()
  })
  
  output$de_stim_vs_ctrl_um = renderPlot({
    
    FeaturePlot(stim_markers(), features = input$de_genes_ifn, split.by = "stim", max.cutoff = 3,cols = c("grey", "red"))
  })
  
  output$de_stim_vs_ctrl_vp = renderPlot({
    
    VlnPlot(stim_markers(), features = input$de_genes_ifn, split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)  
    
  })
  
  
  output$tb <- renderUI({
    tabsetPanel(tabPanel("UMAP plot", 
                         plotOutput("de_stim_vs_ctrl_um")),
                tabPanel("Violin plot", 
                         plotOutput("de_stim_vs_ctrl_vp"))) 
  })
  
  
}


shinyApp(ui = ui, server = server)



