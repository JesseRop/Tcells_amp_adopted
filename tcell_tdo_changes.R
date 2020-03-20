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
library(data.table)
library(shinycssloaders)

# adapt to your path
setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/")

# tcells.combined.umap.list<-readRDS("tcells.combined.umap.list.rds")
# tcells.combined.clusters.tables<-readRDS("tcells.combined.clusters.tables.rds")
# all_genes_final = row.names(tcells.combined.umap.list[[1]][[1]][[1]])

tcells_combined_umap_list_res<-readRDS("tcells_combined_umap_list_res.rds")
# tcells_combined_clusters_tables_res1<-readRDS("tcells_combined_clusters_tables_res.rds")
# saveRDS(tcells_combined_clusters_tables_res1, "tcells_combined_clusters_tables_res1.rds")
# tcells_combined_clusters_tables_res2_3<-readRDS("tcells_combined_clusters_tables_res2_3.rds")
# tcells_combined_clusters_tables_res4_5<-readRDS("tcells_combined_clusters_tables_res4_5.rds")
# tcells_combined_clusters_tables_res = c(tcells_combined_clusters_tables_res1, tcells_combined_clusters_tables_res2_3, tcells_combined_clusters_tables_res4_5)
# saveRDS(tcells_combined_clusters_tables_res, "tcells_combined_clusters_tables_res.rds")
tcells_combined_clusters_tables_res<-readRDS("tcells_combined_clusters_tables_res.rds")



#all_genes_final = row.names(tcells_combined_umap_list_res[1])
all_genes_common_in_all_groups = readRDS("all_genes_common_in_all_groups.rds")

# uniprot_info_raw = fread("/GCRF_UoG/Vicky_JCR_Shiny/uniprot table/unipro-mouseID")
# uniprot_info_raw$uniprot = paste('<a href="https://www.uniprot.org/uniprot/',uniprot_info_raw$Entry,'" target="_blank">', uniprot_info_raw$Entry,'</a>', sep = "")
# write.table(uniprot_info_raw, "/GCRF_UoG/Vicky_JCR_Shiny/uniprot_info_with_link", row.names = F, sep = "\t", quote = F)
uniprot_info = fread("/GCRF_UoG/Vicky_JCR_Shiny/uniprot_info_with_link", stringsAsFactors = F)


dim=15
# dim1 = 10
# dim2 = 15
# diff_dim = 5

res1 = 0.15
res2 = 0.55
diff_res = 0.1
fav_genes = c("Cdk6", "Gzma", "Ly6c2", "Ccr7", "Il7a")

#if (!(all(exists("tcells.combined.umap.list"), exists("tcells.combined.clusters.tables")))) {
if (!(all(exists("tcells_combined_umap_list_res"), exists("tcells_combined_clusters_tables_res")))) {
  #Merge first the two WT 1 and 3
  WT1.data <- Read10X(data.dir = "WT1")
  WT1 <- CreateSeuratObject(counts = WT1.data, project = "WT1", min.cells = 3)
  WT1$sample <- "WT1"
  WT1$group <- "WT"
  WT1[["percent.mt"]] <- PercentageFeatureSet(object = WT1, pattern = "^mt-")
  plot1 <- FeatureScatter(object = WT1, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = WT1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  WT1 <- subset(WT1, subset = nFeature_RNA > 500 & nFeature_RNA < 2200 & percent.mt < 18)
  
  # store mitochondrial percentage in object meta data
  WT1 <- PercentageFeatureSet(WT1, pattern = "^mt-", col.name = "percent.mt")
  
  # run sctransform
  WT1 <- SCTransform(WT1, vars.to.regress = "percent.mt", verbose = FALSE)
  
  ##capturing all genes at this stage
  all_genes_wt1 = rownames(WT1)
  saveRDS(all_genes_wt1, "all_genes_wt1.rds")
  ##capturing all genes at this stage
  
  ### load the KO data
  KO.data <- Read10X(data.dir = "KO1")
  KO <- CreateSeuratObject(counts = KO.data, project = "KO1", min.cells = 3)
  KO$sample <- "KO1"
  KO$group <- "KO"
  KO[["percent.mt"]] <- PercentageFeatureSet(object = KO, pattern = "^mt-")
  plot1 <- FeatureScatter(object = KO, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = KO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  # store mitochondrial percentage in object meta data
  KO <- PercentageFeatureSet(KO, pattern = "^mt-", col.name = "percent.mt")
  
  # run sctransform
  KO <- SCTransform(KO, vars.to.regress = "percent.mt", verbose = FALSE)
  
  ##capturing all genes at this stage
  all_genes_ko = rownames(KO)
  saveRDS(all_genes_ko, "all_genes_ko.rds")
  ##capturing all genes at this stage
  
  WT3.data <- Read10X(data.dir = "WT3")
  WT3 <- CreateSeuratObject(counts = WT3.data, project = "WT3", min.cells = 3)
  WT3$sample <- "WT2"
  WT3$group <- "WT"
  WT3[["percent.mt"]] <- PercentageFeatureSet(object = WT3, pattern = "^mt-")
  plot1 <- FeatureScatter(object = WT3, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = WT3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  WT3 <- subset(WT3, subset = nFeature_RNA > 500 & nFeature_RNA < 3200 & percent.mt < 18)
  
  
  # store mitochondrial percentage in object meta data
  WT3 <- PercentageFeatureSet(WT3, pattern = "^mt-", col.name = "percent.mt")
  
  # run sctransform
  WT3 <- SCTransform(WT3, vars.to.regress = "percent.mt", verbose = FALSE)
  
  ##capturing all genes at this stage
  all_genes_wt3 = rownames(WT3)
  saveRDS(all_genes_wt3, "all_genes_wt3.rds")
  ##capturing all genes at this stage
  
  WT3.downsample = subset(WT3, cells = sample(Cells(WT3), 4000))
  
  ##capturing all genes at this stage
  all_genes_wt3_ds = rownames(WT3)
  saveRDS(all_genes_wt3_ds, "all_genes_wt3_ds.rds")
  ##capturing all genes at this stage
  
  
  
  ## estimated in our script the best value here
  #WT.combined <- JackStraw(object = WT.combined, num.replicate = 100, dims=30)
  #WT.combined <- ScoreJackStraw(object = WT.combined, dims = 1:30)
  #JackStrawPlot(object = WT.combined, dims = 1:30)
  #ElbowPlot(object = WT.combined, ndims = 40)
  # merge the two WT's
  #dim=15
  
  TC.anchors <- FindIntegrationAnchors(object.list = list(WT1,KO), dims = 1:dim)
  Combined <- IntegrateData(anchorset = TC.anchors, dims = 1:dim)
  DefaultAssay(Combined) <- "integrated"
  Combined <- ScaleData(Combined, verbose = FALSE)
  Combined <- RunPCA(Combined, npcs = dim, verbose = FALSE)
  Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
  Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
  Combined <- FindClusters(Combined, resolution = 0.2)
  p1 <- DimPlot(Combined, reduction = "umap", group.by = "sample")
  p2 <- DimPlot(Combined, reduction = "umap", label = TRUE)
  
  #x11()
  plot_grid(p1, p2)
  Combined[["UMI"]] <-  Combined$nCount_RNA  # Why divided by 100
  Combined[["genes"]] <-  Combined$nFeature_RNA
  FeaturePlot(Combined, features = "UMI")
  
  ##capturing all genes at this stage
  all_genes_combined = rownames(Combined)
  saveRDS(all_genes_combined, "all_genes_combined.rds")
  ##capturing all genes at this stage
  
  
  dim=15
  WT.anchors <- FindIntegrationAnchors(object.list = list(Combined, WT3.downsample), dims = 1:dim)
  Three.combined <- IntegrateData(anchorset = WT.anchors, dims = 1:dim)
  DefaultAssay(Three.combined) <- "integrated"
  Three.combined <- ScaleData(Three.combined, verbose = FALSE)
  Three.combined <- RunPCA(Three.combined, npcs = dim, verbose = FALSE)
  Three.combined <- RunUMAP(Three.combined, reduction = "pca", dims = 1:dim)
  Three.combined <- FindNeighbors(Three.combined, reduction = "pca", dims = 1:dim)
  Three.combined <- FindClusters(Three.combined, resolution = 0.1)
  p1 <- DimPlot(Three.combined, reduction = "umap", group.by = "sample")
  p2 <- DimPlot(Three.combined, reduction = "umap", label = TRUE, label.size = 5)
  plot_grid(p1, p2)
  
  # filtering away other clusters
  Combined.filt <- subset(Three.combined, idents = c("0","1","2"), invert = FALSE)
  DimPlot(Combined.filt, reduction = "umap", split.by = "sample")
  
  # re-clusterting, without the other cell types
  DefaultAssay(object = Combined.filt) <- "integrated"
  # Run the standard workflow for visualization and clustering
  Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
  Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
  
  # t-SNE and Clustering
  Combined.filt <- RunUMAP(object = Combined.filt, reduction = "pca", dims = 1:dim)
  Combined.filt <- FindNeighbors(object = Combined.filt, reduction = "pca", dims = 1:dim)
  #Combined.filt <- FindClusters(Combined.filt, resolution = 0.15)
  
  
  ### decrease the amount of clusters, and split the other one (green, killers)
  # filtering away other clusters
  # Combined.filt2 <- subset(Combined.filt, idents = c("0","1","2","3","5"), invert = FALSE)
  # 
  # Combined.filt<-Combined.filt2
  # # re-clusterting, without the other cell types
  # DefaultAssay(object = Combined.filt) <- "integrated"
  # # Run the standard workflow for visualization and clustering
  # Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
  # Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
  # 
  all_genes_common_in_all_groups = Reduce(intersect,list(all_genes_ko,all_genes_wt1,all_genes_wt3))
  saveRDS(all_genes_common_in_all_groups, "all_genes_common_in_all_groups.rds")
  
  
  ##tSNE and clustering allowing change of PCs
  # tcells.combined.umap.list = lapply(seq(dim1, dim2, by = diff_dim), function(x) {
  #   a = RunUMAP(Combined.filt, reduction = "pca", dims = 1:x)
  #   lapply(seq(dim1, dim2, by = diff_dim), function(x) {
  #     b = FindNeighbors(a, reduction = "pca", dims = 1:x)
  #     lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(b, resolution = x))
  #   })
  # })
  
  tcells_combined_umap_list_res = lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(Combined.filt, resolution = x))
  
  saveRDS(tcells_combined_umap_list_res, "tcells_combined_umap_list_res.rds")
  
  # tcells.combined.clusters.tables = lapply(tcells.combined.umap.list, function(x) { 
  #   lapply(x, function(x) {
  #     lapply(x, function(x){
  #       DefaultAssay(x) = "RNA"
  #       #cluster.markers = lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
  #       cluster.markers = lapply(0:1, function(y) {
  #         FindConservedMarkers(x, ident.1 = y, grouping.var = "sample")
  #         
  #       })
  #     }) 
  #   }) 
  # })
  # 
  # saveRDS(tcells.combined.clusters.tables, "tcells.combined.clusters.tables.rds")
  
  tcells_combined_clusters_tables_res = lapply(tcells_combined_umap_list_res[1], function(x) { 
    DefaultAssay(x) = "RNA"
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      FindConservedMarkers(x, ident.1 = y, grouping.var = "group")
      
    })
  })
  
  saveRDS(tcells_combined_clusters_tables_res, "tcells_combined_clusters_tables_res.rds")
  
  
  
}

#cell_types = c( "Il171", "Ccr7", "Lyc2", "Gzma", "Cdk6")

#######################################IFN

#######################################IFN


ui = fluidPage(theme = shinytheme("lumen"),
               navbarPage("Role of CD18 in γδ T cells",
                          tabPanel("Differential expression",
                            #title = "Effect of Beta 2 integrin in T Cell differentiation",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons(inputId = "de_panel", label = strong("Differential expression between CD18 Knockouts (KO) and Wildtypes (WT1 / WT2)"), choices = c("Gene expression visualization between conditions across cell types", "Dotplot for gene expression between conditions and across cell types", "Differential gene expression between conditions per cell type" )),
                                
                                # Display depending on which graph type is selected
                                
                                conditionalPanel(condition = "input.de_panel == 'Gene expression visualization between conditions across cell types'",
                                                 #Select PC to plot
                                                 selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_common_in_all_groups, multiple = F),
                                                 #actionButton(inputId = "go_heatmap", label = "Run"),
                                                 #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                                ),
                                conditionalPanel(condition = "input.de_panel == 'Dotplot for gene expression between conditions and across cell types'",
                                                 #Select PC to plot
                                                 selectInput(inputId = "select_markers_dotplot", label = strong("Select markers for dotplot:"),choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes)
                                ),
                                conditionalPanel(condition = "input.de_panel == 'Differential gene expression between conditions per cell type'",
                                                 #Select PC to plot
                                                 #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                                 uiOutput("cluster_ids"),
                                                 #radioButtons(inputId = "condition2", label = "WT group to compare to KO", choiceNames = c("Wildtype 1", "Wildtype 2"), choiceValues = (c("WT1","WT2"))),
                                                 conditionalPanel(condition = "input.tabslctd == 'gg'",
                                                                  uiOutput("ct_de_plot_dyn")
                                                                  ),
                                                 conditionalPanel(condition = "input.tabslctd == 'tbl'",
                                                                  sliderInput(inputId = "no_of_de_genes", label = strong("Number of top DE genes:"), value = 10, min = 1, max = 100, step = 1),
                                                                  uiOutput("ct_de_table_dyn"),
                                                                  downloadButton("downloadData", "Download table of significantly DE genes")
                                                 )
                                                 
                                                 
                                )
                                
                              ),
                              
                              # Output: Description, lineplot, and reference
                              mainPanel(  
                                conditionalPanel(
                                  condition = "input.de_panel == 'Gene expression visualization between conditions across cell types'", 
                                  uiOutput("tb") ),
                                conditionalPanel(
                                  condition = "input.de_panel == 'Dotplot for gene expression between conditions and across cell types'", 
                                  plotOutput("marker_dotplot")),
                                conditionalPanel(
                                  condition = "input.de_panel == 'Differential gene expression between conditions per cell type'", 
                                  uiOutput("de_outputs"))
                                # conditionalPanel(
                                #   condition = "input.de_panel == 'Differential gene expression between conditions per cell type'", 
                                #   withSpinner(tableOutput("top_de_genes")))
                                
                                
                              )
                            )
                          ),
                          
                          tabPanel("Cluster adjustment",
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons(inputId = "graph_type", label = strong("Cluster adjustment by altering resolution"), choices = c("UMAP clusters", "Top cluster markers conserved between conditions", "Label populations")),
                              # Display depending on which graph type is selected

                              conditionalPanel(condition = "input.graph_type == 'UMAP clusters'",
                                               # Select PC to plot
                                               # sliderInput(inputId = "umap_dim", label = strong("Number of PCs from 1"), value = dim1, min = dim1, max = dim2, step = diff_dim),
                                               # sliderInput(inputId = "neighbours_dim", label = strong("Number of PCs from 1 (KNN clustering in FindNeighbours function)"), value = dim1, min = dim1, max = dim2, step = diff_dim),
                                               sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F),
                                               #actionButton(inputId = "go", label = "Run")
                                               
                              ),
                              conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between conditions'",
                                               
                                               #selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = , multiple = F),
                                               uiOutput("dyn_clusters"),
                                               selectInput(inputId = "select_markers_umap", label = strong("select marker to visualize in clusters:"), choices = all_genes_common_in_all_groups, multiple = T, selected = fav_genes),
                                               sliderInput(inputId = "marker_genes_no", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1)
                                               #actionButton(inputId = "go_marker", label = "Search"),
                                               
                              ),
                              conditionalPanel(condition = "input.graph_type == 'Label populations'",
                                               #Select PC to plot
                                               uiOutput("cluster_annot"),
                                               #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                              )
                            ),

                            # Output: Description, lineplot, and reference
                            mainPanel(

                              conditionalPanel(
                                condition = "input.graph_type == 'UMAP clusters'", 
                                plotOutput("all_groups"), plotOutput("stim_vs_ctrl")),
                              conditionalPanel(
                                condition = "input.graph_type == 'Top cluster markers conserved between conditions'", 
                                plotOutput("conserved_markers_umap"), tableOutput("top_conserved_genes")), 
                              conditionalPanel(
                                condition = "input.graph_type == 'Label populations'", 
                                plotOutput("labelled_umap"))
                              

                            )
                          )
                          )
                          
               )
)


server = function(input, output) {
  
  # umap_clusters = eventReactive(input$go, {
  #   tcells.combined.umap.list[[((input$umap_dim/5) - 1)]][[((input$umap_dim/5) - 1)]][[round(input$clusters_res/0.25)]]
  # })
  
  #function draw clusters
  #umap_clusters = eventReactive(input$go, {
  umap_clusters = reactive({
      
    tcells_combined_umap_list_res[[(input$clusters_res * 10)-0.5]]
  })
  
  output$all_groups = renderPlot({
    
    
    p1 = DimPlot(umap_clusters(), reduction = "umap", group.by = "sample", label.size = 6)
    p2 = DimPlot(umap_clusters(), reduction = "umap", label = T, label.size = 6)
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
  
  # cluster_markers = eventReactive(input$go_marker, {
  #   head(tcells.combined.clusters.tables[[((input$umap_dim/5) - 1)]][[((input$umap_dim/5) - 1)]][[round(input$clusters_res/0.25)]][[(as.numeric(input$marker_genes_cluster) + 1)]], n = input$marker_genes_no)
  # })
  #cluster_markers = eventReactive(input$go_marker, {
  cluster_markers = reactive({
    head(tcells_combined_clusters_tables_res[[(input$clusters_res * 10)-0.5]][[(as.numeric(input$marker_genes_cluster) + 1)]], n = input$marker_genes_no)
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
  
  # umap_cluster_modified_ren_reo = reactive({
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
  
  # umap_names = reactive({
  #   umap_cluster_modified = umap_cluster_modified_rna()
  #   new_cluster_name = unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
  #     new_cluster_name = input[[paste("labeller",x)]]
  #   }))
  #   #print(new_cluster_name)
  #   
  # })
  
  annotation = reactiveValues(annot = c("Il171", "Ccr7", "Lyc2", "Gzma", "Cdk6"))
  
  observe({
    umap_cluster_modified = umap_cluster_modified_rna()
    
    req(unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    })))
    annotation$annot = unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    }))
  })
  
  output$cluster_ids <- renderUI({
    umap_names = annotation$annot
    selectInput(inputId = "select_cell_type", label = strong("Select cell type to compare gene expression across conditions:"),choices = umap_names, multiple = F)
    
  })
  
  
  
  umap_cluster_modified_ren_reo = reactive({
    umap_cluster_modified1 = umap_cluster_modified_rna()
    umap_names = annotation$annot
    #print(umap_names)
    # old_cluster_name = unlist(as.character(0:(length(unique(umap_cluster_modified$seurat_clusters))-1)))
    # new_cluster_name = unlist(lapply(0:(length(unique(umap_cluster_modified$seurat_clusters))-1), function(x) {
    #   new_cluster_name = input[[paste("labeller",x)]]
    #   }))
    # print(old_cluster_name)
    # print(new_cluster_name)
    
    #umap_cluster_modified <- RenameIdents(umap_cluster_modified, `0` = "Il171",`1` = "Ccr7", `2` = "Lyc2", `3` = "Gzma", `4`  = "Cdk6")
    
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
      names(umap_names) <- levels(umap_cluster_modified1)
      umap_cluster_modified1 <- RenameIdents(umap_cluster_modified1, umap_names)
      
    } else {
      umap_cluster_modified1
    }
    
    #umap_cluster_modified@active.ident <- plyr::mapvalues(umap_cluster_modified@active.ident, from = old_cluster_name, to = new_cluster_name)
    #umap_cluster_modified = RenameIdents(umap_cluster_modified, old.ident.na = x, new.ident.name = new_cluster_name)
    
    #})
  })
  
  #umap_cluster_modified_umap = eventReactive(input$go_labelled_umap, {
  umap_cluster_modified_umap = reactive({
  
    DimPlot(umap_cluster_modified_ren_reo(), label = TRUE, label.size = 6)
  })
  output$labelled_umap = renderPlot({
    umap_cluster_modified_umap()
  })
  
  # umap_cluster_modified_ren_reo = reactive({
  #   umap_cluster_modified_reo = umap_cluster_modified_ren_reo()
  #   umap_names = annotation$annot
  #   Idents(umap_cluster_modified_reo) <- factor(Idents(umap_cluster_modified_reo), levels = umap_names)
  #   umap_cluster_modified_reo
  #   
  # })
  
  output$marker_dotplot = renderPlot({
    DotPlot(umap_cluster_modified_ren_reo(), features = input$select_markers_dotplot, cols = c("blue", "red", "green"), dot.scale = 6, split.by = "sample") + RotatedAxis()
    
  })
  
  stim_markers = reactive({
    
    umap_cluster_modified = umap_cluster_modified_ren_reo()
    print(length(unique(umap_cluster_modified$seurat_clusters)) != length(annotation$annot))
    
    # validate(
    #   need(try(length(unique(umap_cluster_modified$seurat_clusters)) != length(annotation$annot)), "Please label clusters in 'cluster adjustment' -> 'label cluster'")
    # )
    umap_cluster_modified$celltype.group <- paste(Idents(umap_cluster_modified), umap_cluster_modified$group, sep = "_")
    umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
    Idents(umap_cluster_modified) <- "celltype.group"
    umap_cluster_modified
  })
  
  genes_in_de_order = reactive({
    umap_names = annotation$annot
    #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO1", sep = "_"), ident.2 = paste(input$select_cell_type, input$condition2, sep = "_"), verbose = FALSE)
    FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
    
  })
  
  ## Define two panels for UMAP and violin plots
   

  
  output$ct_de_plot_dyn <- renderUI({
    actionButton(inputId = "ct_de", label = paste("Plot differential expression (DE) of all genes in ", input$select_cell_type , sep = ""))
  })
  
  output$ct_de_table_dyn <- renderUI({
    actionButton(inputId = "ct_de_table", label = paste("Display table of top DE genes in ", input$select_cell_type , sep = ""))
  })
  
  top_de_g = eventReactive(input$ct_de_table,{
    
    ##Finding DE genes in clusters between conditions
    head(genes_in_de_order(), n = input$no_of_de_genes) %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
  })
  
  output$top_de_genes = renderTable({
    top_de_g()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      #paste("differentially_expressed_genes_in",input$select_cell_type, "KO1","VS",input$condition2, ".csv", sep = "_")
      paste("differentially_expressed_genes_in",input$select_cell_type, "KO","VS","WT", ".csv", sep = "_")
      
    },
    content = function(file) {
      write.csv(genes_in_de_order() %>% filter(p_val_adj <= 0.05) , file)
    }
  )
  
  cell_type_de = eventReactive(input$ct_de,{
    
    cells_type <- subset(umap_cluster_modified_ren_reo(), idents = input$select_cell_type)
    #Idents(cells_type) <- "sample"
    Idents(cells_type) <- "group"
    avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
    avg.cells$gene <- rownames(avg.cells)
    #avg.cells <- merge(avg.cells, uniprot_info[,"Gene names  (primary )", "Protein names", "uniprot"], by.x = "gene", by.y = "Gene names (primary )", all.x = T)
    avg.cells <- avg.cells %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = "Gene names  (primary )")) %>% dplyr::distinct(., gene, .keep_all = T)%>% select(gene, KO, WT , `Protein names`, uniprot)
    avg.cells
    
  })
  
  cell_type_de_plot = eventReactive(input$ct_de,{
    theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=9, fontface="italic")))
    #ggplot(data=cell_type_de(), aes_string("KO1", input$condition2)) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob)
    ggplot(data=cell_type_de(), aes_string("KO", "WT")) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob)
    
    
  })
  
  output$cell_type_plot = renderPlot({
    theme_set(theme_cowplot())
    cell_type_de_plot()
    #ggplotly(p)
  })
  
  output$de_outputs <- renderUI({
    tabsetPanel(tabPanel("DE ggplot", value = "gg", 
                         #withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))), verbatimTextOutput("click_info")),
                         withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))), dataTableOutput("click_info")),
                tabPanel("DE table", value = "tbl",
                         withSpinner(tableOutput("top_de_genes"))), id = "tabslctd")
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
  
  output$click_info <- renderDataTable({
    req(displayed_text())
    
    #cat("Name\n")
    displayed_text = displayed_text()
    displayed_text 
  }, escape = F)
  
  
  ######################
  ### Differential expressino dynamic part
  ######################
  
  #Functions to updated differentially expressed genes
  output$de_stim_vs_ctrl_um = renderPlot({
    
    FeaturePlot(stim_markers(), features = input$de_genes, split.by = "sample", max.cutoff = 3,cols = c("grey", "red"))
  })
  
  output$de_stim_vs_ctrl_vp = renderPlot({
    
    VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE)  
    
  })
  
  ## Define two panels for UMAP and violin plots
  output$tb <- renderUI({
    tabsetPanel(tabPanel("UMAP plot", 
                         plotOutput("de_stim_vs_ctrl_um")),
                tabPanel("Violin plot", 
                         plotOutput("de_stim_vs_ctrl_vp"))) 
  })
  ######################
  ### END Differential expressino dynamic part
  ######################
  
  
}


shinyApp(ui = ui, server = server)
