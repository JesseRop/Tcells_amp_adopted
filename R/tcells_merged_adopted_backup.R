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
setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/")

tcells.combined.umap.list<-readRDS("tcells.combined.umap.list.rds")
tcells.combined.clusters.tables<-readRDS("tcells.combined.clusters.tables.rds")
all_genes_final = row.names(tcells.combined.umap.list[[1]][[1]][[1]])

dim=15
dim1 = 10
dim2 = 15
diff_dim = 5

res1 = 0.15
res2 = 0.55
diff_res = 0.4

if (!(all(exists("tcells.combined.umap.list"), exists("tcells.combined.clusters.tables")))) {
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
  
  
  WT3.downsample = subset(WT3, cells = sample(Cells(WT3), 4000))
  
  
  
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
  Combined.filt <- FindClusters(Combined.filt, resolution = 0.15)
  
  
  ### decrease the amount of clusters, and split the other one (green, killers)
  # filtering away other clusters
  Combined.filt2 <- subset(Combined.filt, idents = c("0","1","2","3","5"), invert = FALSE)
  
  Combined.filt<-Combined.filt2
  # re-clusterting, without the other cell types
  DefaultAssay(object = Combined.filt) <- "integrated"
  # Run the standard workflow for visualization and clustering
  Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
  Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
  
  all_genes_final = row.names(Combined.filt)
  saveRDS(all_genes_final, "all_genes_final.rds")
  
  
  ##tSNE and clustering
  tcells.combined.umap.list = lapply(seq(dim1, dim2, by = diff_dim), function(x) {
    a = RunUMAP(Combined.filt, reduction = "pca", dims = 1:x)
    lapply(seq(dim1, dim2, by = diff_dim), function(x) {
      b = FindNeighbors(a, reduction = "pca", dims = 1:x)
      lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(b, resolution = x))
    })
  })
  
  saveRDS(tcells.combined.umap.list, "tcells.combined.umap.list.rds")
  
  tcells.combined.clusters.tables = lapply(tcells.combined.umap.list, function(x) { 
    lapply(x, function(x) {
      lapply(x, function(x){
        DefaultAssay(x) = "RNA"
        #cluster.markers = lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
        cluster.markers = lapply(0:1, function(y) {
          FindConservedMarkers(x, ident.1 = y, grouping.var = "sample")
          
        })
      }) 
    }) 
  })
  
  saveRDS(tcells.combined.clusters.tables, "tcells.combined.clusters.tables.rds")
  
  
}

#cell_types = c( "Il171", "Ccr7", "Lyc2", "Gzma", "Cdk6")

#######################################IFN

#######################################IFN


ui = fluidPage(theme = shinytheme("lumen"),
               navbarPage("III scRNA results",
                          tabPanel("About",
                                   
                                   mainPanel(width = 12,
                                             
                                             HTML(
                                               "
<h1>Welcome</h1>
<p>This site provides a platform for deeper interogation of sc-RNA outputs from the
<a href='https://www.gla.ac.uk/researchinstitutes/iii/'>Institute of Infection, Immunity & Inflammation (III), University of Glasgow</a>.
<p>


<h1>Published scRNAseq outputs</h1>
<p>
Please cite the relevant paper if you use our data in your work:
</p>

<button onclick='myFunction()'><div class='citation'>
<h4>Beta2 integrins differentially regulate γδ T cell subset thymic development and peripheral maintenance</h4>
<p class='citation-authors'>C. L. McIntyre, L. Monin, T. D. Otto, C. S. Goodyear, A. C. Hayday, and V. L. Morrison.
</p>
</div></button>



<div id='myDIV'>
<h4>Abstract</h4>
<p>γδ T cells reside predominantly at barrier sites and play essential roles in immune protection
against infection and cancer. Despite recent advances in the development of γδ T cell
immunotherapy, our understanding of the basic biology of these cells, including how their
numbers are regulated in vivo, remains poor. This is particularly true for tissue-resident γδ T cells.
We have identified the β 2 family of integrins as novel regulators of γδ T cells. β 2 integrin-deficient
mice displayed a striking increase in numbers of IL-17-producing Vγ6Vδ1 + γδ T cells in the lungs
and uterus, as well as circulation. Thymic development of this population was normal. However,
single cell RNA sequencing revealed the enrichment of genes associated with T cell survival and
proliferation specifically in β 2 integrin-deficient IL-17 + cells compared to their WT counterparts.
Indeed, β 2 integrin-deficient Vγ6 + cells from the lungs showed enhanced survival ex vivo,
suggesting that increased survival contributes to the accumulation of these cells in β 2 integrin-
deficient tissues. Furthermore, our data revealed an unexpected role for β 2 integrins in promoting
the thymic development of the IFNγ-producing CD27 + Vγ4 + γδ T cell subset. Together, our data
reveal that β 2 integrins play important roles in maintaining γδ T cell homeostasis, particularly in
regulating Vγ6 + cell numbers in mucosal tissues by controlling survival and Vγ4 + cell development
in the thymus. Our study indicates new and unprecedented mechanisms of control for γδ T cell
subsets.
</p>


</div>
<script>
function myFunction() {
  var x = document.getElementById('myDIV');
  if (x.style.display == 'none') {
    x.style.display = 'block';
  } else {
    x.style.display = 'none';
  }
}
</script>

<p></p>
<button onclick='myFunction1()'><div class='citation'>
<h4>Multiplexed droplet single-cell RNA-sequencing using natural genetic variation</h4>
<p class='citation-authors'>Hyun Min Kang, Meena Subramaniam, Sasha Targ, Michelle Nguyen, Lenka Maliskova, Elizabeth McCarthy, Eunice Wan, Simon Wong, Lauren Byrnes, Cristina M Lanata, Rachel E Gate, Sara Mostafavi, Alexander Marson, Noah Zaitlen, Lindsey A Criswell & Chun Jimmie Ye 
</p>
</div></button>



<div id='myDIV1'>
<h4>Abstract</h4>
<p>Droplet single-cell RNA-sequencing (dscRNA-seq) has enabled rapid, massively parallel profiling of transcriptomes. However, assessing differential expression across multiple individuals has been hampered by inefficient sample processing and technical batch effects. Here we describe a computational tool, demuxlet, that harnesses natural genetic variation to determine the sample identity of each droplet containing a single cell (singlet) and detect droplets containing two cells (doublets). These capabilities enable multiplexed dscRNA-seq experiments in which cells from unrelated individuals are pooled and captured at higher throughput than in standard workflows. Using simulated data, we show that 50 single-nucleotide polymorphisms (SNPs) per cell are sufficient to assign 97% of singlets and identify 92% of doublets in pools of up to 64 individuals. Given genotyping data for each of eight pooled samples, demuxlet correctly recovers the sample identity of >99% of singlets and identifies doublets at rates consistent with previous estimates. We apply demuxlet to assess cell-type-specific changes in gene expression in 8 pooled lupus patient samples treated with interferon (IFN)-β and perform eQTL analysis on 23 pooled samples.
</p>


</div>


<script>
function myFunction1() {
  var x = document.getElementById('myDIV1');
  if (x.style.display == 'none') {
    x.style.display = 'block';
  } else {
    x.style.display = 'none';
  }
}
</script>

") 
                                             # <div class='text-center'>
                                             #   <button id='button-data' type='button' class='btn btn-primary btn-lg'>Explore sc-RNA seq results for this work</button>
                                             #   </div>
                                             
                                   )
                          ),
                          
                          
                          tabPanel(
                            title = "Effect of Beta 2 integrin in T Cell differentiation",
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
                          ),
                          
                          tabPanel("Interferron beta stimulation of PBMCs",
                                   # sidebarLayout(
                                   #   sidebarPanel(
                                   #     radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("UMAP clusters", "Top cluster markers conserved between conditions", "Label populations", "Dotplot for gene expression between conditions and across cell types", "Differential gene expression between conditions per cell type", "Gene expression visualization between conditions across cell types")),
                                   #     # Display depending on which graph type is selected
                                   #     
                                   #     conditionalPanel(condition = "input.graph_type == 'UMAP clusters'",
                                   #                      # Select PC to plot
                                   #                      sliderInput(inputId = "umap_dim_ifn", label = strong("Number of PCs from 1"), value = dim1_ifn, min = dim1_ifn, max = dim2_ifn, step = diff_dim_ifn),
                                   #                      sliderInput(inputId = "neighbours_dim_ifn", label = strong("Number of PCs from 1 (KNN clustering in FindNeighbours function)"), value = dim1_ifn, min = dim1_ifn, max = dim2_ifn, step = diff_dim_ifn),
                                   #                      sliderInput(inputId = "clusters_res_ifn", label = strong("Louvain algorithm resolution"), value = res1_ifn, min = res1_ifn, max = res2_ifn, step = diff_res_ifn),
                                   #                      actionButton(inputId = "go_ifn", label = "Run")
                                   #                      
                                   #     ),
                                   #     conditionalPanel(condition = "input.graph_type == 'Top cluster markers conserved between conditions'",
                                   #                      
                                   #                      sliderInput(inputId = "marker_genes_no_ifn", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
                                   #                      #selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = , multiple = F),
                                   #                      uiOutput("dyn_clusters_ifn"),
                                   #                      actionButton(inputId = "go_ifn_marker", label = "Search"),
                                   #                      selectInput(inputId = "select_markers_umap_ifn", label = strong("Select conserved cluster markers:"), choices = all_genes_final, multiple = T)
                                   #                      
                                   #     ),
                                   #     conditionalPanel(condition = "input.graph_type == 'Label populations'",
                                   #                      #Select PC to plot
                                   #                      uiOutput("cluster_annot_ifn"),
                                   #                      actionButton(inputId = "go_ifn_labelled_umap", label = "View labelled clusters"),
                                   #     ),
                                   #     conditionalPanel(condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'",
                                   #                      #Select PC to plot
                                   #                      selectInput(inputId = "select_markers_dotplot_ifn", label = strong("Select markers for dotplot:"),choices = all_genes_final, multiple = T)
                                   #     ),
                                   #     conditionalPanel(condition = "input.graph_type == 'Differential gene expression between conditions per cell type'",
                                   #                      #Select PC to plot
                                   #                      #actionButton(inputId = "go_ifn_labelled_umap", label = "View labelled clusters"),
                                   #                      uiOutput("cluster_ids_ifn"),
                                   #                      uiOutput("ct_de_plot_dyn_ifn"),   
                                   #                      sliderInput(inputId = "no_of_de_genes_ifn", label = strong("Number of top differentially expressed (DE) genes:"), value = 10, min = 1, max = 100, step = 1),
                                   #                      uiOutput("ct_de_table_dyn_ifn"),   
                                   #                      
                                   #     ),
                                   #     conditionalPanel(condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'",
                                   #                      #Select PC to plot
                                   #                      selectInput(inputId = "de_genes_ifn", label = strong("Choose gene:"),choices = all_genes_final, multiple = F),
                                   #                      #actionButton(inputId = "go_ifn_heatmap", label = "Run"),
                                   #                      #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                                   #     )
                                   #   ),
                                   #   
                                   #   # Output: Description, lineplot, and reference
                                   #   mainPanel(  
                                   #     
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'UMAP clusters'", 
                                   #       plotOutput("all_groups"), plotOutput("stim_vs_ctrl")),
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'Top cluster markers conserved between conditions'", 
                                   #       tableOutput("top_conserved_genes"), plotOutput("conserved_markers_umap")), 
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'Label populations'", 
                                   #       plotOutput("labelled_umap")),
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'Dotplot for gene expression between conditions and across cell types'", 
                                   #       plotOutput("marker_dotplot")),
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
                                   #       plotOutput("cell_type_plot", hover = hoverOpts(id ="plot_hover")), verbatimTextOutput("hover_info")),
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'Differential gene expression between conditions per cell type'", 
                                   #       tableOutput("top_de_genes")),
                                   #     conditionalPanel(
                                   #       condition = "input.graph_type == 'Gene expression visualization between conditions across cell types'", uiOutput("tb") )
                                   #     
                                   #     
                                   #   )
                                   # )
                          )
                          
               )
)


server = function(input, output) {
  
  umap_clusters = eventReactive(input$go, {
    tcells.combined.umap.list[[((input$umap_dim/5) - 1)]][[((input$umap_dim/5) - 1)]][[round(input$clusters_res/0.25)]]
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
    DimPlot(umap_cluster_modified_renamed(), label = TRUE, label.size = 6)
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
    DotPlot(umap_cluster_modified_ren_reo(), features = input$select_markers_dotplot, cols = c("blue", "red", "green"), dot.scale = 6, split.by = "sample") + RotatedAxis()
    
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
