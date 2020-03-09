
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
