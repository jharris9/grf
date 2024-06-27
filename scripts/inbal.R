
source("/home/rstudio/nmf/scripts/utils.R")
suppressPackageStartupMessages({
  library(RcppML)
  library(tidyverse)
  library(cowplot)
  library(Matrix)
  library(data.table)
  library(Seurat)
})


base_path="/home/rstudio/nmf/"
data_path=paste0(base_path,"raw_data/inbal/")
figure_path=paste0(base_path,"output/figures/inbal/")

cluster_data = function(obj.cluster, resolution,return.model=FALSE,run_umap=FALSE){
  DefaultAssay(obj.cluster) <- "RNA"
  if(!("data" %in% names(obj.cluster@assays[["RNA"]]@layers))){
    obj.cluster <- NormalizeData(obj.cluster)
  }
  for (l in Layers(obj.cluster, search='data')){
    obj.cluster[["RNA"]][l] <- as(obj.cluster[["RNA"]][l], Class = "dgCMatrix")
  }
  obj.cluster <- FindVariableFeatures(obj.cluster)
  obj.cluster <- ScaleData(obj.cluster,features=Features(obj.cluster))
  #obj.cluster <- SCTransform(obj.cluster)
  gc()
  obj.cluster <- RunPCA(obj.cluster, seed.use = 42, npcs=25)
  obj.cluster <- FindNeighbors(obj.cluster, dims = 1:25)
  obj.cluster <- FindClusters(obj.cluster, resolution = resolution, cluster.name = "unintegrated_clusters", random.seed=42)
  if (run_umap){
    obj.cluster <- RunUMAP(obj.cluster, dims = 1:25, reduction = "pca", return.model=return.model, reduction.name = "umap",seed.use = 42)
  }
  return(obj.cluster)
}

file <- gzfile(paste0(data_path,'GSE199317_ONC-retina.mtx.gz'), 'rt')
mat <- readMM(file)

colnames(mat) = fread(paste0(data_path,"GSE199317_ONC-retina_cell.tsv"),sep="\t",header=F)$V1
row.names(mat) = fread(paste0(data_path,"GSE199317_ONC-retina_gene.tsv"),sep="\t",header=F)$V1
mat <- as(mat, "CsparseMatrix")
obj.inbal <-CreateSeuratObject(counts = mat)
cell_types = fread(paste0(data_path,"GSE199317_ONC-retina_celltype.tsv.gz"),sep="\t",header=TRUE)
obj.inbal$cell_type_major = cell_types$cell_type_major


obj.inbal <- JoinLayers(obj.inbal)
obj.inbal = cluster_data(obj.inbal,2,T,T)

plot=DimPlot(obj.inbal,reduction="umap",group.by = "cell_type_major",label = T)
ggsave(file="inbal_et_al_subtype_umap.jpg", path="/home/rstudio/grf/output/figures/inbal/", device="jpeg",dpi="retina",
       units="in",height=8,width=11, plot=plot)

saveRDS(
  object = obj.inbal,
  file = "/home/rstudio/grf/output/processed_data/inbal/obj.inbal.Rds",
)

obj.inbal <- readRDS("/home/rstudio/grf/output/processed_data/obj.inbal.Rds")
