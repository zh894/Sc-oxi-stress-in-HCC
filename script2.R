
library(Seurat)
library(patchwork)
library(dplyr)
library(SeuratData)
library(patchwork)
library(dplyr)
library(pheatmap)
library(magrittr)
library(Biobase)
library(glmGamPoi)
library(ggsci)
library(SCpubr)
setwd("~/oxi/DR")

obj.list <- lapply(obj.list,function(x) {
  SCTransform(x,method = "glmGamPoi", vars.to.regress = "percent.MT", verbose = FALSE)
})


obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)

obj.combined <- IntegrateData(anchorset = obj.anchors, dims = 1:30)
DefaultAssay(obj.combined) <- "integrated"


all.genes <- rownames(obj.combined[["RNA"]]@data)
length(all.genes)
obj.combined <- ScaleData(obj.combined, features = all.genes)


obj.combined <- RunPCA(obj.combined, npcs = 50, verbose = FALSE)

pdf("PCA-ElbowPlot.pdf",width = 10,height = 8)
ElbowPlot(obj.combined,ndims = 30)
dev.off()
save(obj.combined,file = "pca_obj.combined.rdata")

obj.combined <- RunUMAP(obj.combined, dims = 1:20)
obj.combined <- FindNeighbors(obj.combined, dims = 1:20)
obj.combined <- FindClusters(obj.combined, resolution = 0.25,algorithm= 4)

pdf("UMAP_cluster1.pdf", width=20,height=10)
p1 <- DimPlot(obj.combined , 
                 reduction = "umap",
                 pt.size =0.9,
                 group.by = "group",
                 label = F)+
  scale_color_d3("category20")
p2 <- DimPlot(obj.combined, reduction = "umap",
                 pt.size =0.9,
                 label = F)+
  scale_color_d3("category20")
p1 + p2
dev.off()


pdf("UMAP_cluster_split.pdf",width = 20,height = 10)
p3 <- DimPlot(obj.combined, reduction = "umap",  split.by = "group",label = F,
           pt.size = 0.9)+
  scale_color_d3("category20")
print(p3)
dev.off()

save(obj.combined,file = "umap_obj.combined.rdata")
