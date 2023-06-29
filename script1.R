library(Seurat)
library(dplyr)
library(stringr)

path<-"~/oxi/source data"
files <- list.files(path = path,pattern = 'txt$')

fullpath<-paste(path,files,sep="/")
names(fullpath) <- lapply(files,function(x) {str_extract(x, '[\\w\\d]+')})

rawdata <- lapply(fullpath, function(x) {
  cells.table <- read.table(x, sep = "\t",check.names = F, 
                            header = TRUE, row.names = 1)
})
lapply(rawdata, dim)

obj.list <- lapply(names(rawdata), function(x) {
  CreateSeuratObject(counts = rawdata[[x]],
                     project  = x,
                     min.cells = 5,    #??????????3??ϸ???????????Ļ???     
                     min.features = 200)})  #???????ٱ???200????????ϸ??
names(obj.list) <- names(rawdata)
obj.list


obj.list$control@meta.data$group <- "Control"
obj.list$HCC@meta.data$group <- "HCC"


HB_all <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
for (x in names(obj.list)){
  obj.list[[x]][["percent.MT"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^MT-")
  obj.list[[x]][["percent.Ribo"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^RP[SL]")
  HB_genes <- intersect(HB_all, rownames(obj.list[[x]]))
  obj.list[[x]][["percent.HB"]] <- PercentageFeatureSet(obj.list[[x]], features  = HB_genes)
}

library(ggsci)
library(tidyverse)
library(ggplot2)
setwd("~/oxi/qc")

qc_feature <- c("nFeature_RNA", "nCount_RNA", "percent.HB", "percent.MT",  "percent.Ribo")
for(sample in names(obj.list)){
  pdf(file=paste0("1_",sample,"_quality_control.pdf"),width = 20,height=7)
  print(VlnPlot(obj.list[[sample]], features = qc_feature, ncol = 5, pt.size = 0.5
                ))
  dev.off()
}



for(sample in names(obj.list)){
  pdf(file=paste0("1_",sample,"_feature_relationship.pdf"),width = 12,height=10)
  p1 <- FeatureScatter(obj.list[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- FeatureScatter(obj.list[[sample]], feature1 = "nCount_RNA", feature2 = "percent.HB")
  p3 <- FeatureScatter(obj.list[[sample]], feature1 = "nCount_RNA", feature2 = "percent.MT")
  p4 <- FeatureScatter(obj.list[[sample]], feature1 = "nCount_RNA", feature2 = "percent.Ribo")
  print(p1 + p2 + p3 + p4)
  dev.off()
}

gene.freq <- do.call("cbind", tapply(obj.list$control@meta.data$nFeature_RNA,
                                     obj.list$control@meta.data$orig.ident,quantile,
                                     probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(obj.list$control@meta.data$nCount_RNA,
                                    obj.list$control@meta.data$orig.ident,quantile,
                                    probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(obj.list$control@meta.data$percent.MT,
                                   obj.list$control@meta.data$orig.ident,quantile,
                                   probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"))
write.table(freq.combine,file = "QC-gene_frequency.txt",quote = F,sep = "\t")
rm(gene.freq,rna.freq,mt.freq)
View(freq.combine)

gene.freq2 <- do.call("cbind", tapply(obj.list$HCC@meta.data$nFeature_RNA,
                                      obj.list$HCC@meta.data$orig.ident,quantile,
                                      probs=seq(0,1,0.05)))
rna.freq2<- do.call("cbind", tapply(obj.list$HCC@meta.data$nCount_RNA,
                                    obj.list$HCC@meta.data$orig.ident,quantile,
                                    probs=seq(0,1,0.05)))
mt.freq2 <- do.call("cbind", tapply(obj.list$HCC@meta.data$percent.MT,
                                    obj.list$HCC@meta.data$orig.ident,quantile,
                                    probs=seq(0,1,0.05)))
freq.combine2 <- as.data.frame(cbind(gene.freq2,rna.freq2,mt.freq2))
colnames(freq.combine2) <- c(paste(colnames(gene.freq2),"Gene",sep = "_"),
                             paste(colnames(rna.freq2),"RNA",sep = "_"),
                             paste(colnames(mt.freq2),"MT",sep = "_"))
write.table(freq.combine2,file = "QC-gene_frequency2.txt",quote = F,sep = "\t")
rm(gene.freq2,rna.freq2,mt.freq2)
View(freq.combine2)


obj.list <- lapply(obj.list, function(x) {
  subset(x, subset = nFeature_RNA > 300 & percent.MT <5 & nCount_RNA > 500 & 
           nCount_RNA < 20000)})


