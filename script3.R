setwd("~/oxi/cell type")
logFCfilter=0.25               
adjPvalFilter=0.05         
#鉴定marker基因
sample.markers <- FindAllMarkers(object = obj.combined,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
sig.markers=sample.markers[(abs(as.numeric(as.vector(sample.markers$avg_log2FC)))>=logFCfilter & as.numeric(as.vector(sample.markers$p_val_adj))<=adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
save(sample.markers,file = "sample.markers.rdata")
save(sig.markers,file = "sig.markers.rdata")

top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster",names(top10_table))


n=20
cutoff <-ifelse(min(table(sample.markers$cluster))>n,n,min(table(sample.markers$cluster)))
topN <- sample.markers %>% group_by(cluster) %>% top_n(n = cutoff, wt = avg_log2FC)
write.csv(file = paste0("4_top",cutoff,"_cell_markers.csv"),topN,row.names=F)
topN_table=unstack(topN, gene ~ cluster)
names(topN_table)=gsub("X","cluster",names(topN_table))


obj.combined2 <-obj.combined

obj.combined2 <- RenameIdents(object = obj.combined2,
                                "0" = "Dendritic cells",
                                "1" = "Gamma delta T cells", 
                                "2" = "Mucosal-associated invariant T cells",
                                "3" ="Hepatocytes",
                                "4"="Hepatocytes",
                                "5"="Hepatocytes",
                                "6"="Gamma delta T cells",
                                "7"="T memory cells",
                                "8"= "Hepatocytes",
                                "9"= "Hepatocytes",
                                "10"= "Hepatocytes",
                                "11"= "B cells",
                              "12"="Hepatocytes",
                              "13"="Dendritic cells",
                              "14"="B cells",
                              "15"="T memory cells",
                              "16"="Hepatocytes",
                              "17"="B cells",
                              "18"="Hepatocytes",
                              "19"="Hepatocytes")
save(obj.combined2,file="obj.combined2.rdata")
sample.markers2 <- FindAllMarkers(object = obj.combined2,
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0)


pdf('4_cluster_annotations.pdf',width = 23,height = 10)
DimPlot(obj.combined2, reduction = "umap", label = TRUE, pt.size = 1,split.by = "group")+
  theme(legend.position = "right")+
  scale_color_npg()
dev.off()

pdf('4_cluster_annotations2.pdf',width = 15,height = 10)
DimPlot(obj.combined2, 
        reduction = "umap",
        pt.size = 1,
        group.by = "group",
        label = F)+
  theme(legend.position = "right")+
  scale_color_npg()
dev.off()


col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))
top10 <- sample.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
pdf('4_Heatmap_cluster_name.pdf', width = 15, height = 15)
options(repr.plot.width = 13, repr.plot.height=13)
DoHeatmap(obj.combined2, features = top10$gene, size = 3, angle =0,hjust = 0.2,group.colors = col)+
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()


