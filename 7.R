
bmc4 <- bmc2$Control
save(bmc4,file="bm4.rdata")

DefaultAssay(bmc4) <- "RNA"
Idents(bmc4) <- "celltype"
ISG_genelist <- read.csv("ISG_genelist.csv")
ISG_genelist <- unique(as.vector(ISG_genelist$GeneName))
geneSets <- list(OSs=ISG_genelist)

cells_rankings_N <- AUCell_buildRankings(bmc4@assays$RNA@data, nCores=10,splitByBlocks=TRUE)  #基因排序，使用10个核

cells_AUC_N <- AUCell_calcAUC(geneSets, cells_rankings_N,nCores = 10, 
                            aucMaxRank = nrow(cells_rankings)*0.05) 

pdf("N_cells_assignment.pdf",height = 10, width = 15)
cells_assignment_N <- AUCell_exploreThresholds(cells_AUC_N, plotHist=TRUE, assign=TRUE, nCores=10)
dev.off()

geneSet.name <- "OSs"
AUC_Exp <- as.numeric(getAUC(cells_AUC)[geneSet.name, ])
bmc4$AUC <- AUC_Exp
plot.df_N<- data.frame(bmc4@meta.data, 
                         bmc4@reductions$umap@cell.embeddings)

class_avg <- plot.df_N %>%
  group_by(celltype) %>%      
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )


AUCP_N <- ggplot(plot.df_N, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_text_repel(aes(label = celltype),
                           data = class_avg,
                           size = 4.5,
                           box.padding = 0,
                           point.padding = 0, 
                           segment.color = NA)+   theme(legend.position = "none") + theme_classic()+
  ggtitle("Control") +
  theme(plot.title = element_text(hjust = 0.5))
AUCP_N
ggsave(filename = "Normal_UMAP_AUC_value.pdf", plot = AUCP_N, height = 10, width = 15)
