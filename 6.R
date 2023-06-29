library(AUCell)
setwd("~/oxi/AUCCELL")
class(obj.combined3)

bmc <- subset(obj.combined3, SplitObject())
bmc  <- obj.combined3
Idents(bmc) <- "group"
bmc2 <- SplitObject(bmc,split.by = "group")
bmc3 <- bmc2$HCC

DefaultAssay(bmc3) <- "RNA"
Idents(bmc3) <- "celltype"
ISG_genelist <- read.csv("OS_genelist.csv")
ISG_genelist <- unique(as.vector(ISG_genelist$GeneName))
geneSets <- list(OSs=ISG_genelist)

cells_rankings_HCC <- AUCell_buildRankings(bmc3@assays$RNA@data, nCores=10,splitByBlocks=TRUE)  #基因排序，使用10个核

cells_AUC_HCC <- AUCell_calcAUC(geneSets, cells_rankings_HCC,nCores = 10, 
                            aucMaxRank = nrow(cells_rankings)*0.05) 

pdf("HCC_cells_assignment.pdf",height = 10, width = 15)
cells_assignment_HCC <- AUCell_exploreThresholds(cells_AUC_HCC, plotHist=TRUE, assign=TRUE, nCores=10)
dev.off()

geneSet.name <- "OSs"
AUC_Exp <- as.numeric(getAUC(cells_AUC)[geneSet.name, ])
bmc3$AUC <- AUC_Exp
plot.df_HCC<- data.frame(bmc3@meta.data, 
                     bmc3@reductions$umap@cell.embeddings)

class_avg <- plot.df_HCC %>%
  group_by(celltype) %>%       
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
#

AUCP_HCC <- ggplot(plot.df_HCC, aes(UMAP_1, UMAP_2),title())  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
  ggrepel::geom_text_repel(aes(label = celltype),
                            data = class_avg,
                            size = 4.5,
                            box.padding = 0,
                            point.padding = 0, 
                            segment.color = NA)+   theme(legend.position = "none") + theme_classic()+
  
  ggtitle("HCC") +
  theme(plot.title = element_text(hjust = 0.5))
AUCP_HCC
ggsave(filename = "HCC_UMAP_AUC_value.pdf", plot = AUCP_HCC, height = 10, width = 15)
