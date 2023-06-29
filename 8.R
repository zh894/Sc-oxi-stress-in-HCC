library(AUCell)
setwd("~/oxi/AUCCELL")
DefaultAssay(obj.combined3) <- "RNA"
Idents(obj.combined3) <- "celltype"
ISG_genelist <- read.csv("OS_genelist.csv")
ISG_genelist <- unique(as.vector(ISG_genelist$GeneName))
geneSets <- list(OSs=ISG_genelist)

cells_rankings <- AUCell_buildRankings(obj.combined3@assays$RNA@data, nCores=10,splitByBlocks=TRUE)  #基因排序，使用10个核
# Calculates the 'AUC' for each gene-set in each cell
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores = 10, 
                            aucMaxRank = nrow(cells_rankings)*0.05) 

pdf("cells_assignment.pdf")
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE, nCores=10)
dev.off()

geneSet.name <- "OSs"
AUC_Exp <- as.numeric(getAUC(cells_AUC)[geneSet.name, ])
obj.combined3$AUC <- AUC_Exp
plot.df<- data.frame(obj.combined3@meta.data, 
                     obj.combined3@reductions$umap@cell.embeddings)
p <- ggplot() + 
  geom_point(data=plot.df, aes(x=UMAP_1,y=UMAP_2,colour=AUC), size =0.5) +
  viridis::scale_color_viridis(option="A") +
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text= element_text(colour= 'black',size=14),
        axis.title= element_text(size = 14),
        axis.line= element_line(colour= 'black'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black"), 
        aspect.ratio = 1)
p
ggsave(filename = "UMAP_AUC_value.pdf", plot = p, height = 15, width = 15)

geneSets <- list(ISG_genelist)
Inscore <- AddModuleScore(obj.combined3,
                          features = geneSets,
                          ctrl = 500,
                          name = "geneSets")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[15] <- 'Oxidative Stress Score' 
compaired <-list(c("Dendritic cells_HCC","Dendritic cells_Control"),c("Gamma delta T cells_HCC","Gamma delta T cells_Control"),
                 c("Hepatocytes_HCC","Hepatocytes_Control"),c("Mucosal-associated invariant T cells_HCC","Mucosal-associated invariant T cells_Control"),
                 c("T memory cells_HCC","T memory cells_Control"),c("B cells_HCC","B cells_Control"))

VlnPlot(Inscore,features = 'Oxidative Stress Score', 
        pt.size = 0, adjust = 2,group.by = "celltype.group")+
  geom_signif(comparisons = compaired,test=t.test,map_signif_level = F)

data<- FetchData(Inscore,vars = c("celltype.group","Oxidative Stress Score"))

AUCBOXPLOT <- ggplot(data, aes(x=celltype.group,y=`Oxidative Stress Score`,fill = celltype.group)) +
  theme_bw()+RotatedAxis()+ 
  theme(panel.grid = element_blank(),
       axis.text.x=element_text(angle=0,hjust = 1,vjust=1))+
  #labs(x=NULL,y=NULL,title = "Oxidative Stress Score")+ 
  #geom_jitter(col="#00000033", pch=19,cex=2.5, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0),outlier.shape = 19,colour="grey")+labs(x = NULL)+
  #NoLegend()+theme(plot.title = element_text(hjust = 0.5)) +#也就加上这一行,标题居中
  geom_signif(comparisons = compaired,test=wilcox.test,map_signif_level = F) +
  #scale_color_discrete()+
  #facet_grid(.~"Oxidative Stress Score")+ 
  theme(legend.position = "right")+
  #scale_fill_simpsons()+
coord_flip()


p_top <- ggplot(data, aes(x = `Oxidative Stress Score`, color = celltype.group,fill=celltype.group))+
  geom_density(alpha=0.5,adjust=1.75,lwd=0.5) +
 
  theme_classic() + ylab(NULL) + 
  theme(legend.position = "none", 
        #legend.title = element_blank(),
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_blank(), # 原文不显示纵轴的密度
        #axis.text.y = element_text(size = 12,color = "black"), # 如果要显示采用这一行
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()+labs(x="")
p_top 
