library(gtable)
library(cowplot)
library(gridExtra)
setwd("~/oxi/cell ratio")

cell.num <- table(Idents(obj.combined2))
cell.freq <- round(prop.table(table(Idents(obj.combined2)))*100,2)
cell.combined <- rbind(cell.num, cell.freq)


cell.num.group <- table(Idents(obj.combined2), obj.combined2$group) 
colnames(cell.num.group) <- paste0(colnames(cell.num.group),'_cell_counts')
cell.freq.group <- round(prop.table(table(Idents(obj.combined2), obj.combined2$group), margin = 2) *100,2)
colnames(cell.freq.group) <- paste0(colnames(cell.freq.group),'_cell_Freq')
cell.group <- cbind(cell.num.group, cell.freq.group)


pdf('5_UMAP_multi_samples_split_anno.pdf',width = 20)
p<-DimPlot(obj.combined2, reduction = "umap", split.by = "group",label=T,repel=T) + NoLegend()
tb <- tableGrob(cell.group)
plot_grid(p, tb,ncol=2,rel_widths=c(0.6,0.4))
dev.off()

table(obj.combined2$group)#查看各组细胞数
prop.table(table(Idents(obj.combined2)))
table(Idents(obj.combined2), obj.combined2$group)#各组不同细胞群细胞数
Cellratio2 <- prop.table(table(Idents(obj.combined2), obj.combined2$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio2
Cellratio2 <- as.data.frame(Cellratio2)
colourCount = length(unique(Cellratio2$Var1))
pdf('5_stackgroup.pdf',width = 8,height = 1)
p1 <- ggplot(Cellratio2) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.9,size = 0.5)+ 
  theme_classic() +
  labs(x='Group',y = 'Ratio')+
  #coord_flip()+
  scale_fill_npg()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75))+
  guides(fill=F)
print(p1)
dev.off()

table(obj.combined2$orig.ident)
prop.table(table(Idents(obj.combined2)))
table(Idents(obj.combined2), obj.combined2$orig.ident)
Cellratio3<- prop.table(table(Idents(obj.combined2), obj.combined2$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio3
Cellratio3 <- as.data.frame(Cellratio3)
colourCount = length(unique(Cellratio3$Var1))
pdf('5_stacksample.pdf',width = 7,height = 3)
p2 <- ggplot(Cellratio3) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.9,size = 0.5)+ 
  theme_classic() +
  labs(x='Sample',y = '')+
  scale_fill_npg()+
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(angle = 45,vjust = 0.85,hjust = 0.75))
print(p2)
dev.off()

pdf('5_stackremix.pdf',width = 7,height =10)
plot_grid(p1, p2, align = "h", rel_widths = c(1, 5),rel_heights = c(1,1),axis = "tlbr",ncol = 2)
dev.off()
