library(cowplot)
setwd("~/oxi/diff genes in cells")
obj.combined3 <- obj.combined2


DefaultAssay(obj.combined3) <- "RNA"


obj.combined3$celltype <- Idents(obj.combined3)


obj.combined3$celltype.group <- paste(Idents(obj.combined3), obj.combined3$group, sep = "_")

Idents(obj.combined3) <- "celltype.group"

DC.diff <- FindMarkers(obj.combined3, ident.1 = "Dendritic cells_HCC", ident.2 = "Dendritic cells_Control", verbose = FALSE)
DC.diff  <- data.frame(DC.diff  ,gene=rownames(DC.diff ))

GDT.diff <- FindMarkers(obj.combined3, ident.1 = "Gamma delta T cells_HCC", ident.2 = "Gamma delta T cells_Control", verbose = FALSE)
GDT.diff  <- data.frame(GDT.diff  ,gene=rownames(GDT.diff ))

Hepatocytes.diff <- FindMarkers(obj.combined3, ident.1 = "Hepatocytes_HCC", ident.2 = "Hepatocytes_Control", verbose = FALSE)
Hepatocytes.diff  <- data.frame(Hepatocytes.diff  ,gene=rownames(Hepatocytes.diff ))

MIT.diff <- FindMarkers(obj.combined3, ident.1 = "Mucosal-associated invariant T cells_HCC", ident.2 = "Mucosal-associated invariant T cells_Control", verbose = FALSE)
MIT.diff  <- data.frame(MIT.diff  ,gene=rownames(MIT.diff ))

TM.diff <- FindMarkers(obj.combined3, ident.1 = "T memory cells_HCC", ident.2 = "T memory cells_Control", verbose = FALSE)
TM.diff  <- data.frame(TM.diff  ,gene=rownames(TM.diff ))

B.diff <- FindMarkers(obj.combined3, ident.1 = "B cells_HCC", ident.2 = "B cells_Control", verbose = FALSE)
B.diff  <- data.frame(B.diff  ,gene=rownames(B.diff ))

Idents(obj.combined3) <- "celltype"
theme_set(theme_cowplot())

DC <- subset(obj.combined3, idents = "Dendritic cells")
Idents(DC) <- "group"
avg.DC <- log1p(AverageExpression(DC, verbose = FALSE)$RNA)
avg.DC <- data.frame(avg.DC ,gene=rownames(avg.DC))


GDT<- subset(obj.combined3, idents = "Gamma delta T cells")
Idents(GDT) <- "group"
avg.GDT <- log1p(AverageExpression(GDT, verbose = FALSE)$RNA)
avg.GDT <- data.frame(avg.GDT ,gene=rownames(avg.GDT))


Hepatocytes <- subset(obj.combined3, idents = "Hepatocytes")
Idents(Hepatocytes) <- "group"
avg.Hepatocytes <- log1p(AverageExpression(Hepatocytes, verbose = FALSE)$RNA)
avg.Hepatocytes <- data.frame(avg.Hepatocytes ,gene=rownames(avg.Hepatocytes))


MIT<- subset(obj.combined3, idents = "Mucosal-associated invariant T cells")
Idents(MIT) <- "group"
avg.MIT <- log1p(AverageExpression(MIT, verbose = FALSE)$RNA)
avg.MIT <- data.frame(avg.MIT ,gene=rownames(avg.MIT))


TM <- subset(obj.combined3, idents = "T memory cells")
Idents(TM) <- "group"
avg.TM <- log1p(AverageExpression(TM, verbose = FALSE)$RNA)
avg.TM <- data.frame(avg.TM ,gene=rownames(avg.TM))


B <- subset(obj.combined3, idents = "B cells")
Idents(B) <- "group"
avg.B <- log1p(AverageExpression(B, verbose = FALSE)$RNA)
avg.B <- data.frame(avg.B ,gene=rownames(avg.B))

