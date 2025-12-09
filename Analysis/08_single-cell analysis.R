
library(dplyr)
library(Seurat)
library(patchwork)
#GEO
load('D:/项目/olink/result/MR analysis/单细胞测序/scRNA1.Rdata')
scRNA1 <- RunPCA(scRNAb, features = VariableFeatures(scRNA1)) 
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident") 
plot1
plot2 <- ElbowPlot(scRNA1, ndims=20, reduction="pca") 
plot2 
Loadings(object = scRNA1[["pca"]])  
Embeddings(object = scRNA1[["pca"]])  
Stdev(scRNA1)   
pc.num=1:18
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
table(scRNA1@meta.data$seurat_clusters)   
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
embed_tsne <- Embeddings(scRNA1, 'tsne')
plot1 = DimPlot(scRNA1, reduction = "tsne",label = TRUE)
plot1
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
embed_umap <- Embeddings(scRNA1, 'umap')
plot2 = DimPlot(scRNA1, reduction = "umap",label = TRUE) 
plot2
library(SingleR)
library(celldex)
logFCfilter=0.5
adjPvalFilter=0.05
pbmc.markers <- FindAllMarkers(scRNA1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logFCfilter)
x <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)
sig.markers <- pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(pbmc.markers,'D:/项目/olink/result/MR analysis/UKB生存结局/单细胞测序/pbmc.markers.txt',sep=" ",row.names=F,quote=F)
epithelial<- c('EPCAM', 'KRT8', 'KRT18') 
stromal <- c('COL1A1', 'COL1A2', 'COL6A1', 'COL6A2', 'VWF', 'PLVAP', 'CDH5', 'S100B')
myeloid <- c('CD68', 'XCR1', 'CLEC9A', 'CLEC10A', 'CD1C', 'S100A8', 'S100A9', 'TPSAB1', 'OSM')
T_cells <- c('NKG7', 'KLRC1', 'CCR7', 'FOXP3', 'CTLA4', 'CD8B', 'CXCR6', 'CD3D')
B_cells <- c('MZB1', 'IGHA1', 'SELL', 'CD19', 'AICDA')
markers <- c(epithelial,stromal,myeloid,T_cells,B_cells)
DotPlot(scRNA1, features = unique(markers),group.by = "seurat_clusters")+RotatedAxis()
new.cluster.ids <- c("epithelial cells","B cells","epithelial cells","epithelial cells","epithelial cells",
                     "stromal cells","epithelial cells","myeloid cells","epithelial cells",
                     "stromal cells","T cells","B cells","T cells")  
names(new.cluster.ids) <- levels(scRNA1) 
scRNA1 <- RenameIdents(scRNA1, new.cluster.ids)
scRNA1$celltype <- Idents(scRNA1)
clusterCols <- c( "#8C6D31","#E7969C", "#E6550D", "#3182BD", "#54990F",
                  "#BD9E39", "#843C39", "#31A354", "#E41A1C", "#6BAED6",
                  "#9ECAE1", "#AD494A", "#A1D99B", "#C7E9C0", "#99600F",
                  "#C3BC3F", "#D6616B", "#FF7F00", "#1B9E77", "#FDAE6B", 
                  "#66A61E", "#F1788D", "#E6550D")
plotumap = DimPlot(scRNA1, group.by = 'celltype',reduction = "umap",label = TRUE,pt.size=1)
plotumap
plottsne = DimPlot(scRNA1, group.by = 'celltype',reduction = "tsne",label = T,pt.size=1,,cols=clusterCols)
plottsne
library(tidyverse)
library(harmony)
DotPlot <- DotPlot(scRNA1, group.by = 'celltype',features = 'CD274')
plotCD274 = FeaturePlot(scRNA1,features = "CD274",reduction = "tsne",pt.size = 1,order=T)
pdf('D:/项目/olink/result/MR analysis/单细胞测序/PT_Phase_tsne+CD274.pdf',height=6,width=22)
plottsne+plotCD274+DotPlot
dev.off()


#10X
load('D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/scRNA1.Rdata')
options(future.globals.maxSize = 4000 * 1024^2)
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
ElbowPlot(scRNA1, ndims=20, reduction="pca") 
pc.num = 1:14
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
table(scRNA1@meta.data$seurat_clusters)
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
embed_tsne <- Embeddings(scRNA1, 'tsne')
DimPlot(scRNA1, reduction = "tsne")

scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
embed_umap <- Embeddings(scRNA1, 'umap')
DimPlot(scRNA1, reduction = "umap")
new.cluster.ids <- c("stromal cells","epithelial cells","epithelial cells","myeloid cells","epithelial cells",
                     "epithelial cells","stromal cells","B cells","T cells")
names(new.cluster.ids) <- levels(scRNA1) 
scRNA1 <- RenameIdents(scRNA1, new.cluster.ids)
scRNA1$celltype <- Idents(scRNA1)
clusterCols <- c("#E6550D","#8C6D31" , "#3182BD", "#54990F","#E7969C")
plotumap = DimPlot(scRNA1, group.by = 'celltype',reduction = "umap",label = TRUE,pt.size = 1)
plotumap
plottsne = DimPlot(scRNA1, group.by = 'celltype',reduction = "tsne",label = TRUE,pt.size=1,cols=clusterCols)
plottsne
VlnPlot(scRNA1, features = 'CD274')
Dotplot=DotPlot(scRNA1, features = 'CD274',group.by = "celltype")
plotCD274=FeaturePlot(scRNA1,features = "CD274",reduction = "tsne",pt.size = 1,order = TRUE)
pdf('D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/Phase_noSCT_feature1000_14_tsne_CD274.pdf',height=6,width = 22)
plottsne+plotCD274+Dotplot
dev.off()
