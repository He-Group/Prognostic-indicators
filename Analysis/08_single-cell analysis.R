library(dplyr)
library(Seurat)
library(patchwork)
data <- read.table('D:/项目/olink/result/MR analysis/单细胞测序/TPM_data.txt',header=T,comment.char = "#",sep = "\t",fill=TRUE,na.strings = "",row.names = NULL)
test.seu=CreateSeuratObject(counts = data)
scRNA1 <- FindVariableFeatures(test.seu, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA1), 10) 
plot1 <- VariableFeaturePlot(scRNA1)   
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
plot 
scale.genes <-  rownames(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = scale.genes)
load('D:/项目/olink/result/MR analysis/单细胞测序/scRNA1_PT_Phase.Rdata')

CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1)) 
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAa <- RunPCA(scRNA1, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAb, reduction = "pca", group.by = "Phase")
scRNAb <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA1))
VlnPlot(scRNAa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6) ###周期基因打分的小提琴图

#PCA
scRNA1 <- RunPCA(scRNAb, features = VariableFeatures(scRNAb)) 
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident") 
plot1
plot2 <- ElbowPlot(scRNA1, ndims=20, reduction="pca") 
plot2 
Loadings(object = scRNA1[["pca"]])  
Embeddings(object = scRNA1[["pca"]])  
Stdev(scRNA1)   
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
save(scRNA1,file = 'D:/项目/olink/result/MR analysis/单细胞测序/scRNA1_PT_deletPhase.Rdata')


load('D:/项目/olink/result/MR analysis/单细胞测序/scRNA1_PT_Phase.Rdata')
library(SingleR)
library(celldex)
library(dplyr)
library(Seurat)

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

##CD274
library(tidyverse)
library(harmony)
DotPlot <- DotPlot(scRNA1, group.by = 'celltype',features = 'CD274')
plotCD274 = FeaturePlot(scRNA1,features = "CD274",reduction = "tsne",pt.size = 1,order=T)
pdf('D:/项目/olink/result/MR analysis/单细胞测序/PT_Phase_tsne+CD274.pdf',height=6,width=22)
plottsne+plotCD274+DotPlot
dev.off()


##10X
library(dplyr)
library(Seurat)
library(ggplot2)
counts <- Read10X("D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo")
scRNA <- CreateSeuratObject(counts = counts,min.cells = 3, min.features = 200)
str(scRNA)
meta <- scRNA@meta.data

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
x <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
col.num <- length(levels(scRNA@active.ident))
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, 
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin

scRNA1 <- subset(scRNA, subset =  percent.mt < 10 & percent.HB < 3)
scRNA1 <- NormalizeData(scRNA1, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA1), 10) 
scale.genes <-  rownames(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = scale.genes)

cc.genes
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAa <- RunPCA(scRNA1, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
p

scRNAb <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA1))
scRNAb <- RunPCA(scRNAb, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAb, reduction = "pca", group.by = "Phase")
p
save(scRNA1,file='D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/scRNA1.Rdata')
save(scRNAb,file='D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/scRNAb.Rdata')

load('D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/scRNA1.Rdata')
options(future.globals.maxSize = 4000 * 1024^2)
scRNA1 <- SCTransform(scRNAb, 
                      method = "glmGamPoi", 
                      vars.to.regress = c("percent.mt","S.Score","G2M.Score"), 
                      verbose = T)
scRNA1 <- RunPCA(scRNAb, features = VariableFeatures(scRNAb))
#scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1)) 
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
save(scRNA1,file='D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/scRNA1.Rdata')

epithelial<- c('EPCAM', 'KRT8', 'KRT18') 
stromal <- c('COL1A1', 'COL1A2', 'COL6A1', 'COL6A2', 'VWF', 'PLVAP', 'CDH5', 'S100B')
myeloid <- c('CD68', 'XCR1', 'CLEC9A', 'CLEC10A', 'CD1C', 'S100A8', 'S100A9', 'TPSAB1', 'OSM')
T_cells <- c('NKG7', 'KLRC1', 'CCR7', 'FOXP3', 'CTLA4', 'CD8B', 'CXCR6', 'CD3D')
B_cells <- c('MZB1', 'IGHA1', 'SELL', 'CD19', 'AICDA')
markers <- c(epithelial,stromal,myeloid,T_cells,B_cells)
DotPlot(scRNA1, features = unique(markers),group.by = "seurat_clusters")+RotatedAxis()
###SCT，pc=1:40
new.cluster.ids <- c("epithelial cells","stromal cells","epithelial cells","epithelial cells","myeloid cells",
                     "T cells","epithelial cells", "B cells","stromal cells","B cells","epithelial cells",
                     "B cells","B cells","B cells")
###noSCT,pc=1:16
new.cluster.ids <- c("epithelial cells","stromal cells","epithelial cells","epithelial cells","myeloid cells",
                     "B cells","stromal cells","T cells","epithelial cells","stromal cells")
###noSCT,feature=1000,pc=1:14
new.cluster.ids <- c("stromal cells","epithelial cells","epithelial cells","myeloid cells","epithelial cells",
                     "epithelial cells","stromal cells","B cells","T cells")
###noSCT,deletPhase,pc=1:20
new.cluster.ids <- c("epithelial cells","epithelial cells","epithelial cells","stromal cells","myeloid cells",
                     "B cells","T cells","stromal cells","stromal cells","epithelial cells","stromal cells")

names(new.cluster.ids) <- levels(scRNA1) 
scRNA1 <- RenameIdents(scRNA1, new.cluster.ids)
scRNA1$celltype <- Idents(scRNA1)
clusterCols <- c("#E6550D","#8C6D31" , "#3182BD", "#54990F","#E7969C")
plotumap = DimPlot(scRNA1, group.by = 'celltype',reduction = "umap",label = TRUE,pt.size = 1)
plotumap
plottsne = DimPlot(scRNA1, group.by = 'celltype',reduction = "tsne",label = TRUE,pt.size=1,cols=clusterCols)
plottsne

library(tidyverse)
VlnPlot(scRNA1, features = 'CD274')
Dotplot=DotPlot(scRNA1, features = 'CD274',group.by = "celltype")
plotCD274=FeaturePlot(scRNA1,features = "CD274",reduction = "tsne",pt.size = 1,order = TRUE)
pdf('D:/项目/olink/result/MR analysis/空转/10X/scRNA/CRC_octo/Phase_noSCT_feature1000_14_tsne_CD274.pdf',height=6,width = 22)
plottsne+plotCD274+Dotplot
dev.off()
