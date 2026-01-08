#10x
load('D:/项目/olink/result/MR analysis/空转/10X/CRC.Rdata')
library(imager)
a <- SpatialFeaturePlot(CRC, features =c("EPCAM"),,pt.size.factor = 2.3, alpha = c(0.1, 1))
b <-  SpatialDimPlot(CRC, label = TRUE, label.size = 3,pt.size.factor = 2.3)
tumor <- subset(CRC, idents = c(1,2,3,4,8,6,11))  
c <- SpatialDimPlot(tumor, crop = FALSE, label = TRUE, pt.size.factor = 2.8, label.size = 3)
tumor1 <- RenameIdents(tumor,
                       "1"= "tumor",
                       "2"="tumor", 
                       "3"="tumor", 
                       "4"= "tumor", 
                       "8"= "tumor", 
                       "6"= "tumor",
                       "11"= "tumor",
                       "0"= "no_tumor",
                       "5"= "no_tumor",
                       "9"= "no_tumor",
                       "7"= "no_tumor",
                       "10"= "no_tumor")
c <- SpatialDimPlot(tumor1, crop = FALSE, label = F, pt.size.factor = 2.2, label.size = 2.5,cols = c('tumor'='#1DB2B7'))
c
markers <- FindMarkers(tumor1, ident.1 = 'tumor', ident.2 = 'no_tumor')
save(markers,file='D:/项目/olink/result/MR analysis/空转/1/markers.Rdata')
p3 <- SpatialFeaturePlot(CRC, features =c("CD274"), pt.size.factor = 2.8, alpha = c(0, 1))
p3
x <- c+p3
x
pdf("D:/项目/olink/result/MR analysis/空转/10X/scatter_plot.pdf", width = 6, height = 6)
x
dev.off()
topptx(x,file="D:/项目/olink/result/MR analysis/空转/1/CD274.pptx",width = 8,height = 7)

##WCH
load('D:/项目/olink/result/MR analysis/空转/WCH/st_a.rda')
tumor <- subset(st_a, subset = spot_type == "tumor") 
SpatialFeaturePlot(st_a, features =c("CD274"), pt.size.factor = 6)
pdf('D:/项目/olink/result/MR analysis/空转/WCH/scatter_plot.pdf')
a <- SpatialDimPlot(tumor, group.by = "spot_type", pt.size.factor = 5,cols = c('tumor'='#1DB2B7'))
b <- SpatialFeaturePlot(st_a, features =c("CD274"), pt.size.factor = 6,alpha=c(0,1))
x <- a+b
x
dev.off()

