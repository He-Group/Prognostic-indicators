library(dplyr)
library(readxl)
#WCH
##data cleaning
clinical <- read_xlsx("D://项目//olink//150例olink患者临床信息.xlsx")
olink <- read.csv("D://项目//olink//华西医院_NPX.csv",header=T)
olink <- olink[which(olink$SampleID %in% clinical$SampleID),]
length(unique(olink$SampleID))


olink_1 <- olink[which(olink$LOD > olink$NPX),]
pro <- as.data.frame(table(olink_1$Assay))
depro <- pro[which(pro$Freq>150*0.5),] 
olinkdata <- olinkdata[-which(olinkdata$Assay %in% depro$Var1),]
write.table(depro, file="D:/项目/olink/result/00_depro_0.5.csv",col.names = T,row.names = F,sep=",",quote=FALSE)

##clinical data
clinicaldata <- clinical %>% 
  mutate(Location = case_when(Location=='colon'~0,Location=='rectum'~1),
         sex = case_when(sex=='male'~0,sex=='female'~1),
         T = case_when(T=='T2'~2,T=='T3'~3,T=='T4'~4),
         N = case_when(N=='N0'~0,N=='N1'~1,N=='N2'~2),
         M = case_when(M=='M0'~0,M=='M1'~1),
         Stage = case_when(Stage=='I'~1,Stage=='II'~2,Stage=='III'~3,Stage=='IV'~4))

olinkdata <- olinkdata[,c(1,5,12)]
a <- unique(olinkdata[,1])
b <- unique(olinkdata[,2])
data= matrix(nrow = 150, ncol=339, byrow=TRUE)
rownames(data) <- a
colnames(data) <- b
npx <- as.data.frame(data)
for(i in 1:339){
  protein <- olinkdata[which(olinkdata$Assay==b[i]),]
  npx[,i] <- protein[,3]
}
npx <- tibble::rownames_to_column(npx,"SampleID")

alldata <- full_join(clinicaldata, npx, by = "SampleID")
alldata <- as.data.frame(alldata)
alldata$DFS.time[which(alldata$DFS.time == 'PreMeta')] <- NA
write.table(alldata, file="D:/项目/olink/result/00_alldata_ori.csv",col.names = T,row.names = F,sep=",",quote=FALSE)

