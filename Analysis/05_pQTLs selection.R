#pQTLs selection
library(TwoSampleMR)
library(stringr)
library(readxl)
loc <- as.data.frame(read_xlsx('D:/项目/olink/result/MR analysis/lasso_gene_location_BP19.xlsx'))
reference.1000G.maf<-read.table("D:/项目/olink/result/MR analysis/共定位参考/reference.1000G.txt",header=T,sep=' ') 
colnames(reference.1000G.maf)[3] <- 'BP_hg19'
proteins <- read.table('D:/项目/olink/result/MR analysis/lasso proteins.txt',sep='\t',header=T)
pathS <- "D:/项目/olink/result/MR analysis/pQTL_summary/S/cis_coloc"
fileNamesS <- dir(pathS) 
fileNamesS <- fileNamesS[fileNamesS %in% paste(proteins$Proteins,'.regenie',sep='')]
filePathS <- sapply(fileNamesS, function(x){ 
  paste(pathS,x,sep='/')})
dataS <- lapply(filePathS, function(x){
  read.table(x,header=T,sep='')})  
names(dataS) <- gsub(".regenie", "", names(dataS))

#01
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
pQTL <- data.frame()
for (i in 1:length(fileNamesS)){
  pqtl <- dataS[[i]]
  Assay <- names(dataS)[i]
  a <- str_split_fixed(pqtl$ID, ":", n = 6)
  pqtl$BP_hg19 <- as.numeric(a[,2])
  pqtl$pval.exposure <- 1.1^-pqtl$LOG10P
  colnames(pqtl)[1] <- 'CHR'
  if (pqtl$CHR[1] == 6){
    pqtl <- pqtl[-which(pqtl$BP_hg19 > 25500000 & pqtl$BP_hg19 < 34000000),]
  }
  pqtl <- pqtl[which(pqtl$LOG10P>(-log10(5E-08))),]
  pqtl <- pqtl[pqtl$A1FREQ>0.01,]
  if (!(Assay =='FUT3_FUT5')){
    pqtl <- pqtl[which(pqtl$BP_hg19>loc[loc$Assay==Assay,'gene_start']-1000000 & pqtl$BP_hg19<loc[loc$Assay==Assay,'gene_end']+1000000),]}
  if (Assay=='FUT3_FUT5'){
    pqtl <- pqtl[which(pqtl$GENPOS>5842888-1000000 & pqtl$GENPOS<5870540+1000000),]}
  ref <- reference.1000G.maf[which(reference.1000G.maf$CHR %in% pqtl$CHR),]
  ref <- ref[,-2]
  pqtl <- merge(ref,pqtl,by='BP_hg19')
  if (nrow(pqtl)>0){
    input <- clump_data(pqtl, clump_kb=10000, clump_r2=0.001)
    input$Assay <- Assay
    pQTL <- rbind(pQTL,input)}}
pQTL$pval.exposure <- 10^-pQTL$LOG10P

#02
pathP <- "D:/项目/olink/result/MR analysis/pQTL_summary/P"
fileNamesP <- dir(pathP) 
fileNamesP <- fileNamesP[fileNamesP %in% proteins$Proteins]
filePathP <- sapply(fileNamesP, function(x){ 
  paste(pathP,x,sep='/')})
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

pQTL_P_5e8_0.001 <- data.frame()
pQTL_P_1e10_0.001 <- data.frame()
pQTL_P_1e10_0.01 <- data.frame()
pQTL_P_5e8_0.01 <- data.frame()
for (x in 1:length(filePathP)){
  dataP <- read.table(filePathP[x],header=T,sep='')
  Assay <- fileNamesP[x]
  dataSm <- dataS[[Assay]]
  colnames(dataP)[c(11,18,19)] <- c('pval.exposure','CHR','BP_hg19')
  dataP <- dataP[dataP$CHR %in% dataSm$CHROM,]
  ref <- reference.1000G.maf[which(reference.1000G.maf$CHR %in% dataSm$CHROM),]
  ref <- ref[,-2]
  dataP <- merge(ref,dataP,by='BP_hg19')
  if (dataP$CHR[1]==6){
    dataP <- dataP[-which(dataP$CHR==6 & dataP$BP_hg19 > 25500000 & dataP$BP_hg19 < 34000000),]}
  dataP <- dataP[dataP$Freq1>0.01,]
  dataP <- dataP[which(dataP$BP_hg19>loc[loc$Assay==Assay,'gene_start']-1000000 & dataP$BP_hg19<loc[loc$Assay==Assay,'gene_end']+1000000),]
  pqtl5e8 <- dataP[which(dataP$pval.exposure<(5E-08)),]
  pqtl1e10 <- dataP[which(dataP$pval.exposure<(1E-10)),]
  if (nrow(pqtl5e8)>0){
    input_5e8_0.001 <- clump_data(pqtl5e8, clump_kb=10000, clump_r2=0.001)
    input_5e8_0.001$Assay <- Assay
    input_5e8_0.01 <- clump_data(pqtl5e8, clump_kb=10000, clump_r2=0.01)
    input_5e8_0.01$Assay <- Assay
    pQTL_P_5e8_0.01 <- rbind(pQTL_P_5e8_0.01,input_5e8_0.01)
    pQTL_P_5e8_0.001 <- rbind(pQTL_P_5e8_0.001,input_5e8_0.001)}
  if (nrow(pqtl1e10)>0){
    input_1e10_0.001 <- clump_data(pqtl1e10, clump_kb=10000, clump_r2=0.001)
    input_1e10_0.001$Assay <- Assay
    input_1e10_0.01 <- clump_data(pqtl1e10, clump_kb=10000, clump_r2=0.01)
    input_1e10_0.01$Assay <- Assay
    pQTL_P_1e10_0.01 <- rbind(pQTL_P_1e10_0.01,input_1e10_0.01)
    pQTL_P_1e10_0.001 <- rbind(pQTL_P_1e10_0.001,input_1e10_0.001)}}


#03
pathG <- "D:/项目/olink/result/MR analysis/pQTL_summary/G"
fileNamesG <- dir(pathG) 
filePathG <- sapply(fileNamesG, function(x){ 
  paste(pathG,x,sep='/')})
dataG <- lapply(filePathG, function(x){
  read.table(x,header=T,sep='')})  
names(dataS) <- gsub(".regenie", "", names(dataS))
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

pQTL <- data.frame()
for (i in 1:length(fileNamesS)){
  pqtl <- dataS[[i]]
  Assay <- names(dataS)[i]
  a <- str_split_fixed(pqtl$ID, ":", n = 6)
  pqtl$BP_hg19 <- as.numeric(a[,2])
  pqtl$pval.exposure <- 1.1^-pqtl$LOG10P
  colnames(pqtl)[1] <- 'CHR'
  if (pqtl$CHR[1] == 6){
    pqtl <- pqtl[-which(pqtl$BP_hg19 > 25500000 & pqtl$BP_hg19 < 34000000),]
  }
  pqtl <- pqtl[which(pqtl$LOG10P>(-log10(5E-08))),]
  pqtl <- pqtl[pqtl$A1FREQ>0.01,]
  if (!(Assay =='FUT3_FUT5')){
    pqtl <- pqtl[which(pqtl$BP_hg19>loc[loc$Assay==Assay,'gene_start']-1000000 & pqtl$BP_hg19<loc[loc$Assay==Assay,'gene_end']+1000000),]}
  if (Assay=='FUT3_FUT5'){
    pqtl <- pqtl[which(pqtl$GENPOS>5842888-1000000 & pqtl$GENPOS<5870540+1000000),]}
  ref <- reference.1000G.maf[which(reference.1000G.maf$CHR %in% pqtl$CHR),]
  ref <- ref[,-2]
  pqtl <- merge(ref,pqtl,by='BP_hg19')
  if (nrow(pqtl)>0){
    input <- clump_data(pqtl, clump_kb=10000, clump_r2=0.001)
    input$Assay <- Assay
    pQTL <- rbind(pQTL,input)}}
pQTL$pval.exposure <- 10^-pQTL$LOG10P


#04
pathS1 <- "D:/项目/olink/result/MR analysis/pQTL_summary/S1"
fileNamesS1 <- dir(pathS1) 
fileNamesS1 <- fileNamesS1[fileNamesS1 %in% paste(proteins$Proteins,'.tsv',sep='')]
filePathS1 <- sapply(fileNamesS1, function(x){ 
  paste(pathS1,x,sep='/')})
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
pQTL_Sun1 <- data.frame()
for (x in 1:length(filePathP)){
  dataS1 <- read.table(filePathS1[x],header=T,sep='')
  Assay <- gsub(".tsv", "", fileNamesS1[x])
  dataSm <- dataS[[Assay]]
  dataS1$pval.exposure <- 1.1^dataS1$log.P.
  colnames(dataS1)[c(2,3)] <- c('CHR','BP_hg19')
  dataS1 <- dataS1[dataS1$CHR %in% dataSm$CHROM,]
  ref <- reference.1000G.maf[which(reference.1000G.maf$CHR %in% dataSm$CHROM),]
  ref <- ref[,-2]
  dataS1 <- merge(ref,dataS1,by='BP_hg19')
  if (dataS1$CHR[1]==6){
    dataS1 <- dataS1[-which(dataS1$CHR==6 & dataS1$BP_hg19 > 25500000 & dataS1$BP_hg19 < 34000000),]}
  dataS1 <- dataS1[dataS1$MAF>0.01,]
  dataS1 <- dataS1[which(dataS1$BP_hg19>loc[loc$Assay==Assay,'gene_start']-1000000 & dataS1$BP_hg19<loc[loc$Assay==Assay,'gene_end']+1000000),]
  dataS1 <- dataS1[which(dataS1$p_value<(5E-08)),]
  if (nrow(dataS1)>0){
    input <- clump_data(dataS1, clump_kb=10000, clump_r2=0.001)
    input$Assay <- Assay
    pQTL_Sun1 <- rbind(pQTL_Sun1,input)}}
pQTL_Sun1 <- pQTL_Sun1[,-14]
colnames(pQTL_Sun1)[13] <- 'pval.exposure'


#05
pathK <- "D:/项目/olink/result/MR analysis/pQTL_summary/K"
fileNamesK <- dir(pathK) 
fileNamesK <- fileNamesK[fileNamesK %in% paste(proteins$Proteins,'.tsv',sep='')]
filePathK <- sapply(fileNamesK, function(x){ 
  paste(pathK,x,sep='/')})
dataK <- lapply(filePathK, function(x){
  read.table(x,header=T,sep='')})  
names(dataK) <- gsub(".tsv", "", names(dataK))
pQTL_K <- data.frame()
for (i in 1:length(fileNamesK)){
  pqtl <- dataK[[i]]
  Assay <- names(dataK)[i]
  dataSm <- dataS[[Assay]]
  colnames(pqtl)[c(1,2,8)] <- c('CHR','BP_hg19','pval.exposure')
  pqtl <- pqtl[pqtl$CHR %in% dataSm$CHROM,]
  ref <- reference.1000G.maf[which(reference.1000G.maf$CHR %in% dataSm$CHROM),]
  ref <- ref[,-2]
  pqtl <- merge(ref,pqtl,by='BP_hg19')
  if (pqtl$CHR[1]==6){
    pqtl <- pqtl[-which(pqtl$CHR==6 & pqtl$BP_hg19 > 25500000 & pqtl$BP_hg19 < 34000000),]}
  pqtl <- pqtl[pqtl$MAF>0.01,]
  pqtl <- pqtl[which(pqtl$BP_hg19>loc[loc$Assay==Assay,'gene_start']-1000000 & pqtl$BP_hg19<loc[loc$Assay==Assay,'gene_end']+1000000),]
  pqtl <- pqtl[which(pqtl$pval.exposure<(5E-08)),]
  if (nrow(pqtl)>0){
    input <- clump_data(pqtl, clump_kb=10000, clump_r2=0.001)
    input$Assay <- Assay
    pQTL_K <- rbind(pQTL_K,input)}}


#integrate
pQTL_Sun <- pQTL[,c(1,2,6,9:11,15,16,20,22)]
pQTL_Sun <- pQTL_Sun[,c(10,2,3,1,5,4,6,7:9)]
colnames(pQTL_Sun) <- c('Protein','SNP','Chr','BP_hg19','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure','se.exposure','pval.exposure') 
pQTL_Sun$Study <- 'Sun'
pQTL_Sun$N <- 34557
pQTL_P_5e8_0.001$Study <- 'Pietzner'
pQTL_P_5e8_0.001 <- pQTL_P_5e8_0.001[,c(25,2,23,1,8:10,14:16,26,22)]
colnames(pQTL_P_5e8_0.001) <- c('Protein','SNP','Chr','BP_hg19','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure','se.exposure','pval.exposure','Study','N') 
pQTL_Sun1 <- pQTL_Sun1[,c(15,2,7,1,8,9,3,10,11,13)] #EAF为MAF
colnames(pQTL_Sun1) <- c('Protein','SNP','Chr','BP_hg19','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure','se.exposure','pval.exposure') 
pQTL_Sun1$Study <- 'Sun1'
pQTL_Sun1$N <- 3301
pQTL_K <- pQTL_K[,c(14,2,6,1,7,8,11,9,10,12)]
colnames(pQTL_K) <- c('Protein','SNP','Chr','BP_hg19','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure','se.exposure','pval.exposure') 
pQTL_K$Study <- 'Kalnapenkis'
pQTL_K$N <- 500

pQTL <- rbind(pQTL_Sun,pQTL_P_5e8_0.001,pQTL_Sun1,pQTL_K,pQTL_R)
pQTL$Protein <- toupper(pQTL$Protein)
pQTL$Protein <- gsub('_','',pQTL$Protein)
pQTL$Protein <- gsub('-','',pQTL$Protein)
pQTL$Protein <- gsub(' ','',pQTL$Protein)
pQTL <- pQTL[pQTL$Protein %in% proteins$Proteins,]
pQTL$MAF <- NA
for (i in 1:nrow(pQTL)){
  if (pQTL[i,'eaf.exposure']>0.5){
    pQTL[i,'MAF'] <- 1-pQTL[i,'eaf.exposure']}
  if (pQTL[i,'eaf.exposure']<0.5){
    pQTL[i,'MAF'] <- pQTL[i,'eaf.exposure']}
}
pQTL$R2 <- 2*pQTL$eaf.exposure*(1-pQTL$eaf.exposure)*pQTL$beta.exposure^2
pQTL$F <- pQTL$beta.exposure^2/pQTL$se.exposure^2
pQTL[pQTL$Study=='Sun',11] <- 'Sun2'
write.table(pQTL,file='05_pQTLs.csv',sep=",",row.names=F,quote=F)
