#direct
library(TwoSampleMR)
rt <- read.csv("06_mydata_CSS.csv",header=T)
colnames(rt)[27]<-'samplesize.exposure'
rt$r.exposure <-  rt$R2^0.5
rt$r.outcome <- (2*rt$eaf.outcome*(1-rt$eaf.outcome)*rt$beta.outcome^2)^0.5
rownames(rt) <- NULL
direct <- data.frame()
for(i in 1:nrow(rt)){
  dir <- rt[i,]
  out <- directionality_test(dir)
  direct <- rbind(direct,out)
}
direct <- cbind(rt$SNP,rt$Protein,direct)
colnames(direct)[1] <- 'SNP'
colnames(direct)[2] <- 'Protein'
write.table(direct,file = "07_direct_CSS.csv",col.names = T,row.names = F,sep=",",quote=FALSE)

#coloc
library(coloc)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(stringr)
library(LDlinkR)
library(tidyverse)
library(readxl)
library(locuscomparer)
library(Rmpfr)
proteins <- read_xlsx('D:/项目/olink/result/MR analysis/共定位/pQTLs.xlsx')
allgwas = read.table("D:/项目/olink/result/MR analysis/GWAS/GWAS_summer/UKBsummer_all.txt",header=T,sep=" ")
allgwas$BP_hg19 <- as.numeric(allgwas$BP_hg19)
allgwas <- allgwas[,-c(11:13)]
allgwas$beta.outcome <- log(allgwas$HR_OS)
colnames(allgwas)[4] <- 'pos'
colnames(allgwas)[10] <- 'pval.outcome'
colnames(allgwas)[9] <- 'se.outcome' 
allgwas$id.outcome <- 'UKB colorectal Cancer'
allgwas$outcome <- 'colorectal cancer survival'
colnames(allgwas)[5] <- 'effect_allele.outcome'
colnames(allgwas)[6] <- 'other_allele.outcome' 
colnames(allgwas)[7] <- 'eaf.outcome'

pathS <- "D:/项目/olink/result/MR analysis/pQTL_summary/cis_coloc"
fileNamesS <- dir(pathS) 
fileNamesS <- fileNamesS[fileNamesS %in% paste(proteins$Protein,'.regenie',sep='')]
filePathS <- sapply(fileNamesS, function(x){ 
  paste(pathS,x,sep='/')})
dataS <- lapply(filePathS, function(x){
  read.table(x,header=T,sep='')})  
names(dataS) <- gsub(".regenie", "", names(dataS))
ngwas = 2621 
npqtl = 50636 

pathpdf <- "D:/项目/olink/result/MR analysis/共定位"
snps <- paste(proteins$Protein,'_',proteins$SNP,'.pdf')
filePathpdf <- sapply(snps, function(x){ 
  paste(pathpdf,x,sep='/')})
result <- data.frame()
for (i in 8:10){
  co <- proteins[i,]
  gwas =allgwas[which(allgwas$CHR==co$Chr),]
  chd = gwas[which((co$BP_hg19-500000<=gwas$pos) & (gwas$pos<=co$BP_hg19+500000)),] 
  pqtl = dataS[[co$Protein]]
  a <- str_split_fixed(pqtl$ID, ":", n = 6)
  colnames(pqtl)[c(4:6,10,11)] <- c('other_allele.exposure','effect_allele.exposure','eaf.exposure','beta.exposure','se.exposure')
  pqtl$pval.exposure <- 10^(-pqtl$LOG10P)
  pqtl$pos <- a[,2]
  pqtl$id.exposure <- 'pQTL'
  pqtl$exposure <- 'pQTL meta'
  pqtl$pval.exposure[pqtl$pval.exposure == 0] <- 1E-308
  input <- merge(pqtl,chd,by='pos')
  exp <- input[,c(1:21)]
  mydata <- harmonise_data(
    exposure_dat=exp,
    outcome_dat=chd,
    action= 1
  )  
  mydata <- mydata[which(mydata$mr_keep==TRUE),]
  mydata$MAF <- ifelse(mydata$eaf.exposure>0.5, 1-mydata$eaf.exposure, mydata$eaf.exposure)
  mydata = mydata[!duplicated(mydata$SNP),]
  re = coloc.abf(dataset1 = list(beta=mydata$beta.exposure, varbeta=mydata$se.exposure^2, snp=mydata$SNP, type="quant" ,N=npqtl), 
                 dataset2 = list(beta=mydata$beta.outcome,varbeta=mydata$se.outcome^2,snp=mydata$SNP, type="quant" ,N=ngwas),
                 MAF = mydata$MAF,p1 = 1e-04, p2 = 1e-04, p12 =1e-05)
  result1 <- as.data.frame(re[["summary"]])
  result1 <- as.data.frame(t(result1))
  result <- rbind(cbind(result1,co[,c(1:4)]),result)
  write.table(mydata[,c(1:3,6,8,15,16,28,32,41)],"D:/项目/olink/result/MR analysis/共定位/pqtl_fn.txt",sep="\t",row.names=F,quote=F)
  write.table(mydata[,c(1,4,5,7,9,15,16,18,19,41)],"D:/项目/olink/result/MR analysis/共定位/gwas_fn.txt",sep="\t",row.names=F,quote=F)
  marker_col="SNP"
  pqtl_fn="D:/项目/olink/result/MR analysis/共定位/pqtl_fn.txt"
  gwas_fn="D:/项目/olink/result/MR analysis/共定位/gwas_fn.txt"
  pdf(file=filePathpdf[i], width=10, height=6)
  p1 <- locuscompare(in_fn1=pqtl_fn, in_fn2=gwas_fn, 
                     title1="pQTL", title2="GWAS", marker_col1= marker_col, 
                     pval_col1='pval.exposure', marker_col2=marker_col, 
                     pval_col2='pval.outcome',genome='hg19')
  print(p1)
  dev.off()}

#MRlap
library(archive)
library(purrr)
library(stringr)
path <- 'D:/项目/olink/result/Car修改/04_MRlap_ld_reference/MUC13'
MUC13 <- list.files(path, pattern = "\\.regenie$", full.names = TRUE) %>%
  map_dfr(~ {
    data <- read.table(.x, header = TRUE, sep = "")
    data$filename <- basename(.x)
    data$protein <- gsub("\\.regenie$", "", basename(.x))
    return(data)
  })
alldata <- read.table("D:/项目/olink/result/MR analysis/GWAS/GWAS_summer/UKBsummer_all.txt",header=T,sep=" ")
a <- str_split_fixed(MUC13$ID, ":", n = 6)
MUC13$pos <- a[,2]
CRC <- alldata
colnames(MUC13)[c(1,4,5,6,10,11)] <- c('chr','ref','alt','eaf','beta','se')
colnames(CRC)[c(2,4,5,6)] <- c('chr','pos','alt','ref')
pro <- MUC13[MUC13$chr %in% CRC$chr & MUC13$pos %in% CRC$pos,]
pro <- merge(CRC[c(2:4)],pro,by=c('pos','chr'))
CRC$beta <- log(CRC$HR_CSS)
CRC$se <- CRC$SE_CSS
CRC$N <- 2621

library(MRlapPro)
mydata <- mydata[which(mydata$Protein=='MUC13'),]
exp <- mydata[,c("SNP","CHR","BP_hg19.x","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure",
                 "pval.exposure","N")]
chd <- mydata[,c("SNP","CHR","BP_hg19.x","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome","se.outcome",
                 "pval.outcome","samplesize.outcome")]
colnames(exp) <- c('snp','chr','pos','alt','ref','eaf','beta','se','pval','N')
colnames(chd) <- c('snp','chr','pos','alt','ref','eaf','beta','se','pval','N')
result <- MRlap::MRlap(exposure = pro, 
                       outcome = CRC,
                       do_pruning=F,
                       user_SNPsToKeep=mydata$SNP,
                       ld = system.file("Data/eur_w_ld_chr", package="MRlapPro"),
                       hm3 = system.file("Data/eur_w_ld_chr", "w_hm3.snplist", package="MRlapPro"))
save(result,file='04_malap_CSS.Rdata')
lap <- as.data.frame(result[["MRcorrection"]])
lap$corrected_HR <- exp(lap$corrected_effect)
lap$corrected_HR_LCI<- exp(log(lap$corrected_HR) - lap$corrected_effect_se * 1.96)
lap$corrected_HR_UCI <- exp(log(lap$corrected_HR) + lap$corrected_effect_se * 1.96)

save(lap,file='07_mrlap_OS.Rdata')


