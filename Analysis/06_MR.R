#LDlink
pQTL <- read.table('D:/项目/olink/result/MR analysis/最终结果/05_pQTL_MR.csv',sep=',',header=T)
allgwas = read.table("D:/项目/olink/result/MR analysis/GWAS/GWAS_summer/UKBsummer_all.txt",header=T,sep=" ")
survival <- allgwas[allgwas$SNP %in% pQTL$SNP,]
redu <- c(setdiff(pQTL$SNP,allGWAS$SNP),'rs11260014')
library(LDlinkR)
LDlink <- data.frame() 
for (i in 1:length(redu)) {
  LD <- LDproxy(snp = redu[i],pop = 'EUR',r2d = "r2", token = '94f643a397e0') 
  LD <- LD[which(LD$R2>0.8),]
  LD$SNP <- redu[i]
  LDlink <- rbind(LDlink,LD)}
LDlink <- LDlink[LDlink$RS_Number %in% allgwas$SNP,] 

library('TwoSampleMR')
exp <- read.table('D:/项目/olink/result/MR analysis/最终结果/05_pQTL_MR.csv',header=T,sep=',')
chd <- allgwas
chd <- chd[,-c(8:10)]
chd$beta.outcome <- log(chd$HR_CSS)
colnames(chd)[10] <- 'pval.outcome'
colnames(chd)[9] <- 'se.outcome' 
chd$id.outcome <- 'UKB colorectal Cancer'
chd$outcome <- 'colorectal cancer survival'
colnames(chd)[5] <- 'effect_allele.outcome'
colnames(chd)[6] <- 'other_allele.outcome' 
colnames(chd)[7] <- 'eaf.outcome'  
exp$id.exposure <- 'pQTL'
exp$exposure <- 'pQTL meta'
mydata <- harmonise_data(
  exposure_dat=exp,
  outcome_dat=chd,
  action= 2
)  
mydata$samplesize.outcome <- 2621
IVW <- mydata[duplicated(mydata$Protein) == T,] 
IVW <- IVW[,21]
IVW <- mydata[which(mydata$Protein %in% IVW),] 
IVW <- IVW[order(IVW$Protein),]
rownames(IVW) <- NULL
wald <- mydata[(which(!(mydata$Protein %in% IVW$Protein))),] 
rownames(wald) <- NULL
re_wald <- data.frame()
for(i in rownames(wald[1:nrow(wald),])){
  rati <- mr_wald_ratio(wald[i,6],wald[i,7],wald[i,24],wald[i,18])
  rati <- data.frame(rati)
  re_wald <- rbind(re_wald,rati)
}
re_wald <- cbind(wald$SNP,wald$Protein,re_wald)
colnames(re_wald)[1] <- 'SNP'
colnames(re_wald)[2] <- 'Protein'
re_IVW <- data.frame()
for(i in rownames(IVW[1:nrow(IVW),])){
  factor <- IVW[i,21]
  a <- IVW[which(IVW$Protein == factor),]
  rati <- mr_ivw(a[,6],a[,7],a[,24],a[,18])
  rati <- data.frame(rati)
  re_IVW <- rbind(re_IVW,rati)}  
re_IVW <- cbind(IVW$Protein,re_IVW)
re_IVW <- unique(re_IVW)
re_IVW<- re_IVW[!duplicated(re_IVW$`IVW$Protein`),]
colnames(re_IVW)[1] <- 'Protein'
re_wald <- re_wald[,-1]
re_IVW <- re_IVW[,-c(6,7,8)]
re <- rbind(re_wald,re_IVW)
re$HR <- exp(re$b)
re$lower_CI <- exp(log(re$HR) - re$se * 1.96)
re$upper_CI <- exp(log(re$HR) + re$se * 1.96)

a <- mydata[which(mydata$Protein=='FCER2'),]
x <- mr(a)
x$HR <- exp(x$b)
x$lower_CI <- exp(log(x$HR) - x$se * 1.96)
x$upper_CI <- exp(log(x$HR) + x$se * 1.96)
x <- x[-3,]
write.table(re,'06_re_CSS.csv',col.names = T,row.names = F,sep=",",quote=FALSE)
write.table(mydata,'06_mydata_CSS.csv',col.names = T,row.names = F,sep=",",quote=FALSE)
#OS like it