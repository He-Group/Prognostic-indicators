library(survival)
alldata <- read.csv('D:/项目/olink/result/alldata_ori.csv',header=T,sep=',')
depro <- read.csv('D:/项目/olink/result/depro_0.5.csv',header=T)
a <- alldata[,c(16:380)]
a <- a[,-which(colnames(a) %in% depro$Var1)]
alldata <- cbind(alldata[,1:15],a)
#OS
cox_OS <- data.frame()
for(factor in colnames(alldata[,16:ncol(alldata)])){
  cox=coxph(Surv(alldata$OS.time,alldata$OS.event) ~ alldata[,factor]+age+sex+Stage+Location+chemoradiotherapy, data = alldata) 
  coxSummary = summary(cox)                         
  coxP=coxSummary$coefficients["alldata[, factor]","Pr(>|z|)"]
  cox_OS=rbind(cox_OS,
               cbind(Assay=factor,
                     coef = coxSummary$coefficients["alldata[, factor]","coef"],
                     HR=coxSummary$conf.int["alldata[, factor]","exp(coef)"],
                     HR.95L=coxSummary$conf.int["alldata[, factor]","lower .95"],
                     HR.95H=coxSummary$conf.int["alldata[, factor]","upper .95"],
                     pvalue=coxP) )}
cox_OS[,c(2:6)] <- as.data.frame(lapply(cox_OS[,c(2:6)],as.numeric))
FDR <- p.adjust(cox_OS$pvalue,method = 'fdr')
cox_OS$FDR <- p.adjust(cox_OS$pvalue,method = 'fdr')
cox_OS$outcome <- 'CRC OS'
#DFS
alldata_DFS <- alldata[-which(alldata$Stage==4),]
cox_DFS <- data.frame()
for(factor in colnames(alldata_DFS[,16:ncol(alldata_DFS)])){
  cox=coxph(Surv(alldata_DFS$DFS.time,alldata_DFS$DFS.event) ~ alldata_DFS[,factor]+age+sex+Stage+Location+chemoradiotherapy, data = alldata_DFS) 
  coxSummary = summary(cox)                         
  coxP=coxSummary$coefficients["alldata_DFS[, factor]","Pr(>|z|)"]
  cox_DFS=rbind(cox_DFS,
                cbind(Assay=factor,
                      coef = coxSummary$coefficients["alldata_DFS[, factor]","coef"],
                      HR=coxSummary$conf.int["alldata_DFS[, factor]","exp(coef)"],
                      HR.95L=coxSummary$conf.int["alldata_DFS[, factor]","lower .95"],
                      HR.95H=coxSummary$conf.int["alldata_DFS[, factor]","upper .95"],
                      pvalue=coxP) )}
cox_DFS[,c(2:6)] <- as.data.frame(lapply(cox_DFS[,c(2:6)],as.numeric))
FDR <- p.adjust(cox_DFS$pvalue,method = 'fdr')
cox_DFS$FDR <- p.adjust(cox_DFS$pvalue,method = 'fdr')
cox_DFS$outcome <- 'CRC DFS'

across <- rbind(cox_OS,cox_DFS)
across$FDR_across <- p.adjust(across$pvalue,method = 'fdr')
OS <- cox_OS[which(cox_OS$Assay %in% c('MUC13','TNFRSF10B','STC1','ADAM8','FUT3_FUT5','DRAXIN','PTPRN2','SLIT2')),]
DFS <- cox_DFS[which(cox_DFS$Assay %in% c('CD109','BST1','FCER2','CD274','GGT1','KEL','CA6','CLPS','THBS2','PTS','VSIG4')),]
OS$FDR <- p.adjust(OS$pvalue,method='fdr')
DFS$FDR <- p.adjust(DFS$pvalue,method='fdr')
write.table(across, file="04_cox_result.csv",col.names = T,row.names = F,sep=",",quote=FALSE)

