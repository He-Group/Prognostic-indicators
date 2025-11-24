library(survival)
library(glmnet)
library(survminer)
library(timeROC)
library(boot) 
library(nestedcv)
#01 nested CV lasso
alldata <- read.csv('D:/项目/olink/result/00_alldata_ori.csv',header=T)
depro <- read.csv('D:/项目/olink/result/00_depro_0.5.csv',header=T)
a <- alldata[,c(16:380)]
a <- a[,-which(colnames(a) %in% depro$Var1)]
alldata <- cbind(alldata[,1:15],a)
rt <- alldata[,c(12,13,16:ncol(alldata))]
x <- as.matrix(alldata[,c(16:ncol(alldata))])       
y <- Surv(alldata$OS.time,alldata$OS.event)

library(nestedcv)
resulted <- nestcv.glmnet(
  x = x, 
  y = y,
  family = "cox",                  
  n_outer_folds = 5,
  n_inner_folds = 5,         
  alpha = 1,
  min_1se = 0,
  outer_result = T)
summary(resulted)
Coef_OS <- resulted[["final_coef"]]

lassoPro <- rownames(Coef_OS)
rt <- alldata[,c('OS.time','OS.event',lassoPro)]
score <- rt[,rownames(Coef_OS)]
for(i in (1:nrow(Coef_OS))){       
  score[,i] <- (score[,i]*Coef_OS$coef[i])
}
trainScore = rowSums(score)
rt <- cbind(SampleID=(alldata$SampleID),alldata[,c(2:11)],rt,Riskscore=as.vector(trainScore))
shapiro.test(rt$Riskscore)
rt$risk <- as.vector(ifelse(rt$Riskscore > median(rt$Riskscore),  "high", "low"))
table(rt$risk)
rt_OS <- rt
write.csv(rt_OS,file= '01_WCH_rt_OS.csv',row.names = F,quote = F)
write.csv(Coef_OS,file= '01_Coef_LASSO_OS.csv',row.names = T,quote=FALSE)

#nomo
library(foreign)
library(rms)
library(regplot)
d <- rt[,c('age','sex','Location','Stage','chemoradiotherapy','Riskscore','OS.time','OS.event')]
d$sex <- factor(d$sex,levels = c(0,1),
                labels = c('male','female'))
d$Location <- factor(d$Location,levels = c(0,1),
                     labels = c('colon','rectum'))
d$chemoradiotherapy <- factor(d$chemoradiotherapy,levels = c(0,1),
                              labels = c('chemotherapy','Non-chemotherapy'))
dd<-datadist(d) 
options(datadist = 'dd') 
coxm <- cph(Surv(OS.time,OS.event == 1) ~ age+sex+Location+Stage+chemoradiotherapy+Riskscore,data = d,
            x = T,y = T, surv = T)
surv<- Survival(coxm)
surv1<- function(x) surv(12,x) 
surv2<- function(x) surv(36,x) 
surv3<- function(x) surv(60,x) 
nom1<-nomogram(coxm,fun = list(surv1,surv2,surv3),lp = F,
               funlabel = c("1-Yeas overall survival probability", '3-Year overall survival probability','5-Year overall survival probability'), maxscale = 100, 
               fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1'))
plot(nom1)
pdf(file="01_WCH_OS_nomo.pdf", width=20, height=10)
plot((nom1),xfrac=.3)
dev.off()

#KM curve
library(pheatmap)
library(tibble)
library(stringr)
library(survival)
library(dplyr)
library(survminer)
rt1 <- rt
fit <- survfit(Surv(OS.time, OS.event)~risk, data=rt1)
pdf(file="01_WCH_OS_KM.pdf", width=9, height=9)
KM <- ggsurvplot(
  fit,
  data = rt1,
  size = 1, 
  linetype = 1,
  palette =
    c("#E7B800", "#2E9FDF"),
  risk.table = TRUE,        
  risk.table.col = "strata",
  legend.labs =
    c("Risk:High", "Risk:Low"),   
  pval = T,
  risk.table.height = 0.25, 
  ggtheme = theme_bw(),
  fontsize = 7,
  fontsize_row=8,
  fontsize_col=3,
  conf.int=F,
  xlab = "Time(Months)",
  censor=F)
print(KM)
dev.off()

#ROC
library(survival)
library(survminer)
library(timeROC)
rt <- read.csv('01_WCH_rt_OS.csv')
fit <- coxph(Surv(OS.time,OS.event)~age+sex+Location+Stage+chemoradiotherapy+Riskscore,data=rt)
Coef_fit <- as.data.frame(fit$coefficients)
score <- rt[,c('age','sex','Location','Stage','chemoradiotherapy','Riskscore')]
score <- score[,rownames(Coef_fit)]
for(i in (1:nrow(Coef_fit))){       
  score[,i] <- score[,i]*Coef_fit[i,1]}
modelScore = rowSums(score)
rt <- cbind(rt,modelScore=as.vector(modelScore))
fit_cli <- coxph(Surv(OS.time,OS.event)~age+sex+Location+Stage+chemoradiotherapy,data=rt)
Coef_cli <- as.data.frame(fit_cli$coefficients)
score_cli <- rt[,c('age','sex','Location','Stage','chemoradiotherapy')]
score_cli <- score_cli[,rownames(Coef_cli)]
for(i in (1:nrow(Coef_cli))){       
  score_cli[,i] <- score_cli[,i]*Coef_cli[i,1]}
Score_cli = rowSums(score_cli)
rt <- cbind(rt,Score_cli=as.vector(Score_cli))
ROC_rt=timeROC(T=rt$OS.time,delta=rt$OS.event,
               marker=rt$Score_cli,cause=1,
               weighting = "marginal",
               times=c(12,36,60),ROC=TRUE,iid=T)
CI <- as.data.frame(confint(ROC_rt, level = 0.95)$CI_AUC)
CI$AUC <- paste0(sprintf("%.03f",ROC_rt$AUC)," (",sprintf("%.03f",CI$`2.5%`*0.01),"-",sprintf("%.03f",CI$`97.5%`*0.01),")")
pdf(file="01_WCH_OS_ROC_cli.pdf", width=5, height=5)
plot(ROC_rt,time=12,col="#2874C5",title=FALSE,lwd=2)
plot(ROC_rt,time=36,col="#f87669",add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=60,col="#e6b707",add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',CI$AUC[1]),
         paste0('AUC at 3 years: ',CI$AUC[2]),
         paste0('AUC at 5 years: ',CI$AUC[3])),
       col=c("#2874C5", "#f87669", "#e6b707"), lwd=3,cex=0.75,bty = 'n')
dev.off()

##AUC curve
AUC1=timeROC(T=rt$OS.time,delta=rt$OS.event,
             marker=rt$Riskscore,cause=1,
             weighting='marginal',
             times=c(10:60),ROC=TRUE,iid = TRUE)
AUC2=timeROC(T=rt$OS.time,delta=rt$OS.event,
             marker=rt$Score_cli,cause=1,
             weighting='marginal',
             times=c(10:60),ROC=TRUE,iid = TRUE)
AUC3=timeROC(T=rt$OS.time,delta=rt$OS.event,
             marker=rt$modelScore,cause=1,
             weighting='marginal',
             times=c(10:60),ROC=TRUE,iid = TRUE)

pdf(file="01_WCH_OS_AUC.pdf", width=8, height=8)
plot(AUC1$times, AUC1$AUC, 
     type = "n", 
     xlab = "time(Months)",  
     ylab = "AUC",
     xaxt = "n", 
     ylim = c(0.5, 1.0))
axis(1, at = seq(12, 60, by = 12), labels = seq(12, 60, by = 12))
abline(v = c(12, 36, 60), lty = 2, col = "gray50", lwd = 1)
plotAUCcurve(AUC1, col = "#A51C36", add = TRUE)
plotAUCcurve(AUC2, col = "#3D5488", add = TRUE)
plotAUCcurve(AUC3, col = "#DBB428", add = TRUE)
legend("topright", 
       legend = c("Riskscore", "other variables", "all"),
       col = c("#A51C36", "#3D5488", "#DBB428"),
       lwd = 1.5,
       cex = 0.8,
       bty = "o",
       bg = "white")
dev.off()

#heatmap
library(pheatmap)
rt <- rt[order(rt$Riskscore),]
rt2 <- rt[,c(14:18)]
rt2 <- as.matrix(t(rt2))
annotation=data.frame(type=rt$risk)
colnames(annotation)[1] <- 'Group'
annotation[annotation$Group=='high',1] <- 'Risk: High'
annotation[annotation$Group=='low',1] <- 'Risk: Low'
rownames(annotation)=colnames(rt2)
ann_colors = list(
  Group = c('Risk: Low'="#1B9E77",'Risk: High'="#D95F02"))
pdf(file="01_WCH_OS_heat.pdf", width=8, height=6)
pheatmap(rt2,
         cluster_cols = F,
         cluster_rows = FALSE,
         color = colorRampPalette(c(rep("#333399",3.5), "white", rep("#CC3333",3.5)))(60),
         show_rownames = T,
         show_colnames = F,
         scale="row", 
         #border_color ="NA",
         fontsize = 7,
         fontsize_row=8,
         fontsize_col=3,
         annotation_col=annotation,
         annotation_colors = ann_colors)
dev.off()

#DCA curve
library(rms)
library(ggDCA)
library(survival)
rt <- read.csv('01_WCH_rt_OS.csv')
dc <- rt[,c('age','sex','Location','Stage','chemoradiotherapy','Riskscore','OS.time','OS.event')]
cox_fit1 <- coxph(Surv(OS.time,OS.event) ~ ., data = dc)
cox_fit2 <- coxph(Surv(OS.time,OS.event) ~ age+sex+Location+Stage+chemoradiotherapy, data = dc)
cox_fit3 <- coxph(Surv(OS.time,OS.event) ~ Riskscore, data = dc)
df1 <- ggDCA::dca(cox_fit1,cox_fit2,cox_fit3,
                  times = c(12, 36, 60))
library(ggsci)
pdf(file="01_WCH_OS_DCA.pdf", width=10, height=6)
dca <- ggplot(df1,linetype = F)+
  scale_color_manual(name="Model Type",
                     values = c("#DBB428", "#3D5488", "#A51C36", "#374E55", "#79AF97"),
                     labels=c('All variables','other variables','Riskscore','All','None'))+
  theme_bw(base_size = 14)+
  theme(legend.position.inside = c(0.8,0.75),
        legend.background = element_blank())
print(dca)
dev.off()

##clibration curve
mul_cox_1 <- cph(Surv(OS.time,OS.event == 1) ~ age+sex+Location+Stage+chemoradiotherapy+Riskscore,data = rt,
                 x = T,y = T, surv = T,time.inc = 12)
mul_cox_2 <- cph(Surv(OS.time,OS.event == 1) ~ age+sex+Location+Stage+chemoradiotherapy+Riskscore,data = rt,
                 x = T,y = T, surv = T,time.inc = 36)
mul_cox_3 <- cph(Surv(OS.time,OS.event == 1) ~ age+sex+Location+Stage+chemoradiotherapy+Riskscore,data = rt,
                 x = T,y = T, surv = T,time.inc = 60)
cal1 <- calibrate(mul_cox_1, 
                  cmethod = 'KM',  
                  method = "boot",  
                  u = 12,         
                  m = 50,         
                  B = 1000)      
cal2 <- calibrate(mul_cox_2, 
                  cmethod = 'KM',   
                  method = "boot",u = 36, m = 50, B = 1000)
cal3 <- calibrate(mul_cox_3, 
                  cmethod = 'KM',   
                  method = "boot",u = 58, m = 50, B = 1000)
pdf(file="01_WCH_OS_celi.pdf", width=7, height=6)
par(mar = c(6, 6, 3, 3))
plot(cal1,                        
     lwd=1.5,                      
     lty=3,                       
     conf.int=T,                  
     errbar.col="#2874C5",       
     col="#2874C5",               
     xlim=c(0,1),                
     ylim=c(0,1),                 
     xlab="Nomogram-Predicted Probability of OS",
     ylab="Actual OS (proportion)",               
     subtitles = F)              
plot(cal2,lwd = 1.5,lty = 3,errbar.col = c("#f87669"),
     xlim = c(0,1),ylim= c(0,1),col = c("#f87669"),add = T)
plot(cal3,lwd = 1.5,lty = 3,errbar.col = c("#e6b707"),
     xlim = c(0,1),ylim= c(0,1),col = c("#e6b707"),add = T)
legend("topleft", 
       legend = c("1-year","3-year","5-year"), 
       col =c("#2874C5","#f87669","#e6b707"), 
       lwd = 1.5,
       cex = 0.9,
       bty = "n")
dev.off()
