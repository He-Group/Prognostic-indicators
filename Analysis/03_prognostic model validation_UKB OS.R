#UKB validation developing
rt <- read.csv('01_WCH_rt_OS.csv')
Coef <- read.csv('01_Coef_LASSO_OS.csv')
UKB_pro <- read.csv('03_UKB_pro.csv')
lassoPro <- Coef$X
UKB_pro <- UKB_pro[-which(UKB_pro$fu_time<30),]
#UKB_pro <- UKB_pro[-which(UKB_pro$lag_time>365*3),]  #for 3year
#UKB_pro <- UKB_pro[-which(UKB_pro$lag_time<(-365*3)),]
UKB_pro <- UKB_pro[,c(1:16,which(colnames(UKB_pro) %in% colnames(rt)[14:18]))]

library(survival)
library(mice)
library(VIM)
UKB_aggr = aggr(UKB_pro,
                col=mdc(1:2), 
                numbers=TRUE, 
                sortVars=TRUE, 
                labels=names(UKB_pro), 
                cex.axis=.7, gap=3, 
                ylab=c("Proportion of missingness","Missingness Pattern"))
mice_data <- mice(UKB_pro,
                  meth='pmm', seed=500)
myce_data_df <- complete(mice_data)
UKB_pro <- myce_data_df
UKB_pro <- na.omit(UKB_pro)

model <- coxph(
  Surv(time = UKB_pro$lag_time, time2 = UKB_pro$exit.time, event = UKB_pro$OS_event) ~ ADAM8+FUT3_FUT5+MUC13+STC1+TNFRSF10B,
  data = UKB_pro)
summary(model)
score <- UKB_pro
score <- score[,Coef[,1]]
for(i in (1:nrow(Coef))){       
  score[,i] <- (score[,i]*Coef[i,2])}
riskScore = rowSums(score)
UKB_pro <- cbind(UKB_pro,Riskscore=as.vector(riskScore))
UKB_pro$risk <- as.vector(ifelse(UKB_pro$Riskscore > median(UKB_pro$Riskscore),  "high", "low"))
UKB_pro <- UKB_pro[order(UKB_pro$Riskscore),]
UKB_pro$fu_time <- UKB_pro$fu_time/365.25*12
write.table(UKB_pro,"03_UKB_3year_rt.csv",col.names = T,row.names = F,sep=",",quote=FALSE)

rt <- read.csv('03_UKB_3year_rt.csv')
cox_test <- coxph(Surv(time = lag_time, time2 = exit.time, event = OS_event) ~ age_diag+sex+ethnicity_id+Site+Riskscore, data = rt)
summary(cox_test)

#KM
library(pheatmap)
library(tibble)
library(stringr)
library(survival)
library(dplyr)
library(survminer)
#median follow-up time
UKB_pro <- read.csv('03_UKB_pro.csv')
rt2 <- UKB_pro
rt2$reverse_status <- ifelse(rt2$OS.event == 1, 0, 1)
fit <- survfit(Surv(OS.time, reverse_status) ~ 1, data = rt2)
summary(fit)$table["median"]

rt <- read.csv('03_UKB_all_rt.csv')
rt1 <- rt
fit <- survfit(Surv(fu_time, OS_event)~risk, data=rt1)
pdf(file="03_UKB_3year_KM.pdf", width=9, height=9)
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
  xlim=c(0,120),
  break.x.by=30,
  censor=F)
print(KM)
dev.off()

#ROC
library(survival)
library(survminer)
library(timeROC)
rt <- read.csv('03_UKB_3year_rt.csv')
fit <- coxph(Surv(fu_time, OS_event)~age_diag+sex+ethnicity_id+Site+Riskscore,data=rt)
Coef_fit <- as.data.frame(fit$coefficients)
score <- rt[,c('age_diag','sex','ethnicity_id','Site','Riskscore')]
score <- score[,rownames(Coef_fit)]
for(i in (1:nrow(Coef_fit))){       
  score[,i] <- score[,i]*Coef_fit[i,1]}
modelScore = rowSums(score)
rt <- cbind(rt,modelScore=as.vector(modelScore))
fit_cli <- coxph(Surv(fu_time, OS_event)~age_diag+sex+ethnicity_id+Site,data=rt)
Coef_cli <- as.data.frame(fit_cli$coefficients)
score_cli <- rt[,c('age_diag','sex','ethnicity_id','Site')]
score_cli <- score_cli[,rownames(Coef_cli)]
for(i in (1:nrow(Coef_cli))){       
  score_cli[,i] <- score_cli[,i]*Coef_cli[i,1]}
Score_cli = rowSums(score_cli)
rt <- cbind(rt,Score_cli=as.vector(Score_cli))

ROC_rt=timeROC(T=rt$fu_time,delta=rt$OS_event,
               marker=rt$modelScore,cause=1,
               weighting = "marginal",
               times=c(12,36,60),ROC=TRUE,iid=T)
CI <- as.data.frame(confint(ROC_rt, level = 0.95)$CI_AUC)
CI$AUC <- paste0(sprintf("%.03f",ROC_rt$AUC)," (",sprintf("%.03f",CI$`2.5%`*0.01),"-",sprintf("%.03f",CI$`97.5%`*0.01),")")

pdf(file="03_UKB_3year_ROC_model.pdf", width=5, height=5)
plot(ROC_rt,time=12,col="#2874C5",title=FALSE,lwd=2)
plot(ROC_rt,time=36,col="#f87669",add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=60,col="#e6b707",add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',CI$AUC[1]),
         paste0('AUC at 3 years: ',CI$AUC[2]),
         paste0('AUC at 5 years: ',CI$AUC[3])),
       col=c("#2874C5", "#f87669", "#e6b707"), lwd=3,cex=0.75,bty = 'n')
dev.off()
#AUC
AUC1=timeROC(T=rt$fu_time,delta=rt$OS_event,
             marker=rt$Riskscore,cause=1,
             weighting='marginal',
             times=c(10:60),ROC=TRUE,iid = TRUE)
AUC2=timeROC(T=rt$fu_time,delta=rt$OS_event,
             marker=rt$Score_cli,cause=1,
             weighting='marginal',
             times=c(10:60),ROC=TRUE,iid = TRUE)
AUC3=timeROC(T=rt$fu_time,delta=rt$OS_event,
             marker=rt$modelScore,cause=1,
             weighting='marginal',
             times=c(10:60),ROC=TRUE,iid = TRUE)
pdf(file="03_UKB_3year_AUC.pdf", width=8, height=8)
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
       legend = c("Riskscore", "other variables", "3year"),
       col = c("#A51C36", "#3D5488", "#DBB428"),
       lwd = 1.5,
       cex = 0.8,
       bty = "o",
       bg = "white")
dev.off()

#heatmap
library(pheatmap)
rt <- rt[order(rt$Riskscore),]
rt2 <- rt[,c(17:21)]
rt2 <- as.matrix(t(rt2))
annotation=data.frame(type=rt$risk)
colnames(annotation)[1] <- 'Group'
annotation[annotation$Group=='high',1] <- 'Risk: High'
annotation[annotation$Group=='low',1] <- 'Risk: Low'
rownames(annotation)=colnames(rt2)
ann_colors = list(
  Group = c('Risk: Low'="#1B9E77",'Risk: High'="#D95F02"))
pdf(file="03_UKB_3year_heat.pdf", width=8, height=6)
pheatmap(rt2,
         cluster_cols = F,
         cluster_rows = FALSE,
         color = colorRampPalette(c(rep("#333399",3.5), "white", rep("#CC3333",3.5)))(60),
         show_rownames = T,
         show_colnames = F,
         scale="row",  
         fontsize = 7,
         fontsize_row=8,
         fontsize_col=3,
         annotation_col=annotation,
         annotation_colors = ann_colors)
dev.off()

#DCA
library(rms)
library(ggDCA)
library(survival)
rt <- read.csv('03_UKB_3year_rt.csv')
dc <- rt[,c('age_diag','sex','ethnicity_id','Site','Riskscore','fu_time','OS_event')]
cox_fit1 <- coxph(Surv(fu_time, OS_event) ~ ., data = dc)
cox_fit2 <- coxph(Surv(fu_time, OS_event) ~ age_diag+sex+ethnicity_id+Site, data = dc)
cox_fit3 <- coxph(Surv(fu_time, OS_event) ~ Riskscore, data = dc)
df1 <- ggDCA::dca(cox_fit1,cox_fit2,cox_fit3,
                  times = c(12, 36, 60))
library(ggsci)
pdf(file="03_UKB_3year_DCA.pdf", width=10, height=6)
dca <- ggplot(df1,linetype = F)+
  scale_color_manual(name="Model Type",
                     values = c("#DBB428", "#3D5488", "#A51C36", "#374E55", "#79AF97"),
                     labels=c('All variables','other variables','Riskscore','All','None'))+
  theme_bw(base_size = 14)+
  theme(legend.position.inside = c(0.8,0.75),
        legend.background = element_blank())
print(dca)
dev.off()

#calibration
library(survival)
library(rms)
rt <- read.csv('03_UKB_3year_rt.csv')
mul_cox_1 <- cph(Surv(fu_time, OS_event == 1) ~ age_diag+sex+ethnicity_id+Site+Riskscore+Riskscore,data = rt,
                 x = T,y = T, surv = T,time.inc = 12)
mul_cox_2 <- cph(Surv(fu_time, OS_event == 1) ~ age_diag+sex+ethnicity_id+Site+Riskscore+Riskscore,data = rt,
                 x = T,y = T, surv = T,time.inc = 36)
mul_cox_3 <- cph(Surv(fu_time, OS_event == 1) ~ age_diag+sex+ethnicity_id+Site+Riskscore+Riskscore,data = rt,
                 x = T,y = T, surv = T,time.inc = 59)
cal1 <- calibrate(mul_cox_1, 
                  cmethod = 'KM',   
                  method = "boot", 
                  u = 12,         
                  m = 80,     
                  B = 1000)    
cal2 <- calibrate(mul_cox_2, 
                  cmethod = 'KM',   
                  method = "boot",u = 36, m = 80, B = 1000)
cal3 <- calibrate(mul_cox_3, 
                  cmethod = 'KM',   
                  method = "boot",u = 60, m = 80, B = 1000)
pdf(file="03_UKB_3year_celi.pdf", width=7, height=6)
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

