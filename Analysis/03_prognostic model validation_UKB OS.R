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

# Other analyses are similar to those in the above files

