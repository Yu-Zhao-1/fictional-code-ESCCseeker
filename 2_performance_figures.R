
library(tidyverse)
library(glmnet)
library(caret)
library(ROCR)
library(reportROC)
library(pROC)
library(readxl)
library(randomForest)
library(bestglm)
library(tidymodels)
library(bonsai)
library(lightgbm)
library(e1071)
library(readxl)
require(xgboost)
require(Matrix)
require(data.table)
library(survival)
library(dplyr)
library(caret)
library(tidyverse)
#library(MASS)
library(ROCR)
library(reportROC)
library(shapviz)
library(reshape2)
library(ggprism)
library(rms)
library(rmda)
agresti_coull_ci <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level)/2)
  n_tilde <- n + z^2
  p_tilde <- (x + z^2/2) / n_tilde
  se <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  lower <- p_tilde - z * se
  upper <- p_tilde + z * se
  p=x/n
  i=paste0(x,"/",n)
  
  return(c(p,lower, upper))
}
library(PRROC)
library(boot)
##1.Train,Test,Verification(yp1,yp2,yp3)####
load("Performance.RData")
##1.ROC
yp=yp3#yp1,yp2,yp3
yp$class=ifelse(yp$group==0,0,1)
##Control vs ESCC/HGIN
roc0<-roc(class~Methy9,yp,direction = "<")
roc1<-roc(class~ESCCseeker,yp,direction = "<")
roc2<-roc(class~CEA,yp,direction = "<")
roc3<-roc(class~SCC,yp,direction = "<")
roc4<-roc(class~CYR,yp,direction = "<")
ci0=ci(roc0)
ci1=ci(roc1)
ci2=ci(roc2)
ci3=ci(roc3)
ci4=ci(roc4)
plot(roc1, col = "#EF312C", main = "ESCC/HGIN vs Control",lwd=1,lty=1,cex.lab = 1.2) 
lines(roc0, col = "#E39611",lwd=1,lty=1) 
lines(roc2, col = "#3462AA",lwd=1,lty=1) 
lines(roc3, col = "#9B3E88",lwd=1,lty=1) 
lines(roc4, col = "#84C85A",lwd=1,lty=1) 
legend("bottomright",
       legend = c(paste0("ESCCseeker:", round(ci1[2], 3), paste0(";(", round(ci1[1],3), "-", round(ci1[3],3),")")), 
                  paste0("9-Methy:", round(ci0[2], 3), paste0(";(", round(ci0[1],3), "-", round(ci0[3],3),")")),
                  paste0("CEA:", round(ci2[2], 3), paste0(";(", round(ci2[1],3), "-", round(ci2[3],3),")")),
                  paste0("SCC:", round(ci3[2], 3), paste0(";(", round(ci3[1],3), "-", round(ci3[3],3),")")), 
                  paste0("CYFRA21-1:", round(ci4[2], 3), paste0(";(", round(ci4[1],3), "-", round(ci4[3],3),")"))
       ), 
       col = c("#EF312C","#E39611","#3462AA","#9B3E88","#84C85A"), 
       lty = 1, lwd=1,cex = 0.95,bty = "n")
##Control vs HGIN
roc0<-roc(group~Methy9,yp[yp$group!=2,],direction = "<")
roc1<-roc(group~ESCCseeker,yp[yp$group!=2,],direction = "<")
roc2<-roc(group~CEA,yp[yp$group!=2,],direction = "<")
roc3<-roc(group~SCC,yp[yp$group!=2,],direction = "<")
roc4<-roc(group~CYR,yp[yp$group!=2,],direction = "<")
ci0=ci(roc0)
ci1=ci(roc1)
ci2=ci(roc2)
ci3=ci(roc3)
ci4=ci(roc4)
plot(roc1, col = "#EF312C", main = "HGIN vs Control",lwd=1,lty=1,cex.lab = 1.2)
lines(roc0, col = "#E39611",lwd=1,lty=1) 
lines(roc2, col = "#3462AA",lwd=1,lty=1) 
lines(roc3, col = "#9B3E88",lwd=1,lty=1) 
lines(roc4, col = "#84C85A",lwd=1,lty=1) 
legend("bottomright",
       legend = c(paste0("ESCCseeker:", round(ci1[2], 3), paste0(";(", round(ci1[1],3), "-", round(ci1[3],3),")")), 
                  paste0("9-Methy:", round(ci0[2], 3), paste0(";(", round(ci0[1],3), "-", round(ci0[3],3),")")),
                  paste0("CEA:", round(ci2[2], 3), paste0(";(", round(ci2[1],3), "-", round(ci2[3],3),")")),
                  paste0("SCC:", round(ci3[2], 3), paste0(";(", round(ci3[1],3), "-", round(ci3[3],3),")")), 
                  paste0("CYFRA21-1:", round(ci4[2], 3), paste0(";(", round(ci4[1],3), "-", round(ci4[3],3),")"))
       ), 
       col = c("#EF312C","#E39611","#3462AA","#9B3E88","#84C85A"), 
       lty = 1, lwd=1,cex = 0.95,bty = "n")
##Control vs ESCC
roc0<-roc(group~Methy9,yp[yp$group!=1,],direction = "<")
roc1<-roc(group~ESCCseeker,yp[yp$group!=1,],direction = "<")
roc2<-roc(group~CEA,yp[yp$group!=1,],direction = "<")
roc3<-roc(group~SCC,yp[yp$group!=1,],direction = "<")
roc4<-roc(group~CYR,yp[yp$group!=1,],direction = "<")
ci0=ci(roc0)
ci1=ci(roc1)
ci2=ci(roc2)
ci3=ci(roc3)
ci4=ci(roc4)
plot(roc1, col = "#EF312C", main = "ESCC vs Control",lwd=1,lty=1,cex.lab = 1.2) 
lines(roc0, col = "#E39611",lwd=1,lty=1) 
lines(roc2, col = "#3462AA",lwd=1,lty=1) 
lines(roc3, col = "#9B3E88",lwd=1,lty=1) 
lines(roc4, col = "#84C85A",lwd=1,lty=1) 
legend("bottomright",
       legend = c(paste0("ESCCseeker:", round(ci1[2], 3), paste0(";(", round(ci1[1],3), "-", round(ci1[3],3),")")), 
                  paste0("9-Methy:", round(ci0[2], 3), paste0(";(", round(ci0[1],3), "-", round(ci0[3],3),")")),
                  paste0("CEA:", round(ci2[2], 3), paste0(";(", round(ci2[1],3), "-", round(ci2[3],3),")")),
                  paste0("SCC:", round(ci3[2], 3), paste0(";(", round(ci3[1],3), "-", round(ci3[3],3),")")), 
                  paste0("CYFRA21-1:", round(ci4[2], 3), paste0(";(", round(ci4[1],3), "-", round(ci4[3],3),")"))
       ), 
       col = c("#EF312C","#E39611","#3462AA","#9B3E88","#84C85A"), 
       lty = 1, lwd=1,cex = 0.95,bty = "n")
###.Table S4####
data=rbind(yp1,yp2,yp3)
data$class=ifelse(data$group==0,"Control",ifelse(data$group==1,"HGIN","ESCC"))
data$cohort=rep(c("Training set","Test set","Verification set"),times=c(nrow(yp1),nrow(yp2),nrow(yp3)))
data$class=factor(data$class,levels = c("Control","HGIN","ESCC"))
data$group=as.numeric(data$class)-1
yp=data.frame(cohort=data$cohort,class=data$class,group=data$group,ESCCseeker=data$ESCCseeker,Methy9=data$Methy9,CEA=data$CEA,SCC=data$SCC,CYR=data$CYR,TMcombine=data$TMcombine)
da=yp[yp$cohort=="Verification set"&yp$class%in%c("ESCC","HGIN"),]
da$class="ESCC/HGIN"
da$group=1
yp=rbind(yp,da)
cutoffs=c(ESCCseeker = 0.5,Methy9=0.43, CEA = 5, SCC=2.7,CYR = 3.3,TMcombine=0)
sen=data.frame()  
spe=data.frame()
to=data.frame()
for (set in c("Training set","Test set","Verification set")) {
  da=yp[yp$cohort==set,]
  for (h in c("ESCC","HGIN","ESCC/HGIN")) {
    if (sum(da$class==h)==0) next
    data=da[da$class%in%c(h,"Control"),]
    for (i in 4:9) {
      data$predictor=data[,i]
      cutoff=cutoffs[i-3]
      ii=names(cutoffs)[i-3]
      set.seed(225)
      pred_class <- ifelse(data$predictor  >  cutoff, 1, 0)  # 预测分类
      tp <- sum(data$group != 0 & pred_class == 1)  # 真阳性
      fp <- sum(data$group == 0 & pred_class == 1)  # 假阳性
      fn <- sum(data$group != 0 & pred_class == 0)  # 假阴性
      tn <- sum(data$group == 0 & pred_class == 0)  # 真阴性
      a=tp + fn
      current_sens <-pmax(0, pmin(1, agresti_coull_ci(tp , a)))  # 灵敏度
      a=tn + fp
      current_spec <- pmax(0, pmin(1, agresti_coull_ci(tn , a)))  # 特异度
      a=tn + tp
      current_accu <- pmax(0, pmin(1, agresti_coull_ci(a , nrow(data)))) 
      roc1 <- roc(group ~ predictor, data,direction = "<")
      AUC <- as.numeric(round(ci(roc1),3))
      AUC=paste0(AUC[2],"(",AUC[1],"-",AUC[3],")")
      sen <- rbind(sen, cbind(cohort=set,group=h,Sensitivity=names(cutoffs)[i-3], lower = current_sens[2],upper = current_sens[3],current = current_sens[1]))
      spe <- rbind(spe, cbind(cohort=set,group=h,Specificity=names(cutoffs)[i-3], lower = current_spec[2],upper = current_spec[3],current = current_spec[1]))
      auc <- rbind(auc, cbind(cohort=set,group=h,AUC=names(cutoffs)[i-3],lower = AUC[1],upper = AUC[3],current = AUC[2]))
      Sen =paste0(round(current_sens[1],4)*100,"(",round(current_sens[2],4)*100,"-",round(current_sens[3],4)*100,")")
      Spe =paste0(round(current_spec[1],4)*100,"(",round(current_spec[2],4)*100,"-",round(current_spec[3],4)*100,")")
      to=rbind(to,cbind(set=set,type=h,panel=ii,AUC = AUC,Sen=Sen,Spe=Spe))
    }
  }
}
##2.Other cancer####
df=yp4
df$class=factor(df$class,levels = c("ESCC","BC","LC","HCC"))
df$class=as.numeric(df$class)-1
df$group=ifelse(df$class!=0,1,0)
yp=df
roc1<-roc(group~ESCCseeker,yp)
roc2<-roc(class~ESCCseeker,yp[yp$class%in%c(0,1),])
roc3<-roc(class~ESCCseeker,yp[yp$class%in%c(0,2),])
roc4<-roc(class~ESCCseeker,yp[yp$class%in%c(0,3),])
ci1=ci(roc1)
ci2=ci(roc2)
ci3=ci(roc3)
ci4=ci(roc4)
plot(roc1, col = "#B90019", main = "OC",lwd=1,lty=1)
lines(roc2, col = "#EEA210",lwd=1,lty=1) 
lines(roc3, col ="#7EC325",lwd=1,lty=1) 
lines(roc4, col = "#9565A5",lwd=1,lty=1) 
legend("bottomright",
       legend = c(
         paste0("ESCC vs Other tumors (AUC:", round(ci1[2], 3), paste0(",", round(ci1[1],3), "-", round(ci1[3],3),")")), 
       paste0("ESCC vs BC (AUC:", round(ci2[2], 3), paste0(",", round(ci2[1],3), "-", round(ci2[3],3),")")),
       paste0("ESCC vs LC (AUC:", round(ci3[2], 3), paste0(",", round(ci3[1],3), "-", round(ci3[3],3),")")), 
       paste0("ESCC vs HCC (AUC:", round(ci4[2], 3), paste0(",", round(ci4[1],3), "-", round(ci4[3],3),")"))
       ),col = c("#B90019","#EEA210", "#7EC325","#9565A5"),
       lty = 1,lwd=1,
       bty="n",
       cex = 0.75) 
##3.Subtype####
data1$gg=as.factor("All")
cutoffs=c(ESCCseeker = 0.5, CEA = 5, SCC=2.7,CYR = 3.3,TMcombine=0)
sen=data.frame()  
spe=data.frame()
auc=data.frame()
cutoff=yd_cutoff
df_long <- pivot_longer(data1, cols = c(6:11), names_to = "key", values_to = "value")
df_sorted <- df_long %>%
  arrange(key, value)
df_clean <- df_sorted[complete.cases(df_sorted), ]
df_clean$key <- factor(df_clean$key, levels = c("class" ,"stage","LN","MI","position","size"), ordered = TRUE)
df_clean$value <- factor(df_clean$value, levels = c("Control","ESCC" ,"1","2","3","4","LN(-)",  "LN(+)", "MI(-)" , "MI(+)","Upper", "Middle" ,"Lower","≤3cm",">3cm" ), ordered = TRUE)
to=data.frame()
for (ii in 1:4) {
  da=df_clean[,c(1,(ii+1),8:9)]
  colnames(da)[2]="predictor"
  cutoff=cutoffs[ii]
  ii=names(cutoffs)[[ii]]
  for (h in levels(da$gg)) {
    df=da[da$gg==h,-1]
    for (i in levels(df$value)[-(1)]) {
      if (!any(df$value == i)) next
      data=df[df$value%in%c(i,"Control"),]
      colnames(data)=c("predictor","type","class")
      data$class=ifelse(data$class=="Control",0,1)
      set.seed(225)
      pred_class <- ifelse(data$predictor  >  cutoff, 1, 0)  # 预测分类
      tp <- sum(data$class == 1 & pred_class == 1)  # 真阳性
      fp <- sum(data$class == 0 & pred_class == 1)  # 假阳性
      fn <- sum(data$class == 1 & pred_class == 0)  # 假阴性
      tn <- sum(data$class == 0 & pred_class == 0)  # 真阴性
      a=tp + fn
      current_sens <-pmax(0, pmin(1, agresti_coull_ci(tp , a)))  # 灵敏度
      a=tn + fp
      current_spec <- pmax(0, pmin(1, agresti_coull_ci(tn , a)))  # 特异度
      a=tn + tp
      current_accu <- pmax(0, pmin(1, agresti_coull_ci(a , nrow(data))))  # 特异度
      roc <- roc(class ~ predictor, data,direction = "<")
      AUC=as.numeric(round(ci(roc),3))
      data=arrange(data,desc(type))
      sen <- rbind(sen, cbind(gg=h,type=data[1,2],class=i, lower = current_sens[2],upper = current_sens[3],current = current_sens[1]))
      spe <- rbind(spe, cbind(gg=h,type=data[1,2],class=i, lower = current_spec[2],upper = current_spec[3],current = current_spec[1]))
      auc <- rbind(auc, cbind(gg=h,type=data[1,2],class=i,lower = AUC[1],upper = AUC[3],current = AUC[2]))
      AUC=paste0(AUC[2],"(",AUC[1],"-",AUC[3],")")
      Sen =paste0(round(current_sens[1],4)*100,"(",round(current_sens[2],4)*100,"-",round(current_sens[3],4)*100,")")
      Spe =paste0(round(current_spec[1],4)*100,"(",round(current_spec[2],4)*100,"-",round(current_spec[3],4)*100,")")
      n=sum(data$class==1)
      to=rbind(to,cbind(gg=h,panel=ii,type=data[1,2],i=i,n=n,AUC = AUC,Sen=Sen,Spe=Spe))
    }
  }
}
colnames(to)=c("cohort","Panel","type","class","n","AUC(95% CI)","Sensitivity(95% CI),%","Specificity(95% CI),%")
to=pivot_wider(to,names_from = c(1),values_from = c(5:8))
to=dplyr::arrange(to,type,class)
##.sex&age==========================
yp=pivot_longer(data1,cols = 12:13,names_to = "type",values_to = "value")
yp$value=factor(yp$value,levels = c("<60", ">60", "Female", "Male"))
yp$class=ifelse(yp$class=="ESCC",1,0)
sen=data.frame()  
spe=data.frame()
auc=data.frame()
off=c(ESCCseeker=0.5,CEA=5,SCC=2.7,CYR=3.3)
to=data.frame()
for (ii in 1:4) {
  df=yp[,c(1,ii+1,13,6)]
  cutoff=off[ii]
  panel=names(off)[ii]
  for (h in levels(df$gg)) {
    da=df[df$gg==h,-1]
    for (i in levels(da$value)) {
      data=da[da$value==i,]
      colnames(data)=c("predictor","type","class")
      data$class=ifelse(data$class==0,0,1)
      set.seed(225)
      cutoff=cutoff
      pred_class <- ifelse(data$predictor  >=  cutoff, 1, 0)  # 预测分类
      tp <- sum(data$class == 1 & pred_class == 1)  # 真阳性
      fp <- sum(data$class == 0 & pred_class == 1)  # 假阳性
      fn <- sum(data$class == 1 & pred_class == 0)  # 假阴性
      tn <- sum(data$class == 0 & pred_class == 0)  # 真阴性
      a=tp + fn
      current_sens <-pmax(0, pmin(1, agresti_coull_ci(tp , a)))  # 灵敏度
      a=tn + fp
      current_spec <- pmax(0, pmin(1, agresti_coull_ci(tn , a)))  # 特异度
      a=tn + tp
      current_accu <- pmax(0, pmin(1, agresti_coull_ci(a , nrow(data))))  # 特异度
      roc <- roc(class ~ predictor, data,direction = "<")
      AUC=as.numeric(round(ci(roc),3))
      sen <- rbind(sen, cbind(type=data[1,2],class=i, lower = current_sens[2],upper = current_sens[3],current = current_sens[1],panel=panel))
      spe <- rbind(spe, cbind(type=data[1,2],class=i, lower = current_spec[2],upper = current_spec[3],current = current_spec[1],panel=panel))
      
      auc <- rbind(auc, cbind(type=data[1,2],class=i,lower = AUC[1],upper = AUC[3],current = AUC[2],panel=panel))
      AUC=paste0(AUC[2],"(",AUC[1],"-",AUC[3],")")
      Sen =paste0(round(current_sens[1],4)*100,"(",round(current_sens[2],4)*100,"-",round(current_sens[3],4)*100,")")
      Spe =paste0(round(current_spec[1],4)*100,"(",round(current_spec[2],4)*100,"-",round(current_spec[3],4)*100,")")
      
      c=sum(data$class==1)
      n=sum(data$class==0)
      to=rbind(to,cbind(cohort=h,i=i,cancer=c,normal=n,AUC = AUC,Sen=Sen,Spe=Spe,panel=panel))
    }
  }
}
colnames(to)=c("cohort","type","cancer.n","normal.n","AUC(95% CI)","Sensitivity(95% CI),%","Specificity(95% CI),%","Panel")

#Phase Ⅱ####
##Fig.5B####
library(ggplot2)
library(ggpubr)
library(dplyr)
data <- data.frame(group=yp5$class,score=as.numeric(yp5$ESCCseeker))
data <- data[order(-data$score),]
data <- data%>%mutate(level=ifelse(score>=0.5,"High-risk","Low-risk"))
colnames(data)=c("status","predicted_risk_score","risk_group")
data$i=1:nrow(data)
points_data <- data.frame(
  x = max(which(data$risk_group=="High-risk"))-0.5, 
  y = 0.5
)
p3 <- ggplot(data, aes(x = (i-1), y = predicted_risk_score)) +
  geom_area(aes(y = predicted_risk_score), fill = "#8FACDD", alpha = 0.5) +  # 填充曲线下的区域
  geom_line(color = "#8FACDD") +#,linewidth = 16B94D1   878788  #geom_line(color = "#6B94D1", size = 1) +#6B94D1   878788
  # 添加从点向左和向上的线段
  geom_segment(data = points_data, 
               aes(x = x, xend = 0, 
                   y = y, yend = y), 
               linetype = "dashed", color = "gray") +  # ,linewidth=1向左的水平线
  geom_segment(data = points_data, 
               aes(x = x, xend = x, 
                   y = y, yend = 1), 
               linetype = "dashed", color = "gray") +  #,linewidth=1 向上的垂直线
  geom_point(data = points_data, aes(x = x, y = y), 
             color = "black", size = 2.5) +
  labs(y = "Predicted risk score", x = "") +
  theme_classic() +  # 使用经典主题，没有网格线
  theme(
    panel.border = element_blank(),  # 去掉面板边框
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(color = "black"),  #size = 1 , linewidth = 1设置坐标轴线的粗细
    axis.ticks.y = element_line(color = "black"),#, linewidth = 1
    
    panel.grid = element_blank()
  )+
  scale_x_continuous(expand = c(0, 0)) +  # 设置 x 轴从 0 开始
  scale_y_continuous(expand = c(0, 0))    # 设置 y 轴从 0 开始
p3

data$status=factor(ifelse(data$status%in%c("ESCC","HGIN"),"Yes","No"),levels = c("No","Yes"))
data=pivot_longer(data,cols=c("status","risk_group"),names_to = "status",values_to = "value")
p1 <- ggplot(data, aes(x = (i-1), fill = value)) +
  geom_tile(aes(y = status, height = 0.8, width = 1),linewidth=2.5) +
  scale_fill_manual(values = c("No" = "#DBDBDB", "Yes" ="#8FACDD","High-risk" = "#FFE799",  "Low-risk" = "#A8D18C")) +#"#B20D35"
  labs(y = "", x = "") +
  theme_classic() + 
  theme(
    panel.border = element_blank(), 
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )
p1
##Fig.5D####
library(circlize)
color_mapping <- c("0" = "grey70", "1" = "#FFE799")#"grey70","#006BB9","#A8D18C"，"0" = "#0273c2","#8FACDD", "1" = "#efc001"，"#FCC100"
yp=yp5#QZ
yp$ESCCseeker=ifelse(yp$ESCCseeker > 0.5, "1", "0")
yp=yp[yp$class%in%c("ESCC","HGIN"),]
m=yp[,2:3]
circos.par(gap.after = c(90),start.degree = 90)
circos.heatmap(
  m,
  col = color_mapping,  # 使用离散颜色映射
  track.height = 0.6,
  bg.border = "white",#"black",
  cluster = TRUE,
  #cell_width = rep(2, nrow(mat)), # 设置格子宽度
  cell.border ="white", #边缘颜色
  dend.track.height = 0.1,
  show.sector.labels = FALSE
)
circos.clear()
#

matrix<- table(yp$class,yp$ESCCseeker)  #混淆矩阵
E_HGIN <- pmax(0, pmin(1, agresti_coull_ci(matrix[2, 2], sum(matrix[2, ]))))
E_ESCC <- pmax(0, pmin(1, agresti_coull_ci(matrix[1, 2], sum(matrix[1, ]))))

matrix<- table(yp$class,yp$TMcombine)  #混淆矩阵
TM_HGIN <- pmax(0, pmin(1, agresti_coull_ci(matrix[2, 2], sum(matrix[2, ]))))
TM_ESCC <- pmax(0, pmin(1, agresti_coull_ci(matrix[1, 2], sum(matrix[1, ]))))

data=data.frame(rbind(E_HGIN=E_HGIN,E_ESCC=E_ESCC,TM_HGIN=TM_HGIN,TM_ESCC=TM_ESCC))
colnames(data)=c("current","lower","upper")
data$group=factor(c("HGIN","ESCC","HGIN","ESCC"),levels = c("HGIN","ESCC"))
data$panel=factor(c("ESCCseeker","ESCCseeker","TMcombine","TMcombine"),levels = c("ESCCseeker","TMcombine"))
data[,1:3]=apply(data[,1:3],2,as.numeric)
data$label=rownames(data)
ggplot(data, aes(x = panel, y = current, fill = panel)) +
  geom_col(width = 1,color=NA) +#,alpha=0.7
  labs(title = "",y = "", x = "") +
  facet_wrap(~group,strip.position="bottom",nrow = 1) +
  #geom_point(size = 3,aes(color = index)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "black",width = 0.2) +  # 添加垂直误差线
  scale_y_continuous(limits = c(0.0, 1.0), breaks = seq(0.0, 1.0, 0.2), labels = c("0.00", "0.2", "0.40", "0.60", "0.80",  "1.00"), expand = c(0, 0)) +
  theme_classic() +  # 使用经典主题，没有网格线
  theme(
    #panel.border = element_blank(),  # 去掉面板边框
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.placement  = "outside", 
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.spacing = unit(0.1, "cm")
  )+
  scale_fill_manual(values = c("ESCCseeker"="#FFDE79","TMcombine"="#8FACDD"))
####Fig.5E-F####
library(dplyr)
library(ggalluvial)
library(ggplot2)
df=data2
plot_data <- df %>%
  dplyr::count(TMcombine, ESCCseeker, Histopathology)
plot_data$TMcombine=factor(plot_data$TMcombine,levels = c("TM(+)","TM(-)"))
plot_data$ESCCseeker=factor(plot_data$ESCCseeker,levels = c("High_risk","Low_risk"))
plot_data$Histopathology=factor(plot_data$Histopathology,levels = c( "High-grade lesions","LGIN","Benign","Healthy"))

colors=c(
  "High_risk" = "#FFE799", 
  "TM(-)" = "#DCD6E8",
  "Low_risk" = "#A8D18C",
  "TM(+)"="#AFA7D1",
  "High-grade lesions"="#6E99CF",
  "Healthy"="#DDDFE0",
  "Benign"="#DDDFE0",
  "HGIN"="#93B2DC",
  "LGIN"="#BDD4E2",
  "ESCC"="#305B90"
)

ggplot(plot_data,
       aes(y = n, axis1 = TMcombine, axis2 = ESCCseeker, 
           axis3 = Histopathology)) +
  geom_alluvium(aes(fill = Histopathology), width = 1/4) +
  geom_stratum(
    aes(fill = after_stat(stratum)), 
    width = 1/4,
    color = NA,
    size = 0.5
  ) +
  scale_fill_manual(
    values = colors
  )+
  geom_label(stat = "stratum", 
             aes(label = paste0(round(after_stat(count)/nrow(df)*100, 1), "%")),
             label.size = 0,
             fill=NA,
             size = 3)+
  scale_x_discrete(limits = c("TMcombine", "ESCCseeker", "Histopathology"),
                   expand = c(0.05, 0.05)) +
  theme_classic() +
  theme(legend.position = "top",
        legend.key.size = unit(0.35, "cm"),        
        legend.text = element_text(size = 7),    
        legend.title = element_text(size = 7)) +
  labs(title = "",
       y = "Number of Patients")
da <- df %>%
  dplyr::count(TMcombine,histopathology)
da1 <- df %>%
  dplyr::count(ESCCseeker, histopathology)
colnames(da1)[1]="group"
colnames(da)[1]="group"
da=rbind(da,da1)
da$panel=rep(c("TMcombine","ESCCseeker"),time=c(10,10))
da$panel=factor(da$panel,levels = c("TMcombine","ESCCseeker"))
da$group=factor(da$group,levels = c("TM(-)","TM(+)","Low_risk","High_risk"))
da$histopathology=factor(da$histopathology,levels = c("ESCC","HGIN","LGIN","Benign","Healthy"))
da=arrange(da,histopathology)
label <- da %>% 
  group_by(panel, group) %>% 
  dplyr::summarise(
    labels = paste(histopathology, n, sep = "=", collapse = "\n"),
    .groups = "drop"
  )
da_total <- da %>% 
  group_by(panel, group) %>% 
  dplyr::summarise(total = sum(n), .groups = "drop")
label <- left_join(label, da_total, by = c("panel", "group"))
ggplot(da, aes(x = group, y = n, fill = histopathology)) +
  geom_col(width = 0.7,colour = NA)  +
  facet_wrap(~panel, scales = "free_x", strip.position="bottom")+ # 按指标分面
  labs(title = "",y = "Number of patients", x = "") +
  theme_classic() +
  scale_y_continuous(limits = c(0.0, 1100),expand = c(0, 0)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(color = "black"),  # 保留坐标轴线
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    strip.placement = "outside",  # 将分面标签放在面板外部
    strip.background = element_blank(),  # 去掉分面标签的背景
    strip.text.x = element_text(size = 12, face = "bold"),  # 设置分面标签的样式
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    legend.spacing = unit(0.1, "cm")
  )+
  scale_fill_manual(values = c("Healthy"="#E2E4E5","Benign"="#DADCDD","HGIN"="#3768A5","LGIN"="#93B2DC",
                               "ESCC"="#2D5F88"))+
  geom_text(
    data = label,
    aes(x = group, y = total +20, label = labels),
    inherit.aes = FALSE,
    vjust = 0, size = 2, fontface = "bold"
  )
####Fig.6A####
df=data2
da_screen <- df %>%
  group_by(ESCCseeker, histopathology)%>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(ESCCseeker)%>%
  dplyr::mutate(total = sum(n))%>%
  ungroup()
h <-data.frame(ESCCseeker=c("High_risk","Low_risk"),histopathology=c("ESCC/HGIN","ESCC/HGIN"),
               n=c((da_screen[2,3]+da_screen[3,3])[[1]],(da_screen[7,3]+da_screen[8,3])[[1]]),total=c(da_screen[2,4][[1]],da_screen[7,4][[1]]))
da_noscreen <- df %>%
  group_by(histopathology)%>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::mutate(total = sum(n))%>%
  add_row(histopathology="ESCC/HGIN",n=1,total=1)%>%
  ungroup()
da_noscreen[6,2]=(da_noscreen[2,2][[1]]+da_noscreen[3,2][[1]])
da_noscreen[6,3]=da_noscreen[2,3][[1]]
da_noscreen$ESCCseeker="noscreen"
da=rbind(da_screen,h,da_noscreen)
da=da[da$histopathology%in%c("ESCC/HGIN","ESCC","HGIN","LGIN"),]
da=da[da$histopathology%in%c("ESCC/HGIN","ESCC","HGIN"),]

da1=data.frame()
for (i in 1:nrow(da)) {
  histopathology=as.character(da$histopathology[i])
  strategy=as.character(da$ESCCseeker[i])
  data_summary=as.matrix(da[i,3:4])
  set.seed(225)
  ci=pmax(0, pmin(1, agresti_coull_ci(data_summary[1,1], data_summary[1,2]))) 
  probability=ci[1]
  lower=ci[2]
  upper=ci[3]
  da1=rbind(da1,cbind(histopathology=histopathology,strategy=strategy,probability= probability,lower=lower,upper=upper))
}

da1$histopathology=factor(da1$histopathology,levels = rev(c("ESCC/HGIN","ESCC","HGIN")))
da1[,3:5]=apply(da1[,3:5],2,as.numeric)
ggplot(da1, aes(x = probability, y = histopathology, color = strategy)) +
  geom_point(size=3.5) +  # 添加点
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +#width = 0.25,linewidth=1.35  # 添加垂直误差线
  scale_color_manual(values = c("Low_risk" = "#FDB218", "noscreen" = "black", "High_risk" ="#6B94D1")) +  # 设置颜色映射
  scale_x_continuous(limits = c(0, 0.3))+
  #geom_hline(yintercept = 39/894, linetype = "dashed", color = "grey") +#
  theme(
    panel.border = element_rect(color = "black", fill=NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position="none"
  )