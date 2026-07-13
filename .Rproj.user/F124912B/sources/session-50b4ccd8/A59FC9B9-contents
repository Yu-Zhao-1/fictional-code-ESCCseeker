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
library(ROCR)
library(reportROC)
library(shapviz)
library(reshape2)
library(ggprism)
library(rms)
library(rmda)
library(PRROC)
library(boot)
library(tidyverse)
library(ggpubr)
###1. 9-Methy model####
load("Model_construction.RData")
train$group=ifelse(train$group%in%c("ESCC","HGIN"),1,0)
trainy=train$group
trainx=as.matrix(train[,gene0])
trainy=ifelse(trainy==1,"T","N")
set.seed(120)
cntrl = trainControl(
  method = "cv",
  number = 5,
  verboseIter = F,
  search = "random"
)
set.seed(120)
tuned_model = train(
  x = trainx,
  y = trainy,
  trControl = cntrl,
  method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) 
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
Methy9 <- xgboost(data = trainx, label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
y_predicted.train <- predict(Methy9, trainx)
yp1=data.frame(class=trainy,Methy9=y_predicted.train)
testx=as.matrix(test[,gene0])
testy=ifelse(test$group=="ESCC",2,ifelse(test$group=="HGIN",1,0))
verificationy=ifelse(verification$group=="ESCC",2,ifelse(verification$group=="HGIN",1,0))
verificationx=as.matrix(verification[,gene0])
y_predicted.test <- predict(Methy9, testx)
yp2=data.frame(class=testy,Methy9=y_predicted.test)
y_predicted.verification <- predict(Methy9, verificationx)
yp3=data.frame(class=verificationy,Methy9=y_predicted.verification)
###2. ESCCseeker ####
##Machine learning####
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
library(xgboost)
## 1.LR###
data=train[,c("group",gene1)]
data$group=factor(data$group)
trainx=as.matrix(data[,gene1])
trainy=as.factor(data[,1])
LR<-glm(group~., family=binomial, data=data)
y_predicted.train=predict(LR,train,type = "response")
ML=data.frame(class=trainy,LR=y_predicted.train)
## 2.LightGBM##
bt_wf <- workflow() %>%
  add_model(
    boost_tree(
      mode = "classification", 
      mtry = tune(),
      tree_depth = tune(),
      trees = tune(),
      min_n = tune()
    ) %>%
      set_engine("lightgbm", objective = "binary") 
  ) %>%
  add_formula(group ~ .)
set.seed(133)
tree_grid <- grid_space_filling(mtry(range = c(as.integer(ncol(trainx)/4), as.integer(2*ncol(trainx)/3))),
                                tree_depth(range = c(1L,3L)),
                                trees(range = c(20,50)), 
                                min_n(range = c(20,50)),
                                size= 10)
set.seed(133)
bt_folds=vfold_cv(data,v=5)
bt_tune  <- bt_wf %>%
  tune_grid(resamples = bt_folds,
            grid = tree_grid,
            metrics = metric_set(yardstick::accuracy, 
                                 yardstick::roc_auc, 
                                 yardstick::pr_auc),
            control = control_grid(save_pred = T, verbose = F))
bt_best <-select_best(bt_tune,metric="roc_auc")
lgbm <- bt_wf %>%
  finalize_workflow(bt_best) %>%
  fit(data)
ML$LGBM <- predict(lgbm,new_data=data,type="prob")$.pred_1
##3.SVM###
trainx=as.matrix(data[,gene1])
trainy=factor(data$group)
param_grid <- expand.grid(
  sigma = c(0.001,0.01,0.1, 1),
  C = c(0.01,0.1,1,10))
ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
set.seed(1)
tuned_model <- train(
  x = trainx,
  y = trainy,
  method = "svmRadial",
  tuneGrid = param_grid,
  trControl = ctrl)
sigma=as.numeric(tuned_model$bestTune[1])
C=as.numeric(tuned_model$bestTune[2])
svm <- svm(x = trainx, y = trainy,sigma=sigma,C=C,probability = T)
ML$SVM <- attr(predict(svm, trainx,probability = TRUE),"probabilities" )[,1]
## 4.XG###
trainx=as.matrix(data[,gene1])
trainy=data$group
trainy=ifelse(trainy==1,"T","N")
set.seed(16)
cntrl = trainControl(
  method = "cv",
  number = 5,
  verboseIter = F,
  search = "random"
)
set.seed(16)
tuned_model = train(
  x = trainx,
  y = trainy,
  trControl = cntrl,
  method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) 
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
xgb <- xgboost(data = trainx, label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
ML$Xgboost <- predict(xgb, trainx)
yp=ML#
roc3<-roc(class~LR,yp)
roc6<-roc(class~LGBM,yp)
roc5<-roc(class~SVM,yp)
roc4<-roc(class~Xgboost,yp)
ci3=ci(roc3)
ci4=ci(roc4)
ci5=ci(roc5)
ci6=ci(roc6)
plot(roc3, col = "#EEA210", main = "ML ROC",lwd=1,lty=1,cex.lab = 1.2) 
lines(roc4, col = "#B90019",lwd=1,lty=1) 
lines(roc5, col = "#84C85A",lwd=1,lty=1) 
lines(roc6, col = "#0055A5",lwd=1,lty=1) 
legend("bottomright",
       legend = c(
         paste0("LR:", round(ci3[2], 3), paste0("; 95% CI:", round(ci3[1],3), "-", round(ci3[3],3))), 
         paste0("XGboost:", round(ci4[2], 3), paste0("; 95% CI:", round(ci4[1],3), "-", round(ci4[3],3))), 
         paste0("SVM:", round(ci5[2], 3), paste0("; 95% CI:", round(ci5[1],3), "-", round(ci5[3],3))),
         paste0("LightGBM:", round(ci6[2], 3), paste0("; 95% CI:", round(ci6[1],3), "-", round(ci4[3],3)))), 
       col = c("#EEA210","#B90019","#84C85A","#0055A5"), 
       lty = 1,lwd=1,cex = 0.9,bty = "n")
