library("ggplot2")
library("reshape2")
library("ROCR")
setwd("/Users/alex/Desktop/Columbia/Shen_Lab/chd_ndd")

annof <- data.frame(read.csv(file="annos.certain.csv", header=TRUE))
ntrain <- round(length(annof$Subject)*0.8)
tindex <- sample(length(annof$Subject),ntrain)
train <- annof[tindex,]
test <- annof[-tindex,]
feats <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
cols <- names(train)[feats]
fit <- glm(diagnosis ~ LoF_notch+LoF_tgfb+LoF_wnt+LoF_hmg+LoF_txf+LoF_hhe+LoF_hbe+LoF_asd+DMS_notch+DMS_tgfb+DMS_wnt+DMS_hmg+DMS_txf+DMS_hhe+DMS_hbe+DMS_asd+extracardiac, family = binomial(logit), data = train)
#fit <- glm(train$diagnosis ~ ., family = binomial(logit), data = train[,cols])

## Warning message: glm.fit: fitted probabilities numerically 0 or 1 occurred 
## Depends on training/testing split

# logit function
pred <- function(annos, weights, int) {
  t <- sum(annos*weights) + int
  mu <- 1/(1+exp(-t))
  return(mu)
}

# evaluate performance on test set
testout <- apply(test[,cols], 1, pred, weights = coef(fit)[feats], int = coef(fit)[1])
names(testout) <- test[,c("Subject")]
testlabs <- test[,c("diagnosis")]
perf <- function (cutoff, testout, labels) {
  tp <- sum(testout > cutoff & labels > cutoff)
  fp <- sum(testout > cutoff & labels < cutoff)
  tn <- sum(testout < cutoff & labels < cutoff)
  fn <- sum(testout < cutoff & labels > cutoff)
  cp <- tp + fn # condition positives
  cn <- fp + tn # condition negatives
  tpr <- tp/cp # sensitivity or TPR
  #tnr <- tn/cn # specificity or TNR
  fpr <- fp/cn # FPR
  #fnr <- fn/cp # FNR
  acc <- (tp+tn)/(cp+cn) # accuracy
  ppv <- tp/(tp+fp) # ppv
  return(c(tpr, fpr, acc, ppv)) # output list of performance metric
  #return(c(tp,fp,tn,fn))
}

fitpreds1 <- predict(fit, newdata=test, type="response")
fitpred1 <- prediction(fitpreds1, test$diagnosis)
fitperf1 <- performance(fitpred1, "tpr", "fpr")
auc1 <- performance(fitpred1, "auc")
auc1 <- unlist(slot(auc1, "y.values"))

## Generate points for ROC curve
cuts <- seq(from = 0.01, to = 0.99, by = 0.01) # vector of cutoff values from 0 to 1.0 by 0.1
all_tpr <- numeric()
all_fpr <- numeric()
all_acc <- numeric()
all_ppv <- numeric()
tmp1 <- numeric()
for (c in cuts) {
  tmp1 <- perf(c, testout, testlabs) # get list of performance metrics
  all_tpr <- c(all_tpr, tmp1[1])
  all_fpr <- c(all_fpr, tmp1[2])
  all_acc <- c(all_acc, tmp1[3])
}
df1 <- data.frame(all_tpr, all_fpr)
p1 <- ggplot(df1, aes(all_fpr, y = all_tpr, color = variable)) + geom_line(aes(y=all_tpr, col="all_tpr")) + geom_line(aes(y=all_fpr, col="gray")) + theme(legend.position="none") + ggtitle("ROC all features") + annotate("text", x = 0.9, y = 0.1,label = auc1)


eval <- perf(0.5, testout, testlabs)
names(eval) <- c("tpr/recall", "fpr", "acc", "ppv/precision")
eval
