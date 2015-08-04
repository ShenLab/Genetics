library("ggplot2")
#library("ROCR")
#library("AUC")
library("Epi")
setwd("/Users/alex/Dropbox/chd_ndd/del_pred/")
#e <- data.frame(read.csv(file="eval_domains_testing1_multianno.parsed2.csv", header=TRUE, na.strings='.',stringsAsFactors=FALSE))
#e <- data.frame(read.csv(file="eval_domains_testing2_multianno.parsed2.csv", header=TRUE, na.strings='.',stringsAsFactors=FALSE))
#e <- data.frame(read.csv(file="eval_domains_testing3_multianno.parsed2.csv", header=TRUE, na.strings='.',stringsAsFactors=FALSE))
#e <- data.frame(read.csv(file="non_domains_testing1_multianno.parsed2.csv", header=TRUE, na.strings='.',stringsAsFactors=FALSE))
#e <- data.frame(read.csv(file="non_domains_testing2_multianno.parsed2.csv", header=TRUE, na.strings='.',stringsAsFactors=FALSE))
e <- data.frame(read.csv(file="non_domains_testing3_multianno.parsed2.csv", header=TRUE, na.strings='.',stringsAsFactors=FALSE))
f <- e[complete.cases(e),]
perf <- function (cutoff, score, labels) {
  tp <- sum(score >= cutoff & labels == 1, na.rm=TRUE)
  fp <- sum(score >= cutoff & labels == 0, na.rm=TRUE)
  tn <- sum(score <= cutoff & labels == 0, na.rm=TRUE)
  fn <- sum(score <= cutoff & labels == 1, na.rm=TRUE)
  cp <- tp + fn # condition positives
  cn <- fp + tn # condition negatives
  tpr <- tp/cp # sensitivity or TPR
  tnr <- tn/cn # specificity or TNR
  fpr <- fp/cn # FPR
  fnr <- fn/cp # FNR
  acc <- (tp+tn)/(cp+cn) # accuracy
  ppv <- tp/(tp+fp) # ppv
  sens <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  return(c(tpr, fpr, acc, ppv, sens, spec)) # output list of performance metric
}

sift <- f[,6]
pp2hdiv <- f[,7]
pp2hvar <- f[,8]
radialsvm <- f[,9]
lr <- f[,10]
cadd <- f[,11]
lab <- f[,12]

#sum(sift>0.01 & lab>0.01, na.rm=TRUE)
#perf(0.01, sift, lab)

makecuts <- function(v){
  min <- min(na.omit(v))
  max <- max(na.omit(v))
  return(seq(from=min, to=max, by=(max-min)/100))
}

## Generate points for ROC curve
#cuts <- seq(from = 0.01, to = 0.99, by = 0.01) # vector of cutoff values from 0.01 to 0.99 by 0.01
# sift
sift_tpr <- numeric()
sift_fpr <- numeric()
tmpsift <- numeric()
siftcuts <- makecuts(sift)
for (c in siftcuts) {
  tmpsift <- perf(c, sift, lab) # get list of performance metrics
  sift_tpr <- c(sift_tpr, tmpsift[1])
  sift_fpr <- c(sift_fpr, tmpsift[2])
}
dfsift <- data.frame(sift_tpr, sift_fpr)
siftp <- ggplot(dfsift, aes(sift_fpr, y = sift_tpr, color = variable)) + 
      geom_line(aes(y=sift_tpr, col="sift_tpr")) + 
      geom_line(aes(y=sift_fpr, col="grey")) + 
      theme(legend.position="none") + 
      xlim(0,1) + 
      ylim(0,1) + 
      ggtitle("ROC SIFT")

# pp2hdiv
pp2hdiv_tpr <- numeric()
pp2hdiv_fpr <- numeric()
tmppp2hdiv <- numeric()
pp2hdivcuts <- makecuts(pp2hdiv)
for (c in pp2hdivcuts) {
  tmppp2hdiv <- perf(c, pp2hdiv, lab) # get list of performance metrics
  pp2hdiv_tpr <- c(pp2hdiv_tpr, tmppp2hdiv[1])
  pp2hdiv_fpr <- c(pp2hdiv_fpr, tmppp2hdiv[2])
}
dfpp2hdiv <- data.frame(pp2hdiv_tpr, pp2hdiv_fpr)
pp2hdivp <- ggplot(dfpp2hdiv, aes(pp2hdiv_fpr, y = pp2hdiv_tpr, color = variable)) + 
  geom_line(aes(y=pp2hdiv_tpr, col="pp2hdiv_tpr")) + 
  geom_line(aes(y=pp2hdiv_fpr, col="grey")) + 
  theme(legend.position="none") + 
  xlim(0,1) + 
  ylim(0,1) + 
  ggtitle("ROC pp2hdiv")

# pp2hvar
pp2hvar_tpr <- numeric()
pp2hvar_fpr <- numeric()
tmppp2hvar <- numeric()
pp2hvarcuts <- makecuts(pp2hvar)
for (c in pp2hvarcuts) {
  tmppp2hvar <- perf(c, pp2hvar, lab) # get list of performance metrics
  pp2hvar_tpr <- c(pp2hvar_tpr, tmppp2hvar[1])
  pp2hvar_fpr <- c(pp2hvar_fpr, tmppp2hvar[2])
}
dfpp2hvar <- data.frame(pp2hvar_tpr, pp2hvar_fpr)
pp2hvarp <- ggplot(dfpp2hvar, aes(pp2hvar_fpr, y = pp2hvar_tpr, color = variable)) + 
  geom_line(aes(y=pp2hvar_tpr, col="pp2hvar_tpr")) + 
  geom_line(aes(y=pp2hvar_fpr, col="grey")) + 
  theme(legend.position="none") + 
  xlim(0,1) + 
  ylim(0,1) + 
  ggtitle("ROC pp2hvar")

# radialsvm
radialsvm_tpr <- numeric()
radialsvm_fpr <- numeric()
tmpradialsvm <- numeric()
radialsvmcuts <- makecuts(radialsvm)
for (c in radialsvmcuts) {
  tmpradialsvm <- perf(c, radialsvm, lab) # get list of performance metrics
  radialsvm_tpr <- c(radialsvm_tpr, tmpradialsvm[1])
  radialsvm_fpr <- c(radialsvm_fpr, tmpradialsvm[2])
}
dfradialsvm <- data.frame(radialsvm_tpr, radialsvm_fpr)
radialsvmp <- ggplot(dfradialsvm, aes(radialsvm_fpr, y = radialsvm_tpr, color = variable)) + 
  geom_line(aes(y=radialsvm_tpr, col="radialsvm_tpr")) + 
  geom_line(aes(y=radialsvm_fpr, col="grey")) + 
  theme(legend.position="none") + 
  xlim(0,1) + 
  ylim(0,1) + 
  ggtitle("ROC radialsvm")

# lr
lr_tpr <- numeric()
lr_fpr <- numeric()
tmplr <- numeric()
lrcuts <- makecuts(lr)
for (c in lrcuts) {
  tmplr <- perf(c, lr, lab) # get list of performance metrics
  lr_tpr <- c(lr_tpr, tmplr[1])
  lr_fpr <- c(lr_fpr, tmplr[2])
}
dflr <- data.frame(lr_tpr, lr_fpr)
lrp <- ggplot(dflr, aes(lr_fpr, y = lr_tpr, color = variable)) + 
  geom_line(aes(y=lr_tpr, col="lr_tpr")) + 
  geom_line(aes(y=lr_fpr, col="grey")) + 
  theme(legend.position="none") + 
  xlim(0,1) + 
  ylim(0,1) + 
  ggtitle("ROC lr")

# cadd
cadd_tpr <- numeric()
cadd_fpr <- numeric()
tmpcadd <- numeric()
#caddscaled <- cadd/40 # fix it so cadd scores on same scale as cuts
caddcuts <- makecuts(cadd)
for (c in caddcuts) {
  #tmpcadd <- perf(c, caddscaled, lab) # get list of performance metrics
  tmpcadd <- perf(c, cadd, lab)
  cadd_tpr <- c(cadd_tpr, tmpcadd[1])
  cadd_fpr <- c(cadd_fpr, tmpcadd[2])
}
dfcadd <- data.frame(cadd_tpr, cadd_fpr)
caddp <- ggplot(dfcadd, aes(cadd_fpr, y = cadd_tpr, color = variable)) + 
  geom_line(aes(y=cadd_tpr, col="cadd_tpr")) + 
  geom_line(aes(y=cadd_fpr, col="grey")) + 
  theme(legend.position="none") + 
  xlim(0,1) + 
  ylim(0,1) + 
  ggtitle("ROC cadd")

#calculate AUC for all tools
aucs <- c(ROC(sift, lab)$AUC, 
          ROC(pp2hdiv, lab)$AUC, 
          ROC(pp2hvar, lab)$AUC, 
          ROC(radialsvm, lab)$AUC,
          ROC(lr, lab)$AUC,
          ROC(cadd, lab)$AUC)
names(aucs) <- c("sift", "pp2hdiv", "pp2hvar", "radialsvm", "lr", "cadd")
aucs

#plot of all curves
ggplot() +
  geom_line(aes(x=cadd_fpr, y=cadd_tpr, col="cadd_tpr")) +
  geom_line(aes(x=lr_fpr, y=lr_tpr, col="lr_tpr")) +
  geom_line(aes(x=radialsvm_fpr, y=radialsvm_tpr, col="radialsvm_tpr")) +
  geom_line(aes(x=pp2hvar_fpr, y=pp2hvar_tpr, col="pp2hvar_tpr")) + 
  geom_line(aes(x=pp2hdiv_fpr, y=pp2hdiv_tpr, col="pp2hdiv_tpr")) +
  geom_line(aes(x=sift_fpr, y=sift_tpr, col="sift_tpr")) +
  geom_line(aes(x=cadd_fpr, y=cadd_fpr, col="y=x")) + 
  theme(legend.position=c(0.8, 0.3)) +
  xlim(0,1) + 
  ylim(0,1) +
  ggtitle("ROC all tools")

length(f[,1])
