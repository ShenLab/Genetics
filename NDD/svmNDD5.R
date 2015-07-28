library("kernlab")
library("ROCR")
setwd("/Users/alex/Desktop/Columbia/Shen_Lab/chd_ndd")

annof <- data.frame(read.csv(file="annos.certain.csv", header=TRUE))

svm.loocv <- function(data, cval){
  err <- c()
  pred_out <- c()
  # select features to use
  #feats <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18) # all features
  feats <- c(5, 9, 13, 17, 18) # LoF_hmg, LoF_asd, D.Mis3_hmg, D.Mis3_asd, extracardiac
  #feats <- c(9, 17, 18)
  #feats <- c(seq(from=10, to=18, by = 1),34)
  # iterate through every row of dataset
  for(i in 1:length(annof[,1])){
    print(i) # stand-in status message
    test <- data[i,,drop=FALSE] # leave one out as test
    train <- data[-i,] # use rest of data as training
    cols <- names(train)[feats] # store feature column names
    svp <- ksvm(train$diagnosis ~ .,data=train[,cols],type="C-svc",kernel='rbfdot', kpar='automatic',C=cval, scale=FALSE) # train svm
    pred <- predict(svp,test[cols]) # predict using svm
    pred_out <- c(pred_out, pred)
    res <- test$diagnosis - pred # get error
    err <- c(err, res) # append to error vector
  }
  mse <- sum(err^2)/length(data[,1]) # mean squared error
  rmse <- sqrt(mse) # RMSE
  acc <- table(err==0)[2]/(table(err==0)[2]+table(err==0)[1])
  #return(c(rmse,acc))
  return(c(acc,rmse))
}
svm.loocv(annof, 0.05)

# test c values
cs <- seq(from=0.01, to=0.1, by=0.01)
out <- sapply(cs,function(x) svm.loocv(annof,as.numeric(x))) # save accuracy results from each c value
best.c <- cs[which.max(out)]

