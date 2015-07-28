library("randomForest")
library("ROCR")
setwd("/Users/alex/Desktop/Columbia/Shen_Lab/chd_ndd")

annof <- data.frame(read.csv(file="annos.certain.csv", header=TRUE))

rf.loocv <- function(data,nt){
  set.seed(1)
  err <- c()
  pred_out <- c()
  feats <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18) # all features
  #feats <- c(seq(from=10, to=18, by=1),34)
  # iterate through every row of feature vector
  for(i in 1:length(data[,1])){
    print(i) # status variable
    test <- data[i,] # hold out one entry for testing
    train <- data[-i,] # use rest of data as training
    cols <- names(train)[feats]
    fit <- randomForest(as.factor(train$diagnosis) ~ . , data = train[,cols], ntree=nt, nodesize=5, mtry=4, importance=TRUE)
    pred <- as.numeric(levels(predict(fit,test[cols])))[predict(fit,test[cols])] # convert from factor to int
    pred_out <- c(pred_out, pred)
    res <- test$diagnosis - pred # get error --> is this correct?
    err <- c(err, res) # append to error vector
  }
  mse <- sum(err^2)/length(data[,1]) # mean squared error
  rmse <- sqrt(mse) # RMSE
  acc <- table(err==0)[2]/(table(err==0)[2]+table(err==0)[1])
  return(c(acc,rmse))
}

rf.loocv(annof,23) 
nts <- seq(from=50, to=100, by=1)
out <- sapply(nts,function(x) rf.loocv(annof,as.numeric(x))) # save accuracy results from each ntree value
best.nt <- nts[which.max(out)]

# Variable Importance Plot
cols <- names(annof)[c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)]
tmp <- randomForest(as.factor(annof$diagnosis) ~ . , data = annof[,cols], ntree=21, nodesize=5, mtry=4, importance=TRUE)
varImpPlot(tmp)