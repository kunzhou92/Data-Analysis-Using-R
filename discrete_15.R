introducedummy<-function(x,b){
  class <- as.character(sort(unique(x[,b])))
  length <- length(class)
  n<-nrow(x)
  for (k in 1: length){ 
    p<-colnames(x)[b]
    p2<-paste(p,class[k],sep="_")
    new1<-numeric(n)
    state<-x[,b]
    interest<-class[k]
    for(i in 1:n){
      if(state[i]==interest){
        new1[i]=1
      }
      else{
        new1[i]=0
      }
    }
    x$p2<-new1
    names(x)[ncol(x)] <- p2
  }
  return(x)
}
data_15 <- read.table("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\project\\german.data.txt")
colnames(data_15) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","Output")
head(data_15)
data_15 <- data_15[,-c(4,9,10,15,17)]
data_15 <- introducedummy(data_15,1)
data_15 <- introducedummy(data_15,3)
data_15 <- introducedummy(data_15,5)
data_15 <- introducedummy(data_15,6)
data_15 <- introducedummy(data_15,9)
data_15 <- introducedummy(data_15,11)
data_15 <- introducedummy(data_15,14)
data_15 <- introducedummy(data_15,15)
data_15 <- data_15[,-c(1,3,5,6,9,11,14,15)]
data_15 <- cbind(data_15[,c(9:12,1,13:17,2,18:27,3,4,28:31,5,32:34,6,7,35:38,8)])
data_15$Output[data_15$Output==2] <- 0 # 1 is good 0 is bad
data_15$Output <- as.factor(data_15$Output)

set.seed(21)
indexes = sample(1:nrow(data_15), size=0.7*nrow(data_15))
data_15.train = data_15[indexes,]
data_15.test = data_15[-indexes,]
for (j in c(5,11,22,23,28,32,33))
  data_15.test[,j]=(data_15.test[,j]-mean(data_15.train[,j]))/sd(data_15.train[,j])
for (j in c(5,11,22,23,28,32,33))
  data_15.train[,j] = scale(data_15.train[,j])
X_15_train <- data_15.train[,-38]
Y_15_train <- data_15.train[,38]
X_15_test <- data_15.test[,-38]
Y_15_test <- data_15.test[,38]

X_15_train_mat <- as.matrix(X_15_train)
Y_15_train_mat <- as.matrix(as.numeric(as.character(Y_15_train)))
X_15_test_mat <- as.matrix(X_15_test)
Y_15_test_mat <- as.matrix(as.numeric(as.character(Y_15_test)))

library(Rcpp)
sourceCpp("kernel.cpp")

KLR <- function(par, lambda, K, X, Y, iters){
  for(i in 1:iters){
    p <- 1/(1 + exp(-K %*% par))
    w <- p * (1-p)  
    W <- diag(as.vector(w),nrow(X),nrow(X))
    W_s <- solve(W)
    S <- K %*% par + W_s %*% (Y-p)
    par <- solve(K + lambda*W_s, S) 
  }
  return (par)
}

KLR.cv <-function(par, K, X, Y, iters, lam.beg, lam.end, lam.step, k)
{
  n = nrow(X)
  lambda = seq(from = lam.beg, to = lam.end, by = lam.step)
  cv = rep(0, length(lambda))
  #divide data into k folds
  select.num = n %/% k
  index = 1:n
  fold = list()
  for(i in 1:(k-1))
  {
    num.within = sample(index, select.num, replace = F)
    fold = c(fold, list(num.within))
    index = setdiff(index, num.within)
  }
  fold = c(fold, list(index))
  for(i in (1:length(lambda)))
  {
    for(j in 1:k)
    {
      par.train = par[-fold[[j]], , drop = F]
      K.train = K[-fold[[j]], -fold[[j]], drop = F]
      X.train = X[-fold[[j]], ,drop = F]
      Y.train = Y[-fold[[j]], ,drop = F]
      alpha_j = KLR(par.train, lambda[i], K.train, X.train, Y.train, iters)
      p.fitted = 1/ (1+exp(-K[fold[[j]], -fold[[j]], drop = F] %*% alpha_j))
      err.sq = sum((Y[fold[[j]],1] - as.numeric(p.fitted>0.5))^2)
      cv[i] = cv[i] + err.sq
    }
  }
  cv = cv / n
  return(list(lam = lambda, c = cv))
}
library("caret")
iters=30
K_15 <- GaussianKernel(X_15_train_mat, sqrt(37))
alpha = matrix(rep(0.001, nrow(X_15_train_mat)))
result_15 = KLR.cv(alpha, K_15, X_15_train_mat, Y_15_train_mat, iters, 0.05, 1, 0.025, 5)
#result_15_1 = KLR.cv(alpha, K_15, X_15_train_mat, Y_15_train_mat, iters, 0.05, 0.3, 0.025, 5)
lambda= 0.8
alpha = KLR(alpha, lambda, K_15, X_15_train_mat, Y_15_train_mat, iters)
K1 = GaussianKernel(rbind(X_15_train_mat, X_15_test_mat), sqrt(37))
p.fit = 1 / (1+exp(-K1[701:1000,1:700]%*%alpha))
y.fit = as.matrix(as.numeric(p.fit>0.5))
confusionMatrix(y.fit,Y_15_test_mat)
plot(result_15$lam, result_15$c, xlab="lambda", ylab="MSE", 
     main="cross validation")


#logist
glmmodel <- glm(Output~.-1,data=data_15.train,family=binomial(link = "logit"))
glmpredit <- predict(glmmodel,X_15_test,na.action=na.pass)
glmresult <- confusionMatrix(y.fit,Y_15_test_mat)
par[which(is.na(par))] = 0




library("e1071")
svm_model <- svm("vote"~.,data=wordMatrix)
tune.out.radial <- tune(svm, Output~.,data=data_15.train, scale=FALSE)
best.svm.radial <- tune.out.radial$best.model
class.test.radial=predict(best.svm.radial, X_15_test)
svmresult_15 <- confusionMatrix(class.test.radial,Y_15_test)
svmresult_15
