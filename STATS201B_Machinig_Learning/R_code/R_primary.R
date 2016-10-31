library(Rcpp)
sourceCpp("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\project\\kernel.cpp")

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

robust.loss<-function(x, m, epsilon)
{
  S1.index = which(x > (epsilon + m))
  S2.index = which(x > epsilon & x <= (epsilon + m))
  S3.index = which(x >= -epsilon & x <= epsilon )
  S4.index = which(x >= -(epsilon + m) & x < -epsilon)
  S5.index = which(x < -(epsilon + m))
  y = rep(0, length(x))
  y[S1.index] = m*(2*(x[S1.index]-epsilon) -m)
  y[S2.index] = x[S2.index]^2
  y[S4.index] = x[S4.index]^2
  y[S5.index] = -m*(2*(x[S5.index]+epsilon)+m)
  return(y)
}
x = seq(-5, 5, 0.1)
y = robust.loss(x, 2, 0)
plot(x, y, type="l", main="Huber Loss Function", ylim = c(0, 20))
points(x, x^2, type="l", col="red")

RRKR <- function(par, lambda, K, Y, iters, epsilon, m, step)
{
  for(i in 1:iters)
  {
    n = nrow(Y)
    alpha = par[1:n,,drop = F]
    b = par[n+1,]
    u = Y - K %*% alpha - b 
    S1.index = which(u > (epsilon + m))
    S2.index = which(u > epsilon & u <= (epsilon + m))
    S3.index = which(u >= -epsilon & u <= epsilon )
    S4.index = which(u >= -(epsilon + m) & u < -epsilon)
    S5.index = which(u < -(epsilon + m))
    #I I0 vec1
    temp = rep(1,n)
    I = diag(temp)
    temp[S1.index] = 0
    temp[S3.index] = 0
    temp[S5.index] = 0
    I0 = diag(temp)
    vec1 = matrix(rep(1,n))
    #e
    e = rep(0, n)
    e[S1.index] = m
    e[S2.index] = epsilon
    e[S4.index] = -epsilon
    e[S5.index] = -m
    e = matrix(e)
    
    sol.mat = rbind(cbind(lambda * I + I0 %*% K, I0 %*% vec1), cbind(t(vec1), 0))
    q = I0 %*% (b * vec1 - Y) + m * e
    par = rbind((1-step)*alpha, b) - step *  solve(sol.mat, rbind(q, 0))
  }
  return(par)
}






RRKR.cv <- function(par, K, Y, 
                    iters, epsilon, m, step, lam.beg, lam.end, lam.step, k)
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
      Y.train = Y[-fold[[j]], ,drop = F]
      par_j = RRKR(par.train, lambda[i], K.train, 
                   Y.train, iters, epsilon, m, step)
      alpha_j = par_j[1:nrow(X.train),,drop=F]
      b_j = par_j[nrow(X.train)+1,]
      y.fitted =  b_j + K[fold[[j]], -fold[[j]], drop=F] %*% alpha_j
      err.sq = sum((Y[fold[[j]],1] - y.fitted)^2)
      cv[i] = cv[i] + err.sq
    }
  }
  cv = cv / n
  return(list(lam = lambda, c = cv))
}

#########################################

pred.ker.log <-function(Y.test, K, alpha)
{
  p.fit = 1 / (1+exp(-K%*%alpha))
  accuracy = sum((p.fit>0.5) == Y.test)
  return(list(p.fit = p.fit, accuracy = accuracy))
}

###########################################




############X_15########################
lambda = 0.2
par.train = matrix(0.001, nrow = nrow(X_15_train_matrix), ncol = 1)
K.train = GaussianKernel(X_15_train_matrix, 6.08)
iters = 30
alpha = KLR(par.train, lambda, K.train, X_15_train_matrix, Y_15_train_matrix, iters)
K1 = GaussianKernel(rbind(X_15_train_matrix, as.matrix(X_15_test)), 6.08)
Y_15_test_matrix = as.matrix(as.numeric(as.character(Y_15_test)))
pred.ker.log(Y_15_test_matrix , K1[701:1000,1:700], alpha)
###########################################

lambda = 0.2
par.train = matrix(0.001, nrow = nrow(X_train_matrix), ncol = 1)
K.train = GaussianKernel(X_train_matrix, 7.87)
iters = 30
alpha = KLR(par.train, lambda, K.train, X_train_matrix, Y_train_matrix, iters)
K1 = GaussianKernel(rbind(X_train_matrix, as.matrix(X_test)), 7.87)
Y_test_matrix = as.matrix(as.numeric(as.character(Y_test)))
temp = pred.ker.log(Y_test_matrix , K1[701:1000,1:700], alpha)
y.fit = as.numeric(temp$p.fit > 0.5)
y = as.numeric(Y_test_matrix)
table(y, y.fit)
###########################################

lambda = 0.2
par.train = matrix(0.001, nrow = nrow(X_10_train_matrix), ncol = 1)
K.train = GaussianKernel(X_10_train_matrix, 3.87)
iters = 30
alpha = KLR(par.train, lambda, K.train, X_10_train_matrix, Y_10_train_matrix, iters)
K1 = GaussianKernel(rbind(X_10_train_matrix, as.matrix(X_10_test)), 3.87)
Y_10_test_matrix = as.matrix(as.numeric(as.character(Y_10_test)))
pred.ker.log(Y_10_test_matrix , K1[701:1000,1:700], alpha)

###########################################

X = matrix(rnorm(200 * 2), ncol = 2)
X = cbind(X, rbinom(200, size = 1, prob = 0.7))
beta = matrix(c(1, -2, 3))
epsilon = matrix(rnorm(200, sd = 0.1), ncol = 1)
p = 1 / (1 + exp(-(X %*% beta + epsilon)))
Y =  as.matrix(as.numeric((p > 0.5)))
index.0 =  which(Y==0)[1:30]
index.1 = which(Y==1)
X1 = X[c(index.0, index.1),]
Y1 = Y[c(index.0, index.1),, drop = F]
K = GaussianKernel(X1, 1)
par = as.matrix(rnorm(nrow(X1), sd = 0.1))
alpha = KLR(par, 1, K, X1, Y1, 50)




Y2 =as.factor(Y1)
data = data.frame(X1, Y2)


X.test = matrix(rnorm(2*200), ncol = 2)
X.test = cbind(X.test, rbinom(200, size = 1, prob = 0.3))
epsilon.test = matrix(rnorm(200, sd = 0.1), ncol = 1)
p.test = 1 / (1 + exp(-(X %*% beta + epsilon)))
Y.test =  as.numeric((p > 0.5))

K2 = GaussianKernel(rbind(X1, X.test), 1)
p.fit = 1 / (1 + exp(-K2[-(1:nrow(X1)), 1:nrow(X1)] %*% alpha))
y.fit = as.numeric(p.fit > 0.5)
table(Y.test, y.fit)


library("e1071")
svm_model <- svm(Y2~., data=data, scale=FALSE)
svmpred <- predict(svm_model,X.test)
Y.fit = as.numeric(as.character(svmpred))
table(Y.test, Y.fit)





############################################################

K <- GaussianKernel(X, 3)
par = matrix(rep(0.01, nrow(X)))
lambda = 20
a = KLR(alpha, lambda, K, X, Y, iters)
p.fitted <- 1/(1 + exp(-K %*% a))
plot(p, p.fitted)
cor(p, p.fitted)


X = cbind(rnorm(50), matrix(rnorm(50)))
p = 1 / (1 + exp(X %*% matrix(c(1,2), ncol=1)))
Y = matrix(0, nrow=nrow(p), ncol=1)
for(i in 1:nrow(p))
  Y[i,1] = rbinom(1, 1, prob=p[i,1])


iters = 1000

X = cbind(rnorm(50), matrix(rnorm(50)))
Y =1 + X %*% matrix(c(1,2), ncol=1) + matrix(rnorm(nrow(X), sd=0.1))
K <- GaussianKernel(X, 3)
par = matrix(rnorm(nrow(X)+1))
lambda = 2
K <- GaussianKernel(X, 2)
iters = 1000
epsilon = 0.1
m = 10
step = 0.1
par = RRKR(par, lambda, K, X, Y, iters, epsilon, m, step)
alpha = par[1:nrow(X),,drop=F]
b = par[1+nrow(X),]
Y.fitted = K %*% alpha + b
plot(Y, Y.fitted)
