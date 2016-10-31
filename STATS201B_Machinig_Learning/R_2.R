RRKR.cv(par, K, Y, iters, epsilon, m, step, lam.beg, lam.end, lam.step, k)

X = Xc_15_train_mat 
Y = Yc_15_train_mat
par = matrix(c(rep(0.005, nrow(X)), -0.01))
par = matrix(c(rep(0.005, nrow(X))))
K <- GaussianKernel(X, sqrt(29))
m = 4
epsilon = 0
lambda = 0.2
step = 0.5
lam.beg = 1
lam.end = 1
lam.step =1
k = 5
iters = 25

load("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\project\\condata_15.RData")
t = RRKR(par, lambda, K, Y, iters, epsilon, m, step)

n = nrow(X)
plot(Y, K %*% t[1:n,,drop = F] + t[n+1,])

alpha = par[1:n,,drop = F]
b = par[n+1,]

s=  RKR(K, Y, lambda)
plot(Y, K %*% s)
################################
RKR<-function(K, Y, lambda) {
  return(solve(K + lambda * diag(nrow(Y)), Y))
}
  
RKR.cv <- function(K, Y, lam.beg, lam.end, lam.step, k)
{
  n = nrow(K)
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
      Y.train = Y[-fold[[j]], ,drop = F]
      K.train = K[-fold[[j]], -fold[[j]], drop = F]
      alpha_j = RKR(K.train, Y.train, lambda[i])
      y.fitted =  K[fold[[j]], -fold[[j]], drop=F] %*% alpha_j
      err.sq = sum((Y[fold[[j]],1] - y.fitted)^2)
      cv[i] = cv[i] + err.sq
    }
  }
  cv = cv / n
  return(list(lam = lambda, c = cv))
}

K.all = GaussianKernel(rbind(Xc_15_train_mat, Xc_15_test_mat), sqrt(29))

lambda= 0.2
alpha = RKR(K.all[1:700, 1:700], Yc_15_train_mat, lambda)
y.rkr =  K.all[701:1000, 1:700, drop=F] %*% alpha
plot(Yc_15_test_mat, y.rkr, xlab="Y", ylab="fitted Y", 
     main=expression(paste("Kernel Ridge Regression ",lambda,"=0.2")), ylim = c(-1.5, 2.8))
t.test(Yc_15_test_mat-y.rkr,mu=0)


par = matrix(c(rep(0.005, nrow(X)), -0.01))
m = 3.5
epsilon = 0.05
lambda = 0.2
step = 0.5
iters = 25
t = RRKR(par, lambda, K.all[1:700, 1:700], Yc_15_train_mat, iters, epsilon, m, step)
y.rrkr = K.all[701:1000, 1:700, drop=F] %*% t[1:(nrow(t)-1),,drop = F] +t[nrow(t),]
plot(Yc_15_test_mat, y.rrkr, xlab="Y", ylab="fitted Y", 
     main=expression(paste("Robust Kernel Regression m=3.5, ", epsilon,"=0.05, ",lambda,"=0.2")), ylim = c(-1.5, 2.8))
t.test(Yc_15_test_mat - y.rrkr,mu=0)





par.train.2 = par.train[-length(par.train),,drop=F]
s = RRKR2(par.train.2, lambda, K.train, Y.train, iters, epsilon, m, step)
y.fitted.s = K[index, -(index), drop=F] %*% s
plot(Y[(index), ,drop = F], y.fitted.s)


err.sq = sum((Y[index,1] - y.fitted)^2)
cv[i] = cv[i] + err.sq