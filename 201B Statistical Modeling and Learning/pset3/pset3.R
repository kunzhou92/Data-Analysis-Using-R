brader = read.csv("E:\\brader.csv")
#logit link function
logit<-function(X) #X is a matrix
{
   return(1 / (1 + exp(-X)))
}
#log-likelihood function
l.ordlogit.4<-function(param, X, Y) 
{
  X = as.matrix(X)
  parNum = ncol(X)
  Y = matrix(Y, ncol = 1)
  param = matrix(param, ncol = 1)
  beta = param[1:parNum,,drop = FALSE]
  psi = param[(parNum+1):nrow(param),,drop = FALSE]
  p4 = log(1-logit(psi[3]-X[Y[,1]==4,]%*%beta))
  p3 = log(logit(psi[3]-X[Y[,1]==3,]%*%beta) - logit(psi[2]-X[Y[,1]==3,]%*%beta))
  p2 = log(logit(psi[2]-X[Y[,1]==2,]%*%beta) - logit(psi[1]-X[Y[,1]==2,]%*%beta))  
  p1 = log(logit(psi[1]-X[Y[,1]==1,]%*%beta))
  penalty = 10^6*as.numeric(any(psi[2]>psi[3], psi[1]>psi[2]))
  
  l = sum(p1) + sum(p2)+ sum(p3)+sum(p4) - penalty
  return(l)
}
#process data
X = as.matrix(brader[,-1])
Y = as.matrix(brader[,1])
param.1 = c(0.74, 0.17 , 0.010, 0.004,-1.6, 0.131, 1.3)
#calculate optimum value
orlogit = optim(par=param.1, fn=l.ordlogit.4, X=X, Y=Y, method="BFGS",
                control=list(fnscale=-1), hessian=TRUE)
orlogit$par
diag(-solve(orlogit$hessian))
brader.df = data.frame(immigr=as.factor(brader[,1]), tone=brader[,2], 
                       eth=brader[,3], ppage=brader[,4], ppincimp=brader[,5])
library(MASS)
polr(immigr~tone+eth+ppage+ppincimp,data = brader.df, 
     method = "logistic", Hess = TRUE)


l.ordlogit.4.beta<-function(beta, X, Y, psi) 
{
  beta = matrix(beta, ncol=1)
  X = as.matrix(X)
  Y = matrix(Y, ncol = 1)
  beta.1 = beta[1:4,,drop = FALSE]
  beta.2 = beta[5:8,,drop = FALSE]
  beta.3 = beta[9:12,,drop = FALSE]
  p4 = log(1-logit(psi[3]-X[Y[,1]==4,]%*%beta.3))
  p3 = log(logit(psi[3]-X[Y[,1]==3,]%*%beta.3) 
           - logit(psi[2]-X[Y[,1]==3,]%*%beta.2))
  p2 = log(logit(psi[2]-X[Y[,1]==2,]%*%beta.2) 
           - logit(psi[1]-X[Y[,1]==2,]%*%beta.1))  
  p1 = log(logit(psi[1]-X[Y[,1]==1,]%*%beta.1))
  penalty = 10^6*as.numeric(any(psi[2]>psi[3], psi[1]>psi[2]))
  
  l = sum(p1) + sum(p2)+ sum(p3)+sum(p4) - penalty
  return(l)
}
l.ordlogit.4.psi<-function(psi, X, Y, beta) 
{
  X = as.matrix(X)
  Y = matrix(Y, ncol = 1)
  beta.1 = beta[1:4,,drop = FALSE]
  beta.2 = beta[5:8,,drop = FALSE]
  beta.3 = beta[9:12,,drop = FALSE]
  p4 = log(1-logit(psi[3]-X[Y[,1]==4,]%*%beta.3))
  p3 = log(logit(psi[3]-X[Y[,1]==3,]%*%beta.3) 
           - logit(psi[2]-X[Y[,1]==3,]%*%beta.2))
  p2 = log(logit(psi[2]-X[Y[,1]==2,]%*%beta.2) 
           - logit(psi[1]-X[Y[,1]==2,]%*%beta.1))  
  p1 = log(logit(psi[1]-X[Y[,1]==1,]%*%beta.1))
  penalty = 10^6*as.numeric(any(psi[2]>psi[3], psi[1]>psi[2]))
  
  l = sum(p1) + sum(p2)+ sum(p3)+sum(p4) - penalty
  return(l)
}

beta = c(0.892194514, -0.092005636,  0.011226191,  0.004547339,  0.815988147,
        0.112129601, 0.006075800,  0.020328998,  0.643340882,  0.455857262,
        0.009901077, -0.006459265)
psi = c(-1.6657782,  0.1511321,  1.3543657)
for(i in 1:100)
{

  beta.optim = optim(par=beta, fn=l.ordlogit.4.beta, X=X, Y=Y, psi=psi,
                     method="BFGS", control=list(fnscale=-1), hessian = T)
  beta = beta.optim$par
  psi.optim = optim(par=psi, fn=l.ordlogit.4.beta, X=X, Y=Y, beta=beta,
                    method="BFGS", control=list(fnscale=-1), hessian = T)
  psi = psi.optim$par
}
beta.optim$par
diag(-solve(beta.optim$hessian))
psi.optim$par
diag(-solve(psi.optim$hessian))

###########################
nielsenaid = read.csv("E:\\nielsenaid.csv")
l.tobit<-function(param, X, Y)
{
  X = as.matrix(X)
  Y = as.matrix(Y)
  param = as.matrix(param)
  beta = param[-nrow(param),,drop=FALSE]
  std = param[nrow(param),,drop=TRUE]
  l = sum(log(dnorm(X[Y[,1]>0, ]%*%beta/std))) + 
    sum(log(1-pnorm(X[Y[,1]==0,]%*%beta/std)))
  return(l)
}
X = as.matrix(nielsenaid[,3:8])
Y = as.matrix(nielsenaid[,2])
param = c(rnorm(6), 10)
op =  optim(par=param, fn=l.tobit, X=X, Y=Y, method="Nelder-Mead",
      control=list(fnscale=-1), hessian=TRUE)


itera = 200
mean.pd = rep(NaN, itera)
X.pd = matrix(c(colMeans(X)[-6], 0), nrow=1)
i=0
N = nrow(X)
while(i <= itera)
{
  index = sample(N, N, replace = T)
  X.resample = X[index,]
  Y.resample = Y[index,]
  param = c(rnorm(6), 10)
  par.re = optim(par=param, fn=l.tobit, X=X, Y=Y, method="BFGS",
                 control=list(fnscale=-1))$par
  beta.re = as.matrix(par.re[1:6])
  std.re = par.re[7]
  mean.pd[i] = pnorm(X.pd%*%beta.re/std.re) * (X.pd%*%beta.re) + 
    std.re*dnorm(X.pd%*%beta.re/std.re)
  i=i+1   
}
mean(mean.pd)
dens = density(mean.pd, bw = 0.225)
plot(dens, main = "bootstrap")
q_25 = quantile(mean.pd, 0.025)
q_75 = quantile(mean.pd, 0.975)
x1 = min(which(dens$x >= q_25))
x2 = max(which(dens$x < q_75))
with(dens, polygon(x=c(x[c(x1, x1:x2, x2)]), y=c(0, y[x1:x2], 0), col = "gray"))

#######################

# X is a n*1 matrix. 
poly.cv <-function(X, Y, degree, k)
{
  #construct the polynomial matirx
  X.matrix = matrix(1, nrow = nrow(X), ncol = 1)
  if(degree > 0)
  {
    for(i in 1:degree)
      X.matrix = cbind(X.matrix, X^i)
  }
  n = nrow(X)
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
  error = rep(NaN, k)
  #calculate the test error for each fold
  for(i in 1:k)
  {
    X.train = X.matrix[-fold[[i]], , drop = F]
    Y.train = Y[-fold[[i]], , drop = F]
    beta = solve(crossprod(X.train), crossprod(X.train, Y.train))
    X.test = X.matrix[fold[[i]], , drop = F]
    Y.test = Y[fold[[i]], , drop = F]
    res = Y.test - X.test %*% beta
    sse = t(res) %*% res
    error[i] = sse[1,1] / (nrow(X.test))
  }
  return(mean(error))
}
#######

#calculate the in-sample mse.
poly.mse <-function(X, Y, degree)
{
  X.matrix = matrix(1, nrow = nrow(X), ncol = 1)
  if(degree > 0)
  {
    for(i in 1:degree)
      X.matrix = cbind(X.matrix, X^i)
  }
  beta = solve(crossprod(X.matrix), crossprod(X.matrix, Y))
  res = Y - X.matrix %*% beta
  sse = t(res) %*% res / (nrow(X) - degree -1)
  return(sse[1,1])
}

X = as.matrix(runif(100, min=-4, max=4))
epsilon = as.matrix(rnorm(100))
Y.1 = -2*(X<(-3)) + 2.55*(X>(-2)) - 2*(X>0) + 4*(X>2) - (X>3) + epsilon
Y.2 = 6 + 0.4*X - 0.36*X^2 + 0.005*X^3 + epsilon 
Y.3 = 2.83 * sin(pi/2*X) + epsilon
Y.4 = 4 * sin(3*pi*X)*(X>0) + epsilon

#calcuate cv and mse for DGP1
cv.1 = rep(NaN, 11)
mse.1 = rep(NaN, 11)
k = 10
for(degree in 0:10)
{
  cv.1[degree+1] = poly.cv(X, Y.1, degree, k)
  mse.1[degree+1] = poly.mse(X, Y.1, degree)
}
plot(0:10, cv.1, pch = 1, ylim = c(0, 5), ylab="error", xlab="d")
points(0:10, mse.1, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  

cv.2 = rep(NaN, 11)
mse.2 = rep(NaN, 11)
for(degree in 0:10)
{
  cv.2[degree+1] = poly.cv(X, Y.2, degree, k)
  mse.2[degree+1] = poly.mse(X, Y.2, degree)
}
plot(0:10, cv.2, pch = 1, ylim = c(0, 6), ylab="error", xlab="d")
points(0:10, mse.2, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  

cv.3 = rep(NaN, 11)
mse.3 = rep(NaN, 11)
for(degree in 0:10)
{
  cv.3[degree+1] = poly.cv(X, Y.3, degree, k)
  mse.3[degree+1] = poly.mse(X, Y.3, degree)
}
plot(0:10, cv.3, pch = 1, ylim = c(0, 6), ylab="error", xlab="d")
points(0:10, mse.3, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  

cv.4 = rep(NaN, 11)
mse.4 = rep(NaN, 11)
for(degree in 0:10)
{
  cv.4[degree+1] = poly.cv(X, Y.4, degree, k)
  mse.4[degree+1] = poly.mse(X, Y.4, degree)
}
plot(0:10, cv.4, pch = 1, ylim = c(0, 6), ylab="error", xlab="d")
points(0:10, mse.4, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  
###################
poly.plot <- function(X, Y, degree)
{
  X.matrix = matrix(1, nrow = nrow(X), ncol = 1)
  if(degree > 0)
  {
    for(i in 1:degree)
      X.matrix = cbind(X.matrix, X^i)
  }
  beta = solve(crossprod(X.matrix), crossprod(X.matrix, Y))
  plot(X, Y)
  fit.Y = X.matrix%*%beta
  points(X[order(X),],fit.Y[order(X),], type ="l")
}
poly.plot(X, Y.1, 5)
poly.plot(X, Y.2, 2)
poly.plot(X, Y.3, 5)
poly.plot(X, Y.4, 0)

#####
X = as.matrix(runif(10000, min=-4, max=4))
epsilon = as.matrix(rnorm(10000))
Y.1 = -2*(X<(-3)) + 2.55*(X>(-2)) - 2*(X>0) + 4*(X>2) - (X>3) + epsilon
Y.2 = 6 + 0.4*X - 0.36*X^2 + 0.005*X^3 + epsilon 
Y.3 = 2.83 * sin(pi/2*X) + epsilon
Y.4 = 4 * sin(3*pi*X)*(X>0) + epsilon

cv.1 = rep(NaN, 11)
mse.1 = rep(NaN, 11)
k = 10
for(degree in 0:10)
{
  cv.1[degree+1] = poly.cv(X, Y.1, degree, k)
  mse.1[degree+1] = poly.mse(X, Y.1, degree)
}
plot(0:10, cv.1, pch = 1, ylim = c(0, 5), ylab="error", xlab="d")
points(0:10, mse.1, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  

cv.2 = rep(NaN, 11)
mse.2 = rep(NaN, 11)
for(degree in 0:10)
{
  cv.2[degree+1] = poly.cv(X, Y.2, degree, k)
  mse.2[degree+1] = poly.mse(X, Y.2, degree)
}
plot(0:10, cv.2, pch = 1, ylim = c(0, 6), ylab="error", xlab="d")
points(0:10, mse.2, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  

cv.3 = rep(NaN, 11)
mse.3 = rep(NaN, 11)
for(degree in 0:10)
{
  cv.3[degree+1] = poly.cv(X, Y.3, degree, k)
  mse.3[degree+1] = poly.mse(X, Y.3, degree)
}
plot(0:10, cv.3, pch = 1, ylim = c(0, 6), ylab="error", xlab="d")
points(0:10, mse.3, pch = 2)
legend(8.5, 5, pch=c(1,2), legend=c("cv", "mse"))  

cv.4 = rep(NaN, 11)
mse.4 = rep(NaN, 11)
for(degree in 0:10)
{
  cv.4[degree+1] = poly.cv(X, Y.4, degree, k)
  mse.4[degree+1] = poly.mse(X, Y.4, degree)
}
plot(0:10, cv.4, pch = 1, ylim = c(3, 6), ylab="error", xlab="d")
points(0:10, mse.4, pch = 2)
legend(8.5, 6, pch=c(1,2), legend=c("cv", "mse"))  

