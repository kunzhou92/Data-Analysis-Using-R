library(Rcpp)
sourceCpp("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\project\\kernel.cpp")
KLR <- function(par, lambda, K, Y, iters){
  for(i in 1:iters){
    p <- 1/(1 + exp(-K %*% par))
    w <- p * (1-p)  
    W <- diag(as.vector(w),nrow(K),nrow(K))
    W_s <- solve(W)
    S <- K %*% par + W_s %*% (Y-p)
    par <- solve(K + lambda*W_s, S) 
  }
  return (par)
}

KLR.cv <-function(par, K, Y, iters, lam.beg, lam.end, lam.step, k)
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
      par.train = par[-fold[[j]], , drop = F]
      K.train = K[-fold[[j]], -fold[[j]], drop = F]
      Y.train = Y[-fold[[j]], ,drop = F]
      alpha_j = KLR(par.train, lambda[i], K.train, Y.train, iters)
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
  y[S2.index] = (x[S2.index] - epsilon)^2
  y[S4.index] = (x[S4.index] + epsilon)^2
  y[S5.index] = -m*(2*(x[S5.index]+epsilon)+m)
  return(y)
}

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
      par.train = par[-fold[[j]], , drop = F]
      K.train = K[-fold[[j]], -fold[[j]], drop = F]
      Y.train = Y[-fold[[j]], ,drop = F]
      par_j = RRKR(par.train, lambda[i], K.train, 
                   Y.train, iters, epsilon, m, step)
      alpha_j = par_j[1:nrow(K.train),,drop=F]
      b_j = par_j[nrow(K.train)+1,]
      y.fitted =  b_j + K[fold[[j]], -fold[[j]], drop=F] %*% alpha_j
      err.sq = sum((Y[fold[[j]],1] - y.fitted)^2)
      cv[i] = cv[i] + err.sq
    }
  }
  cv = cv / n
  return(list(lam = lambda, c = cv))
}


KRR<-function(K, Y, lambda) {
  return(solve(K + lambda * diag(nrow(Y)), Y))
}

KRR.cv <- function(K, Y, lam.beg, lam.end, lam.step, k)
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

