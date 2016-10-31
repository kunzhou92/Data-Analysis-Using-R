RRKR2 <- function(par, lambda, K, Y, iters, epsilon, m, step)
{
  for(i in 1:iters)
  {
    n = nrow(Y)
    u = Y - K %*% par
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
    
    sol.mat = lambda * I + I0 %*% K
    q = -I0 %*% Y + m * e
    par = (1-step)*par - solve(sol.mat, q)*step
  }
  return(par)
}


RRKR2.cv <- function(par, K, Y, 
                    iters, epsilon, m, step, lam.beg, lam.end, lam.step, k)
{
  n = nrow(Y)
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
      par_j = RRKR2(par.train, lambda[i], K.train, Y.train, iters, 
                    epsilon, m, step)
      y.fitted =  K[fold[[j]], -fold[[j]], drop=F] %*% par_j
      err.sq = sum((Y[fold[[j]],1] - y.fitted)^2)
      cv[i] = cv[i] + err.sq
    }
  }
  cv = cv / n
  return(list(lam = lambda, c = cv))
}
