EM <- function(n, alpha, mu1, mu2, iter)
{
  i = 1:length(n) - 1
  alpha.m = alpha
  mu1.m = mu1
  mu2.m = mu2
  z.m = alpha.m * exp(-mu1.m) * mu1.m^i/(alpha.m * exp(-mu1.m) * mu1.m^i + 
                                          (1-alpha.m) * exp(-mu2.m) * mu2.m^i)
  for(j in 1:iter)
  {
    alpha.m = sum(z.m*n) / sum(n)
    mu1.m = sum(n*i*z.m) / sum(n*z.m)
    mu2.m = sum(n*i*(1-z.m)) / sum(n*(1-z.m))
    z.m = alpha.m * exp(-mu1.m) * mu1.m^i/
      (alpha.m * exp(-mu1.m) * mu1.m^i + (1-alpha.m) * exp(-mu2.m) * mu2.m^i)
  }
  return(list(alpha=alpha.m, mu1=mu1.m, mu2=mu2.m))
}
n = c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
a = seq(from = 0, to = 8, by = 0.2)
iters = floor(exp(a))
errors = c()
for(i in iters)
{
  result = EM(n, 0.3, 1, 2.5, i) 
  err = (result[[1]] - 0.3599)^2 + (result[[2]] - 1.2561)^2 
  + (result[[3]] - 2.6634)^2
  errors = c(errors, err)
}
plot(a, errors, xlab = "log iters", ylab = "squared error")

EM.acc <- function(n, alpha, mu1, mu2, lambda)
{
  i = 1:length(n) - 1
  alpha.m = alpha
  mu1.m = mu1
  mu2.m = mu2
  z.m = alpha.m * exp(-mu1.m) * mu1.m^i/(alpha.m * exp(-mu1.m) * mu1.m^i + 
                                           (1-alpha.m) * exp(-mu2.m) * mu2.m^i)
  step = 0
  while(err > 1e-6)
  {
    alpha.n = sum(z.m*n) / sum(n)
    mu1.n = sum(n*i*z.m) / sum(n*z.m)
    mu2.n = sum(n*i*(1-z.m)) / sum(n*(1-z.m))
    alpha.m = (alpha.n - alpha.m) * lambda + alpha.m
    mu1.m = (mu1.n - mu1.m) * lambda + mu1.m
    mu2.m = (mu2.n - mu2.m) * lambda + mu2.m
    z.m = alpha.m * exp(-mu1.m) * mu1.m^i/
      (alpha.m * exp(-mu1.m) * mu1.m^i + (1-alpha.m) * exp(-mu2.m) * mu2.m^i)
    err = (alpha.m - 0.3599)^2 + (mu1.m - 1.2561)^2 
    + (mu2.m - 2.6634)^2
    step = step + 1
  }
  return(list(alpha=alpha.m, mu1=mu1.m, mu2=mu2.m, step = step))
}

lambda = seq(from = 1.5, to = 1, by = -0.05)
step = c()
for(i in lambda)
  step = c(step,  EM.acc(n, 0.3, 1, 2.5, i)[[4]])
plot(lambda, step)
