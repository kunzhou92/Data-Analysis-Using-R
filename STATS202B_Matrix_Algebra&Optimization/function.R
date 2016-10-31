tls<-function(x,y,kappa=1,lambda=1) {
  # X is an n x m matrix with m <= n
  # We wish to find the argmon over z and b of
  # kappa*|x-z|^2 + lambda*|y-z %*% b|^2
  # This can be restated as on page 20 of 202Blecture3.pdf
  #	
  n<-nrow(x)
  m<-ncol(x)
  u<-cbind(kappa*x,lambda*y)
  # Compute the SVD of X (or at least the first m vectors and singular values)
  s<-svd(u,nu=m,nv=m)
  # Form z recognizing the row and column restrictions.
  z<-(s$u)%*%diag(s$d[1:m]) #; c<-s$v
  c1<-t(s$v)[,1:m]; c2<-t(s$v)[,-(1:m)]
  m<-kappa*solve(c1)
  b<-m%*%c2/lambda
  z<-z%*%c1/kappa
  # Compute the sums-of-squares
  ssq1<-sum((x-z)^2); ssq2<-sum((y-z%*%b)^2)
  ssqt<-kappa*ssq1+lambda*ssq2
  ssq <- c(ssq1,ssq2,ssqt)
  names(ssq) <- c("X","Y","Total")
  return(list(b=b,z=z,ssq=ssq))
}