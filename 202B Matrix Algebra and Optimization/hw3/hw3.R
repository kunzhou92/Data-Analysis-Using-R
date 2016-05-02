library(MMST)
data(pet)
help(pet)
head(pet)
data(tobacco)
help(tobacco)
head(tobacco)


#1)
tls<-function(x,y,kappa=1,lambda=1) {
  m<-ncol(x)
  u<-cbind(kappa*x,lambda*y)
  s<-svd(u,nu=m,nv=m)
  z<-(s$u)%*%diag((s$d)[1:m]); c<-s$v
  c1<-t(s$v)[,1:m]; c2<-t(s$v)[,-(1:m)]
  m<-kappa*solve(c1); b<-m%*%c2/lambda
  z<-z%*%c1/kappa
  ssq1<-sum((x-z)^2); ssq2<-sum((y-z%*%b)^2)
  ssqt<-kappa*ssq1+lambda*ssq2
  return(list(b=b,z=z,ssq=c(ssq1,ssq2,ssqt)))
}

X.toba = as.matrix(tobacco[,4:9])
Y.toba = as.matrix(tobacco[,1:3])
beta_tls = tls(X.toba, Y.toba)$b
fitted_Y_tls = tls(X.toba, Y.toba)$z
beta_ols = solve(crossprod(X.toba),crossprod(X.toba,Y.toba))
fitted_Y_ols = X.toba %*% beta_ols
plot(Y.toba[,1], fitted_Y_tls[,1], pch=0, ylim = c(1, 3), main = "burn rate")
points(Y.toba[,1], fitted_Y_ols[,1], pch=1)
legend(x = 2.0, y = 3.0, legend = c("tls", "ols"), pch= c(0, 1))
plot(Y.toba[,2], fitted_Y_tls[,2], pch=0,  ylim = c(1, 20),main = "sugar")
points(Y.toba[,2], fitted_Y_ols[,2], pch=1)
legend(x = 19, y = 15, legend = c("tls", "ols"), pch= c(0, 1))
plot(Y.toba[,3], fitted_Y_tls[,3], pch=0,  ylim = c(1, 4),main = "nicotine")
points(Y.toba[,3], fitted_Y_ols[,3], pch=1)
legend(x = 3, y = 3, legend = c("tls", "ols"), pch= c(0, 1))

#3)
ridge<-function(x,y,kappa=0) {
  b<-solve(crossprod(x)
           +kappa*diag(ncol(x)),crossprod(x,y))
  ssq<-sum((y-x%*%b)^2)
  return(list(b=b,ssq=ssq))
}
ridgeTrace<-function(x,y,kmax,nmax) {
  kk<-seq(0,kmax,length=nmax+1)[-1]
  bb<-array(0,c(ncol(x),ncol(y),nmax))
  for (i in 1:nmax)
    bb[,,i]<-ridge(x,y,kappa=kk[i])$b
  plot(kk,seq(min(bb),max(bb),length=nmax),
       type="n",xlab="kappa",ylab="B")
  for (i in 1:60) for (j in 1:ncol(y))
    lines(kk,bb[i,j,], col = i)
}

X = as.matrix(pet[,1:268])
Y = as.matrix(pet[,269, drop = FALSE])
ridgeTrace(X, Y, 1, 60)

#4)
ridgeParameter<-function(x, y, kmax, nmax) #x, y are matrices
{
  kk<-seq(0,kmax,length=nmax+1)[-1]
  crossv = rep(NaN, nmax)
  for(i in 1:nmax)
  {
    sse = matrix(NaN, nrow = nrow(y), ncol = ncol(y))
    for(j in 1:nrow(x))
    {
      xTemp = x[-j, ,drop = FALSE]
      yTemp = y[-j, ,drop = FALSE]
      sse[j,] = y[j,,drop = FALSE] - cbind(1,x[j,,drop = FALSE]) %*% 
        ridge(xTemp, yTemp, kk[i])$b
    }
    crossv[i] = mean(sse^2)
  }
  return(crossv)
}
X.scale = scale(X)
Y.scale = scale(Y)
nmax = 100
kmax = 0.1
kk = seq(0,kmax,length=nmax+1)[-1]
crossv = ridgeParameter(X.scale, Y.scale, kmax, nmax)
minIndex = match(min(crossv),crossv)[1]
kk[minIndex]
plot(kk, crossv)



#5)
data_bank = read.table("http://www.stat.ucla.edu/~handcock/
                       202B/datasets/SwissBankNotes.txt",
           skip=22,header=TRUE)
data_bank = as.matrix(data_bank)
genuine = princomp(data_bank[1:100,])
counterfeit = princomp(data_bank[101:200,])
together = princomp(data_bank)
plot(genuine$scores[,1], genuine$scores[,2], pch = 0, 
     ylim = c(-2, 2), xlim = c(-3, 3.5))
points(counterfeit$scores[,1], counterfeit$scores[,2], pch = 1)
points(together$scores[,1], together$scores[,2], pch = 2)
legend(x = 2, y = 2, legend = c("genuine", "counterfeit", "together"),
       pch= c(0, 1, 2))


z = data_bank 
m<-colSums(z)/ nrow(z)
z<-z-outer(rep(0,nrow(z)),m,"+")
s<-svd(z,nu=2,nv=2)
u<-s$u; v<-s$v; d<-s$d; w<-v*outer(rep(1,ncol(z)),d[1:2])
f<-u%*%t(w)

plot(as.matrix(z),f,xlab="data",ylab="fit")
abline(0,1)


plot(u,type="n",xlab="dim1",ylab="dim2")
text(u,as.character(1:20))
for (i in 1:ncol(z))
  abline(0,w[i,2]/w[i,1])

