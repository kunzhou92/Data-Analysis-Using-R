#plot Huber loss function
x1 = seq(-8, 8, 0.1)
y = robust.loss(x1, 2, 0.3)
plot(x1, y, type="l", main="Huber Loss Function", xlim=c(-9, 9), ylim = c(0, 70),
     xlab = "u", ylab="V(u)")
x2 = seq(-8, 8, 0.1)
points(x2, x2^2, type="l", col="red")
legend(0, 70, col = c("black", "red"), lty=c(1,1), 
       legend=c("Huber loss", "Squared loss"))

plot(density(Yc_15_test_mat), xlab="scaled credit amount", ylab="density",
     main="density function")

K.all = GaussianKernel(rbind(Xc_15_train_mat, Xc_15_test_mat), sqrt(29))
#Kernel Ridge Regression
cv = KRR.cv(K.all[1:700, 1:700], Yc_15_train_mat, 0.1, 2, 0.025, 5)
plot(cv$lam, cv$c, xlab="lambda", ylab="cv")
lambda= 0.5
alpha = KRR(K.all[1:700, 1:700], Yc_15_train_mat, lambda)
y.krr =  K.all[701:1000, 1:700, drop=F] %*% alpha
plot(Yc_15_test_mat, y.krr, xlab="Y", ylab="fitted Y", 
     main=expression(paste("Kernel Ridge Regression ",lambda,"=0.15")), ylim = c(-1.5, 2.8))
t.test(Yc_15_test_mat-y.krr,mu=0)



#Robust Regularized Kernel Regression
para = matrix(c(rep(0.005, nrow(Yc_15_train_mat)), -0.01))
m = 3.5
epsilon = 0.05

step = 0.5
iters = 25
cv = RRKR.cv(para, K.all[1:700, 1:700], Yc_15_train_mat,
             iters, epsilon, m, step, 1, 2, 0.1, 3)
plot(cv$lam,cv$c, xlab="lambda", ylab="cv")
lambda = 1.2
t = RRKR(para, lambda, K.all[1:700, 1:700], Yc_15_train_mat, iters, epsilon, m, step)
y.rrkr = K.all[701:1000, 1:700, drop=F] %*% t[1:(nrow(t)-1),,drop = F] +t[nrow(t),]
plot(Yc_15_test_mat, y.rrkr, xlab="Y", ylab="fitted Y", 
     main=expression(paste("Robust Kernel Regression m=3.5, ", epsilon,"=0.05, ",lambda,"=1.2")), ylim = c(-1.5, 2.8))
t.test(Yc_15_test_mat - y.rrkr,mu=0)


index = which(Yc_15_test_mat<Inf)
library(MASS)
t.krr = wilcox.test(Yc_15_test_mat[index], y.krr[index], paired=TRUE) 
t.rrkr = wilcox.test(Yc_15_test_mat[index], y.rrkr[index], paired=TRUE) 

##################################3

load("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\project\\kangXin\\condata_10.RData")

K.all = GaussianKernel(rbind(Xc_10_train_mat, Xc_10_test_mat), sqrt(29))
#Kernel Ridge Regression
cv = KRR.cv(K.all[1:700, 1:700], Yc_10_train_mat, 0.1, 2, 0.025, 5)
plot(cv$lam, cv$c, xlab="lambda", ylab="cv")
lambda= 0.5
alpha = KRR(K.all[1:700, 1:700], Yc_10_train_mat, lambda)
y.krr =  K.all[701:1000, 1:700, drop=F] %*% alpha
plot(Yc_10_test_mat, y.krr, xlab="Y", ylab="fitted Y", 
     main=expression(paste("Kernel Ridge Regression ",lambda,"=0.15")), ylim = c(-1.5, 2.8))


#Robust Regularized Kernel Regression
para = matrix(c(rep(0.005, nrow(Yc_10_train_mat)), -0.01))
m = 3.5
epsilon = 0.05
step = 0.5
iters = 25
cv = RRKR.cv(para, K.all[1:700, 1:700], Yc_10_train_mat,
             iters, epsilon, m, step, 1, 10, 1, 3)
plot(cv$lam,cv$c, xlab="lambda", ylab="cv")
lambda = 1.2
t = RRKR(para, lambda, K.all[1:700, 1:700], Yc_10_train_mat, iters, epsilon, m, step)
y.rrkr = K.all[701:1000, 1:700, drop=F] %*% t[1:(nrow(t)-1),,drop = F] +t[nrow(t),]
plot(Yc_10_test_mat, y.rrkr, xlab="Y", ylab="fitted Y", 
     main=expression(paste("Robust Kernel Regression m=3.5, ", epsilon,"=0.05, ",lambda,"=1.2")), ylim = c(-1.5, 2.8))



index = which(Yc_10_test_mat<Inf)
library(MASS)
t.krr = wilcox.test(Yc_10_test_mat[index], y.krr[index], paired=TRUE) 
t.rrkr = wilcox.test(Yc_10_test_mat[index], y.rrkr[index], paired=TRUE) 


