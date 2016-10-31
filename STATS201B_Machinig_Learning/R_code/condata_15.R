condata_15 <- read.table("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\project\\german.data.txt")
colnames(condata_15) <- c("V1","V2","V3","V4","Y.ca","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21")
condata_15 <- model.matrix(~V1 + V2 + V3 + V6 + V7 + V8 + V11 + V12 + V13 + V14 + V16 + V18 + V19 + V20 + V21 + Y.ca, condata_15)
condata_15 <- condata_15[,-1]
condata_15 <- as.data.frame(condata_15)
condata_15$V21[condata_15$V21==2] <- 0 

set.seed(31)
indexes = sample(1:nrow(condata_15), size=0.7*nrow(condata_15))
condata_15.train = condata_15[indexes,]
condata_15.test = condata_15[-indexes,]
for (j in c(4,17,18,22,25,26,30))
  condata_15.test[,j]=(condata_15.test[,j]-mean(condata_15.train[,j]))/sd(condata_15.train[,j])
for (j in c(4,17,18,22,25,26,30))
  condata_15.train[,j] = (condata_15.train[,j]-mean(condata_15.train[,j]))/sd(condata_15.train[,j])


Xc_15_train <- condata_15.train[,-30]
Yc_15_train <- condata_15.train[,30]
Xc_15_test <- condata_15.test[,-30]
Yc_15_test <- condata_15.test[,30]

Xc_15_train_mat <- as.matrix(Xc_15_train)
Yc_15_train_mat <- as.matrix(Yc_15_train)
Xc_15_test_mat <- as.matrix(Xc_15_test)
Yc_15_test_mat <- as.matrix(Yc_15_test)

