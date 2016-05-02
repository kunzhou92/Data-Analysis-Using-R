library(randomForest)
data = read.csv(file="E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\pset5\\Heart.csv", header = T)
#Divide dataset into training set and testing set
Index = sample(nrow(data), 150)
TrainSet = data[Index,]
TestSet = data[-Index,]
Training = randomForest(AHD~., data=TrainSet, mtry=3,         #Training
                        ntree=1000, na.action = na.omit)
Testing = predict(Training, newdata=TestSet, type="response", 
                               na.action = na.omit)           #Testing


logic = Testing==TestSet[,"AHD"]
True = length(which(logic == T))
False = length(which(logic == F))
accuracy = True / (True+False)
Error = 1 - accuracy

ActualValue = data[as.numeric(names(Training$predicted)),"AHD"]
oobError = sum(ActualValue!=Training$predicted) / length(Training$predicted)


RFplot <- function(data, ntry)
{
  Index = sample(nrow(data), 150)
  TrainSet = data[Index,]
  TestSet = data[-Index,]
  Test.error = c()
  oob.error = c()
  for(m in 1:ntry)
  {
    Training = randomForest(AHD~., data=TrainSet, mtry=3,     
                            ntree=500, na.action = na.omit)
    Testing = predict(Training, newdata=TestSet, type="response", 
                      na.action = na.omit)       
    
    logic = Testing==TestSet[,"AHD"]
    True = length(which(logic == T))
    False = length(which(logic == F))
    accuracy = True / (True+False)
    Error = 1 - accuracy
    Test.error = c(Test.error, Error)
    ActualValue = data[as.numeric(names(Training$predicted)),"AHD"]
    oobError = sum(ActualValue!=Training$predicted) / length(Training$predicted)
    oob.error = c(oob.error, oobError)
  }
  plot(1:ntry, Test.error, type="l", col=1, ylim=c(0,0.4), ylab="errors",
       xlab="mtry")
  points(1:ntry, oob.error, type="l", col=2)
  legend(x= 10, y=0.4, lty=c(1,1), col=c(1,2), 
         legend=c("Testing error", "oob error"))
}
RFplot(data, ncol(data))


Rdata = data
for(i in 1:14)
  Rdata = Rdata[which(!is.na(Rdata[,i])),]
Index = sample(nrow(Rdata), 150)
TrainSet = Rdata[Index,]
TestSet = Rdata[-Index,]
library(glmnet)
dataMat = model.matrix(AHD~.-1, data = TrainSet)
testMat = model.matrix(AHD~.-1, data = TestSet)
Respond = rep(1,nrow(TrainSet))
Respond[which(as.matrix(TrainSet[,"AHD"]) == "No")] = 0
tRespond = rep(1,nrow(TestSet))
tRespond[which(as.matrix(TestSet[,"AHD"]) == "No")] = 0
cv = cv.glmnet(x=dataMat, y=Respond,family="binomial")
lambda = cv$lambda.min
Training = glmnet(x=dataMat, y=Respond, family = "binomial", lambda = lambda)
predict = as.numeric(predict(Training, newx = testMat, type="class"))
error = 1 - mean(predict == tRespond)



load("E:\\Papers\\Statistics\\201B Statistical Modeling and Learning\\pset5\\pset5words.RData")
wordDat = as.data.frame(wordMatrix)
wordDat$vote = as.factor(wordDat$vote)
index = sample(nrow(wordDat), floor(nrow(wordDat))*0.7)
train = wordDat[index,]
test = wordDat[-index,]
trainMat = wordMatrix[index,]
testMat = wordMatrix[-index,]
library(e1071)
#SVM
#linear knernel
tune.out.linear=tune(svm, vote~., data=train ,kernel="linear")
best.svm.linear=tune.out.linear$best.model
class.test.linear=predict(best.svm.linear, test)
table(class.test.linear, test$vote)

#Gaussian kernel
tune.out.radial=tune(svm, vote~., data=train,kernel="radial")
best.svm.radial=tune.out.radial$best.model
class.test.radial=predict(best.svm.radial, test)
table(class.test.radial, test$vote)

#polynomial kernel
tune.out.polynomial=tune(svm, vote~., data=train,kernel="polynomial")
best.svm.polynomial=tune.out.polynomial$best.model
class.test.polynomial=predict(best.svm.polynomial, test)
table(class.test.polynomial, test$vote)

#LASSo logit
library(glmnet)
vote.index = match("vote", colnames(trainMat))
cv = cv.glmnet(x=trainMat[,-vote.index], y=trainMat[,vote.index, drop=F],
               family="binomial")
lambda = cv$lambda.min
Training = glmnet(x=trainMat[,-vote.index], y=trainMat[,vote.index, drop=F],
                  family = "binomial", lambda = lambda)
predict = as.numeric(predict(Training, newx = testMat[,-vote.index],
                             type="class"))
table(predict, test$vote)


#naive Bayes
library(e1071)
tune.NB = tune(naiveBayes, vote~., data=train)
best.NB=tune.NB$best.model
class.test.NB =predict(best.NB, test)
table(class.test.NB, test$vote)


#random forests
Training = randomForest(vote~., data=train[,-(1:4000)], mtry=2500, ntree=1000)
Testing = predict(Training, newdata=test, type="response")   
table(Testing, test$vote)

#ensemble
ensemble = (as.numeric(class.test.linear) + as.numeric(class.test.radial) +
              (predict+1) + as.numeric(Testing)) / 4  
ensemble[which(ensemble<=1.4)] = 0
ensemble[which(ensemble>=1.6)] = 1
ensemble[which(ensemble==1.5)] = rbinom(1,1,0.5)
ensemble = as.factor(ensemble)
table(ensemble, test$vote)
