# Since we use BFGS method in optim(), there may be errors due to gradient 
#calculation. In that case, reuse probit(), because I use randomized intial
#paramenters in optim() to guarantee the program runs well for most initial
#values 
#
#y is a n-length vector, X is a n*p matrix
probit <-function(y, X) 
{
  #switch y to a n*1 matrix 
  Y = matrix(y, ncol=1)  
  #check if the dimensions of y and X are compatible
  if(nrow(Y) != nrow(X))   
    return("the number of responds and the number of 
           observations are not consistent")
  #set the initial coefficients
  initial_beta = matrix(rnorm(ncol(X), 0, 0.1),ncol=1) 
  #check if X contains intercept.  If not, add it to X and modify the dimension
  # of initial parameter. 
  if(prod(apply(X, 2, var))!= 0)
  {
    X = cbind(1, X)
    initial_beta = rbind(rnorm(1, 0, 0.1), initial_beta)
  }
  
  #establish the targetd log-likelyhood function 
  loglik <- function(beta) 
  {
    psi_value = pnorm(X %*% beta)  #use matrix multiplication for conveniece 
    return( t(Y)%*%log(psi_value) + t(1-Y) %*% log(1-psi_value)) 
  }
  #optimization
  probit_optim = optim(par = initial_beta, fn=loglik, 
                       method = "BFGS", control = list(fnscale = -1), 
                       hessian = TRUE)  
  #result
  result = list("coefficients"=probit_optim$par, "log-likelihood value"=
                  probit_optim$value, "variance"=
                  solve(-probit_optim$hessian)) 
  return(result)
}


#generate matrix X
X = matrix(NaN,nrow=200, ncol=3)
for(i in 1:3)
  X[,i] = rnorm(nrow(X))
#choose beta as (1, 2, 1)
beta = matrix(c(1, 2, 1), ncol=1)
#generate y
prob = pnorm(X%*%beta)
y = matrix(NaN, nrow=nrow(X), ncol=1)
for(i in (1:nrow(X)))
  y[i,] = rbinom(1, 1, prob[i])
#calculate the result
probit(y, X)  

result = probit(y, X)



#generate beta from the asymptotic distribution
#
library(MASS)
#generate beta from the asymptotic distribution
beta_a = mvrnorm(n = 8000, mu = result$coefficients[-1], 
        Sigma = result$variance[-1, -1])


#generate beta gotten from bootstrapping
#
Itera = 8000
#the set of all the estimators
beta_b = matrix(NaN, nrow = Itera, ncol = 3)
i=1
N = nrow(X)
while(i <= Itera)
{
  #There may be errors due to method BFGS, so use tryCatch
  tryCatch(
    {
      #bootstrap sample from X, y
      index = sample(1:N, N, replace = TRUE)
      X_new = X[index,]
      y_new = y[index, , drop=FALSE]
      #calculate the new beta
      beta_b[i, ] = probit(y_new, X_new)$coefficients[-1]
      i = i + 1
    },
    error = function(e){}
  )
}



#Distribution gotten from drawing entirely new samples
#
Itera2 = 8000
beta_c = matrix(NaN, nrow = Itera2, ncol = 3)
i = 1

while(i <= Itera2)
{
  tryCatch(
    {
      #generate X, Y
      new_X = matrix(NaN,nrow=200, ncol=3)
      for(j in 1:3)
        new_X[,j] = rnorm(nrow(new_X))
      #choose beta as (1, 2, 1)
      new_beta = matrix(c(1, 2, 1), ncol=1)
      #generate y
      prob = pnorm(new_X%*%beta)
      new_y = matrix(NaN, nrow=nrow(new_X), ncol=1)
      for(j in (1:nrow(new_X)))
        new_y[j,] = rbinom(1, 1, prob[j])
      beta_c[i, ] = probit(new_y, new_X)$coefficients[-1]
      i = i + 1
    },
    error = function(e){}
  )
}
#plot beta_1
plot(density(beta_a[,1], , bw = 0.05), col = 1, lty = 1, 
     xlim=c(0.25, 1.75), ylim=c(0, 2.5), main = "beta_1")
lines(density(beta_b[,1], , bw = 0.05), col = 2, lty = 2)
lines(density(beta_c[,1], , bw = 0.05), col = 4, lty = 4)
legend(x = 1.3, y = 2.5, legend = c("method a", "method b", "method c"),
       col = c(1, 2, 4), lty = c(1, 2, 4))

#plot beta_2
plot(density(beta_a[,2], , bw = 0.05), col = 1, lty = 1, 
     xlim=c(1.25, 3.5), ylim=c(0, 1.5), main = "beta_2")
lines(density(beta_b[,2], , bw = 0.05), col = 2, lty = 2)
lines(density(beta_c[,2], , bw = 0.05), col = 4, lty = 4)
legend(x = 2.8, y = 1.5, legend = c("method a", "method b", "method c"),
       col = c(1, 2, 4), lty = c(1, 2, 4))

#plot beta_3
plot(density(beta_a[,3], , bw = 0.05), col = 1, lty = 1, 
     xlim=c(0.4, 2.25), ylim=c(0, 2.2), main = "beta_3")
lines(density(beta_b[,3], , bw = 0.05), col = 2, lty = 2)
lines(density(beta_c[,3], , bw = 0.05), col = 4, lty = 4)
legend(x = 1.7, y = 2.0, legend = c("method a", "method b", "method c"),
       col = c(1, 2, 4), lty = c(1, 2, 4))



#read_data
PS2_data = read.csv(file = "E:\\PS2_data.csv", header = TRUE)

cook_inc_risk = PS2_data$cook_inc_risk
diff = PS2_data$diff
same = PS2_data$same
scaled_comp_challenger = c(scale(PS2_data$comp_challenger)) #standardize the data
scaled_comp_incumbent = c(scale(PS2_data$comp_incumbent)) #standardize the data
inc_tenure = PS2_data$inc_tenure
tenure_squared = inc_tenure^2
inc_age = PS2_data$inc_age
age_squared = inc_age^2
vote = PS2_data$vote_incumbent
#construct data frame
Senate_1 = data.frame(cook_inc_risk, diff, same, scaled_comp_challenger,
                    scaled_comp_incumbent, inc_tenure,  
                    inc_age, vote)
glm_Senate_1 = glm(vote~., data = Senate_1, 
                   family = binomial(link = "probit"))
summary(glm_Senate_1)



library(sandwich)
library(zoo)
library(lmtest)
#regenerate a data frame about all the factors.

state = PS2_data$state
Senate_2 = data.frame(cook_inc_risk, diff, same, scaled_comp_challenger,
                            scaled_comp_incumbent, inc_tenure,  
                            inc_age, vote, state)
#perform glm with binomial and link function, probit 
glm_Senate_2 = glm(vote~.-state, data = Senate_2, 
                 family = binomial(link = "probit"))
unique_state = unique(state)
m = length(unique_state)
p = ncol(model.matrix(glm_Senate_2))
#calculate score
scores=estfun(glm_Senate_2)
#clustered score
Senate_clust=matrix(NA,nrow=m,ncol=p)
for(j in 1:p)
  Senate_clust[,j]= tapply(scores[,j], state, sum)
meat.cl = (m/(m-1)) * t(Senate_clust) %*% Senate_clust
vcov.cl = vcov(glm_Senate_2)
vcov.cl = vcov.cl %*% meat.cl %*% vcov.cl
#retest.
round(coeftest(glm_Senate_2, vcov=vcov.cl),4)

###############
#simulate the facial competency, which is 3 standard deviations below the mean
comp_3_sigma = pnorm(-3)
r_comp_challenger =  qnorm(runif(length(scaled_comp_challenger)) * comp_3_sigma)
Senate_3 = cbind(1, cook_inc_risk, diff, same, 
                      scaled_comp_challenger = r_comp_challenger,
                      scaled_comp_incumbent, inc_tenure,  
                      inc_age)
#simulate the vote
prob = pnorm(Senate_3 %*% matrix(glm_Senate_2$coefficients, ncol =1))
y = rbinom(length(prob), 1, prob)
mean(y)
############

itera3 = 5000
expect_share = rep(NaN, itera3)
mu_3 = matrix(glm_Senate_1$coefficients, ncol = 1)
sigma_3 = vcov.cl
for(i in 1:itera3)
{
  #generate 
  coefficients_3 = mvrnorm(1, mu = mu_3, Sigma = sigma_3)
  #the number of simulation for the politician
  N = 500 
  n = nrow(Senate_2)
  #gengerate joint distribution for same and diff
  p_s0_d0 = sum(same == 0 & diff == 0) / n
  p_s1_d0 = sum(same == 1 & diff == 0) / n
  p_s0_d1 = sum(same == 0 & diff == 1) / n
  sample_prob = c(p_s0_d0, p_s1_d0, p_s0_d1)
  sample_value = matrix(c(0, 0, 1, 0, 0, 1), nrow = 3,byrow = TRUE)
  #generate random sample
  s_d_sample = sample(1:3, size = N,
                      prob = sample_prob, replace = TRUE)
  r_same = sample_value[s_d_sample, 1]
  r_diff = sample_value[s_d_sample, 2]
  r_cook_inc_risk = sample(cook_inc_risk, N, replace = TRUE)
  comp_3_sigma = pnorm(-3)
  r_comp_challenger =  qnorm(runif(N) * comp_3_sigma)
  r_comp_incumbent = rnorm(N)
  r_inc_tenure = sample(inc_tenure, N, replace = TRUE)
  r_inc_age = rnorm(N, mean = mean(inc_age), sd = sd(inc_age))
  #calculate the expectation.
  Senate_3 = cbind(1, r_cook_inc_risk, r_diff, r_same, r_comp_challenger,
                   r_comp_incumbent, r_inc_tenure,  
                   r_inc_age)
  expect_share[i] = mean(pnorm(Senate_3 %*% coefficients_3))
}
mean(expect_share)


###################

itera4 = 8000
Senate_4 = data.frame(cook_inc_risk, diff, same, 
                 scaled_comp_challenger,
                 scaled_comp_incumbent, inc_tenure,  
                 inc_age)
Senate_4_25 = cbind(1, cook_inc_risk, diff, same, qnorm(0.25),
                 scaled_comp_incumbent, inc_tenure,  
                 inc_age)
Senate_4_75 = cbind(1, cook_inc_risk, diff, same, qnorm(0.75),
                    scaled_comp_incumbent, inc_tenure,  
                    inc_age)
N = nrow(Senate_4)
marginal_effect = matrix(rep(NaN, itera4), ncol =1)
for(i in 1:itera4)
{
  #generate 
  index = sample(N, N, replace = TRUE);
  glm_Senate_4 = glm(vote[index]~., data = Senate_4[index,], 
                     family = binomial(link = "probit"))
  coefficients = matrix(glm_Senate_4$coefficients, ncol =1)
  prob_25 =  pnorm(Senate_4_25%*%coefficients)
  y_25 = rbinom(length(prob_25), 1, prob_25)
  prob_75 =  pnorm(Senate_4_75%*%coefficients)
  y_75 = rbinom(length(prob_75), 1, prob_75)
  marginal_effect[i,] = mean(prob_75) - mean(prob_25)
}

dens = density(marginal_effect)
plot(dens, main = "bootstrap")
q_25 = quantile(marginal_effect, 0.025)
q_75 = quantile(marginal_effect, 0.975)
x1 = min(which(dens$x >= q_25))
x2 = max(which(dens$x < q_75))
with(dens, polygon(x=c(x[c(x1, x1:x2, x2)]), y=c(0, y[x1:x2], 0), col = "gray"))


marginal_effect2 = matrix(rep(NaN, itera4), ncol =1)
for (i in 1:itera4){
  beta_new = mvrnorm(1,glm_Senate_2$coefficients, vcov.cl)
  prob_25=pnorm(Senate_4_25%*%beta_new)
  prob_75=pnorm(Senate_4_75%*%beta_new)
  Y_25=rbinom(length(prob_25),1,prob_25)
  Y_75=rbinom(length(prob_75),1,prob_75)
  marginal_effect2[i] = mean(Y_75-Y_25)
}

dens2 = density(marginal_effect2)
plot(dens2, main = "Simulation Approximation")
q_25 = quantile(marginal_effect2, 0.025)
q_75 = quantile(marginal_effect2, 0.975)
x1 = min(which(dens$x >= q_25))
x2 = max(which(dens$x < q_75))
with(dens, polygon(x=c(x[c(x1, x1:x2, x2)]), y=c(0, y[x1:x2], 0), col = "gray"))


