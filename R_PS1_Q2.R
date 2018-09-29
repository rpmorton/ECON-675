#Created by RM on 2018.09.14 for ECON 675 PS 1, Q 2
#Subpart 4

#load package
library(MASS)
library(plyr)
library(data.table)
library(matlib)
library(xtable)

########
######Q2, Part 4: Implementation of OLS
#######
rm(list=ls())

######Generate Data
nobs <- 1000
x_1 <- runif(nobs,0,100)
x_2 <- runif(nobs,min=30,max=50)
intercept <- rep.int(1,nobs)
x_1_2 <- c(intercept,x_1,x_2)
X <- matrix(x_1_2, nrow = nobs, ncol = 3)

epsilon <- rnorm(nobs, mean=0, sd=10)

beta <- matrix( c(0,1,5), nrow = 3, ncol = 1)

y <- X %*% beta + epsilon

W <-  diag(x=1, nobs, nobs)


#####i: Symmetric Inverse

####i)a:Estimate Beta and Variance
forinv <- diag(x=1,ncol(X),ncol(X))
XWX_inv <- solve(t(X) %*% W  %*% X,forinv)
beta_hat <- XWX_inv %*% ( t(X)  %*% W  %*% y)

h_0_hat <- nrow(X)^(-1) * (t(X) %*% X)
epsilon_hat <- (y - X %*% beta_hat)
s2 <- (nrow(X)-nrow(beta_hat))^(-1)*(t(epsilon_hat) %*% epsilon_hat)

v_0_hat <- (as.integer(s2)* nrow(X)^(-1)) * (t(X) %*% X)
h_0_hat_inv <- solve(h_0_hat,diag(x=1,ncol(X)))
V_hat <- h_0_hat_inv %*% v_0_hat %*% h_0_hat_inv

####i)b: T-test by element of beta against H_0 of 0
t_numerator <- beta_hat - rep.int(0,nrow(beta_hat))
asyvar <- diag(V_hat)
t_denominator <- sqrt(asyvar * nrow(X)^(-1) )

t_stat <- t_numerator / t_denominator

####i)c: p-value
pval <- 2*pt(-abs(t_stat),df=nrow(X)-nrow(beta_hat))

####i)d: Confidence Interal
alpha <- .90

lb <- beta_hat + qt((1-alpha)/2,df)*t_denominator
ub <-  beta_hat + qt(alpha+(1-alpha)/2,df)*t_denominator

CI <- data.frame(lb,ub,alpha,beta_hat,beta)


#####ii: Choleskey Inverse

####ii)a:Estimate Beta and Variance with Cholesky
XWX_chol <- chol(t(X) %*% W  %*% X)
XWX_cholinv <- chol2inv(XWX_chol)
beta_hat_chol <- XWX_cholinv %*% ( t(X)  %*% W  %*% y)

h_0_hat <- nrow(X)^(-1) * (t(X) %*% X)
epsilon_hat_chol <- (y - X %*% beta_hat_chol)
s2_chol <- (nrow(X)-nrow(beta_hat_chol))^(-1)*(t(epsilon_hat_chol) %*% epsilon_hat_chol)

v_0_hat_chol <- (as.integer(s2_chol)* nrow(X)^(-1)) * (t(X) %*% X)
h_0_hat_chol <- chol(h_0_hat)
h_0_hat_chol_inv <- chol2inv(h_0_hat_chol)
V_hat_chol <- h_0_hat_chol_inv %*% v_0_hat_chol %*% h_0_hat_chol_inv

####ii)b: T-test by element of beta against H_0 of 0 with Cholesky
t_numerator_chol <- beta_hat_chol - rep.int(0,nrow(beta_hat))
asyvar_chol <- diag(V_hat_chol)
t_denominator_chol <- sqrt(asyvar_chol * nrow(X)^(-1) )

t_stat_chol <- t_numerator_chol / t_denominator_chol

####ii)c: p-value
pval_chol <- 2*pt(-abs(t_stat_chol),df=nrow(X)-nrow(beta_hat_chol))

####i)d: Confidence Interal
alpha <- .90

lb_chol <- beta_hat_chol + qt((1-alpha)/2,df)*t_denominator_chol
ub_chol <-  beta_hat_chol + qt(alpha+(1-alpha)/2,df)*t_denominator_chol

CI_chol <- data.frame(lb_chol,ub_chol,alpha,beta_hat_chol,beta)


###Up to at least 7 digits, the Symmetric and Cholesky Confidence Intervals
###(which are based on beta and variance estimators) are identical!

########
######Q2, Part 5: Lalonde Data
######

#####a: Use Estimator Coded Above

###Estimate Beta and Variance Matrix
rm(list=ls())

lalonde <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/LaLonde_1986.csv")
 
lalonde$educ2 <- lalonde$educ * lalonde$educ
lalonde$black_earn74 <- lalonde$black * lalonde$earn74
lalonde$intercept <- rep.int(1,nrow(lalonde))

x_vars <- c("intercept","treat","black","age","educ","educ2","earn74","black_earn74","u74","u75")
X_lalonde <- as.matrix(subset(lalonde,select=x_vars))
W_lalonde <-  diag(x=1, nrow(X_lalonde), nrow(X_lalonde))
forinv_lalonde <- diag(x=1,ncol(X_lalonde),ncol(X_lalonde))
XWX_inv_lalonde <- solve(t(X_lalonde) %*% W_lalonde  %*% X_lalonde,forinv_lalonde)
beta_hat_lalonde <- XWX_inv_lalonde %*% ( t(X_lalonde)  %*% W_lalonde  %*% lalonde$earn78)

h_0_hat_lalonde <- nrow(X_lalonde)^(-1) * (t(X_lalonde) %*% X_lalonde)
epsilon_hat_lalonde <- (lalonde$earn78 - X_lalonde %*% beta_hat_lalonde)
s2_lalonde <- (nrow(X_lalonde)-nrow(beta_hat_lalonde))^(-1)*(t(epsilon_hat_lalonde) %*% epsilon_hat_lalonde)

v_0_hat_lalonde <- (as.integer(s2_lalonde)* nrow(X_lalonde)^(-1)) * (t(X_lalonde) %*% X_lalonde)
h_0_hat_inv_lalonde <- solve(h_0_hat_lalonde,diag(x=1,ncol(X_lalonde)))
V_hat_lalonde <- h_0_hat_inv_lalonde %*% v_0_hat_lalonde %*% h_0_hat_inv_lalonde

##Find t-stat, p-value and CI (alpha = .95) against null of zero
t_numerator_lalonde <- beta_hat_lalonde - rep.int(0,nrow(beta_hat_lalonde))
asyvar_lalonde <- diag(V_hat_lalonde)
t_denominator_lalonde <- sqrt(asyvar_lalonde * nrow(X_lalonde)^(-1) )

t_stat_lalonde <- t_numerator_lalonde / t_denominator_lalonde

###Make Confidence Interval
alpha <- .95

df <- as.integer(nrow(lalonde) - ncol(lalonde))

lb_lalonde <- beta_hat_lalonde + qt((1-alpha)/2,df)*t_denominator_lalonde
ub_lalonde <-  beta_hat_lalonde + qt(alpha+(1-alpha)/2,df)*t_denominator_lalonde

CI_lalonde <- data.frame(lb_lalonde,ub_lalonde,rep.int(alpha,nrow(beta_hat_lalonde)),beta_hat_lalonde)


####a: Use Built in Function
beta_0_lalonde <- lm(earn78 ~ treat + black + age + educ + educ2 + earn74 + black_earn74 + u74 + u75, lalonde)
summary(beta_0_lalonde)

###The Different Estimates All Match

