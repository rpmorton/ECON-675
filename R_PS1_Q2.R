
#Created by RM on 2018.09.14 for ECON 675 PS 1, Q 2
#Subpart 4

#load package
library(MASS)
library(plyr)
library(data.table)


########
######Q2, Part 4: Implementation of OLS
#######

######Generate Data

x_1 <- runif(1000,0,100)
x_2 <- runif(1000,min=30,max=50)
x_1_2 <- c(x_1,x_2)
X <- matrix(x_1_2, nrow = 1000, ncol = 2)

epsilon <- rnorm(1000, mean=0, sd=10)

beta <- matrix( c(1,5), nrow = 2, ncol = 1)

y <- X %*% beta + epsilon

W <-  diag(x=1, 1000, 1000)


#####i: Symmetric Inverse

####i)a:Estimate Beta and Variance
X_W_X_inv <- ginv(t(X) %*% W  %*% X)
beta_hat <- X_W_X_inv %*% ( t(X)  %*% W  %*% y)

h_0_hat <- nrow(X)^(-1) * (t(X) %*% X)
epsilon_hat <- (y - X %*% beta_hat)
x_1_epsilon <- epsilon_hat * x_1
X_esp <- cbind(x_1_epsilon,x_2)
#Note: not squaring since will be squared in matrix multiplication
v_0_hat <- nrow(X)^(-1) * (t(X_esp) %*% X_esp)

V_hat <- ginv(h_0_hat) %*% v_0_hat %*% ginv(h_0_hat)


####i)b: T-test by element of beta
t_numerator <- beta_hat - beta
asyvar <- diag(V_hat)
t_denominator <- sqrt(asyvar * nrow(X)^(-1) )

t_stat <- t_numerator / t_denominator

####i)c: p-value
pval <- 2*pt(-abs(t_stat),df=nrow(X)-nrow(beta_hat))

####i)d: Confidence Interal
alpha <- .90

lb <- beta_hat + qnorm((1-alpha)/2,0,1)*t_denominator
ub <-  beta_hat + qnorm(alpha+(1-alpha)/2,0,1)*t_denominator

CI <- data.frame(lb,ub,alpha,beta_hat,beta)


#####ii: Choleskey Inverse

####ii)a:Estimate Beta and Variance with Cholesky
X_W_X_chol <- chol(t(X) %*% W  %*% X)
X_W_X_cholinv <- chol2inv(X_W_X_chol)
beta_hat_chol <- X_W_X_cholinv %*% ( t(X)  %*% W  %*% y)

h_0_hat <- nrow(X)^(-1) * (t(X) %*% X)
epsilon_hat_chol <- (y - X %*% beta_hat_chol)
x_1_epsilon_chol <- epsilon_hat_chol * x_1
X_esp_chol <- cbind(x_1_epsilon_chol,x_2)
#Note: not squaring since will be squared in matrix multiplication
v_0_hat_chol <- nrow(X)^(-1) * (t(X_esp_chol) %*% X_esp_chol)

V_hat_chol <- ginv(h_0_hat) %*% v_0_hat_chol %*% ginv(h_0_hat)


####ii)b: T-test by element of beta with Cholesky
t_numerator_chol <- beta_hat_chol - beta
asyvar_chol <- diag(V_hat_chol)
t_denominator_chol <- sqrt(asyvar_chol * nrow(X)^(-1) )

t_stat_chol <- t_numerator_chol / t_denominator_chol

####ii)c: p-value
pval_chol <- 2*pt(-abs(t_stat_chol),df=nrow(X)-nrow(beta_hat_chol))

####i)d: Confidence Interal
alpha <- .90

lb_chol <- beta_hat_chol + qnorm((1-alpha)/2,0,1)*t_denominator_chol
ub_chol <-  beta_hat_chol + qnorm(alpha+(1-alpha)/2,0,1)*t_denominator_chol

CI_chol <- data.frame(lb_chol,ub_chol,alpha,beta_hat_chol,beta)


###Up to at least 7 digits, the Symmetric and Cholesky Confidence Intervals
###(which are based on beta and variance estimators) are identical!

########
######Q2, Part 5: Lalonde Data
######

####a: Use Estimator Above

lalonde <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/LaLonde_1986.csv")
 
lalonde$educ2 <- lalonde$educ * lalonde$educ
lalonde$black_earn74 <- lalonde$black * lalonde$earn74
lalonde$intercept <- rep.int(1,nrow(lalonde))

x_vars <- c("intercept","treat","black","age","educ","educ2","earn74","black_earn74","u74","u75")
X_lalonde <- as.matrix(subset(lalonde,select=x_vars))
#X_lalonde <- ta(lalonde$treat,lalonde$black,lalonde$age,)
W_lalonde <-  diag(x=1, nrow(X_lalonde), nrow(X_lalonde))
XWX_inv_lalonde <- ginv(t(X_lalonde) %*% W_lalonde  %*% X_lalonde)
beta_hat_lalonde <- XWX_inv_lalonde %*% ( t(X_lalonde)  %*% W_lalonde  %*% lalonde$earn78)
#test <- ginv(t(X_lalonde) %*% X_lalonde) %*% ( t(X_lalonde) %*% lalonde$earn78)

h_0_hat_lalonde <- nrow(X_lalonde)^(-1) * (t(X_lalonde) %*% X_lalonde)
epsilon_hat_lalonde <- (lalonde$earn78 - X_lalonde %*% beta_hat_lalonde)
X_lalonde_table <- data.table(X_lalonde)
X_lalonde_table$age <- epsilon_hat_lalonde * X_lalonde_table$age
#X_esp <- cbind(x_1_epsilon,x_2)
#Note: not squaring since will be squared in matrix multiplication
v_0_hat <- nrow(X)^(-1) * (t(X_esp) %*% X_esp)

V_hat <- ginv(h_0_hat) %*% v_0_hat %*% ginv(h_0_hat)



beta_0_lalonde <- lm(earn78 ~ treat + black + age + educ + educ2 + earn74 + black_earn74 + u74 + u75, lalonde)

summary(beta_0_lalonde)


