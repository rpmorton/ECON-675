#####
#Created by RM on 2018.10.28
#Econ 675: PS 3, Q 1
#####


library(MASS)
library(plyr)
library(data.table)
library(matlib)
library(randomizr)
library(tidyverse)
library(boot)
library(ggplot2)
library(gmm)
library(foreign)
library(xtable)
library(Hmisc)


export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"


reps <- 200
obs <- 50

r_beta_2 <- rep(-50,reps)
r_zscore_z <- rep(-50,reps)
r_beta_3 <- rep(-50,reps)
r_beta_select <- rep(-50,reps)
results_data <- data.frame(r_beta_2,r_zscore_z,r_beta_3,r_beta_select)

for (j in 1:reps) {

  #Create Data
  epsilon <- rnorm(obs)
  x <- rnorm(obs)
  secondnorm <- rnorm(obs)
  z <- .85 * x + sqrt(1 - (.85 * .85) ) * secondnorm
  y <- 1 + .5 * x + z + epsilon
  seq <- (1:50)
  data <- data.frame(y,x,z,seq)
  
  #Run OLS
  ols <- lm(y ~  x + z ,data=data)
  olsresults <- coef(summary(ols))
  beta_2 <- olsresults[2,1]
  beta_z <- olsresults[3,1]
  se_z <- olsresults[3,2]
  zscore_z <- abs(beta_z/se_z)
  
  results_data[j,1] <- beta_2
  results_data[j,2] <- zscore_z
  
  olsalt <- lm(y ~  x ,data=data)
  olsaltresults <- coef(summary(olsalt))
  beta_3 <- olsaltresults[2,1]
  
  results_data[j,3] <- beta_3

}
  
results_data$r_beta_select <- ifelse(results_data$r_zscore_z >= 1.96,
                                 results_data$r_beta_2, results_data$r_beta_3)











































