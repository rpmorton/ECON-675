#Created by RM on 2018.09.14 for ECON 675 PS 1, Q 3
#Subpart 1b

library(MASS)
library(plyr)
library(data.table)
library(matlib)

rm(list=ls())

lalonde <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/LaLonde_1986.csv")

lalonde_treated <- subset(lalonde, treat > 0)
lalonde_control <- subset(lalonde, treat < 1)

diff_means <- function(x = NULL, y = NULL) {
  mean_earn78 <- mean(y)
  pre_nobs_data <- nrow(as.data.frame(y))
  numobs <- as.integer(pre_nobs_data)
  
  deviations <- y - mean_earn78
  deviation_data <- data.frame(deviations)
  deviation_matrix <- as.matrix(deviation_data)
  sum_squared_devs <- t(deviation_matrix) %*% deviation_matrix 
  pre_s2 <- sum_squared_devs[1,1]
  s2 <- pre_s2 * (1/(numobs-1))
 
  out_list <- list() 
  out_list[["nobs"]] <- numobs
  out_list[["mean"]] <- mean_earn78
  out_list[["s2"]] <- s2
  
  return(out_list)
  
}

earn78_treated <- lalonde_treated$earn78
treated_list <- diff_means(lalonde_treated,earn78_treated) 

earn78_control <- lalonde_control$earn78
control_list <- diff_means(lalonde_control,earn78_control)
  
treated_diff_mean <- treated_list[["mean"]] - control_list[["mean"]]
variance <- (treated_list[["s2"]]/treated_list[["nobs"]]) + (control_list[["s2"]]/control_list[["nobs"]])
  
alpha <- .95

df <- as.integer(nrow(lalonde) - 2)

lb_diff_mean_lalonde <- treated_diff_mean + qt((1-alpha)/2,df)*sqrt(variance)
ub_diff_mean_lalonde <-  treated_diff_mean + qt(alpha+(1-alpha)/2,df)*sqrt(variance)




  