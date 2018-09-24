#Created by RM on 2018.09.14 for ECON 675 PS 1, Q 3
#Subpart 1b

library(MASS)
library(plyr)
library(data.table)
library(matlib)
library(randomizr)

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
ub_diff_mean_lalonde <- treated_diff_mean + qt(alpha+(1-alpha)/2,df)*sqrt(variance)

##Question 2

#First Generate Random Vectors
rm(list=ls())
lalonde <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/LaLonde_1986.csv")

count_treated <- nrow(subset(lalonde, treat > 0))
count_control <- nrow(subset(lalonde, treat < 1))
count_total <- nrow(lalonde)

#
rand_vec <- complete_ra(N = count_total, m = count_treated)
#rand_vec2 <-  complete_ra(N = count_total, m = count_treated)
randomvectors <- data.frame(rand_vec)

loops <- 500
loopcounter <- 1
loopsused <-1 

while (loopcounter < loops) {
  addrandvec <- complete_ra(N = count_total, m = count_treated)
  loopsused <- loopsused + 1
  
  total_match <- 0
  for (i in 1:loopcounter) {
    pre_total_match <- sum(addrandvec == randomvectors[,i] ) 
    total_match <- max(total_match,pre_total_match)
  }
  
  if (total_match < count_total) {
    randomvectors <- cbind(addrandvec,randomvectors)
    loopcounter <- loopcounter + 1
  }
  
}

lalonde$obs <- as.integer(1:nrow(randomvectors))

fishertreat_mean <- rep.int(-10,loops)
fisherctrl_mean <- rep.int(-10,loops)
fisher_ks <- rep.int(-10,loops)

fisherresults <- data.frame(fishertreat_mean,fisherctrl_mean,fisher_ks)
fisherresults$iteration <- as.integer(1:nrow(fisherresults))

for (j in 1:loops) {
  
  #Merge on randvecuse
  randvecuse <-as.data.frame(randomvectors[,j])
  randvecuse$obs <- as.integer(1:nrow(randomvectors))
  lalonde_rand <- NULL
  lalonde_rand <- merge(lalonde,randvecuse)
  
  #diff in means
  means <- ddply(lalonde_rand,.(randomvectors[, j]),summarize, mean = mean(earn78))
  mean_fish_treat <- means[which(means[1,] > .75),]
  mean_fish_treat <- na.omit(mean_fish_treat)
  mean_fish_ctrl <- means[which(means[1,] < .5),]
  mean_fish_ctrl <- na.omit(mean_fish_ctrl)
  
  fisherresults[j,1] <- as.numeric(mean_fish_treat$mean)
  fisherresults[j,2] <- as.numeric(mean_fish_ctrl$mean)
  
  #ks
  outcomes_fish_treat <- na.omit(subset(lalonde_rand,randomvectors[, j] > .5))
  outcomes_fish_ctrl <- na.omit(subset(lalonde_rand,randomvectors[, j] < .5))
  ks_fisher <- ks.test(outcomes_fish_treat$earn78,outcomes_fish_ctrl$earn78)
  
  fisherresults[j,3] <- ks_fisher[["statistic"]]
  
}


act_treat_means <- ddply(lalonde_rand,.(treat),summarize, mean = mean(earn78))
mean_treat_act <-  na.omit(act_treat_means[which(act_treat_means[1,] > .75),])
mean_ctrl_act <-  na.omit(act_treat_means[which(act_treat_means[1,] < .75),])
diff_means_act <- mean_treat_act$mean - mean_ctrl_act$mean



  