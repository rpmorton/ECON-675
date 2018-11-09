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


set.seed(123)

export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

reps <- 1000
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
results_data$bygroupnull <- rep(1,reps)
ci_low <- round(reps * .025)
ci_high <- round(reps * .975)

sorted_beta2 <- sort(results_data$r_beta_2)
ci_low_beta2 <- sorted_beta2[ci_low]
ci_high_beta2 <- sorted_beta2[ci_high]
ci_data <- data.frame(ci_low_beta2,ci_high_beta2)
sorted_beta3 <- sort(results_data$r_beta_3)
ci_data$ci_low_beta3 <- sorted_beta3[ci_low]
ci_data$ci_high_beta3 <- sorted_beta3[ci_high]
sorted_beta_select <- sort(results_data$r_beta_select)
ci_data$ci_low_beta_select <- sorted_beta_select[ci_low]
ci_data$ci_high_beta_select <- sorted_beta_select[ci_high]

#Now Gather Summary Stats by Group
sim_results <-  ddply(results_data, .(bygroupnull), mutate
                    , mean_beta_2 = mean(r_beta_2)
                    , sd_beta_2 = sd(r_beta_2)
                    , mean_beta_3 = mean(r_beta_3)
                    , sd_beta_3 = sd(r_beta_3)
                    , mean_beta_select = mean(r_beta_select)
                    , sd_beta_select = sd(r_beta_select)
                    )[1,6:11]

plot_beta_2 <- ggplot(results_data,aes(x=r_beta_2)) + geom_density() 
plot_beta_2 <- plot_beta_2 + labs(title = expression(paste("PS4 Q3: Density of ", hat(beta))), y = "Density", x = "")
plot_beta_2

outpath <- paste(export,"/R_density_PS4_Q3_Beta_2.pdf", sep = "")
ggsave(outpath)

plot_beta_3 <- ggplot(results_data,aes(x=r_beta_3)) + geom_density() 
plot_beta_3 <- plot_beta_3 + labs(title = expression(paste("PS4 Q3: Density of ", tilde(beta))), y = "Density", x = "")
plot_beta_3

outpath <- paste(export,"/R_density_PS4_Q3_Beta_3.pdf", sep = "")
ggsave(outpath)

plot_beta_select <- ggplot(results_data,aes(x=r_beta_select)) + geom_density() 
plot_beta_select <- plot_beta_select + labs(title = expression(paste("PS4 Q3: Density of Selected ", beta)), y = "Density", x = "")
plot_beta_select

outpath <- paste(export,"/R_density_PS4_Q3_Beta_Select.pdf", sep = "")
ggsave(outpath)



































