#####
#Created by RM on 2018.11.14
#Econ 675: PS 5, Q 2
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
library(sem)
library(AER)

##
#parameters
#sims <- 5000
sims <- 50
obs <- 200

cov <- .99 

gamma2 <-(1/obs) *  c(0, .25, 9, 99)
beta <- 0

#Now Fill Empty Results Table
beta_OLS <- rep(-100,sims*length(gamma2))
se_beta_OLS <- rep(-100,sims*length(gamma2))
beta_OLS_sig <- rep(-100,sims*length(gamma2))
beta_2SLS <- rep(-100,sims*length(gamma2))
se_beta_2SLS <- rep(-100,sims*length(gamma2))
beta_2SLS_sig <- rep(-100,sims*length(gamma2))
F_2SLS <- rep(-100,sims*length(gamma2))
gamma <- rep(-100,sims * length(gamma2))
simcount <- rep(-100,sims*length(gamma2))

resultsdata <- data.frame(beta_OLS, se_beta_OLS, beta_OLS_sig, beta_2SLS, se_beta_2SLS,
                          beta_2SLS_sig, F_2SLS, gamma, simcount)

len_gam2 <- length(gamma2)*1 

for(i in 1:sims) {
  for(g in 1:len_gam2) {
    z <- rnorm(obs,0,1)
    u <- rnorm(obs,0,1)
    indep_norm <- rnorm(obs,0,1)
    weight_indep_norm <- sqrt(1- (cov * cov) )
    v <- cov * u + weight_indep_norm * indep_norm
  
    gamma <- sqrt(gamma2[g])
    
    x <- gamma * z + v
    y <- beta * x + u
    regdata <- data.frame(x,y,z)
    
    row <- 4 * (i-1) + g
    
    olsmodel <- lm(y ~  x ,data=regdata)
    olsresults <- coef(summary(olsmodel))
    resultsdata[row,1] <- olsresults[2,1]
    resultsdata[row,2] <- olsresults[2,2]
    resultsdata[row,3] <- ifelse(abs(olsresults[2,1])/olsresults[2,2] >= 1.96, 1, 0)
    
    tslsmodel <- ivreg(y ~ x | z )
    tslsresults <- coef(summary(tslsmodel))
    resultsdata[row,4] <- tslsresults[2,1]
    resultsdata[row,5] <- tslsresults[2,2]
    resultsdata[row,6] <- ifelse(abs(tslsresults[2,1])/tslsresults[2,2] >= 1.96, 1, 0)
    
    
    resultsdata[row,8] <- gamma
    resultsdata[row,9] <- row
    
    #tslsmodel_alt <- tsls(y ~ x, ~ z )
    #tsls_test <- tslsmodel_alt$V
  }
}

#bydummy <- rep(1,200)


#Check Data Generated Correctly
#checkdata <- data.frame(z,u,v,bydummy)
#data_summary_mean <- ddply(checkdata, .(bydummy), summarise, 
#                      mean_z = mean(z), mean_u = mean(u), mean_v = mean(v) )  
#data_summary_sd <- ddply(checkdata, .(bydummy), summarise, 
#                      sd_z = sd(z), sd_u = sd(u), sd_v = sd(v) )
#data_summary_cov <-  ddply(checkdata, .(bydummy), summarise, 
#                           cov_z_u = cov(z,u), cov_z_v = cov(z,v), cov_v_u = cov(v,u) )


          
