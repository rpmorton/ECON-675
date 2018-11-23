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
sims <- 5000
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
    regdata$dummy <- 1
    
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
    
    #Do F Test

    instr_included <- lm(x ~ z, data = regdata)
    F_val <- summary(instr_included)$fstatistic[1]

    resultsdata[row,7] <- F_val
    
    resultsdata[row,8] <- gamma * gamma * obs
    resultsdata[row,9] <- row

  }
}


summary <- ddply(resultsdata, .(gamma), summarise,
                mean_beta_OLS = mean(beta_OLS), mean_se_beta_OLS = mean(se_beta_OLS),
                mean_beta_OLS_sig = mean(beta_OLS_sig),
                mean_beta_2SLS = mean(beta_2SLS), mean_se_beta_2SLS = mean(se_beta_2SLS),
                mean_beta_2SLS_sig = mean(beta_2SLS_sig), mean_F_2SLS = mean(F_2SLS),
                
                sd_beta_OLS = sd(beta_OLS), sd_se_beta_OLS = sd(se_beta_OLS),
                sd_beta_OLS_sig = sd(beta_OLS_sig),
                sd_beta_2SLS = sd(beta_2SLS), sd_se_beta_2SLS = sd(se_beta_2SLS),
                sd_beta_2SLS_sig = sd(beta_2SLS_sig), sd_F_2SLS = sd(F_2SLS),
                
                quant_beta_OLS = quantile(beta_OLS,.1 ) 
)
              
##Summarize Data; Then Rotate To Output

#Summary

resultstosummarize <- resultsdata
resultstosummarize$simcount <- NULL

results_mean <- resultstosummarize %>% group_by(gamma) %>% summarise_all(funs(mean))
tomerge_results_mean <- melt(results_mean, id.vars = "gamma", variable.name = "variable", value.name = "mean")

results_sd <- resultstosummarize %>% group_by(gamma) %>% summarise_all(funs(sd ))
tomerge_results_sd <- melt(results_sd, id.vars = "gamma", variable.name = "variable", value.name = "sd")

results_q10 <- resultstosummarize %>% group_by(gamma) %>% summarise_all(funs(quantile), probs = .1   )
tomerge_results_q10 <- melt(results_q10, id.vars = "gamma", variable.name = "variable", value.name = "q10")

results_q50 <- resultstosummarize %>% group_by(gamma) %>% summarise_all(funs(quantile), probs = .5   )
tomerge_results_q50 <- melt(results_q50, id.vars = "gamma", variable.name = "variable", value.name = "q50")

results_q90 <- resultstosummarize %>% group_by(gamma) %>% summarise_all(funs(quantile), probs = .9   )
tomerge_results_q90 <- melt(results_q90, id.vars = "gamma", variable.name = "variable", value.name = "q90")

results_export <- merge(tomerge_results_mean, merge(tomerge_results_sd,
                                                merge(tomerge_results_q10
                                                      , merge(tomerge_results_q50, tomerge_results_q90)) ) )

variables <- data.frame(unique(select(results_export, variable)))
variables$sort_order <- seq(1:nrow(variables))
variables$sort_order <- ifelse(variables$variable == "beta_OLS", 1, variables$sort_order)
variables$sort_order <- ifelse(variables$variable == "se_beta_OLS", 2, variables$sort_order)
variables$sort_order <- ifelse(variables$variable == "beta_OLS_sig", 3, variables$sort_order)
variables$sort_order <- ifelse(variables$variable == "beta_2SLS", 4, variables$sort_order)
variables$sort_order <- ifelse(variables$variable == "se_beta_2SLS", 5, variables$sort_order)
variables$sort_order <- ifelse(variables$variable == "beta_2SLS_sig", 6, variables$sort_order)
variables$sort_order <- ifelse(variables$variable == "F_2SLS", 7, variables$sort_order)

merge_sort_order <- merge(results_export, variables)

resultssorted <- merge_sort_order[ order(merge_sort_order$gamma, merge_sort_order$sort_order), ]


