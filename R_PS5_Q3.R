#####
#Created by RM on 2018.11.23
#Econ 675: PS 5, Q 3
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
library(sqldf)

rm(list=ls())

angkreug <-read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Angrist_Krueger.csv")
export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

##3.1: Angrist Kreuger ##


##Create Dummy Variables
angkreug$region.f <- factor(angkreug$region)
angkreug$YoB_ld.f <- factor(angkreug$YoB_ld)

##OLS 1 Estimates

olsmodel1 <- lm(l_w_wage ~  educ + non_white + married + SMSA + region.f + YoB_ld.f,data=angkreug)
ols1results_all <- summary(olsmodel1)
ols1results <- coef(ols1results_all)

##OLS 2 Estimates

olsmodel2 <- lm(l_w_wage ~  educ + age + age_q + non_white + married + SMSA + region.f + YoB_ld.f,data=angkreug)
ols2results_all <- summary(olsmodel2)
ols2results <- coef(ols2results_all)

##2SLS 1 Estimates
angkreug$QoB.f <- as.factor(angkreug$QoB)
angkreug$QoB_YoB_ld.f <- with(angkreug, interaction(QoB.f,  YoB_ld.f), drop = FALSE )

tslsmodel1 <- ivreg(l_w_wage ~  educ + non_white + married + SMSA + region.f + YoB_ld.f | QoB_YoB_ld.f + non_white + married + SMSA + region.f + YoB_ld.f , data = angkreug)
tsls1results <-  as.data.frame(coef(summary(tslsmodel1)))
#colnames(tsls1results) <- paste( colnames(tsls1results),"1", sep = "_")
tsls1results$model <- "Model1"
tsls1results$vars <- rownames(tsls1results)
tsls1results_subset <- subset(tsls1results, tsls1results$vars %in% c("(Intercept)", "educ", "non_white", "married", "SMSA", "age", "age_q") )

##2SLS 2 Estimates

tslsmodel2 <- ivreg(l_w_wage ~  educ + age + age_q + non_white + married + SMSA + region.f + YoB_ld.f | QoB_YoB_ld.f + age + age_q + non_white + married + SMSA + region.f + YoB_ld.f , data = angkreug)
tsls2results <- as.data.frame(coef(summary(tslsmodel2)))
#colnames(tsls2results) <- paste( colnames(tsls2results),"2", sep = "_")
tsls2results$model <- "Model2"
tsls2results$vars <- rownames(tsls2results)
tsls2results_subset <- subset(tsls2results, tsls2results$vars %in% c("(Intercept)", "educ", "non_white", "married", "SMSA", "age", "age_q") )

tsls_results <- rbind(tsls1results_subset, tsls2results_subset )

outpath <- paste(export,"/R_PS5_Q3.1.csv", sep = "")
write.csv(tsls_results,outpath)

##3.2: Bound et al. ##

reps <- 500

##2SLS 1: Permute

permute.2SLS1 <- function(x, ind) {
  permsample  <- angkreug[ind,] # allows boot to select sample 
  permuted_date <- select(permsample, c("QoB.f",  "YoB_ld.f") )
  angkreug_perm <- angkreug
  angkreug_perm$perm_QoB_YoB_ld.f <-  with(permuted_date, interaction(QoB.f,  YoB_ld.f), drop = TRUE )
  
  tslsmodel1_perm <- ivreg(l_w_wage ~  educ + non_white + married + SMSA + region.f + YoB_ld.f | perm_QoB_YoB_ld.f + non_white + married + SMSA + region.f + YoB_ld.f , data = angkreug_perm)
  tsls1results_perm <- as.data.frame(coef(summary(tslsmodel1_perm)))
  tsls1results_perm$model <- "Model1"
  tsls1results_perm$vars <- rownames(tsls1results_perm)
  tsls1results_perm <- subset(tsls1results_perm, tsls1results_perm$vars %in% c("educ") )
  
  tslsmodel2_perm <- ivreg(l_w_wage ~  educ + + age + age_q + non_white + married + SMSA + region.f + YoB_ld.f | perm_QoB_YoB_ld.f + age + age_q + non_white + married + SMSA + region.f + YoB_ld.f , data = angkreug_perm)
  tsls2results_perm <- as.data.frame(coef(summary(tslsmodel2_perm)))
  tsls2results_perm$model <- "Model2"
  tsls2results_perm$vars <- rownames(tsls2results_perm)
  tsls2results_perm <- subset(tsls2results_perm, tsls2results_perm$vars %in% c("educ") )
  
  tsls_perm_results <- rbind(tsls1results_perm, tsls2results_perm)
  tsls_perm_est <- tsls_perm_results$Estimate
  temp <- matrix(NaN, ncol=1, nrow=2)
  temp[1,1] <- tsls_perm_est[1]
  temp[2,1] <- tsls_perm_est[2]
  return(temp)
  
}

permute_2sls <- boot(data = angkreug, R = reps, statistic = permute.2SLS1, sim = "permutation", stype = "i")

permute_2sls_t <- permute_2sls$t

mean_2sls_model1 <- mean(permute_2sls_t[,1])
sd_2sls_model1 <- sd(permute_2sls_t[,1])

mean_2sls_model2 <- mean(permute_2sls_t[,2])
sd_2sls_model2 <- sd(permute_2sls_t[,2])

means <- c(mean_2sls_model1, mean_2sls_model2)
sigmas <- c(sd_2sls_model1, sd_2sls_model2)
models <- c("Model1","Model2")

bound_et_al_export <- data.frame(models,means,sigmas)

outpath <- paste(export,"/R_PS5_Q3.2.xlsx", sep = "")
write.csv(bound_et_al_export,outpath)





