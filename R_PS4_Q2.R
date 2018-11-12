#Created by RM on 2018.11.07
#For ECON 675, PS 4, Q 3
##############################  

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
library("CasualGAM")
library("robustbase")
library("lmtest")
library("sandwich")
library(MatchIt)
library(survey)

###Prepare Data

rm(list=ls())
export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

alpha <- .05

##Add Extra Vars
lalonde <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/LaLonde_all.csv")

lalonde$ln_re74 <- log(lalonde$re74 + 1)
lalonde$ln_re75 <- log(lalonde$re75 + 1)
lalonde$age2 <- lalonde$age * lalonde$age
lalonde$educ2 <- lalonde$educ * lalonde$educ
lalonde$age3 <- lalonde$age2 * lalonde$age
lalonde$black_u74 <- lalonde$black * lalonde$u74
lalonde$educ_ln_re74 <- lalonde$educ * lalonde$ln_re74

lalonde$missing_re78 <- is.na(lalonde$re78) 

lalonde_nopsid <- subset(lalonde,treat < 2, missing_re78 = FALSE)

OLSresults <- function(lm_model) {
  df <- nrow(coef(summary(lm_model)))
  pre_n <- as.data.frame(lm_model$fitted.values)
  n <- nrow(pre_n)
  results_hc1 <- coeftest(lm_model, vcov = vcovHC(lm_model, type="HC1"))
  results_output <- as.data.frame(results_hc1[2,1:4])
  results_output$beta_treat <- results_output[1,1]
  results_output$se <- results_output[2,1]
  results_output$low_ci <- results_output$beta_treat + qt(alpha/2,n - df) *  results_output$se
  results_output$high_ci <- results_output$beta_treat + qt(1 - (alpha/2),n - df) *  results_output$se
  return(results_output[1,2:5])
}

diff_means_model <- lm(re78~ treat, lalonde_nopsid)
export_diff_means_results <- OLSresults(diff_means_model)
export_diff_means_results$type <- "DiffMeans"
export_diff_means_results$spec <- ""


#OLS
  
#speca = "age educ black hisp married nodeg ln_earn74 ln_earn75"
#specb = "age2 educ2 u74 u75"
#specc ="age3 black_u74 educ_ln_earn74"

OLS_speca_model <- lm(re78~ treat + age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 
                       , lalonde_nopsid)
OLS_speca_results <- OLSresults(OLS_speca_model)

OLS_specb_model <- lm(re78~ treat + age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 
                      + age2 + educ2 + u74 + u75, lalonde_nopsid)
OLS_specb_results <- OLSresults(OLS_specb_model)

OLS_specc_model <- lm(re78~ treat + age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 
                      + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74
                      , lalonde_nopsid)
OLS_specc_results <- OLSresults(OLS_specc_model)

export_OLS_output <- rbind(OLS_speca_results,OLS_specb_results,OLS_specc_results)
export_OLS_output$type <- "OLS"
export_OLS_output$spec <- c("A","B","C")

#Reg Impute, IPW, and Doubly Robust Estimated using the CausalGAM package

outputGAM <- function(ATEresults) {
  type <- "RegImpute"
  beta <- ATEresults$ATE.reg.hat
  se <- ATEresults$ATE.reg.asymp.SE
  reg_impute <- data.frame(type,beta,se)
  reg_impute$low_ci <- reg_impute$beta + qnorm(alpha/2) *  reg_impute$se
  reg_impute$high_ci <- reg_impute$beta + qnorm(1 - (alpha/2)) *  reg_impute$se
  
  type <- "IPW"
  beta <- ATEresults$ATE.IPW.hat
  se <- ATEresults$ATE.IPW.asymp.SE
  IPW <- data.frame(type,beta,se)
  IPW$low_ci <- IPW$beta + qnorm(alpha/2) *  IPW$se
  IPW$high_ci <- IPW$beta + qnorm(1 - (alpha/2)) *  IPW$se
  
  type <- "DoublyRobust"
  beta <- ATEresults$ATE.AIPW.hat
  se <- ATEresults$ATE.AIPW.asymp.SE
  DR <- data.frame(type,beta,se)
  DR$low_ci <- DR$beta + qnorm(alpha/2) *  DR$se
  DR$high_ci <- DR$beta + qnorm(1 - (alpha/2)) *  DR$se
  
  output <- rbind(reg_impute,IPW,DR)
  return(output)
  
}

speca <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 ,
                     pscore.family = binomial(probit),
                     outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 ,
                     outcome.formula.c =  re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                     outcome.family = gaussian, 
                     treatment.var = "treat",
                     data = lalonde_nopsid, 
                     divby0.action = c("fail", "truncate", "discard"), 
                     divby0.tol = 0, nboot = 0, 
                     var.gam.plot = FALSE)
ate_speca <- outputGAM(speca)
ate_speca$spec <- "A"

specb <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 ,
                      pscore.family = binomial(probit),
                      outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 ,
                      outcome.formula.c =  re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75  + age2 + educ2 + u74 + u75,
                      outcome.family = gaussian, 
                      treatment.var = "treat",
                      data = lalonde_nopsid, 
                      divby0.action = c("fail", "truncate", "discard"), 
                      divby0.tol = 0, nboot = 0, 
                      var.gam.plot = FALSE)
ate_specb <- outputGAM(specb)
ate_specb$spec <- "B"


specc <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74 ,
                      pscore.family = binomial(probit),
                      outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                      outcome.formula.c =  re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75  + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                      outcome.family = gaussian, 
                      treatment.var = "treat",
                      data = lalonde_nopsid, 
                      divby0.action = c("fail", "truncate", "discard"), 
                      divby0.tol = 0, nboot = 0, 
                      var.gam.plot = FALSE)
ate_specc <- outputGAM(specc)
ate_specc$spec <- "C"

export_RegImpute_IPW_DR_output <- rbind(ate_speca, ate_specb, ate_specc)
colnames(export_RegImpute_IPW_DR_output)[2] <- "beta_treat"

##Matching Estimators

lalonde_nopsid$ctrl <- 1 - lalonde_nopsid$treat

lalonde_rownum <- lalonde_nopsid
lalonde_rownum$rownum <- seq(1:nrow(lalonde_rownum))
mergeon_re78 <- select(lalonde_rownum, re78, rownum)
colnames(mergeon_re78)[1] <- "matched_re78"

matchestimators <- function(treat_match,ctrl_match) {
  matched_data_treat <- match.data(treat_match, 'treat')
  matched_data_treat$mergeon <- treat_match$match.matrix
  matched_treat_mergeon <- merge(matched_data_treat, mergeon_re78, 
                                 by.x = "mergeon", by.y = "rownum")
  
  matched_data_ctrl <- match.data(ctrl_match, 'treat')
  matched_data_ctrl$mergeon <- ctrl_match$match.matrix
  matched_ctrl_mergeon <- merge(matched_data_ctrl, mergeon_re78, 
                                 by.x = "mergeon", by.y = "rownum")
  
  matched_data_all <- rbind(matched_treat_mergeon, matched_ctrl_mergeon)
  matched_data_all$te <- matched_data_all$treat*(matched_data_all$re78-matched_data_all$matched_re78) 
  matched_data_all$te_all <-  matched_data_all$te + (1-matched_data_all$treat) * (matched_data_all$matched_re78- matched_data_all$re78)
  
  matched_output <- select(matched_data_all, treat, te, te_all, re78, matched_re78)
  #return(matched_output)
  
  #Now Summarize for Output
  matched_output$dummy <- 1
  ate_model <- lm(te_all~ dummy - 1, matched_output)
  df <- nrow(coef(summary(ate_model)))
  ate_hc1 <- coeftest(ate_model, vcov = vcovHC(ate_model, type="HC1"))
  ate_output <- as.data.frame(ate_hc1[1,1:4])
  ate_output$beta_treat <- ate_output[1,1]
  ate_output$se <- ate_output[2,1]
  ate_output$low_ci <- ate_output$beta_treat + qt(alpha/2,nrow(matched_output) - df) *  ate_output$se
  ate_output$high_ci <- ate_output$beta_treat + qt(1 - (alpha/2),nrow(matched_output) - df) *  ate_output$se
  ate_output$esttype <- "ATE"
  
  att_model <- lm(te_all~ dummy - 1, subset(matched_output, treat > 0) )
  df <- nrow(coef(summary(att_model)))
  att_hc1 <- coeftest(att_model, vcov = vcovHC(att_model, type="HC1"))
  att_output <- as.data.frame(att_hc1[1,1:4])
  att_output$beta_treat <- att_output[1,1]
  att_output$se <- att_output[2,1]
  att_output$low_ci <- att_output$beta_treat + qt(alpha/2,nrow( subset(matched_output, treat > 0)) - df) *  att_output$se
  att_output$high_ci <- att_output$beta_treat + qt(1 - (alpha/2),nrow( subset(matched_output, treat > 0)) - df) *  att_output$se
  att_output$esttype <- "ATT"
  
  ate_att_output <- rbind(ate_output[1,2:6],att_output[1,2:6])
  return(ate_att_output)
  
}

#Nearest Neighbor Matching

match_treat_nn_a <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                     method = "nearest", distance = "mahalanobis", 
                     mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75"),
                     replace = T, data = lalonde_nopsid)
match_ctrl_nn_a <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                           method = "nearest", distance = "mahalanobis", 
                           mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75"),
                           replace = T, data = lalonde_nopsid)
matched_nn_a_output <- matchestimators(match_treat_nn_a,match_ctrl_nn_a)
matched_nn_a_output$spec <- "A"
matched_nn_a_output$mtype <- "NNeighbor_"
matched_nn_a_output$type <- paste0(matched_nn_a_output$mtype,matched_nn_a_output$esttype)

match_treat_nn_b <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75"),
                            replace = T, data = lalonde_nopsid)
match_ctrl_nn_b <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75"),
                            replace = T, data = lalonde_nopsid)
matched_nn_b_output <- matchestimators(match_treat_nn_b,match_ctrl_nn_b)
matched_nn_b_output$spec <- "B"
matched_nn_b_output$mtype <- "NNeighbor_"
matched_nn_b_output$type <- paste0(matched_nn_b_output$mtype,matched_nn_b_output$esttype)

match_treat_nn_c <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75","age3","black_u74","educ_ln_re74"),
                            replace = T, data = lalonde_nopsid)
match_ctrl_nn_c <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75","age3","black_u74","educ_ln_re74"),
                            replace = T, data = lalonde_nopsid)
matched_nn_c_output <- matchestimators(match_treat_nn_c,match_ctrl_nn_c)
matched_nn_c_output$spec <- "C"
matched_nn_c_output$mtype <- "NNeighbor_"
matched_nn_c_output$type <- paste0(matched_nn_c_output$mtype,matched_nn_c_output$esttype)

#Propensity Score Matching

match_treat_ps_a <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_nopsid)
match_ctrl_ps_a <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_nopsid)
matched_ps_a_output <- matchestimators(match_treat_ps_a,match_ctrl_ps_a)
matched_ps_a_output$spec <- "A"
matched_ps_a_output$mtype <- "PS_"
matched_ps_a_output$type <- paste0(matched_ps_a_output$mtype,matched_ps_a_output$esttype)

match_treat_ps_b <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_nopsid)
match_ctrl_ps_b <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_nopsid)
matched_ps_b_output <- matchestimators(match_treat_ps_b,match_ctrl_ps_b)
matched_ps_b_output$spec <- "B"
matched_ps_b_output$mtype <- "PS_"
matched_ps_b_output$type <- paste0(matched_ps_b_output$mtype,matched_ps_b_output$esttype)

match_treat_ps_c <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_nopsid)
match_ctrl_ps_c <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_nopsid)
matched_ps_c_output <- matchestimators(match_treat_ps_c,match_ctrl_ps_c)
matched_ps_c_output$spec <- "C"
matched_ps_c_output$mtype <- "PS_"
matched_ps_c_output$type <- paste0(matched_ps_c_output$mtype,matched_ps_c_output$esttype)

export_match_results <- rbind(matched_nn_a_output,matched_nn_b_output,matched_nn_c_output,
                      matched_ps_a_output,matched_ps_b_output,matched_ps_c_output)
export_match_results$esttype <- NULL
export_match_results$mtype <- NULL

#Combine All For Output to TeX

export_lalonde_output <- rbind(export_diff_means_results,export_OLS_output,
                               export_RegImpute_IPW_DR_output,export_match_results)
export_lalonde_output$beta_treat <- round(export_lalonde_output$beta_treat, digits = 2)
export_lalonde_output$se <- round(export_lalonde_output$se, digits = 2)
export_lalonde_output$low_ci <- round(export_lalonde_output$low_ci, digits = 2)
export_lalonde_output$high_ci <- round(export_lalonde_output$high_ci, digits = 2)

export_lalonde_nopsid <- export_lalonde_output

#outpath <- paste(export,"/R_PS4_Q2_LaLonde_no_psid.tex", sep = "")

#latex(export_lalonde_output,
#      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="LaLonde No PSID",
#      colheads=c("$\\hat\\tau$", "s.e", "95\\% C.I.", "95\\% C.I.", "Type", "Spec")
#)

#####
#LaLonde PSID Control
#####

lalonde_psid <- subset(lalonde, treat > 0, missing_re78 = FALSE)
lalonde_psid$treat_orig <- lalonde_psid$treat
lalonde_psid$treat <- ifelse(lalonde_psid$treat_orig > 1, 0, 1)

diff_means_model <- lm(re78~ treat, lalonde_psid)
export_diff_means_results <- OLSresults(diff_means_model)
export_diff_means_results$type <- "DiffMeans"
export_diff_means_results$spec <- ""

#OLS

#speca = "age educ black hisp married nodeg ln_earn74 ln_earn75"
#specb = "age2 educ2 u74 u75"
#specc ="age3 black_u74 educ_ln_earn74"

OLS_speca_model <- lm(re78~ treat + age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 
                      , lalonde_psid)
OLS_speca_results <- OLSresults(OLS_speca_model)

OLS_specb_model <- lm(re78~ treat + age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 
                      + age2 + educ2 + u74 + u75, lalonde_psid)
OLS_specb_results <- OLSresults(OLS_specb_model)

OLS_specc_model <- lm(re78~ treat + age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 
                      + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74
                      , lalonde_psid)
OLS_specc_results <- OLSresults(OLS_specc_model)

export_OLS_output <- rbind(OLS_speca_results,OLS_specb_results,OLS_specc_results)
export_OLS_output$type <- "OLS"
export_OLS_output$spec <- c("A","B","C")

#Reg Impute, IPW, and Doubly Robust Estimated using the CausalGAM package

speca <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 ,
                      pscore.family = binomial(probit),
                      outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 ,
                      outcome.formula.c =  re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                      outcome.family = gaussian, 
                      treatment.var = "treat",
                      data = lalonde_psid, 
                      divby0.action = c("discard"), 
                      divby0.tol = .01, nboot = 0, 
                      var.gam.plot = FALSE)
ate_speca <- outputGAM(speca)
ate_speca$spec <- "A"

specb <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 ,
                      pscore.family = binomial(probit),
                      outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 ,
                      outcome.formula.c =  re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75  + age2 + educ2 + u74 + u75,
                      outcome.family = gaussian, 
                      treatment.var = "treat",
                      data = lalonde_psid, 
                      divby0.action = c("discard"), 
                      divby0.tol = .05, nboot = 0, 
                      var.gam.plot = FALSE)
ate_specb <- outputGAM(specb)
ate_specb$spec <- "B"

specc <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74 ,
                      pscore.family = binomial(probit),
                      outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                      outcome.formula.c =  re78 ~ age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75  + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                      outcome.family = gaussian, 
                      treatment.var = "treat",
                      data = lalonde_psid, 
                      divby0.action = c("discard"), 
                      divby0.tol = .05, nboot = 0, 
                      var.gam.plot = FALSE)
ate_specc <- outputGAM(specc)
ate_specc$spec <- "C"

export_RegImpute_IPW_DR_output <- rbind(ate_speca, ate_specb, ate_specc)
colnames(export_RegImpute_IPW_DR_output)[2] <- "beta_treat"

##Matching Estimators

lalonde_psid$ctrl <- 1 - lalonde_psid$treat

lalonde_rownum <- lalonde_psid
lalonde_rownum$rownum <- seq(1:nrow(lalonde_rownum))
mergeon_re78 <- select(lalonde_rownum, re78, rownum)
colnames(mergeon_re78)[1] <- "matched_re78"

#Nearest Neighbor Matching

match_treat_nn_a <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75"),
                            replace = T, data = lalonde_psid)
match_ctrl_nn_a <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75"),
                            replace = T, data = lalonde_psid)
matched_nn_a_output <- matchestimators(match_treat_nn_a,match_ctrl_nn_a)
matched_nn_a_output$spec <- "A"
matched_nn_a_output$mtype <- "NNeighbor_"
matched_nn_a_output$type <- paste0(matched_nn_a_output$mtype,matched_nn_a_output$esttype)

match_treat_nn_b <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75"),
                            replace = T, data = lalonde_psid)
match_ctrl_nn_b <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75"),
                            replace = T, data = lalonde_psid)
matched_nn_b_output <- matchestimators(match_treat_nn_b,match_ctrl_nn_b)
matched_nn_b_output$spec <- "B"
matched_nn_b_output$mtype <- "NNeighbor_"
matched_nn_b_output$type <- paste0(matched_nn_b_output$mtype,matched_nn_b_output$esttype)

match_treat_nn_c <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75","age3","black_u74","educ_ln_re74"),
                            replace = T, data = lalonde_psid)
match_ctrl_nn_c <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "mahalanobis", 
                            mahvars=c("age","educ","black","hisp","married","nodegr","ln_re75","ln_re75","age2","educ2","u74","u75","age3","black_u74","educ_ln_re74"),
                            replace = T, data = lalonde_psid)
matched_nn_c_output <- matchestimators(match_treat_nn_c,match_ctrl_nn_c)
matched_nn_c_output$spec <- "C"
matched_nn_c_output$mtype <- "NNeighbor_"
matched_nn_c_output$type <- paste0(matched_nn_c_output$mtype,matched_nn_c_output$esttype)

#Propensity Score Matching

match_treat_ps_a <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_psid)
match_ctrl_ps_a <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_psid)
matched_ps_a_output <- matchestimators(match_treat_ps_a,match_ctrl_ps_a)
matched_ps_a_output$spec <- "A"
matched_ps_a_output$mtype <- "PS_"
matched_ps_a_output$type <- paste0(matched_ps_a_output$mtype,matched_ps_a_output$esttype)

match_treat_ps_b <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_psid)
match_ctrl_ps_b <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_psid)
matched_ps_b_output <- matchestimators(match_treat_ps_b,match_ctrl_ps_b)
matched_ps_b_output$spec <- "B"
matched_ps_b_output$mtype <- "PS_"
matched_ps_b_output$type <- paste0(matched_ps_b_output$mtype,matched_ps_b_output$esttype)

match_treat_ps_c <- matchit(treat ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_psid)
match_ctrl_ps_c <- matchit( ctrl ~  age + educ + black + hisp + married + nodegr + ln_re74 + ln_re75 + age2 + educ2 + u74 + u75 + age3 + black_u74 + educ_ln_re74,
                            method = "nearest", distance = "probit", 
                            replace = T, data = lalonde_psid)
matched_ps_c_output <- matchestimators(match_treat_ps_c,match_ctrl_ps_c)
matched_ps_c_output$spec <- "C"
matched_ps_c_output$mtype <- "PS_"
matched_ps_c_output$type <- paste0(matched_ps_c_output$mtype,matched_ps_c_output$esttype)

export_match_results <- rbind(matched_nn_a_output,matched_nn_b_output,matched_nn_c_output,
                              matched_ps_a_output,matched_ps_b_output,matched_ps_c_output)
export_match_results$esttype <- NULL
export_match_results$mtype <- NULL

#Combine All For Output to TeX

export_lalonde_output <- rbind(export_diff_means_results,export_OLS_output,
                               export_RegImpute_IPW_DR_output,export_match_results)
export_lalonde_output$beta_treat <- round(export_lalonde_output$beta_treat, digits = 2)
export_lalonde_output$se <- round(export_lalonde_output$se, digits = 2)
export_lalonde_output$low_ci <- round(export_lalonde_output$low_ci, digits = 2)
export_lalonde_output$high_ci <- round(export_lalonde_output$high_ci, digits = 2)

export_lalonde_psid <- export_lalonde_output

#outpath <- paste(export,"/R_PS4_Q2_LaLonde_psid.tex", sep = "")

#latex(export_lalonde_output,
#      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="LaLonde PSID",
#      colheads=c("$\\hat\\tau$", "s.e", "95\\% C.I.", "95\\% C.I.", "Type", "Spec")
#)


##Combine Results for Output 

colnames(export_lalonde_psid)[1] <- "psid_beta_treat"
colnames(export_lalonde_psid)[2] <- "psid_se"
colnames(export_lalonde_psid)[3] <- "psid_low_ci"
colnames(export_lalonde_psid)[4] <- "psid_high_ci"

colnames(export_lalonde_nopsid)[1] <- "nopsid_beta_treat"
colnames(export_lalonde_nopsid)[2] <- "nopsid_se"
colnames(export_lalonde_nopsid)[3] <- "nopsid_low_ci"
colnames(export_lalonde_nopsid)[4] <- "nopsid_high_ci"

all_lalonde_output <- merge(export_lalonde_nopsid, export_lalonde_psid, by = c("type", "spec") )

outpath <- paste(export,"/R_PS4_Q2_LaLonde.tex", sep = "")

latex(all_lalonde_output,
      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="LaLonde Results",
      n.cgroup(2,4,4), cgroup = c("Details","No PSID","PSID"), 
      colheads=c("Type", "Spec", "$\\hat\\tau$", "s.e", "95\\% C.I.", "95\\% C.I.", "$\\hat\\tau$", "s.e", "95\\% C.I.", "95\\% C.I.")
)

outpathcsv <- paste0(export,"/R_PS4_Q2_LaLonde.csv")
write.csv(all_lalonde_output, file = outpathcsv)



