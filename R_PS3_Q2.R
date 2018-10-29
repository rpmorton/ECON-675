#####
#Created by RM on 2018.10.28
#Econ 675: PS 3, Q 2
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
library(Hmisc)

rm(list=ls())

pisofirme <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/pisofirme.csv")
set.seed(123)

reps <- 100

export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"


####
#Question 2b
####

notmissing <- 1 - pisofirme$dmissing
s_age <- pisofirme$S_age
s_hhpeople <- pisofirme$S_HHpeople
ln_s_income <- log(pisofirme$S_incomepc + 1)
dpisofirme <- pisofirme$dpisofirme
danemia <- pisofirme$danemia

allrowsdata <- data.frame(danemia,notmissing,s_age,s_hhpeople,ln_s_income,dpisofirme)
mcardrop <- allrowsdata[complete.cases(allrowsdata),]

mcarmoments <- function(theta, data) {
  mom1 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$dpisofirme
  mom2 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$s_age
  mom3 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$s_hhpeople
  mom4 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$ln_s_income
  cbind(mom1, mom2, mom3, mom4)
}

initial_mcar_2b <- gmm(mcarmoments, mcardrop, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef

#Now Bootstrap

boot.GMM_mcar <- function(x, ind) {
  d <-mcardrop[ind,] # allows boot to select sample 
  temp <- gmm(mcarmoments, d, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef
  return(temp)
}

bootresult_2b <- boot(data = mcardrop, R = reps, statistic = boot.GMM_mcar)

#Now find the CI from here
boot_results_2b <- bootresult_2b$t
meltboot_2b <- melt(boot_results_2b)
meltboot_sorted_2b <- meltboot_2b[order(meltboot_2b$Var2,meltboot_2b$value),]
meltboot_sorted_2b$seq <- seq(1: nrow(meltboot_sorted_2b))
meltboot_sorted_2b$seq_adj <- meltboot_sorted_2b$seq - (reps * (meltboot_sorted_2b$Var2 - 1) )

lower_ci = round(reps * .05)
upper_ci = round(reps * .95)

lower_ci_2b <-  meltboot_sorted_2b[meltboot_sorted_2b$seq_adj == lower_ci,]
upper_ci_2b <- meltboot_sorted_2b[meltboot_sorted_2b$seq_adj == upper_ci,]

varnames <- c("dpisofirme","S_age","S_HHpeople", "log(S_incomepc+1)")

outtable <- data.frame(varnames,round(initial_mcar_2b,4), round(lower_ci_2b$value,4),  round(upper_ci_2b$value,4))
outpath <- paste(export,"/R_table_PS3_Q2_2b.tex", sep = "")
latex(outtable,
      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(2, 2),
      cgroup=c("Estimate", "Inference"),
      colheads=c("Variable", "$\\hat\\beta$", "95\\% C.I.", "95\\% C.I.")
)


####
#Question 3c
####

allrows_subset1 <- subset(allrowsdata,!is.na(allrowsdata$dpisofirme))
allrows_subset2 <- subset(allrowsdata,!is.na(allrowsdata$ln_s_income))


logitmodel_missing <- glm(notmissing ~ dpisofirme + ln_s_income + s_hhpeople + s_age
                          ,family=binomial(link='logit'),data=allrows_subset2)
allrows_predval <- allrows_subset2
allrows_predval$logitpredval <- predict(logitmodel_missing, type = "response")

marmoments <- function(theta, data) {
  mom1 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$dpisofirme * data$notmissing / data$logitpredval
  mom2 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$s_age * data$notmissing / data$logitpredval
  mom3 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$s_hhpeople * data$notmissing / data$logitpredval
  mom4 <- (data$danemia - plogis(theta[1]*data$dpisofirme + theta[2]*data$s_age + theta[3]*data$s_hhpeople + theta[4]*data$ln_s_income)) * 
    data$ln_s_income * data$notmissing / data$logitpredval
  cbind(mom1, mom2, mom3, mom4)
}

forgmm <- allrows_predval[complete.cases(allrows_predval),]
inital_mar_3c <- gmm(marmoments, forgmm, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef
#inital_mar_3c <- gmm(marmoments, allrows_predval, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef

boot.GMM_mar_3c <- function(x, ind) {
  d <<- allrows_subset2[ind,] # allows boot to select sample 
  logitmodel_missing_bstrap <<- glm(notmissing ~ dpisofirme + ln_s_income + s_hhpeople + s_age
                            ,family=binomial(link='logit'),data=d)
  
  logitpredval_bstrap <<- predict(logitmodel_missing_bstrap, type = "response")
  d$logitpredval <- logitpredval_bstrap
  forgmmboot <- d[complete.cases(d),]
  temp <<- gmm(marmoments, forgmmboot, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef
  #return(temp)
}

bootresult_3c <- boot(data = allrows_subset2, R = reps, statistic = boot.GMM_mar_3c)

#Now find the CI from here
boot_results_3c <- bootresult_3c$t
meltboot_3c <- melt(boot_results_3c)
meltboot_sorted_3c <- meltboot_3c[order(meltboot_3c$Var2,meltboot_3c$value),]
meltboot_sorted_3c$seq <- seq(1: nrow(meltboot_sorted_3c))
meltboot_sorted_3c$seq_adj <- meltboot_sorted_3c$seq - (reps * (meltboot_sorted_3c$Var2 - 1) )

lower_ci = round(reps * .05)
upper_ci = round(reps * .95)

lower_ci_3c <-  meltboot_sorted_3c[meltboot_sorted_3c$seq_adj == lower_ci,]
upper_ci_3c <- meltboot_sorted_3c[meltboot_sorted_3c$seq_adj == upper_ci,]

outtable <- data.frame(varnames,round(inital_mar_3c,4), round(lower_ci_3c$value,4),  round(upper_ci_3c$value,4))
outpath <- paste(export,"/R_table_PS3_Q2_3c.tex", sep = "")
latex(outtable,
      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(2, 2),
      cgroup=c("Estimate", "Inference"),
      colheads=c("Variable", "$\\hat\\beta$", "95\\% C.I.", "95\\% C.I.")
)


####
#Question 3d
####

forgmm3d <- forgmm[forgmm$logitpredval >= .1, ]
inital_mar_3d <- gmm(marmoments, forgmm3d, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef
#inital_mar_3c <- gmm(marmoments, allrows_predval, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef

boot.GMM_mar_3d <- function(x, ind) {
  d <<- allrows_subset2[ind,] # allows boot to select sample 
  logitmodel_missing_bstrap <<- glm(notmissing ~ dpisofirme + ln_s_income + s_hhpeople + s_age
                                    ,family=binomial(link='logit'),data=d)
  
  logitpredval_bstrap <<- predict(logitmodel_missing_bstrap, type = "response")
  d$logitpredval <- logitpredval_bstrap
  forgmmboot <- d[complete.cases(d),]
  forgmmboot3d <- forgmmboot[forgmmboot$logitpredval >= .1, ]
  temp <<- gmm(marmoments, forgmmboot3d, t0=c(0,0,0,0), wmatrix="ident", vcov="iid")$coef
  #return(temp)
}

bootresult_3d <- boot(data = allrows_subset2, R = reps, statistic = boot.GMM_mar_3d)

#Now find the CI from here
boot_results_3d <- bootresult_3d$t
meltboot_3d <- melt(boot_results_3d)
meltboot_sorted_3d <- meltboot_3d[order(meltboot_3d$Var2,meltboot_3d$value),]
meltboot_sorted_3d$seq <- seq(1: nrow(meltboot_sorted_3d))
meltboot_sorted_3d$seq_adj <- meltboot_sorted_3d$seq - (reps * (meltboot_sorted_3d$Var2 - 1) )

lower_ci = round(reps * .05)
upper_ci = round(reps * .95)

lower_ci_3d <-  meltboot_sorted_3d[meltboot_sorted_3d$seq_adj == lower_ci,]
upper_ci_3d <- meltboot_sorted_3d[meltboot_sorted_3d$seq_adj == upper_ci,]

outtable <- data.frame(varnames,round(inital_mar_3d,4), round(lower_ci_3d$value,4),  round(upper_ci_3d$value,4))
outpath <- paste(export,"/R_table_PS3_Q2_3d.tex", sep = "")
latex(outtable,
      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(2, 2),
      cgroup=c("Estimate", "Inference"),
      colheads=c("Variable", "$\\hat\\beta$", "95\\% C.I.", "95\\% C.I.")
)
























