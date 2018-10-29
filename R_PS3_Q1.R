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

#q9: a

reps <- 100

pisofirme <- read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/pisofirme.csv")

pisofirme$notmissing <-  1 - pisofirme$dmissing
pisofirme$ln_income_pc <- log(pisofirme$S_incomepc + 1)
pisofirmecols <- select(pisofirme, S_age, S_HHpeople, ln_income_pc,  notmissing)
pisofirmecomplete <- pisofirmecols[complete.cases(pisofirmecols),]

logitmodel <- glm(notmissing ~  S_age + S_HHpeople + ln_income_pc 
                                  , family=binomial(link='logit'),data=pisofirmecomplete)
initial_coefs <- as.data.frame(summary(logitmodel)$coefficients)
alpha <- .95
initial_coefs$lb <- initial_coefs$Estimate + qt((1-alpha)/2,logitmodel$df.residual)*(initial_coefs$"Std. Error")
initial_coefs$ub <- initial_coefs$Estimate + qt(alpha+(1-alpha)/2,logitmodel$df.residual)*(initial_coefs$"Std. Error")
initial_coefs <- round(initial_coefs,4)
initial_coefs$varnames <- c("Constant","S_age","S_HHpeople","log(S_incomepc+1)")


outtable <- select(initial_coefs, varnames, Estimate, "Std. Error", "z value", "Pr(>|z|)", lb, ub )

outpath <- paste(export,"/R_table_PS3_Q1_9a.tex", sep = "")

latex(outtable,
      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(4, 3),
      cgroup=c("Estimate", "Inference"),
      colheads=c("Variable", "$\\hat\\beta$", "Std. Error", "z-value", "p-value", "95\\% C.I.", "95\\% C.I.")
)



#q9: b

boot.logit <- function(x, ind) {
  d <-pisofirmecomplete[ind,] 
  logitboot <- glm(notmissing ~ S_age + S_HHpeople + ln_income_pc 
                    ,family=binomial(link='logit'),data=d)
  logitbootcoefs <- summary(logitboot)$coefficients
  temp <- logitbootcoefs[,1]
  return(temp)
}  

bootresult <- boot(data = pisofirmecomplete, R = reps, statistic = boot.logit)
boot_results <- bootresult$t

meltboot <- melt(boot_results)
meltboot_sorted <- meltboot[order(meltboot$Var2,meltboot$value),]

initial_coefs <- data.frame(initial_coefs[,1], seq(1:4))
boot_initial_merged <- merge(meltboot, initial_coefs, by.x = "Var2", by.y = "seq.1.4.")

bstrap_se <-  ddply(boot_initial_merged, .(Var2), mutate
                    , mean_value = mean( value) ) 
bstrap_se <-  ddply(bstrap_se, .(Var2), mutate
                       , sum_squared_errors = sum( (value-mean_value)^2 ) ) 
bstrap_deduped <- distinct(select(bstrap_se,initial_coefs...1., sum_squared_errors
                                      , Var2) )
bstrap_deduped$se <- sqrt(bstrap_deduped$sum_squared_errors * (1/(reps-1)))

alpha <- .95
df <- reps - 1
bstrap_deduped$lb <- bstrap_deduped$initial_coefs...1. + qt((1-alpha)/2,df)*(bstrap_deduped$se)
bstrap_deduped$ub <- bstrap_deduped$initial_coefs...1. + qt(alpha+(1-alpha)/2,df)*(bstrap_deduped$se)
bstrap_deduped$t <- bstrap_deduped$initial_coefs...1. / bstrap_deduped$se
bstrap_deduped$p_val <- pt(-1 *   abs(bstrap_deduped$t), df)
            + 1 - pt( abs(bstrap_deduped$t) / bstrap_deduped$se, df)
bstrap_deduped <- round(bstrap_deduped,4)
bstrap_deduped$varnames <- c("S_age","S_HHpeople","log(S_incomepc+1)","Constant")
outtable <- select(bstrap_deduped,varnames,initial_coefs...1., p_val
                                  , lb, ub) 


outpath <- paste(export,"/R_table_PS3_Q1_9b.tex", sep = "")

latex(outtable,
      file=paste0(outpath),append=FALSE,table.env=FALSE,center="none",title="",
      n.cgroup=c(2, 3),
      cgroup=c("Estimate", "Inference"),
      colheads=c("Variable", "$\\hat\\beta$", "p-value Using Bootstrapped SE", "95\\% C.I.", "95\\% C.I.")
)


###
#3c
###

logitmodel <- glm(notmissing ~  S_age + S_HHpeople + ln_income_pc 
                  ,family=binomial(link='logit'),data=pisofirmecomplete)

summary(logitmodel)
logitpredval <- predict(logitmodel, type = "response")

logitpredval_data <- data.frame(logitpredval)
plot_1_9c <- ggplot(logitpredval_data,aes(x=logitpredval)) + geom_density() 
plot_1_9c <- plot_1_9c + labs(title = "PS3: Q1 9c", y = "Density", x = "x")
plot_1_9c

outpath <- paste(export,"/R_density_PS3_Q1_9c.pdf", sep = "")
ggsave(outpath)


