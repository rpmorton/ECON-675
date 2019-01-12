#####
#Created by RM on 2018.12.04
#Econ 675: PS 6, Q2
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

#load in RD packages

#install.packages('rdrobust')
#install.packages('rdlocrand')
#install.packages('rddensity')

library(rdrobust)
library(rdlocrand)
library(rddensity)
library(lattice)

#Begin Pset

rm(list=ls())
export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

#load in data 
headstart <-read.csv("/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/HeadStart.csv")
headstart$treat <- ifelse(headstart$povrate60 > 0, 1, 0)
  
#2.1: Plots and Falsification
#Plots
outpath <- paste(export,"/R_PS6_Q2_1_imse_even.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart$mort_related_pre, headstart$povrate60, c = 0, p = 2, binselect = "es")
dev.off()

outpath <- paste(export,"/R_PS6_Q2_1_imse_quant.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart$mort_related_pre, headstart$povrate60, c = 0, p = 2, binselect = "qspr")
dev.off()

outpath <- paste(export,"/R_PS6_Q2_1_var_even.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart$mort_related_pre, headstart$povrate60, c = 0, p = 2, binselect = "esmv")
dev.off()

outpath <- paste(export,"/R_PS6_Q2_1_var_quant.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart$mort_related_pre, headstart$povrate60, c = 0, p = 2, binselect = "qsmvpr")
dev.off()

#Density Test
density_test_p2 <- rddensity(headstart$povrate60)
density_test_results <- density_test_p2$test
density_test_pval <- density_test_results$p_jk

#Binomial Test
binomtest <-  rdwinselect (headstart$povrate60 , headstart$mort_related_pre, wmin = 1 , nwindows = 1)

#Histogram
zero <- 0
zero_data <- data.frame(zero)
plot_hist <- ggplot(headstart,aes(x=povrate60)) + geom_histogram(data = subset(headstart,povrate60 > 0), breaks = seq(-50,25,by = 2.5), fill = "dark red", alpha = .6)
plot_hist <- plot_hist + geom_histogram(data = subset(headstart, povrate60 < 0), breaks = seq(-50,25,by = 2.5), fill = "dark blue", alpha = .6)
plot_hist <- plot_hist + geom_vline(data = zero_data, aes(xintercept = zero), linetype = "dashed" ) +   theme(legend.position="top")
plot_hist <- plot_hist + labs(title = "PS6: Q2.1, Histogram", y = "Observation Count",
                x = "Normalized Poverty Rate")
plot_hist
outpath <- paste(export,"/R_PS6_Q2_1_Hist.pdf", sep = "")
ggsave(outpath)

#2.2: Parametric Methods
treat_coeff <-rep(-100,2 * 4)
treat_se <-rep(-100,2 * 4)
parametric_results <- data.frame(treat_coeff, treat_se)
parametric_results$type <- "test"
parametric_results$poly_order <- -1

headstart_2_2 <- headstart
headstart_2_2$povrate60_2 <- headstart_2_2$povrate ^2
headstart_2_2$povrate60_3 <- 1
headstart_2_2$povrate60_4 <- 1
headstart_2_2$povrate60_5 <- 1
headstart_2_2$povrate60_6 <- 1

coefs_const <- rep(-100,4)
stderr_const <- rep(-100,4)
coefs_het <- rep(-100,4)
stderr_het <- rep(-100,4)
results_2_2 <- data.frame(coefs_const,stderr_const,coefs_het,stderr_het)
results_2_2$poly <- NaN

for(p in 3:6) {
  
  headstart_2_2[,4 + p] <- headstart_2_2[,3 + p] * headstart_2_2$povrate60
  
  #Constant Treatment
  constant_lm <- lm(mort_related_post ~  treat + povrate60 + povrate60_2 + povrate60_3 + povrate60_4 + povrate60_5 + povrate60_6,data=headstart_2_2)
  constant_ols <- coef(summary(constant_lm))
  constant_ols_data <- as.data.frame(constant_ols)
  constant_ols_data$vars <- rownames(constant_ols)
  constant_treat <- subset(constant_ols_data,vars=="treat")
  p_less_two <- p - 2
  results_2_2[p_less_two,1] <- constant_treat$Estimate
  results_2_2[p_less_two,2] <- constant_treat$`Std. Error`
  
  #Het Treatment
  het_lm <- lm(mort_related_post ~  treat + (povrate60*treat) + (povrate60_2*treat)
               + (povrate60_3*treat) + (povrate60_4*treat) + (povrate60_5*treat) + (povrate60_6*treat)
               + povrate60 + povrate60_2 + povrate60_3 + povrate60_4 + povrate60_5 + povrate60_6,data=headstart_2_2)
  het_ols <- coef(summary(het_lm))
  het_ols_data <- as.data.frame(het_ols)
  het_ols_data$vars <- rownames(het_ols)
  het_treat <- subset(het_ols_data,vars=="treat")
  results_2_2[p_less_two,3] <- het_treat$Estimate
  results_2_2[p_less_two,4] <- het_treat$`Std. Error` 

  results_2_2[p_less_two,5] <- p
} 

results_2_2 <- round(results_2_2, digits = 4)

outpath <- paste(export,"/R_PS6_Q2_2ConstHetTreat.csv", sep = "")
write.csv(results_2_2,outpath)

##Plot the Fitted Values 
#Constant
graph_constant <- stack(data.frame(Observed = headstart$mort_related_post, Predicted = fitted(constant_lm)))
graph_constant <- cbind(graph_constant, x = rep(headstart$povrate60, 2) )
graph_constant_subset <- subset(graph_constant, values < 10)
plot_const <- ggplot(graph_constant_subset,aes(x=x, y = values, color=ind)) + geom_point()
plot_const <- plot_const + labs(title = "PS6: Q2.2, Constant Treatment", y = "Mort Related Post",
                                    x = "Normalized Poverty Rate", caption = "Note: Observations with Mort Related Post over 10 not shown.")
plot_const <- plot_const + scale_color_manual(values = c("dark red", "dark blue")
                                          , labels = c("Observed", "6th Order Poly Fit")
                                          , name = NULL)
plot_const
outpath <- paste(export,"/R_PS6_Q2_2_Const.pdf", sep = "")
ggsave(outpath)

#Het
graph_het <- stack(data.frame(Observed = headstart$mort_related_post, Predicted = fitted(het_lm)))
graph_het <- cbind(graph_het, x = rep(headstart$povrate60, 2) )
graph_het_subset <- subset(graph_het, values < 10)
plot_het <- ggplot(graph_het_subset,aes(x=x, y = values, color=ind)) + geom_point()
plot_het <- plot_het + labs(title = "PS6: Q2.2, Heterogeneous Treatment", y = "Mort Related Post",
                                x = "Normalized Poverty Rate", caption = "Note: Observations with Mort Related Post over 10 not shown.")
plot_het <- plot_het + scale_color_manual(values = c("dark red", "dark blue")
                                              , labels = c("Observed", "6th Order Poly Fit")
                                              , name = NULL)
plot_het
outpath <- paste(export,"/R_PS6_Q2_2_Het.pdf", sep = "")
ggsave(outpath)

#Bandwidth Local
bandwidths <- c(1,5,9,18) 
num_bandwidths <- length(bandwidths)

bandwidth <- rep(-100,num_bandwidths * 3)
order <- rep(-100,num_bandwidths * 3)
coef <- rep(-100,num_bandwidths * 3)
stderr <- rep(-100,num_bandwidths * 3)

results_local_param <- data.frame(bandwidth,order,coef,stderr)

for(b in 1:num_bandwidths) {
  h <- bandwidths[b]
  
  row <- (b-1)*3 + 1
  results_local_param[row,1] <- h
  
  reg_subset <- subset(headstart_2_2, abs(povrate60) <= h)
  
  local_poly_0_result <-  coef(summary(lm(mort_related_post ~  treat, data = reg_subset)))
  local_poly_0_data <- as.data.frame(local_poly_0_result)
  local_poly_0_data$vars <- rownames(local_poly_0_result)
  local_poly_0_treat <- subset(local_poly_0_data,vars=="treat")
  results_local_param[row,2] <- 0
  results_local_param[row,3] <- local_poly_0_treat$Estimate
  results_local_param[row,4] <- local_poly_0_treat$`Std. Error`
  
  rplusone <- row + 1
  results_local_param[rplusone,1] <- h
  local_poly_1_result <-  coef(summary(lm(mort_related_post ~  treat + povrate60, data = reg_subset)))
  local_poly_1_data <- as.data.frame(local_poly_1_result)
  local_poly_1_data$vars <- rownames(local_poly_1_result)
  local_poly_1_treat <- subset(local_poly_1_data,vars=="treat")
  results_local_param[rplusone,2] <- 1
  results_local_param[rplusone,3] <- local_poly_1_treat$Estimate
  results_local_param[rplusone,4] <- local_poly_1_treat$`Std. Error`
  
  rplustwo <- row + 2
  results_local_param[rplustwo,1] <- h
  local_poly_2_result <-  coef(summary(lm(mort_related_post ~  treat + povrate60 + povrate60_2, data = reg_subset)))
  local_poly_2_data <- as.data.frame(local_poly_2_result)
  local_poly_2_data$vars <- rownames(local_poly_2_result)
  local_poly_2_treat <- subset(local_poly_2_data,vars=="treat")
  results_local_param[rplustwo,2] <- 2
  results_local_param[rplustwo,3] <- local_poly_2_treat$Estimate
  results_local_param[rplustwo,4] <- local_poly_2_treat$`Std. Error`
}

results_local_param <- round(results_local_param, digits = 4)

#Graph: h = 18
local_poly_2_result <-  lm(mort_related_post ~  treat + povrate60 + povrate60_2, data = reg_subset)
graph_local <- stack(data.frame(Observed = reg_subset$mort_related_post, Predicted = fitted(local_poly_2_result)))
graph_local <- cbind(graph_local, x = rep(reg_subset$povrate60, 2) )
graph_local_subset <- subset(graph_local, values < 10)
plot_local <- ggplot(graph_local_subset,aes(x=x, y = values, color=ind)) + geom_point()
plot_local <- plot_local + labs(title = "PS6: Q2.2, Local Parametric Model", y = "Mort Related Post",
                            x = "Normalized Poverty Rate", caption = "Note: A bandwidth of 18 is used.  Observations with Mort Related Post over 10 not shown.")
plot_local <- plot_local + scale_color_manual(values = c("dark red", "dark blue")
                                          , labels = c("Observed", "2nd Order Poly Fit")
                                          , name = NULL)
plot_local

outpath <- paste(export,"/R_PS6_Q2_2_Local.pdf", sep = "")
ggsave(outpath)

##2.3: Robust Local Polynomial Methods

estimate <- rep(-100,9)
stderr <- rep(-100,9)
ci_l <- rep(-100,9)
ci_h <- rep(-100,9)
results_2_3 <- data.frame(estimate,stderr,ci_l,ci_h)
results_2_3$type <- "test"
results_2_3$poly_order <- -100
p <- 1
for (p in 0:2) {
  row <- p * 3 + 1
  
  local_poly_0 <- rdrobust(headstart$mort_related_post, headstart$povrate60, c = 0, p = p, all = TRUE)
  local_poly_res <- as.data.frame(local_poly_0$Estimate)
  local_poly_ci <- as.data.frame(local_poly_0$ci)
  
  results_2_3[row,5] <- "Conventional"
  results_2_3[row,1] <- local_poly_res$tau.us
  results_2_3[row,2] <- local_poly_res$se.us
  results_2_3[row,3] <- local_poly_ci[1,1]
  results_2_3[row,4] <- local_poly_ci[1,2]
  results_2_3[row,6] <- p
  
  row <- row + 1
  results_2_3[row,5] <- "Bias-Corrected"
  results_2_3[row,1] <- local_poly_res$tau.bc
  results_2_3[row,2] <- local_poly_res$se.us
  results_2_3[row,3] <- local_poly_ci[2,1]
  results_2_3[row,4] <- local_poly_ci[2,2]
  results_2_3[row,6] <- p
  
  row <- row + 1
  results_2_3[row,5] <- "Robust"
  results_2_3[row,1] <- local_poly_res$tau.bc
  results_2_3[row,2] <- local_poly_res$se.rb
  results_2_3[row,3] <- local_poly_ci[3,1]
  results_2_3[row,4] <- local_poly_ci[3,2]
  results_2_3[row,6] <- p

}

for(i in 1:4) {
  results_2_3[,i] <- round(results_2_3[,i],digits = 4)
}

outpath <- paste(export,"/R_PS6_Q2_3.csv", sep = "")
write.csv(results_2_3,outpath)

##2.3: Robustness Checks

#Bandwidths

bandwidths_rob <- seq(1:10)
size_band_rob <- length(bandwidths_rob)

bandwidth <- rep(-100,size_band_rob*3)
bias_corrected_est <- rep(-100,size_band_rob*3)

results_2_3_rob_band <- data.frame(bandwidth,bias_corrected_est)
results_2_3_rob_band$kernel <- "test"

for (i in 1:size_band_rob) {
  row <- (i-1)* 3 + 1
  
  h <- bandwidths_rob[i]
  
  local_poly_rob_band_tri <- rdrobust(headstart$mort_related_post, headstart$povrate60, c = 0, p = 1, h = h, all = TRUE, kernel = "triangular")
  local_poly_rob_band_tri_res <- as.data.frame(local_poly_rob_band_tri$Estimate)
  results_2_3_rob_band[row,1] <- h
  results_2_3_rob_band[row,2] <- local_poly_rob_band_tri_res$tau.bc
  results_2_3_rob_band[row,3] <- "tri"
  
  row <- row + 1
  local_poly_rob_band_uni <- rdrobust(headstart$mort_related_post, headstart$povrate60, c = 0, p = 1, h = h, all = TRUE, kernel = "uniform")
  local_poly_rob_band_uni_res <- as.data.frame(local_poly_rob_band_uni$Estimate)
  results_2_3_rob_band[row,1] <- h
  results_2_3_rob_band[row,2] <- local_poly_rob_band_uni_res$tau.bc
  results_2_3_rob_band[row,3] <- "uni"
 
  row <- row + 1
  local_poly_rob_band_epa <- rdrobust(headstart$mort_related_post, headstart$povrate60, c = 0, p = 1, h = h, all = TRUE, kernel = "epanechnikov")
  local_poly_rob_band_epa_res <- as.data.frame(local_poly_rob_band_epa$Estimate)
  results_2_3_rob_band[row,1] <- h
  results_2_3_rob_band[row,2] <- local_poly_rob_band_epa_res$tau.bc
  results_2_3_rob_band[row,3] <- "epa"
}

results_2_3_rob_band$bias_corrected_est <- round(results_2_3_rob_band$bias_corrected_est, digits = 4)
results_2_3_rob_band_wide <- reshape(results_2_3_rob_band, timevar = "bandwidth", idvar = "kernel", direction = "wide" )

outpath <- paste(export,"/R_PS6_Q2_3_rob_band_kern.csv", sep = "")
write.csv(results_2_3_rob_band_wide,outpath)

#Donut Hole
headstart_sorted <- headstart[order(headstart$povrate60),] 
nrow_headstart <- nrow(headstart)
headstart_sorted$seq <- seq(1:nrow_headstart)
max_lt_cutoff <- max(subset(subset(headstart_sorted,povrate60 < 0),select = seq ))
headstart_sorted$donut <- ifelse(headstart_sorted$povrate60 < 0,
                                 -1 * (headstart_sorted$seq  - max_lt_cutoff ) + 1,
                                 headstart_sorted$seq - max_lt_cutoff)

results_donut <-rep(-100,10) 
donut_size <- rep(-100,10)
results_2_3_rob_donut <- data.frame(results_donut,donut_size)

for (i in 1:10) {
  headstart_subset <- subset(headstart_sorted,donut > i)
  local_poly_rob_donut <- rdrobust(headstart_subset$mort_related_post, headstart_subset$povrate60, c = 0, p = 1, all = TRUE, kernel = "triangular")
  local_poly_rob_donut_res <- as.data.frame(local_poly_rob_donut$Estimate)
  results_2_3_rob_donut[i,1] <- local_poly_rob_donut_res$tau.bc
  results_2_3_rob_donut[i,2] <- i

}

results_2_3_rob_donut$results_donut <- round(results_2_3_rob_donut$results_donut, digits = 4)

outpath <- paste(export,"/R_PS6_Q2_3_rob_donut.csv", sep = "")
write.csv(results_2_3_rob_donut,outpath)

##Cutoff

cutoffs <- seq(-10, 10, 2 )
cutoffs_len <- length(cutoffs)

est_cutoff <- rep(-100,cutoffs_len)
pval_cutoff <- rep(-100,cutoffs_len)
cutoff_point <- rep(-100,cutoffs_len)
results_2_3_rob_cutoff <- data.frame(est_cutoff,stderr_cutoff,cutoff_point)

for (i in 1:cutoffs_len ) {
  c <- cutoffs[i]
  local_poly_rob_cutoff <- rdrobust(headstart$mort_related_post, headstart$povrate60, c = c, p = 1, all = TRUE, kernel = "triangular")
  local_poly_rob_cutoff_res <- as.data.frame(local_poly_rob_cutoff$Estimate)
  local_poly_rob_cutoff_res_pv <- as.data.frame(local_poly_rob_cutoff$pv)
  local_poly_rob_cutoff_res_pv$type <- rownames(local_poly_rob_cutoff_res_pv)
  local_poly_rob_cutoff_res_pv_bc <- subset(local_poly_rob_cutoff_res_pv, type == "Bias-Corrected")
  results_2_3_rob_cutoff[i,1] <- local_poly_rob_cutoff_res$tau.bc
  results_2_3_rob_cutoff[i,2] <- local_poly_rob_cutoff_res_pv_bc$`P>|z|`
  results_2_3_rob_cutoff[i,3] <- c
}

results_2_3_rob_cutoff <- round(results_2_3_rob_cutoff, digits = 4)

outpath <- paste(export,"/R_PS6_Q2_3_rob_cutoff.csv", sep = "")
write.csv(results_2_3_rob_cutoff,outpath)

###2.4: Local Randomization Methods

win_local <- rdwinselect(headstart$povrate60,headstart$mort_injury_post)
win_local_act <- as.vector(win_local$window)
headstart_local_rand <- subset(headstart, abs(povrate60) < abs(win_local_act[1]) )
test_win <- lm(headstart_local_rand$mort_related_pre ~ headstart_local_rand$treat)
summary(test_win)
  
outpath <- paste(export,"/R_PS6_Q2_4_injury.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart_local_rand$mort_injury_post, headstart_local_rand$povrate60, c = 0, p = 0, binselect = "es")
dev.off()

outpath <- paste(export,"/R_PS6_Q2_4_premort.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart_local_rand$mort_related_pre, headstart_local_rand$povrate60, c = 0, p = 0, binselect = "es")
dev.off()

outpath <- paste(export,"/R_PS6_Q2_4_postmort.pdf", sep = "")
pdf(outpath,width=7,height=5)
rdplot(headstart_local_rand$mort_related_post, headstart_local_rand$povrate60, c = 0, p = 0, binselect = "es")
dev.off()

#Bootstrap Approach

permute.diff_mean <- function(data, ind) {
  count_treat <- sum(data$povrate60 > 0)
  permsample  <- data[ind,] # allows boot to select sample 
  total_rows <- nrow(permsample)
  permsample$seq <- seq(1:total_rows)
  permsample$treat <- ifelse( permsample$seq <= count_treat, 1, 0)
  means_by_group <<- ddply(permsample, .(treat), summarise, mean_outcome = mean(mort_related_post) )
  diff_means <- as.numeric(subset(subset(means_by_group, treat == 1), select = mean_outcome) - subset(subset(means_by_group, treat == 0), select = mean_outcome))
  ret <- as.numeric(diff_means)
  return(ret)
  
}

reps <- 1000
perm_diff_mean <- boot(data = headstart_local_rand, R = reps, statistic = permute.diff_mean, sim = "permutation", stype = "i")
perm_diff_mean_res <- as.data.frame(perm_diff_mean$t)
perm_diff_mean_res <- perm_diff_mean_res[order(perm_diff_mean_res$V1),] 
low_ci_obs <- round(.025*reps)
low_ci <- perm_diff_mean_res[low_ci_obs]
high_ci_obs <- round(.975*reps)
high_ci <- perm_diff_mean_res[high_ci_obs]

actual_means_by_group <- ddply(headstart_local_rand, .(treat), summarise, mean_mort = mean(mort_related_post) )
actual_diff_means <-as.numeric(subset(subset(actual_means_by_group, treat == 1), select = mean_mort) - subset(subset(actual_means_by_group, treat == 0), select = mean_mort) )

results_2_4_fisher <- data.frame(actual_diff_means,low_ci,high_ci)

##Sensitivity: Change Windows:

window <- seq(.8,2.6,.2)

sensitivty_local_rand <- function(window_sense) {
  mat_ret <- matrix(NaN,nrow = 4, ncol = 1)
  
  subset_headstart <- subset(headstart,abs(povrate60) < window_sense)
  subset_headstart$treat <- ifelse( subset_headstart$povrate60 >0 , 1, 0)
  actual_means_by_group <- ddply(subset_headstart, .(treat), summarise, mean_mort = mean(mort_related_post) )
  actual_diff_means <- as.numeric(subset(subset(actual_means_by_group, treat == 1), select = mean_mort) - subset(subset(actual_means_by_group, treat == 0), select = mean_mort) )
  mat_ret[1,1] <- actual_diff_means
  
  perm_diff_mean <- boot(data = subset_headstart, R = reps, statistic = permute.diff_mean, sim = "permutation", stype = "i")
  perm_diff_mean_res <- as.data.frame(perm_diff_mean$t)
  perm_diff_mean_res <- data.frame(perm_diff_mean_res)
  std_error <- sd(perm_diff_mean_res$V1)
  mat_ret[2,1] <- std_error
  
  perm_diff_mean_res$act <- rep(actual_diff_means,nrow(perm_diff_mean_res))
  num_over <- sum(abs(perm_diff_mean_res$V1) > abs(perm_diff_mean_res$act) )
  p_val <- num_over / reps
  mat_ret[3,1] <- p_val
  
  mat_ret[4,1] <- window_sense
  
  return(mat_ret)
}

ani_way <- sapply(1:length(window),function(i) sensitivty_local_rand(window[i]))

outpath <- paste(export,"/R_PS6_Q2_4_sensitivity.csv", sep = "")
write.csv(ani_way,outpath)


 
