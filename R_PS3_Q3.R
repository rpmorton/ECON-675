#####
#Created by RM on 2018.10.28
#Econ 675: PS 3, Q 3
#####

library(MASS)
library(plyr)
library(data.table)
library(matlib)
library(randomizr)
library(tidyverse)
library(boot)
library(ggplot2)

export <- "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

####Part 1

set.seed(123)

obs <- 1000
reps <- 599

rand_uni <- runif(obs)
counter <- seq(1:obs)

rand_data <- data.frame(rand_uni,counter)

boot.T <- function(x, ind) {
  d <-rand_data[ind,] 
  temp <- max(d$rand_uni) 
  return(temp)
}

bootresult <- boot(data = rand_data, R = reps, statistic = boot.T)
bstrap_max <- bootresult$t

max_rand <- max(rand_uni)
max_rand_rep <- rep(max_rand,reps)

bstrap_stat <- obs * (max_rand_rep - bstrap_max)

rand_exp <- rexp(reps, rate = 1)

density_data <- data.frame(bstrap_stat, rand_exp)
meltdensitynonpar <- melt(density_data)

plot_3_1 <- ggplot(meltdensitynonpar,aes(x=value, color=variable)) + geom_density() 
plot_3_1 <- plot_3_1 + labs(title = "PS3: Q3 1", y = "Density", x = "x")
plot_3_1 <- plot_3_1 + scale_color_manual(values = c("dark red", "dark blue")
          , breaks = c("bstrap_stat", "rand_exp")
          , labels = c("Nonparametric Bootstrap Statistic", "Simulated Exponential Distribution") 
          , name = NULL)
plot_3_1


outpath <- paste(export,"/R_density_PS3_Q3_1.pdf", sep = "")
ggsave(outpath)

####Part 2

set.seed(123)

obs <- 1000
reps <- 599

rand_uni <- runif(obs)
max_rand <- max(rand_uni)

rand_data <- data.frame(rep(max_rand,reps),rep(-10,reps))


for (j in 1:reps) {
      bstrap_rand <- runif(obs)*max_rand
      max_bstrap <- max(bstrap_rand)
      rand_data[j,2] <- max_bstrap
}

para_bstrap_stat <- obs * (rand_data[,1] - rand_data[,2])
rand_exp <- rexp(reps, rate = 1)


density_data <- data.frame(para_bstrap_stat, rand_exp)
meltdensitypar <- melt(density_data)

plot_3_2 <- ggplot(meltdensitypar,aes(x=value, color=variable)) + geom_density() 
plot_3_2 <- plot_3_2 + labs(title = "PS3: Q3 2", y = "Density", x = "x")
plot_3_2 <- plot_3_2 + scale_color_manual(values = c("dark red", "dark blue")
            , breaks = c("para_bstrap_stat", "rand_exp")
            , labels = c("Parametric Bootstrap Statistic", "Simulated Exponential Distribution")
            , name = NULL)
plot_3_2

outpath <- paste(export,"/R_density_PS3_Q3_2.pdf", sep = "")
ggsave(outpath)