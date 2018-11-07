***Created by RM on 2018.11.04
***For ECON 675, PS 4, Q 2
*********************

set more off

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"

global alpha = .95

/* First, ATE using LaLonde Data */

import delim using "$data/LaLonde_1986.csv", delim(",") clear

*Diff in Means

g treat_earn = earn78 if treat == 1
g ctrl_earn = earn78 if treat == 0

su treat_earn, d
g treat_avg_earn =  `r(mean)'
g treat_sd_earn =  `r(sd)'
g treat_count_earn = `r(N)'

su ctrl_earn, d
g ctrl_avg_earn =  `r(mean)'
g ctrl_sd_earn =  `r(sd)'
g ctrl_count_earn = `r(N)'

g se_diff_means = sqrt(  ( (treat_sd_earn * treat_sd_earn ) / treat_count_earn ) + ( (ctrl_sd_earn * ctrl_sd_earn ) / ctrl_count_earn ) )

g diff_means_earn = treat_avg_earn - ctrl_avg_earn
g diff_means_ci_low = diff_means_earn +  invttail(treat_count_earn - 1 + ctrl_count_earn - 1, 1 - ($alpha / 2) ) * se_diff_means
g diff_means_ci_high = diff_means_earn +  invttail(treat_count_earn - 1 + ctrl_count_earn - 1, $alpha / 2 ) * se_diff_means

*OLS

g ln_earn74 = ln(earn74 + 1)
g ln_earn75 = ln(earn75 + 1)
global speca = "age educ black hisp married nodeg ln_earn74 ln_earn75"

g age2 = age^2
g educ2 = educ^2
global specb = "age2 educ2 u74 u75"

g age3 = age^3
g black_u74 = black * u74
g educ_ln_earn74 = educ * ln_earn74
global specc = "age3 black_u74 educ_ln_earn74"

reg earn78 treat $speca , r
reg earn78 treat $speca $specb , r
reg earn78 treat $speca $specb $specc , r

*Reg Impute

teffects ra (earn78 $speca) (treat), ate
teffects ra (earn78 $speca $specb) (treat), ate
teffects ra (earn78 $speca $specb $specc) (treat), ate

* IPW

teffects ipw (earn78) (treat $speca , probit) , ate
teffects ipw (earn78) (treat $speca $specb , probit) , ate
teffects ipw (earn78) (treat $speca $specb $specc , probit) , ate 

*DR

teffects ipwra (earn78 $speca) (treat $speca , probit) , ate
*DOES NOT CONVERGE, SO MODEL DIFFERENTLY
probit treat $speca $specb 
predict predval
teffects ipwra (earn78 $speca $specb) (treat $speca $specb , probit) , ate
teffects ipwra (earn78 $speca $specb $specc) (treat $speca $specb $specc , probit) , ate 


*N1 Match

teffects nnmatch (earn78 $speca) (treat) , ate nneighbor(1) metric(maha)
teffects nnmatch (earn78 $speca $specb) (treat) , ate nneighbor(1) metric(maha)
teffects nnmatch (earn78 $speca $specb $specc) (treat) , ate nneighbor(1) metric(maha)

*P Match


teffects psmatch (earn78) (treat $speca ) , ate
teffects psmatch (earn78) (treat $speca $specb ) , ate 
teffects psmatch (earn78) (treat $speca $specb $specc ) , ate 


/* Second, ATT using LaLonde Data */

import delim using "$data/LaLonde_1986.csv", delim(",") clear


*Reg Impute

teffects ra (earn78 $speca) (treat), atet
teffects ra (earn78 $speca $specb) (treat), atet
teffects ra (earn78 $speca $specb $specc) (treat), atet

*IPW

teffects ipw (earn78) (treat $speca , probit) , atet
teffects ipw (earn78) (treat $speca $specb , probit) , atet
teffects ipw (earn78) (treat $speca $specb $specc , probit) , atet

*DR

teffects ipwra (earn78 $speca ) (treat $speca , probit) , atet
teffects ipwra (earn78 $speca $specb ) (treat $speca $specb , probit) , atet
teffects ipwra (earn78 $speca $specb $specc ) (treat $speca $specb $specc , probit) , atet

*N1 Match

teffects nnmatch (earn78 $speca ) (treat) , atet nneighbor(1) metric(maha)
teffects nnmatch (earn78 $speca $specb ) (treat) , atet nneighbor(1) metric(maha)
teffects nnmatch (earn78 $speca $specb $specc ) (treat) , atet nneighbor(1) metric(maha)

*P Match

teffects psmatch (earn78) (treat $speca ) , atet
teffects psmatch (earn78) (treat $speca $specb ) , atet 
teffects psmatch (earn78) (treat $speca $specb $specc ) , atet 

/* Third, ATE using PSID Data */

import delim using "$data/LaLonde_all.csv", delim(",") clear

*Prep Data

drop if treat == 0
g treat_alt = treat == 1
rename re78 earn78
rename re75 earn75
rename re74 earn74

* Diff Means

g treat_earn = earn78 if treat_alt == 1
g ctrl_earn = earn78 if treat_alt == 0

su treat_earn, d
g treat_avg_earn =  `r(mean)'
g treat_sd_earn =  `r(sd)'
g treat_count_earn = `r(N)'

su ctrl_earn, d
g ctrl_avg_earn =  `r(mean)'
g ctrl_sd_earn =  `r(sd)'
g ctrl_count_earn = `r(N)'

g se_diff_means = sqrt(  ( (treat_sd_earn * treat_sd_earn ) / treat_count_earn ) + ( (ctrl_sd_earn * ctrl_sd_earn ) / ctrl_count_earn ) )

g diff_means_earn = treat_avg_earn - ctrl_avg_earn
g diff_means_ci_low = diff_means_earn +  invttail(treat_count_earn - 1 + ctrl_count_earn - 1, 1 - ($alpha / 2) ) * se_diff_means
g diff_means_ci_high = diff_means_earn +  invttail(treat_count_earn - 1 + ctrl_count_earn - 1, $alpha / 2 ) * se_diff_means

*OLS

g ln_earn74 = ln(earn74 + 1)
g ln_earn75 = ln(earn75 + 1)
global speca = "age educ black hisp married nodeg ln_earn74 ln_earn75"

g age2 = age^2
g educ2 = educ^2
global specb = "age2 educ2 u74 u75"

g age3 = age^3
g black_u74 = black * u74
g educ_ln_earn74 = educ * ln_earn74
global specc = "age3 black_u74 educ_ln_earn74"

reg earn78 treat_alt $speca , r
reg earn78 treat_alt $speca $specb , r
reg earn78 treat_alt $speca $specb $specc , r

*Reg Impute

teffects ra (earn78 $speca) (treat_alt), ate
teffects ra (earn78 $speca $specb) (treat_alt), ate
teffects ra (earn78 $speca $specb $specc) (treat_alt), ate

*IPW

*teffects ipw (earn78) (treat_alt $speca , probit) , ate

probit treat_alt $speca
preserve
predict predval
drop if predval < .01 | predval > .99
teffects ipw (earn78) (treat_alt $speca , probit) , ate
restore

*teffects ipw (earn78) (treat_alt $speca $specb , probit) , ate

probit treat_alt $speca $specb
preserve
predict predval
drop if predval < .01 | predval > .99
teffects ipw (earn78) (treat_alt $speca $specb , probit) , ate
restore

teffects ipw (earn78) (treat_alt $speca $specb $specc , probit) , ate



