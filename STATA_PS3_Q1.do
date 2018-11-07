***Created by RM on 2018.10.19
***For ECON 675, PS 3, Q 1
*********************
 
set more off

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"
global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"


/************
***Question 1: 9
***********/

import delim using "$data/pisofirme", delim(",") clear

g notmissing = 1 - dmissing

g ln_s_incomepc = ln(s_incomepc + 1)

*Parametric Version
logit notmissing s_age s_hhpeople ln_s_incomepc, vce(robust)

*Bootstrap
logit notmissing s_age s_hhpeople ln_s_incomepc, vce(bootstrap, reps(100)) 

*Predicted value
predict pred_val

twoway (kdensity pred_val, xtitle("x") ytitle("Estimated Density") title("Density of Propensity Scores" ) ) 
graph export "$out/STATA_PS3_Q1_9c.pdf", as (pdf) replace

