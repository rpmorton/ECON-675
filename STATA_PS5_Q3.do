***Created by RM on 2018.11.25
***For ECON 665, PS 5, Q3
*********************

clear
set more off

global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"
global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"
global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"

/* 3.1: Angrist Krueger */

global controls "non_white married SMSA i.region i.YoB_ld"

use "$data/Angrist_Krueger", clear

**OLS 1

reg l_w_wage educ $controls
regsave educ using "$temp/PS5_Q3_1_results.dta", tstat addlabel(model,"OLS_1") replace

**OLS 2

reg l_w_wage educ age_q age_sq $controls
regsave educ age_q age_sq using "$temp/PS5_Q3_1_results.dta", tstat addlabel(model,"OLS_2") append

**2SLS 1

ivregress 2sls l_w_wage $controls (educ = i.QoB##i.YoB_ld $controls )
regsave educ using "$temp/PS5_Q3_1_results.dta", tstat addlabel(model,"2SLS_1") append

**2SLS 2

ivregress 2sls l_w_wage age_q age_sq $controls (educ = i.QoB##i.YoB_ld $controls age_q age_sq )
regsave educ age_q age_sq using "$temp/PS5_Q3_1_results.dta", tstat addlabel(model,"2SLS_2") append

use "$temp/PS5_Q3_1_results.dta", clear 

export excel using "$out/STATA_PS5_Q3_1.xlsx", replace first(variables)


/* 3.2 : Bound Et Al. */

capture program drop IV_quick_1
program define IV_quick_1, rclass
    *syntax varlist(max=1)
    *local x "‘varlist’"
        capture drop educ_hat
        qui reg educ non_white married SMSA i.region i.YoB_ld i.YoB_ld##i.QoB
        predict educ_hat
        qui reg l_w_wage educ_hat non_white married SMSA i.region i.YoB_ld
        return scalar beta = _b[educ_hat]
end

capture program drop IV_quick_2
program define IV_quick_2, rclass
    *if (‘model’ == 2) 
        capture drop educ_hat
        qui reg educ non_white married SMSA age_q age_sq i.region i.YoB_ld i.YoB_ld##i.QoB
        predict educ_hat
        qui reg l_w_wage educ_hat non_white married SMSA age_q age_sq i.region i.YoB_ld
        return scalar beta = _b[educ_hat]
 end


use "$data/Angrist_Krueger", clear

* Now Permute 

IV_quick_1
g actual_2sls_1 = r(beta)
local actual_2sls_1 = actual_2sls_1

capture erase "$temp/PS5_Q3_2_2SLS_1_results.dta"
permute QoB beta_2sls_1 = r(beta), reps(500) sa("$temp/PS5_Q3_2_2SLS_1_results.dta") : IV_quick_1


IV_quick_2
g actual_2sls_2 = r(beta)
local actual_2sls_2 = actual_2sls_2

capture erase "$temp/PS5_Q3_2_2SLS_2_results.dta"
permute QoB beta_2sls_2 = r(beta), reps(500) sa("$temp/PS5_Q3_2_2SLS_2_results.dta") : IV_quick_2

* Summarize Data

forv i = 1(1)2 {
	use "$temp/PS5_Q3_2_2SLS_`i'_results.dta", clear
	su beta_2sls_`i', d
	g mean_beta_2sls_`i' = round(`r(mean)',.0001)
	local mean_beta_2sls_`i' = mean_beta_2sls_`i'
	g stddev_beta_2sls_`i' = round(`r(sd)',.0001) 
	local stddev_beta_2sls_`i' = stddev_beta_2sls_`i'
	g p05_beta_2sls_`i' =  round(`r(p05)',.0001) 
	local p05_beta_2sls_`i' = p05_beta_2sls_`i'	
	g p95_beta_2sls_`i' = round(`r(p95)',.0001) 
	local p95_beta_2sls_`i' = p95_beta_2sls_`i'	


}

clear
set obs 2
g obs = [_n]

g actual_beta_2sls = -100
g mean_beta_2sls = -100
g stddev_beta_2sls = -100
g p05_beta_2sls = -100
g p95_beta_2sls = -100

g model_2sls = "Model 1: No Age" if obs == 1
replace model_2sls = "Model 2: With Age" if obs == 2

forv i = 1(1)2 {

	replace actual_beta_2sls = `actual_2sls_`i'' if obs == `i'
	replace mean_beta_2sls = `mean_beta_2sls_`i'' if obs == `i'
	replace stddev_beta_2sls = `stddev_beta_2sls_`i'' if obs == `i'
	replace p05_beta_2sls = `p05_beta_2sls_`i'' if obs == `i'
	replace p95_beta_2sls = `p95_beta_2sls_`i'' if obs == `i'

}






