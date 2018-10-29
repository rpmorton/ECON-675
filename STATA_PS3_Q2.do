***Created by RM on 2018.10.27
***For ECON 675, PS 3, Q 2
*********************
 
set more off

set seed 42338992

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"
global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"
global output "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

capture program drop prepexport
program define prepexport

g problem = "`1'"
g type = "`2'"
g point_estimate = .
g ci_lower_95 = .
g ci_upper_95 = .

	forv v = 1(1)4{
		replace point_estimate = initial_b`v' if obs == `v'
		replace ci_lower_95 = ci_lower_`v' if obs == `v'
		replace ci_upper_95 = ci_upper_`v' if obs == `v'

	}


sort obs
keep if obs <= 4
keep problem type point_estimate ci_lower_95 ci_upper_95 

export excel using "$output/`2'_`1'.xlsx", replace
	
end

	
/************
***Question 2: 2b
***********/

global bstrapreps = 100
global cialpha = .05

import delim using "$data/pisofirme", delim(",") clear

/* Booststrap Estimator: MCAR --drop obs with missing ovservations*/
drop if danemia == . | dpisofirme == . | s_hhpeople == . | s_age == . | s_incomepc == .
keep danemia dpisofirme s_hhpeople s_age s_incomepc

g ln_s_incomepc = ln(s_incomepc + 1)

*gmm (crime - policepc*{b1} - legalwage*{b2} - {b3}), ///	. gmm (crime - policepc*{b1} - legalwage*{b2} - {b3}), ///

gmm (danemia - (1 + exp( (-1) * ( dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4} ) ) )^(-1) ), ///
	instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep
	
matrix beta = e(b)
	
g initial_b1 = beta[1,1]
g initial_b2 = beta[1,2]
g initial_b3 = beta[1,3]
g initial_b4 = beta[1,4]

/* Now bootstrap: skip to STATA predone code for it */


forv i = 1(1)4 {
	g bootstrap_`i' = .
}


g obs = [_n]

egen obscount = max(obs)

local obscount = obscount

save "$temp/preboostrap", replace

preserve

forv b = 1(1)$bstrapreps {

	bsample `obscount'
	di " b is `b'"

	quietly: gmm (danemia - (1 + exp( (-1) * ( dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4} ) ) )^(-1) ), ///
		instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep
		
	*gmm (danemia - (1 + exp( dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4}) )^(-1) ), ///
		instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep	
		
	matrix bstrapbeta = e(b)

	forv i = 1(1)4 {
		g beta`i' = bstrapbeta[1,`i']
		local beta`i' = beta`i'
	}

	restore 

	forv j = 1(1)4 {
		replace bootstrap_`j' = `beta`j'' if obs == `b'
	} 
	
	preserve
}

restore

/* Find CI */

g ci_lower_obs = round($bstrapreps * (($cialpha )/2))
g ci_upper_obs = round($bstrapreps *  (1-($cialpha /2)))

	forv i = 1(1)4 {
		sort bootstrap_`i'
		
		g obs_bstrap_`i' = [_n]
		g pre_ci_lower_`i' = bootstrap_`i' if obs_bstrap_`i' == ci_lower_obs
		g pre_ci_upper_`i' = bootstrap_`i' if obs_bstrap_`i' == ci_upper_obs
		
		egen ci_lower_`i' = max(pre_ci_lower_`i')
		egen ci_upper_`i' = max(pre_ci_upper_`i')
		
		drop obs_bstrap_`i' pre*

	}
su ci_lower_1 ci_upper_1 ci_lower_2 ci_upper_2 ci_lower_3 ci_upper_3 ci_lower_4 ci_upper_4



prepexport mcar 2b

/*

 bootstrap, reps($bstrapreps): 	gmm (danemia - (1 + exp( dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4}) )^(-1) ), ///
		instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep
		
*/
/************
***Question 2: 3c
***********/	


import delim using "$data/pisofirme", delim(",") clear


/* Predict First Stage */

g notmissing = 1 - dmissing

g ln_s_incomepc = ln(s_incomepc + 1)

keep danemia dpisofirme s_hhpeople s_age ln_s_incomepc notmissing

logit notmissing dpisofirme s_hhpeople s_age ln_s_incomepc

predict pred_val

g si_over_po = 1 / pred_val 

gmm (si_over_po * (danemia - (1 + exp( (-1) * (dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4} ) ) )^(-1) )), ///
	instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep

matrix beta = e(b)
	
g initial_b1 = beta[1,1]
g initial_b2 = beta[1,2]
g initial_b3 = beta[1,3]
g initial_b4 = beta[1,4]

/* Now bootstrap: skip to STATA predone code for it */


forv i = 1(1)4 {
	g bootstrap_`i' = .
}


g obs = [_n]

egen obscount = max(obs)

local obscount = obscount

save "$temp/preboostrap", replace

preserve

forv b = 1(1)$bstrapreps {

	bsample `obscount'
	di " b is `b'"
	
	logit notmissing dpisofirme s_hhpeople s_age ln_s_incomepc

	predict pred_val_bstrap
	
	g si_over_po_bstrap = 1 / pred_val_bstrap 

	quietly: gmm (si_over_po_bstrap*(danemia - (1 + exp( (-1) * (dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4} ) ) )^(-1) ) ), ///
		instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep

	matrix bstrapbeta = e(b)

	forv i = 1(1)4 {
		g beta`i' = bstrapbeta[1,`i']
		local beta`i' = beta`i'
	}

	restore 

	forv j = 1(1)4 {
		replace bootstrap_`j' = `beta`j'' if obs == `b'
	} 
	
	preserve
}

restore

/* Find CI */

g ci_lower_obs = round($bstrapreps * (($cialpha )/2))
g ci_upper_obs = round($bstrapreps *  (1-($cialpha /2)))

	forv i = 1(1)4 {
		sort bootstrap_`i'
		
		g obs_bstrap_`i' = [_n]
		g pre_ci_lower_`i' = bootstrap_`i' if obs_bstrap_`i' == ci_lower_obs
		g pre_ci_upper_`i' = bootstrap_`i' if obs_bstrap_`i' == ci_upper_obs
		
		egen ci_lower_`i' = max(pre_ci_lower_`i')
		egen ci_upper_`i' = max(pre_ci_upper_`i')
		
		drop obs_bstrap_`i' pre*

	}
su ci_lower_1 ci_upper_1 ci_lower_2 ci_upper_2 ci_lower_3 ci_upper_3 ci_lower_4 ci_upper_4


prepexport mar 3c



/************
***Question 2: 3d
***********/	



import delim using "$data/pisofirme", delim(",") clear


/* Predict First Stage */

g notmissing = 1 - dmissing

g ln_s_incomepc = ln(s_incomepc + 1)

keep danemia dpisofirme s_hhpeople s_age ln_s_incomepc notmissing

logit notmissing dpisofirme s_hhpeople s_age ln_s_incomepc

predict pred_val


g si_over_po = 1 / pred_val 

preserve

drop if pred_val < .1

gmm (si_over_po * (danemia - (1 + exp( (-1) * ( dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4} ) ) )^(-1) )), ///
	instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep

restore
	
matrix beta = e(b)
	
g initial_b1 = beta[1,1]
g initial_b2 = beta[1,2]
g initial_b3 = beta[1,3]
g initial_b4 = beta[1,4]

/* Now bootstrap: skip to STATA predone code for it */


forv i = 1(1)4 {
	g bootstrap_`i' = .
}


g obs = [_n]

egen obscount = max(obs)

local obscount = obscount

save "$temp/preboostrap", replace

preserve

forv b = 1(1)$bstrapreps {

	bsample `obscount'
	di " b is `b'"
	
	logit notmissing dpisofirme s_hhpeople s_age ln_s_incomepc

	predict pred_val_bstrap
	
	g si_over_po_bstrap = 1 / pred_val_bstrap 

	drop if pred_val_bstrap < .1

	quietly: gmm (si_over_po_bstrap*(danemia - (1 + exp( (-1) * ( dpisofirme * {b1} + s_age * {b2} + s_hhpeople * {b3} + ln_s_incomepc * {b4} ) ) )^(-1) ) ), ///
		instruments( dpisofirme s_age s_hhpeople ln_s_incomepc) nolog onestep

	matrix bstrapbeta = e(b)

	forv i = 1(1)4 {
		g beta`i' = bstrapbeta[1,`i']
		local beta`i' = beta`i'
	}

	restore 

	forv j = 1(1)4 {
		replace bootstrap_`j' = `beta`j'' if obs == `b'
	} 
	
	preserve
}

restore

/* Find CI */

g ci_lower_obs = round($bstrapreps * (($cialpha )/2))
g ci_upper_obs = round($bstrapreps *  (1-($cialpha /2)))

	forv i = 1(1)4 {
		sort bootstrap_`i'
		
		g obs_bstrap_`i' = [_n]
		g pre_ci_lower_`i' = bootstrap_`i' if obs_bstrap_`i' == ci_lower_obs
		g pre_ci_upper_`i' = bootstrap_`i' if obs_bstrap_`i' == ci_upper_obs
		
		egen ci_lower_`i' = max(pre_ci_lower_`i')
		egen ci_upper_`i' = max(pre_ci_upper_`i')
		
		drop obs_bstrap_`i' pre*

	}
su ci_lower_1 ci_upper_1 ci_lower_2 ci_upper_2 ci_lower_3 ci_upper_3 ci_lower_4 ci_upper_4


prepexport martrim 3d




	
