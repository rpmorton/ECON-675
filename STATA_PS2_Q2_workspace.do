***Created by RM on 2018.09.30
***For ECON 675, PS 2, Q 2
*********************

/************
***Question 2: 5a and 5b
***********/

clear
set more off

global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"
global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

set matadebug off


clear all


**Generate Data

global sims = 1000
global obs = 1000


forv s = 1(1)$sims {
	di "s is `s'"

	clear

	set obs $obs
	quietly g x = runiform(-1,1)


	forv i = 1(1)5 {
		g rnorm`i'_sq = (rnormal())^2
	}

	quietly g epsilon = x * x * (rnorm1_sq + rnorm2_sq + rnorm3_sq + rnorm4_sq + rnorm5_sq - 5)

	quietly g y = exp(-.1 * (4*(x^2)-1)^2) * sin(5*x) + epsilon

	quietly g intercept = 1
	
	forv k = 1(1)20 {
		g x_`k' = 0
	}
	
	drop rnorm*
	*pause
	quietly g obs = [_n]
	
	mata: mata clear
	mata: X = st_data(.,("intercept"))
	
	forv k = 1(1)20 {
		replace x_`k' = x^`k'
		
		quietly reg y x_*
		quietly predict resid_`k', resid
		
		*capture eret clear
		
		quietly g x_add = x_`k'
		mata: add = st_data(.,("x_add"))
		mata: X = (X,add)
		mata:st_matrix("X",X)
		*matrix list X
		
		mata: Xinv = luinv(X' * X)
		mata: wmat = X * Xinv * X'
		mata: wmatdiag = diagonal(wmat)
		mata: st_matrix("wmatdiag", wmatdiag)
		*matrix list wmatdiag
		
		quietly g weight_i = wmatdiag[obs,1]
		quietly g overn_cv_`k'= (1 / $obs) * (resid_`k' / (1- weight_i) )^2
		egen cv_`k' = sum(overn_cv_`k' )
		
		drop overn_cv_`k' weight_i x_add
		
		}

	keep if obs == 1
	keep cv_*
	quietly g sim = `s'
	*su *	
		
	if `s' < 1.5 {
		quietly save "$temp/cross val results", replace
	}
		
	if `s' > 1.5 {
		quietly append using "$temp/cross val results"
		quietly save "$temp/cross val results", replace
	}
			

}


use "$temp/cross val results", clear

g k = sim if sim <= 20
g avg_cv_k = .

	forv k = 1(1)20 {
	egen avg_cv_`k' = mean(cv_`k')
	replace avg_cv_k = avg_cv_`k' if k == `k'
	}
	
sort sim

keep if avg_cv_k != .

twoway (connected avg_cv_k k, xtitle("K") ytitle("Average CV(K)") title("PS2 Q2 5b: Average CV(K) as a Function of K"))
graph export "$out/PS2 Q2_5b_CVK.pdf", as (pdf) replace

egen min_avg_cv_k = min(avg_cv_k)
g pre_k_min = k if min_avg_cv_k == avg_cv_k
egen k_min = max(pre_k_min)
drop pre_k_min

global kmin = k_min

/************
***Question 2: 5c
***********/


**FIRST DO SIMULATIONS OF FUNCTION OF X TO GET AVERAGE RELATIONSHIP


**Generate Data

forv s = 1(1)$sims {

clear
set obs $obs

quietly g x = runiform(-1,1)

forv i = 1(1)5 {
		g rnorm`i'_sq = (rnormal())^2
	}

quietly g epsilon = x * x * (rnorm1_sq + rnorm2_sq + rnorm3_sq + rnorm4_sq + rnorm5_sq - 5)

quietly g y = exp(-.1 * (4*(x^2)-1)^2) * sin(5*x) + epsilon

keep x y

g obs = [_n]

expand = 11 if obs == 1, g(expanded)
drop obs
g obs = [_n]

replace y = . if expanded == 1

g expanded_obs = obs - $obs if obs > $obs

replace x = (expanded_obs - 1)  / ( (11 - 2) / 2) - 1 if expanded_obs != .

forv k = 1(1)$kmin {
g x_`k' = x^`k'
}

quietly reg y x_*

predict fitted, xb

g x_0 = 1

mata: mata clear
mata: X = st_data(.,("x_0"))
forv k = 1(1)$kmin {
	capture drop x_add
	quietly g x_add = x_`k'
	mata: add = st_data(.,("x_add"))
	mata: X = (X,add)
}
mata:st_matrix("X",X)
		*matrix list X
		
mata: Xinv = luinv(X' * X)
mata: wmat = X * Xinv * X'
mata: wmatdiag = diagonal(wmat)
mata: st_matrix("wmatdiag", wmatdiag)
*matrix list wmatdiag
		
quietly g weight_i = wmatdiag[obs,1]

*di "finishes mata piece"

g resid = abs(y - fitted)
egen sigma_sq_hat = sum( (1 / ($obs - 1 ) ) * resid^2)
egen var = sum( (1 / $obs ) * weight_i^2 * sigma_sq_hat)

g ub_1 = fitted + 1.96 * sqrt(var)
g ub_2 = fitted + 1.96 * sqrt(var / $obs)
g lb_1 = fitted - 1.96 * sqrt(var)
g lb_2 = fitted - 1.96 * sqrt(var / $obs)

keep if expanded == 1

g sim = `s'
keep x fitted sim ub* lb*


if `s' < 1.5 {
	
	quietly save "$temp/sim 5c", replace
	
	}

if `s' > 1.5 {
	
	quietly append using "$temp/sim 5c"
	quietly save "$temp/sim 5c", replace

	}

}

use "$temp/sim 5c", clear

bys x: egen avg_reg_fct = mean(fitted)
bys x: egen avg_ub_1 = mean(ub_1)
bys x: egen avg_ub_2 = mean(ub_2)
bys x: egen avg_lb_1 = mean(lb_1)
bys x: egen avg_lb_2 = mean(lb_2)

egen tagx = tag(x)
keep if tagx == 1

keep x avg*

save "$temp/sim 5c for merge", replace


clear

global obsc = 1002
	
set obs $obsc
g obs_n = [_n]
g x = (obs_n - 1)  / ( ($obs - 2) / 2) - 1
drop if x > 1
	
g true_gdp = exp(-.1*(4*x-1)^2)*sin(5*x)

forv k = 1(1)$kmin {
	g x_`k' = x^`k'
	}
	
reg true_gdp x_*
predict truereg, xb

merge 1:1 x using "$temp/sim 5c for merge"

sort x
	
twoway (connected avg_reg_fct x, msize(vsmall)) (scatter avg_ub_1 x, msize(small)) (scatter avg_lb_1 x, msize(small))  (scatter truereg x, msize(vsmall) xtitle("x") )

	
	
	
