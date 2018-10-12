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

forv s = 1(1)$sims {
	di "s is `s'"

	clear

	global obs = 1000

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
	
