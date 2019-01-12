***Created by RM on 2018.11.24
***For ECON 665, PS 5, Q2
*********************

clear
set more off

global sims = 5000
*global sims = 50
global obs = 200
global cov = .99
global beta = 0

global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"
global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"

local first = 1
forv s = 1(1)$sims {

clear
set obs $obs
g z = rnormal(0,1)
g u = rnormal(0,1)

g indep_norm = rnormal(0,1)
g weight_indep_norm = sqrt(1- ($cov * $cov) )
g v = $cov * u + weight_indep_norm * indep_norm

local gamma1 = sqrt( (1/$obs) * 0 )
local gamma2 = sqrt( (1/$obs) * 0.25 )
local gamma3 = sqrt( (1/$obs) * 9 )
local gamma4 = sqrt( (1/$obs) * 99 )

forv i = 1(1)4 {

	g x`i' = `gamma`i'' * z + v
	g y`i' = $beta * x`i' + u
	qui: reg y`i' x`i'
	local stderr = _se[x`i']
	
	if `first' < 1 {
		regsave x`i' using "$temp/PS5_Q2_results.dta", tstat addlabel(gamma,`gamma`i'',model,"OLS",sim,`s') append
	}
	
	if `first' > 0 {
		regsave x`i' using "$temp/PS5_Q2_results.dta", tstat addlabel(gamma,`gamma`i'',model,"OLS",sim,`s') replace
		local first = 0
	}
	
	qui: reg x`i' z
	local fstat = `e(F)'
	
	qui: ivregress 2sls y`i' (x`i' = z)
	regsave x`i' using "$temp/PS5_Q2_results.dta", tstat addlabel(gamma,`gamma`i'',model,"2SLS",sim,`s',Fstat,`fstat') append 

	
	
}
}

/* Clean Results to Make Table */
	
use "$temp/PS5_Q2_results.dta", clear
g gamma_actual = round(gamma^2 * $obs , .01)

g sig_beta = abs(tstat) > 1.96 

*Now Collapse to make one row per table *
bys gamma_actual sim: egen pre_beta_ols = max(coef) if model == "OLS"
bys gamma_actual sim: egen pre_beta_ols_se = max(stderr) if model == "OLS"
bys gamma_actual sim: egen pre_sig_beta_ols = max(sig_beta) if model == "OLS"
bys gamma_actual sim: egen pre_beta_2sls = max(coef) if model == "2SLS"
bys gamma_actual sim: egen pre_beta_2sls_se = max(stderr) if model == "2SLS"
bys gamma_actual sim: egen pre_sig_beta_2sls = max(sig_beta) if model == "2SLS"
bys gamma_actual sim: egen pre_F_2sls = max(Fstat) if model == "2SLS"

local varlist "beta_ols beta_ols_se sig_beta_ols beta_2sls beta_2sls_se sig_beta_2sls F_2sls"

foreach var of local varlist {
	by gamma_actual sim: egen `var' = max(pre_`var')
}

egen tag = tag(gamma_actual sim)
keep if tag == 1


levelsof gamma, local(gams)

foreach var of local varlist {
	g mean_`var' =   -10000
	g sd_`var' = -10000
	g p10_`var' =  -10000
	g p50_`var' =  -10000
	g p90_`var' =  -10000
	
	foreach g of local gams {
		*di "g is `g' and var is `var'"
		su `var' if gamma == `g', d
		
		replace mean_`var' =   round(`r(mean)',.0001) if gamma == `g'
		replace sd_`var' = round(`r(sd)',.0001) if gamma == `g'
		replace p10_`var' =  round(`r(p10)',.0001) if gamma == `g'
		replace p50_`var' =  round(`r(p50)',.0001) if gamma == `g'
		replace p90_`var' =  round(`r(p90)',.0001) if gamma == `g'
	}
	
}
	
keep gamma_actual mean* sd* p10* p50* p90*
duplicates drop
	
reshape long mean_ sd_ p10_ p50_ p90_ , i(gamma_actual) j(parameter) string

g param_sort = 1 if parameter == "beta_ols"
replace param_sort = 2 if parameter == "beta_ols_se"
replace param_sort = 3 if parameter == "sig_beta_ols"
replace param_sort = 4 if parameter == "beta_2sls"
replace param_sort = 5 if parameter == "beta_2sls_se"
replace param_sort = 6 if parameter == "sig_beta_2sls"
replace param_sort = 7 if parameter == "F_2sls"

sort gamma param_sort
drop param_sort

texsave using "$out/STATA_PS5_Q2.tex", replace

g amper = "&"
egen tex = concat(mean_ amper sd_ amper p10_ amper p50_ amper p90_)
drop amper

export excel using "$out/STATA_PS5_Q2_excel.xlsx", replace firstrow(variables)

