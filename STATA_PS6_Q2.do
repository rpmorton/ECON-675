/**Created by RM on 2018.12.05
***For ECON 675, PS 6, Q2
*********************/

set more off

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"
global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"
global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"


/* Load Packages */
/*
net install rdrobust, from(https://sites.google.com/site/rdpackages/rdrobust/stata) replace
net install rdlocrand, from(https://sites.google.com/site/rdpackages/rdlocrand/stata) replace
net install rddensity, from(https://sites.google.com/site/rdpackages/rddensity/stata) replace
*/

/*** 2.1: Plots and Falsification ****/

use "$data/HeadStart", clear

**Plots

/* IMSE Evenly Spaced: */

rdplot mort_related_pre povrate60,  binselect(es) p(2)
graph export "$out/STATA_PS6_Q2_1_imse_even.pdf", as (pdf) replace

/* IMSE Quantile Binning: */

rdplot mort_related_pre povrate60,  binselect(qspr) p(2)
graph export "$out/STATA_PS6_Q2_1_imse_quant.pdf", as (pdf) replace

/* Var Quantile Binning: */

rdplot mort_related_pre povrate60,  binselect(esmv) p(2)
graph export "$out/STATA_PS6_Q2_1_var_even.pdf", as (pdf) replace

/* Var Quantile Binning: */

rdplot mort_related_pre povrate60,  binselect(qsmvpr) p(2)
graph export "$out/STATA_PS6_Q2_1_var_quant.pdf", as (pdf) replace


**Falsification

*Histogram
twoway (histogram povrate60 if povrate60 < 0, bcolor(navy) w(2.5) start(-50) freq) ///
       (histogram povrate60 if povrate60 >=0, bcolor(maroon) w(2.5) start(-50) freq)
	   
graph export "$out/STATA_PS6_Q2_1_histplot.pdf", as (pdf) replace

g treat = povrate60 >= 0	   
tab treat if abs(povrate60) <= 1
bitesti 64 35 .5 

/*** 2.2: Parametric Methods ***/

local first = 5

g treat_povrate = povrate60 * treat
g povrate60_2 = povrate60 * povrate60
g treat_povrate_2 = treat * povrate60_2

forv p = 3(1)6 {

	local pminus = `p' - 1
	g povrate60_`p' = povrate60 * povrate60_`pminus'
	g treat_povrate_`p' = treat * povrate60_`p'
	
	reg mort_related_post treat povrate60*
	*su treat2
		
	if `first' < 1 {
		regsave treat using "$temp/PS6_Q2_2_results.dta", tstat addlabel(model,"Constant",Order,`p') append
	}
		
	if `first' > 0 {		
		regsave treat using "$temp/PS6_Q2_2_results.dta", tstat addlabel(model,"Constant",Order,`p') replace
		local first = 0 
	}
		
	reg mort_related_post treat povrate60* treat_povrate*
	regsave treat using "$temp/PS6_Q2_2_results.dta", tstat addlabel(model,"Het",Order,`p') append

}

drop povrate60_*
g povrate60_1 = povrate60
g povrate60_0 = 0

foreach h in 1 5 9 18 {
	
	g h = `h'
	
	forv p = 0(1)2 {

		if `p' < 1 {

		reg mort_related_post treat if -h <= povrate60 & povrate60 <= h
		regsave treat using "$temp/PS6_Q2_2_results.dta", tstat addlabel(model,"Local",Order,`p',Band,`h') append			
		
		}
		
		if `p' > 0 {
		
		capture drop povrate60_2
		local pminus = `p' - 1
		capture g povrate60_`p' = povrate60  * povrate60_`pminus'

		reg mort_related_post treat povrate60* if -h <= povrate60 & povrate60 <= h
		regsave treat using "$temp/PS6_Q2_2_results.dta", tstat addlabel(model,"Local",Order,`p',Band,`h') append			
		
		}
	}
	
	drop h
}

preserve
		
	use "$temp/PS6_Q2_2_results.dta", clear
	export excel using "$out/STATA_2_2.xlsx", replace first(var)
	
restore

g povrate60_3 = povrate60_2 * povrate60
g povrate60_4 = povrate60_3 * povrate60
g povrate60_5 = povrate60_4 * povrate60
g povrate60_6 = povrate60_5 * povrate60 
reg mort_related_post treat povrate* 
predict mort_post_const, xb

twoway (scatter mort_post_const povrate60) (scatter mort_related_post povrate60 if mort_related_post < 10, xline(0) ///
	lcolor(black) ytitle("Mort Related Post") xtitle("Normalized Poverty Rate")  ///
	legend( label(1 "6th Order Poly Fit") label(2 "Observed") ) title("PS6: Q2.2, Constant Treatment") ///
	cap("Note: Observations with Mort Related Post over 10 not shown. ") )
graph export "$out/STATA_PS6_Q2_2_Constant.pdf", as (pdf) replace

reg mort_related_post treat treat_po* povrate* 
predict mort_post_het, xb
twoway (scatter mort_post_het povrate60) (scatter mort_related_post povrate60 if mort_related_post < 10, xline(0) ///
	lcolor(black) ytitle("Mort Related Post") xtitle("Normalized Poverty Rate")  ///
	legend( label(1 "6th Order Poly Fit") label(2 "Observed") ) title("PS6: Q2.2, Heterogeneous Treatment") ///
	cap("Note: Observations with Mort Related Post over 10 not shown. ") )
graph export "$out/STATA_PS6_Q2_2_Het.pdf", as (pdf) replace

reg mort_related_post treat povrate60 povrate60_2  if -18 <= povrate60 & povrate60 <= 18
predict mort_post_local, xb
twoway (scatter mort_post_local povrate60) (scatter mort_related_post povrate60 if mort_related_post < 10, xline(0) ///
	lcolor(black) ytitle("Mort Related Post") xtitle("Normalized Poverty Rate")  ///
	legend( label(1 "2nd Order Poly Fit") label(2 "Observed") ) title("PS6: Q2.2, Local Parametric Method") ///
	cap("Note: A bandwidth of 18 is used. Observations with Mort Related Post over 10 not shown. ") )
graph export "$out/STATA_PS6_Q2_2_Local.pdf", as (pdf) replace


/*** 2.3: Robust Local Polynomial Methods ***/

g row = [_n]
g p = -1
g type = "test"
g point = -100
g std_err = -1
g ci_l = -100000
g ci_u = 100000

forv p = 0(1)2 {
	local row = (3 * `p') + 1

	rdrobust mort_related_post povrate60, all p(`p')
	
	replace point =  e(tau_cl) if row == `row'
	replace point =  e(tau_bc) if row == `row' + 1
	replace point =  e(tau_bc) if row == `row' + 2

	replace std_err =  e(se_tau_cl) if row == `row'
	replace std_err =  e(se_tau_cl) if row == `row' + 1
	replace std_err =  e(se_tau_rb) if row == `row' + 2
	
	replace ci_l =  e(ci_l_cl) if row == `row'
	*replace ci_l =  point + invt(e(N),.025)*std_err   if row == `row' + 1
	replace ci_l =  point + invnormal(.025)*std_err   if row == `row' + 1
	replace ci_l =  e(ci_l_rb) if row == `row' + 2
		
	replace ci_u =  e(ci_r_cl) if row == `row'
	*replace ci_u =  point + invt(e(N),.975)*std_err  if row == `row' + 1
	replace ci_u =  point + invnormal(.975)*std_err   if row == `row' + 1
	replace ci_u =  e(ci_r_rb) if row == `row' + 2
	
	replace p = `p' if row == `row'
	replace p = `p' if row == `row' + 1
	replace p = `p' if row == `row' + 2

	replace type = "Conventional" if row == `row'
	replace type = "Bias-Corrected" if row == `row' + 1
	replace type = "Robust" if row == `row' + 2
	
	
}

preserve

keep point std_err ci* p type
keep if type != "test"

save "$out/PS6_Q2_3_results", replace 
export excel using "$out/STATA_2_3_MSEOptimal.xlsx", replace first(var)

restore

drop point std_err ci* p type

**Different Bandwidths and Kernels

g h = -1
g type = "test"
g point = -100

forv h = 1(1) 10 {
	local row = 3 * (`h' -1) + 1
	
	rdrobust mort_related_post povrate60, all p(1) h(`h') kernel(tri)
	
	replace point = e(tau_bc) if row == `row'
	replace type = "Tri" if row == `row'
	replace h = `h' if row == `row'
	
	local row = `row' + 1
	
	rdrobust mort_related_post povrate60, all p(1) h(`h') kernel(uni)
	
	replace point = e(tau_bc) if row == `row'
	replace type = "Uni" if row == `row'
	replace h = `h' if row == `row'
		
	local row = `row' + 1

	rdrobust mort_related_post povrate60, all p(1) h(`h') kernel(epa)
	
	replace point = e(tau_bc) if row == `row'
	replace type = "Epa" if row == `row'
	replace h = `h' if row == `row'
		
}
	

preserve

	keep point type h
	keep if type != "test"

	save "$out/PS6_Q2_3_results_BandKernel", replace 
	export excel using "$out/STATA_2_3_DiffBandKernel.xlsx", replace first(var)

restore	
	
**Donut Hole	
/* How Many at Cutoff?? */
replace point = -100
replace type = "test"
replace h = -100

count if povrate60 == 0
sort povrate60

g povrate60_order = [_n]
egen pre_max_ctrl = max(povrate60_order) if povrate60 < 0
egen max_ctrl = max(pre_max_ctrl)

g donut_obs = povrate60_order - max_ctrl if povrate60 > 0
replace donut_obs = max_ctrl + 1 - povrate60_order if povrate60 < 0

forv d = 1(1)10 {
		local row = `d'
		rdrobust mort_related_post povrate60 if donut_obs > `d', all p(1) 
		
		replace point = e(tau_bc) if row == `row'
		replace type = "Donut" if row == `row'
		replace h = `d'	 if row == `row'
		
}

preserve

	keep point type h
	keep if type != "test"

	save "$out/PS6_Q2_3_results_DonutHole", replace 
	export excel using "$out/STATA_2_3_DonutHole.xlsx", replace first(var)

restore	

**Placecbo Cutoff
local row = 0
replace point = -100
replace type = "test"
replace h = -100
g pval = -100

foreach c in -10 -8 -6 -4 -2 2 4 6 8 10 {
		local row = `row' + 1
		rdrobust mort_related_post povrate60, all p(1) c(`c')
		
		replace point = e(tau_bc) if row == `row'
		replace pval = 1 - normal(abs(e(tau_bc)) / e(se_tau_cl))  + normal(-abs(e(tau_bc)) / e(se_tau_cl))  if row == `row'
		replace type = "Placebo Cutoff" if row == `row'
		replace h = `c' if row == `row'

}

preserve

	keep  type h point pval
	keep if type != "test"

	save "$out/PS6_Q2_3_results_PlaceboCutoff", replace 
	export excel using "$out/STATA_2_3_PlaceboCutoff.xlsx", replace first(var)

restore	

/*** 2.4: Local Randomization Methods ***/


/* Use the mort_injury_post theoretically unaffected by the program */

rdwinselect povrate60 mort_injury_post 

/* Verify that this works for the pre period as well for the pre-intervention covariate */

reg mort_related_pre treat if abs(povrate60) < 1.44

rdplot mort_injury_post povrate60 if abs(povrate60) < 1.44,  binselect(es) p(0)
graph export "$out/STATA_PS6_Q2_4_imse_injury.pdf", as (pdf) replace

rdplot mort_related_pre povrate60 if abs(povrate60) < 1.44,  binselect(es) p(0)
graph export "$out/STATA_PS6_Q2_4_imse_pre_mort.pdf", as (pdf) replace

rdplot mort_related_post povrate60 if abs(povrate60) < 1.44,  binselect(es) p(0)
graph export "$out/STATA_PS6_Q2_4_imse_treat_mort_post.pdf", as (pdf) replace

*rdrandinf mort_related_post povrate60, wr(1.44) wl(-1.44) reps(999) 
permute treat diffmean=(r(mu_2)-r(mu_1)), reps(999) nowarn: ttest mort_related_post if abs(povrate60)<=1.44, by(treat) 

/* Verify command:
permute treat diffmean=(r(mu_2)-r(mu_1)), reps(999) nowarn: ttest mort_related_post if abs(povrate60)<=5, by(treat) 
*/

/* Use permute command (better run time) to estimate: */

g w = -1
g point_rand = -1000
g std_err_rand = -1000
g pval_rand = -1000

local row = 0

forv w = .8(.2)2.6 {
	di "row is `row'"
	local row = `row' + 1
	replace w = `w' if row == `row'
	
	local window = `w'
	
	permute treat diffmean=(r(mu_2)-r(mu_1)), reps(999) nowarn: ttest mort_related_post if abs(povrate60)<=`window', by(treat)
	
	matrix se = r(se)
	replace std_err_rand = se[1,1] if row == `row'
	matrix pval = r(p)
	replace pval_rand = pval[1,1] if row == `row'

	preserve
	keep if abs(povrate60)<=`window'
	egen count = count(povrate60)
	egen count_treat = sum(treat)
	egen ctrl_mean = sum(mort_related_post * (1-treat)) 
	egen treat_mean = sum(mort_related_post * treat)  
	g diff_mean = ( treat_mean / count_treat) - ( ctrl_mean / (count - count_treat) )
	local diff_mean = diff_mean
	restore
	
	replace point_rand = `diff_mean' if row == `row'

	
}


sort row

preserve

	keep w point_rand std_err_rand pval_rand
	keep if pval_rand != -1000
	
	save "$out/PS6_Q2_4_windowsize", replace 
	export excel using "$out/STATA_2_4_Windows.xlsx", replace first(var)

restore	


