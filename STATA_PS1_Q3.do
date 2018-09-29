***Created by RM on 2018.09.33
***For ECON 675, PS 1, Q 3
*********************

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"
global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

clear
set more off

import delim "$data/LaLonde_1986", delim(",")

*Question 3:1b
g ones = 1

bys treat: egen count_by_treat = sum(ones)
bys treat: egen mean_earn78_by_treat = mean(earn78)

g dev_from_mean_by_treat = (earn78 - mean_earn78_by_treat) *  (earn78 - mean_earn78_by_treat)

bys treat: egen sum_dev_from_mean_by_treat =  sum(dev_from_mean_by_treat)
g s2_by_treat = (1/(count_by_treat-1)) * sum_dev_from_mean_by_treat

forv i = 0(1)1 {

	egen pre_s2_treat_`i' = max(s2_by_treat) if treat == `i'
	egen s2_treat_`i' = max(pre_s2_treat_`i')
	drop pre_s2_treat_`i'
	
	
	egen pre_mean_treat_`i' = max(mean_earn78_by_treat) if treat == `i'
	egen mean_treat_`i' = max(pre_mean_treat_`i')
	drop pre_mean_treat_`i'
	
	egen pre_count_treat_`i' = max(count_by_treat) if treat == `i'
	egen count_treat_`i' = max(pre_count_treat_`i')
	drop pre_count_treat_`i'	
	
}

g treat_dm = mean_treat_1 - mean_treat_0
g conserv_var = (s2_treat_1/count_treat_1) + (s2_treat_0/count_treat_0)

*g df = 100 * ((count_treat_1 + count_treat_0) - 4)
g df = (count_treat_1 + count_treat_0) - 2

local alpha = .95
g lb_treat_lalonde = treat_dm + invttail(df,`alpha'+(1-`alpha')/2)*sqrt(conserv_var)
g ub_treat_lalonde =  treat_dm + invttail(df,(1-`alpha')/2)*sqrt(conserv_var)

local lb = lb_treat_lalonde
local ub = ub_treat_lalonde

di "the CI is `lb' to `ub'"
ttest earn78, by(treat)

/************
***Question 3: 2a and 2b
***********/

clear
set more off

import delim "$data/LaLonde_1986", delim(",")

global loops = 999

egen count_treat = sum(treat)
egen count_total = count(treat)

local treatcount = count_treat
local totalcount = count_total

**Create Random Vectors
clear

set obs `totalcount'

g ones = 1
g obs = [_n]

g rand1 = runiform()
egen pre_rank1 = rank(rand1), field

bys rand1: egen rank1 = max(pre_rank1)

g fisher1 = rank1 <= `treatcount'

drop rand1 pre_rank1

di "breaks -1"

local counter = 2

local loopsplus = `loops' + 1

while `counter' < `loopsplus' {

	di "counter is `counter'"
	*di "enters first loop"
	
	g rand`counter' = runiform()
	
	egen pre_rank`counter' = rank(rand`counter'), field
	
	*di "breaks .5"

	bys rand`counter': egen rank`counter' = max(pre_rank`counter')

	g fisher`counter' = rank`counter' <= `treatcount'

	*di "breaks 0"
	
	local counterminus = `counter' - 1
	
	local matchcounter = 0
	local matchon = 0
	
	*di "breaks 1"
	
	drop rand`counter' pre_rank`counter' rank`counter'
	
	while `matchon'  < 1 {
	
	local matchcounter = `matchcounter' + 1
	
	*di "matchcounter is `matchcounter'"
	
		g pre_compare_rand = fisher`matchcounter' == fisher`counter'
		egen compare_rand = min(pre_compare_rand)
		
		g perfectmatch = compare_rand == 1
		local match = perfectmatch
		
		drop pre_compare_rand compare_rand perfectmatch
		
		*di "breaks 2"
		
		if `match' > 0 {
		
			drop fisher`matchcounter'
			local matchon = 10
			local matchcounter = 1
		
			}
		
		*di "breaks 3"

		if `matchcounter' == `counterminus' {
		
			local matchon = 10
			local counter = `counter' + 1
			local matchcounter = 1
		
		}
		
		g matchcounter = `matchcounter'
		g matchnext = matchcounter + 1
		*su matchcounter
		*su matchnext
		
		*di "matchcounter is `matchcounter'"
		macro drop matchcounter
		local `matchcounter' = matchnext
		local j = matchnext
		
		*di "matchcounter is `matchcounter' and matchnext is `j'"
		
		drop matchcounter matchnext
	
		*di "breaks 4"

	}

}


save "$data/random vectors", replace

*begin to do analysis

clear
set more off

import delim "$data/LaLonde_1986", delim(",")

g obs = [_n]

merge 1:1 obs using "$data/random vectors"

forv i = 1(1)$loops {
*forv i = 1(1)`test' {

	*difference in means
	bys fisher`i': egen mean_effect_`i' = mean(earn78)
	
	egen pre_control_`i' = min(mean_effect_`i') if fisher`i' == 0
	egen control_`i' = max(pre_control_`i')
	
	egen pre_treatment_`i' = min(mean_effect_`i') if fisher`i' == 1
	egen treatment_`i' = max(pre_treatment_`i')	
	
	g treat_diff_means_`i' = treatment_`i' - control_`i'
	
	drop mean_effect_`i' pre_control_`i' control_`i' pre_treatment_`i' treatment_`i'
	
	*ks test
	ksmirnov earn78, by(fisher`i')
		
	g ks_stat_`i' = r(D)

}

*Actual Treatment Mean

bys treat: egen mean_effect = mean(earn78)
	
egen pre_control = min(mean_effect) if treat == 0
egen control = max(pre_control)
	
egen pre_treatment = min(mean_effect) if treat == 1
egen treatment = max(pre_treatment)	
	
g treat_diff_means_actual = treatment - control

global treatdiffmeansact = treat_diff_means_actual
	
drop mean_effect pre_control control pre_treatment treatment

*Actual Treatment KS

ksmirnov earn78, by(treat)
g ks_stat_act = r(D)

global ksstatact = ks_stat_act

save "$data/fisher treatment effect distribution", replace

keep if obs == 1

keep treat_diff* ks*

drop treat_diff_means_actual ks_stat_act

save "$data/fisher results", replace

**Now get  p value for effect and confidence interval for diff means

use "$data/fisher results", clear

keep treat_diff*

xpose, clear
sort v1
g treat_diff_means_actual = $treatdiffmeansact
g abs_treat_diff = abs(v1)

g diff_gt_actual = abs_treat_diff > treat_diff_means_actual
egen count_gt_acutal = sum(diff_gt_actual)

g pval_diff_means = count_gt_acutal / $loops
tab pval_diff_means

g counter = [_n]

g lb_ci_diff_means_flag = counter == round(.025 * $loops)
g pre_lb_ci_diff_means = v1 if lb_ci_diff_means_flag == 1
egen lb_ci_diff_means = max(pre_lb_ci_diff_means)
drop pre_lb_ci_diff_means

g ub_ci_diff_means_flag = counter == round(.975 * $loops)
g pre_ub_ci_diff_means = v1 if ub_ci_diff_means_flag == 1
egen ub_ci_diff_means = max(pre_ub_ci_diff_means)
drop pre_ub_ci_diff_means

tab lb_ci_diff_means
tab ub_ci_diff_means


**Now get  p value for effect and confidence interval for ks test

use "$data/fisher results", clear

keep ks*

xpose, clear
sort v1
g treat_ks_actual = $ksstatact
g abs_ks = abs(v1)

g diff_gt_actual = abs_ks > treat_ks_actual
egen count_gt_acutal = sum(diff_gt_actual)

g pval_ks = count_gt_acutal / $loops
tab pval_ks

g counter = [_n]

g lb_ci_ks_flag = counter == round(.025 * $loops)
g pre_lb_ci_ks = v1 if lb_ci_ks_flag == 1
egen lb_ci_ks = max(pre_lb_ci_ks)
drop pre_lb_ci_ks

g ub_ci_ks_flag = counter == round(.975 * $loops)
g pre_ub_ci_ks = v1 if ub_ci_ks_flag == 1
egen ub_ci_ks = max(pre_ub_ci_ks)
drop pre_ub_ci_ks

tab lb_ci_ks
tab ub_ci_ks


/************
***Question 3: 3a and 3b
***********/

*3a: Power Fct Using Size of .05 

clear
set more off

global size = .05

import delim "$data/LaLonde_1986", delim(",")

g ones = 1

bys treat: egen count_by_treat = sum(ones)
bys treat: egen mean_earn78_by_treat = mean(earn78)

g dev_from_mean_by_treat = (earn78 - mean_earn78_by_treat) *  (earn78 - mean_earn78_by_treat)

bys treat: egen sum_dev_from_mean_by_treat =  sum(dev_from_mean_by_treat)
g s2_by_treat = (1/(count_by_treat-1)) * sum_dev_from_mean_by_treat

forv i = 0(1)1 {

	egen pre_s2_treat_`i' = max(s2_by_treat) if treat == `i'
	egen s2_treat_`i' = max(pre_s2_treat_`i')
	drop pre_s2_treat_`i'
	
	egen pre_mean_treat_`i' = max(mean_earn78_by_treat) if treat == `i'
	egen mean_treat_`i' = max(pre_mean_treat_`i')
	drop pre_mean_treat_`i'
	
	egen pre_count_treat_`i' = max(count_by_treat) if treat == `i'
	egen count_treat_`i' = max(pre_count_treat_`i')
	drop pre_count_treat_`i'	
	
}

g treat_dm = mean_treat_1 - mean_treat_0
g conserv_var = (s2_treat_1/count_treat_1) + (s2_treat_0/count_treat_0)

g df =  (count_treat_1 + count_treat_0) - 2


*global t_hat = treat_dm
global conserv_se = sqrt(conserv_var)
global df = df

clear

set obs 4001
g obs = [_n]
g tau_power = obs - 1 - 2000

g t_input_lt_tau = -(tau_power / $conserv_se) - invttail($df,1 - ($size/2))
g t_prob_lt_tau = ttail($df,t_input_lt_tau)
g t_input_gt_tau = -(tau_power / $conserv_se) + invttail($df,1 - ($size/2))
g t_prob_gt_tau = ttail($df,t_input_gt_tau)
g power_tau = t_prob_lt_tau + 1 - t_prob_gt_tau
*g power_tau = 1 + ttail(-(tau_power / $conserv_se) - invttail($df,1 - (`size'/2))),$df)
g size = $size


twoway (scatter power_tau tau_power, msize(tiny) xtitle("{&tau}") ytitle("{&beta}({&tau})") title("Power Function against H{sub:0}: {&tau} = 0") note("Power function for df = $df with size = $size") yline($size) )
graph export "$out/PS1 Q3a STATA.pdf", as (pdf) replace


*Question 3:3b--GRID SEARCH
clear
set more off

set obs 1
g obs = 1
g tildet = 0

local on = 1
local tildet = 10

while `on' > 0 {

	capture drop power on_update
	g power =  normal(-1.96 - `tildet' ) - normal(1.96 - `tildet') + 1
	g on_update = power >= .8
	local on = on_update
	
	replace tildet = `tildet'
	
	if `on' > 0 {
		local tildet = `tildet' - .01
	}
	
}


	








