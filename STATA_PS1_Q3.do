***Created by RM on 2018.09.33
***For ECON 675, PS 1, Q 3
*********************

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"

clear
set more off

import delim "$data/LaLonde_1986", delim(",")

*Part 1b
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

g df = (count_treat_1 + count_treat_0) - 4

local alpha = .95
g lb_treat_lalonde = treat_dm + invttail(df,`alpha'+(1-`alpha')/2)*sqrt(conserv_var)
g ub_treat_lalonde =  treat_dm + invttail(df,(1-`alpha')/2)*sqrt(conserv_var)
