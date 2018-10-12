***Created by RM on 2018.09.30
***For ECON 675, PS 2, Q 1
*********************

/************
***Question 1: 3a
***********/

set more off

global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"

**Generate Data
clear

global obs = 20000

set obs $obs

g obs = [_n]
g x = (obs - ($obs /2))*(1/4)

/*
\frac{1}{2} \bigg(\frac{1}{\sqrt{3\pi}}exp(-\frac{x^2+3x + 2.25}{3})(-\frac{2x}{3} - 1)^2 - \frac{2}{3\sqrt{3\pi}}exp(-\frac{x^2+3x + 2.25}{3})\\
    &+ \frac{1}{\sqrt{2\pi}}exp(-\frac{x^2-2x + 1}{2})(-x + 1)^2-\frac{1}{\sqrt{2\pi}}exp(-\frac{x^2-2x + 1}{2}) \bigg) 
*/

g pi = _pi
g f_term1 = (1/sqrt(3*pi))*exp(-(x^2 + 3 * x+ 2.25)/3)*((-2*x)/3 - 1)^2
g f_term2 = (2/(3*sqrt(3*pi)))*exp(-(x^2 + 3 * x+ 2.25)/3)
g f_term3 = (1/sqrt(2*pi))*exp(-(x^2 - 2 * x + 1) / 2)*((-1*x)+1)^2
g f_term4 = (1/sqrt(2*pi))*exp(-(x^2 - 2 * x + 1) / 2)

g f_comb = (1/2) * (f_term1 - f_term2 + f_term3 - f_term4)
g f_comb_square = f_comb^2

preserve
keep if obs <= 2
egen max = max(x) 
egen min = min(x)
g unitsize = max - min
local unit = unitsize
restore


g f_comb_square_unit = `unit' * f_comb_square
egen integral_val = sum(f_comb_square_unit)


/************
***Question 1: 3b
***********/

**First form the actual density estimator per value of h

global minh = 5
global maxh = 8

*forv h = 5(1)15 {
forv h = $minh (1)$maxh {

clear

global obs = 1000

set obs $obs

g obs_n = [_n]

g norm1 = rnormal()*1.5-1.5
g norm2 = rnormal()*1 + 1

g x = .5 * norm1 + .5 * norm2
sort x

egen min_x = min(x)
egen max_x = max(x)
g distance = round(max_x,.1)+ .3 - (round(min_x,.1) - .3)

g x_estimator = min_x + (obs_n - 1 ) / distance

g kernel_estimator_`h' = .

g h_`h' = `h' / 10

forv i = 1(1)$obs {
*forv i = 100(1)100 {

	g pre_relevant_obs = x_estimator if obs_n == `i'
	egen relevant_obs = max(pre_relevant_obs)
	drop pre_relevant_obs
	
	g kernel_u = (x - relevant_obs) / h_`h'
	g kernel_rel = abs(kernel_u) <= 1
	g kernel_estimator_val_obs = kernel_rel * .75 * (1 - kernel_u^2)
	egen kernel_estimator_sum = sum(kernel_estimator_val_obs)
	
	replace kernel_estimator_`h' = kernel_estimator_sum / (h_`h' * $obs)  if obs_n == `i'

	drop kernel_u kernel_rel relevant_obs kernel_estimator_val_obs kernel_estimator_sum

}

* twoway (scatter kernel_estimator x)

save "$temp/kernel estimator `h'", replace
save "$temp/kernel estimator for merge `h'", replace

**Now do leave out one estimator

global sims = 100

forv s= 1(1)$sims {

use "$temp/kernel estimator `h'", clear

g norm1loop = rnormal()*1.5-1.5
g norm2loop = rnormal()*1 + 1

g xloop = .5 * norm1loop + .5 * norm2loop
sort xloop

g kernel_one_out_estimator_`s'_`h' = .

*di "enters sims loop.  s is `s'"

drop if obs_n == `s'

g ones = 1

egen total_obs_leave_one_out = sum(ones)
local totaloneout = total_obs_leave_one_out

g obs_counter_leave_one_out = [_n]

	forv j = 1(1)`totaloneout' {
		
	g pre_relevant_obs = x_estimator if obs_counter_leave_one_out == `j'
	egen relevant_obs = max(pre_relevant_obs)
	drop pre_relevant_obs
	
	g kernel_u = (xloop - relevant_obs) / h_`h'
	g kernel_rel = abs(kernel_u) <= 1
	g kernel_estimator_val_obs = kernel_rel * .75 * (1 - kernel_u^2)
	egen kernel_estimator_sum = sum(kernel_estimator_val_obs)
	
	replace kernel_one_out_estimator_`s'_`h' = kernel_estimator_sum / (h_`h' * `totaloneout')  if obs_counter_leave_one_out == `j'
	
	drop kernel_u kernel_rel relevant_obs kernel_estimator_val_obs kernel_estimator_sum

	}
	
*di "right before crash"	
	
*g crash = 
	
keep x_estimator kernel_one_out_estimator_`s'_`h'

merge 1:1 x_estimator using "$temp/kernel estimator for merge `h'"
drop _merge

save "$temp/kernel estimator for merge `h'", replace
	

}

use "$temp/kernel estimator for merge `h'", clear

}

use "$temp/kernel estimator for merge $minh", clear
g minh = $minh
g minhplus = minh + 1
local minhplus =minhplus
keep x_estimator  kernel_estimator* kernel_one*	

save "$temp/combine kernel estimators", replace

forv h = `minhplus'(1)$maxh {
	use "$temp/combine kernel estimators", clear
	merge 1:1 x_estimator using "$temp/kernel estimator for merge `h'"
	
	save "$temp/combine kernel estimators", replace
	
}

g pi = _pi
g density_actual_norm1 = (1 / (sqrt(2*pi * 1.5)) * exp(-1*(x +1.5)^2*(1/(2*1.5)))
g density_actual_norm2 = (1 / (sqrt(2*pi)) * exp(-1*(x - 1)^2*(1/(2)))

g density_acutal = .5 * (density_actual_norm1 + density_actual_norm2)


forv h = `minhplus'(1)$maxh {
	forv s= 1(1)$sims {
		g pre_imse_li_`s'_`h' =  



