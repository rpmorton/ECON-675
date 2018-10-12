***Created by RM on 2018.09.30
***For ECON 675, PS 2, Q 1
*********************

/************
***Question 1: 3a
***********/

set more off

global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"
global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

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

g unit = `unit'
su unit
g f_comb_square_unit = `unit' * f_comb_square
egen integral_val = sum(f_comb_square_unit)


/************
***Question 1: 3b
***********/

**First form the actual density estimator per value of h

global minh = 5
global maxh = 15

*forv h = 5(1)15 {

clear

global obs = 500

set obs $obs

g obs_n = [_n]

g norm1 = rnormal()*1.5-1.5
g norm2 = rnormal()*1 + 1

g x = .5 * norm1 + .5 * norm2
sort x

egen min_x = min(x)
egen max_x = max(x)
g distance = round(max_x,.1)+ .5 - (round(min_x,.1) - .5)

g x_estimator = min_x + (obs_n - 1 ) / distance

keep x_estimator obs_n

forv h = $minh (1)$maxh {
	
	g imse_l1_`h' = 0
	g imse_l0_`h' = 0
	g imse_mle_`h' = 0

}

di "breaks before simulation loops"

global sims = 1000

forv s= 1(1)$sims {

	g norm1 = rnormal()*1.5-1.5
	g norm2 = rnormal()*1 + 1

	g x = .5 * norm1 + .5 * norm2
	sort x
		
	g pi = _pi
	g density_actual_norm1 = (1 / (sqrt(2*pi * 1.5))) * exp(-1*(x +1.5)^2*(1/(2*1.5)))
	g density_actual_norm2 = (1 / (sqrt(2*pi))) * exp(-1*(x - 1)^2*(1/(2)))
	g density_acutal = .5 * (density_actual_norm1 + density_actual_norm2)
	
	
	di "breaks before h loop"
	
	forv h = $minh (1)$maxh {

		capture drop h
		g h = `h' / 10
		local hsmall = h
		drop h


		kdensity(x), bwidth(`hsmall') at(x_estimator) nograph generate(kernel_est_`h'_`s')
		
		g add_imse_l1_`h' = abs(kernel_est_`h'_`s' - density_acutal)^2
		egen sum_imse_l1_`h' = sum(add_imse_l1_`h')
		replace imse_l1_`h' = imse_l1_`h' + sum_imse_l1_`h'

			forv i = 1(1)$obs {
				g x_leaveout = x if obs_n != `i'
				kdensity(x_leaveout), bwidth(`hsmall') at(x_estimator) nograph generate(kernel_no_`i'_est_`h'_`s')
				
				g add_imse_l0_`h' = abs(kernel_no_`i'_est_`h'_`s' - density_acutal)^2 if obs_n == `i'
				egen sum_imse_l0_`h' = sum(add_imse_l0_`h')
				replace imse_l0_`h' = imse_l0_`h' + sum_imse_l0_`h'
				
				drop sum* add* x_leaveout kernel_no_`i'_est_`h'_`s'
				
			}
		
		drop kernel_est_`h'_`s'

		g ones = 1
		egen num_obs = sum(ones)
		egen mu_hat = sum((1 / num_obs) * x)
		local muhat = mu_hat
		egen sigma_hat = sum((1 / num_obs) * (x - mu_hat)^2)
		local sigmahat = sigma_hat
		
		preserve
			clear
			set obs 20000
			g obs = [_n]
			g x = (obs - (20000 /2))*(1/4)
			g pi = _pi
			g f_term1 = (1/sqrt(2*pi*`sigmahat'))*exp(-(x - `muhat')^2/(2*`sigmahat'))*(-1*(x - `muhat')/(`sigmahat'))^2
			g f_term2 = (1/sqrt(2*pi*`sigmahat'))*exp(-(x - `muhat')^2/(2*`sigmahat'))*(1/`sigmahat')
			g f_comb = f_term1 + f_term2 
			*di "error before f comb square"
			*su f_comb
			g f_comb_square = f_comb^2
			g f_comb_square_unit = .25 * f_comb_square
			*su f_comb_square_unit
			egen integral_val = sum(f_comb_square_unit)
			local integralval = integral_val	
			di "error before restore"
		restore
		
		
		g pre_haimse_mle = ((1/ $obs) * (1/5)^2 * `integralval')^(1/5)
		replace imse_mle_`h' = imse_mle_`h' + pre_haimse_mle
		drop pre_haimse_mle ones num_obs mu_hat sigma_hat
		
		}
		
	drop norm1 norm2 x pi 
	drop density_acutal density*
		
	}	
	
	
g h_forgraph = .
g imse_li_h = .
g imse_l0_h = .
g aimse_hat_h = .

local counter = 1

sort obs_n

forv h = $minh (1)$maxh {
	

	g adj_imse_l1_`h' = imse_l1_`h' / ($sims * $obs)
	g adj_imse_l0_`h' = imse_l0_`h' / ($sims * $obs)
	g adj_isme_mle_`h' = imse_mle_`h' / $sims
	
	replace h_forgraph = `h' / 10 if obs_n == `counter'
	
	replace imse_li_h = adj_imse_l1_`h' if obs_n == `counter'
	replace imse_l0_h = adj_imse_l0_`h' if obs_n == `counter'
	replace aimse_hat_h = adj_isme_mle_`h' if obs_n == `counter'
	
	local counter = `counter' + 1
	
}

twoway (scatter imse_l0_h h_forgraph,xtitle("Bandwidth") ytitle("Predicted IMSE L0"))  
graph export "$out/PS2 Q3c_ISME_L0 STATA.pdf", as (pdf) replace

twoway (scatter imse_li_h h_forgraph,xtitle("Bandwidth") ytitle("Predicted IMSE LI"))  
graph export "$out/PS2 Q3c_IMSE_LI STATA.pdf", as (pdf) replace

twoway (scatter aimse_hat_h h_forgraph,xtitle("Bandwidth") ytitle("Average Predicted AIMSE"))  
graph export "$out/PS2 Q3c_AIMSE_Predicted STATA.pdf", as (pdf) replace

