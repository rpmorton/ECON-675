***Created by RM on 2018.11.07
***For ECON 675, PS 4, Q 3
*********************

clear
set more off

global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"

/* Create Data */

global sims = 1000

set seed 42338992

set obs 1
g totalobs = 50 * ($sims + 1)
local totalobs = totalobs

clear

set obs `totalobs'

g obscounter = [_n]

g simgroup = floor( (obscounter - 1) / 50)
egen maxsimgroup = max(simgroup)
local maxsim = maxsimgroup

g epsilon = rnormal()

g x = rnormal()

g secondnorm = rnormal()
*Cov(ax+by,x) = E[(ax+by)x]-E[ax+by]E[x] = aE[x^2]; V[ax+by] = E[a^2x^2 + 2abxy + b^2y^2] = a^2 + b^2 for x indep y, x, y std normal
*as a = .85, 1 - (.85)^2 = b^2 
g b = sqrt(1 - (.85 * .85) )
g z = .85 * x + b * secondnorm
drop secondnorm

g y = 1 + .5 * x + z + epsilon


*Run Models

g beta_2 = .
g beta_3 = .
g beta_z = .
g se_beta_z = .

forv s = 0(1)`maxsim' {

	qui: reg y x z if simgroup == `s'

	replace beta_2 = _b[x] if simgroup == `s'
	replace beta_z = _b[z] if simgroup == `s'
	replace se_beta_z = _se[z] if simgroup == `s'
	
	qui: reg y x if simgroup == `s'
	
	replace beta_3 = _b[x] if simgroup == `s'
	
}

g beta_select = beta_2 if abs(beta_z / se_beta_z) >= 1.96
replace beta_select = beta_3 if abs(beta_z / se_beta_z) < 1.96

egen tag = tag(simgroup)
keep if tag == 1

/* Now Summary Stats (mean, sd, empirical 95\% CI and Kernel Density Estimators) */

local betalist "select 2 3"

g perc_025 = round(.025 * $sims)
g perc_975 = round(.975 * $sims)


foreach l of local betalist {
	di "betalist is `l'"
	
	quietly su beta_`l', d
	g mean_beta_`l' =  `r(mean)'
	g sd_mean_`l' =  `r(sd)'
	
	sort beta_`l'
	g obs_beta_sorted = [_n]
	g pre_ci_low_beta_`l' = beta_`l' if obs_beta_sorted == perc_025
	egen ci_low_beta_`l' = max(pre_ci_low_beta_`l')
	g pre_ci_high_beta_`l' = beta_`l' if obs_beta_sorted == perc_975
	egen ci_high_beta_`l' = max(pre_ci_high_beta_`l')
	
	if "`l'" == "2" {
		g title = "Density of {&beta} from Model (2)"
		local title = title
		drop title
		}
	if "`l'" == "3" {
		g title = "Density of {&beta} from Model (3)"
		local title = title
		drop title
		}
	if "`l'" == "select" {
		g title = "Density of {&beta} from Model Selection Rule"
		local title = title
		drop title
		}	
	
	twoway (kdensity beta_`l', xtitle("Predicted {&beta}") ytitle("Estimated Density") title("PS4 Q3: `title'" ) ) 
	graph export "$out/STATA_PS4_Q3_`l'.pdf", as (pdf) replace

	drop obs_beta_sorted pre*	
	
	
	}


