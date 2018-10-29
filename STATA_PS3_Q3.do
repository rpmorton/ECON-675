***Created by RM on 2018.10.28
***For ECON 675, PS 3, Q 3
*********************

global output "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"


/************
***Question 3:1
***********/

clear
set more off

set seed 123

global obs = 1000

global bstrap = 599

set obs $obs

g obs = [_n]

g rand = runiform()

egen max_rand = max(rand)
g bstrap_max = .

forv b = 1(1)$bstrap {
	preserve

	bsample $obs
	egen bstrap_max_`b' = max(rand)
	
	local bstrapmax = bstrap_max_`b'
	
	restore
	
	replace bstrap_max = `bstrapmax' if obs == `b'
	
}
	

g bstrap_stat = $obs * (max_rand - bstrap_max)

g exp1 = rexponential(1)

 twoway (kdensity exp1, xtitle("x") ytitle("Estimated Density") ) ///
	(kdensity bstrap_stat, xtitle("x") ytitle("Estimated Density") ///
	title("PS3: Q3 1") legend( label(1 "Simulated Exponential Distribution") label(2 "Nonparametric Bootstrap Statistic") ) )

graph export "$out/PS3_Q3_1_STATA.pdf", as (pdf) replace

/************
***Question 3:2
***********/

clear
set more off

set seed 123


set obs $obs

g obs = [_n]

g rand = runiform()

egen max_rand = max(rand)
g bstrap_max = .

forv b = 1(1)$bstrap {
	preserve

	g bstrapsample = runiform() * max_rand
	egen bstrap_max_`b' = max(bstrapsample)
	
	local bstrapmax = bstrap_max_`b'
	
	restore
	
	replace bstrap_max = `bstrapmax' if obs == `b'
	
}


g bstrap_stat = $obs * (max_rand - bstrap_max)

g exp1 = rexponential(1)

*kdensity exp1, addplot(kdensity bstrap_stat)

 twoway (kdensity exp1, xtitle("x") ytitle("Estimated Density") ) ///
	(kdensity bstrap_stat, xtitle("x") ytitle("Estimated Density") ///
	title("PS3: Q3 2") legend( label(1 "Simulated Exponential Distribution") label(2 "Parametric Bootstrap Statistic") ) )

graph export "$out/PS3_Q3_2_STATA.pdf", as (pdf) replace
