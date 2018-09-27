*Created by RM on 2018.09.17
*For ECON 675, PS 1, Q 2

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"

clear
set more off

/* PS 2 Q4: Matrix Implementation of WLS with W = I */

***Generate Data


local obs = 1000
local beta0 = 0
local beta1 = 1
local beta2 = 5
local alpha = .9

set obs `obs'
g x_1 = runiform() * 100
g x_2 = runiform() * 20 + 30
g intercept = 1
g epsilon = rnormal() * 10

g y = x_1 * `beta1' + x_2 * `beta2' + epsilon

**4a: symmetric point estimate of beta
*bring data into mata and estimate

mata: mata clear
capture eret clear

mata:
	
	X = st_data(.,("intercept", "x_1", "x_2"))
	y = st_data(.,("y"))
	Xrows = rows(X)
	Xcols = cols(X)

	betahat = luinv(X' * X) * (X' * y)
		
	epsilonhat = y - X * betahat
	s2 = epsilonhat' * epsilonhat :* (1/(Xrows-Xcols))
	v0_hat = s2 * X' * X :* (1/Xrows)
	
	h0inv = luinv(X'* X :* (1/Xrows))
		
	AsyVar = h0inv * v0_hat * h0inv
	denominator = sqrt(diagonal(AsyVar) * (1/Xrows) )
	t_stats = betahat :/ denominator
	
	st_matrix("betahat", betahat)
	st_matrix("denominator", denominator)
	st_matrix("tstats", t_stats)
	st_matrix("s2",s2)
	st_matrix("v0",v0_hat)
	st_matrix("h0inv",h0inv)
	st_matrix("asyvar",AsyVar)
	st_matrix("countcols",Xcols)
	st_matrix("countrows",Xrows)

end

**4a: now use choleksy inverse


** Now compute the pvalues and confidence intervals
g cols = countcols[1,1]
local cols = cols

g ones = 1
egen obs = sum(ones)
local obs = obs

g df = `obs' - `cols'

forv i = 1(1)`cols' {

	local j = `i' - 1
	g beta_hat_`j' = betahat[`i',1]
	g t_stat_beta_`j' = tstats[`i',1]
	g pval_beta_`j' =  2*ttail(df,t_stat_beta_`j')
	g lb_beta_`j' = beta_hat_`j' + invttail(df,`alpha'+(1-`alpha')/2)*denominator[`i',1]
	g ub_beta_`j' =  beta_hat_`j' + invttail(df,(1-`alpha')/2)*denominator[`i',1]

}


**4a: now use choleksy inverse

mata: mata clear
capture eret clear

mata:
	
	X = st_data(.,("intercept", "x_1", "x_2"))
	y = st_data(.,("y"))
	Xrows = rows(X)
	Xcols = cols(X)

	betahat_chol = cholinv(X' * X) * (X' * y)
		
	epsilonhat_chol = y - X * betahat_chol
	s2_chol = epsilonhat_chol' * epsilonhat_chol :* (1/(Xrows-Xcols))
	v0_hat_chol = s2_chol * X' * X :* (1/Xrows)
	
	h0inv_chol = cholinv(X'* X :* (1/Xrows))
		
	AsyVar_chol = h0inv_chol * v0_hat_chol * h0inv_chol
	denominator_chol = sqrt(diagonal(AsyVar_chol) * (1/Xrows) )
	t_stats_chol = betahat_chol :/ denominator_chol
	
	st_matrix("betahatchol", betahat_chol)
	st_matrix("denominatorchol", denominator_chol)
	st_matrix("tstatschol", t_stats_chol)
	st_matrix("s2chol",s2_chol)
	st_matrix("v0chol",v0_hat_chol)
	st_matrix("h0invchol",h0inv_chol)
	st_matrix("asyvarchol",AsyVar_chol)

end

**4a: now use choleksy inverse


** Now compute the pvalues and confidence intervals

forv i = 1(1)`cols' {

	local j = `i' - 1
	g beta_hat_chol_`j' = betahatchol[`i',1]
	g t_stat_beta_chol_`j' = tstatschol[`i',1]
	g pval_beta_chol_`j' =  2*ttail(df,t_stat_beta_chol_`j')
	g lb_beta_chol_`j' = beta_hat_chol_`j' + invttail(df,`alpha'+(1-`alpha')/2)*denominatorchol[`i',1]
	g ub_beta_chol_`j' =  beta_hat_chol_`j' + invttail(df,(1-`alpha')/2)*denominatorchol[`i',1]

}

**Compare Symmetric and Cholesky Inverse

local compare = "lb_beta ub_beta"

di "cols is `cols'"
	
foreach cibound of local compare {	
	
	forv i = 1(1)`cols' {
	local k = `i' - 1
	
		g diff_`cibound'_`k' = round(`cibound'_`k',.00001) - round(`cibound'_chol_`k',.00001)
		su diff_`cibound'_`k'
		
	}

}


/**************
/* PS 2 Q 5 */
**************/

clear

import delim "$data/LaLonde_1986", delim(",")

local alpha = .95

g educ2 = educ * educ
g black_earn74 = black * earn74
g intercept = 1

mata: mata clear
capture eret clear

mata:
	
	X = st_data(.,("intercept", "treat", "black", "age", "educ", "educ2", "earn74", "black_earn74", "u74", "u75"))
	y = st_data(.,("earn78"))
	Xrows = rows(X)
	Xcols = cols(X)

	betahat = luinv(X' * X) * (X' * y)
		
	epsilonhat = y - X * betahat
	s2 = epsilonhat' * epsilonhat :* (1/(Xrows-Xcols))
	v0_hat = s2 * X' * X :* (1/Xrows)
	
	h0inv = luinv(X'* X :* (1/Xrows))
		
	AsyVar = h0inv * v0_hat * h0inv
	denominator = sqrt(diagonal(AsyVar) * (1/Xrows) )
	t_stats = betahat :/ denominator
	
	st_matrix("betahat", betahat)
	st_matrix("denominator", denominator)
	st_matrix("tstats", t_stats)
	st_matrix("s2",s2)
	st_matrix("v0",v0_hat)
	st_matrix("h0inv",h0inv)
	st_matrix("asyvar",AsyVar)
	st_matrix("countcols",Xcols)
	st_matrix("countrows",Xrows)

end

g cols = countcols[1,1]
local cols = cols

g ones = 1
egen obs = sum(ones)
local obs = obs

g df = `obs' - `cols'

forv i = 1(1)`cols' {

	local j = `i' - 1
	g beta_hat_`j' = betahat[`i',1]
	g t_stat_beta_`j' = tstats[`i',1]
	g pval_beta_`j' =  2*ttail(df,abs(t_stat_beta_`j'))
	g lb_beta_`j' = beta_hat_`j' + invttail(df,`alpha'+(1-`alpha')/2)*denominator[`i',1]
	g ub_beta_`j' =  beta_hat_`j' + invttail(df,(1-`alpha')/2)*denominator[`i',1]

}

reg earn78 treat black age educ educ2 earn74 black_earn74 u74 u75

local colsminus = `cols' - 1

forv i = 1(1)`colsminus' {
	local j = `i' 
	
	local betahat = beta_hat_`j'
	local pval = pval_beta_`j'
	local lb = lb_beta_`j'
	local ub = ub_beta_`j'
	
	di "beta is `betahat'; pval is `pval', lb is `lb', ub is `ub'"
}
	
	local betahat = beta_hat_0
	local pval = pval_beta_0
	local lb = lb_beta_0
	local ub = ub_beta_0
	di "beta is `betahat'; pval is `pval', lb is `lb', ub is `ub'"
