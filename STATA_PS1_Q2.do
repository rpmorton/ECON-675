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

set obs `obs'
g x_1 = runiform() * 100
g x_2 = runiform() * 20 + 30
g intercept = 1
g epsilon = rnormal() * 10

g y = intercept * `beta0' + x_1 * `beta1' + x_2 * `beta2' 

**ia: symmetric point estimate of beta
*bring data into mata

mata: mata clear
capture eret clear



mata:
	
	X = st_data(.,("intercept", "x_1", "x_2"))
	y = st_data(.,("y"))
	Xrows = rows(X)
	Xcols = cols(X)

	function olsbruteforce(real matrix Xmat, real matrix Ymat,| real scalar chol) {
	
		if  (args()==2) {
		
		real colvector betahat	
		betahat = luinv(Xmat' * Xmat) * (Xmat' * Ymat)
		return(betahat)
		
		real colvector epsilon
		epsilon = Ymat - Xmat * betahat
		real scalar s2
		s2 = epsilon' * epsilon *(1/(Xrows-Xcols))
		
		real matrix v0
		v0 = s2 * Xmat' * Xmat
		
		real matrix h0inv
		h0inv = luinv(Xmat'*Xmat)
		
		real colvector AsyVar 
		AsyVar = diagonal(h0inv * v0 * h0inv)
		
		real colvector denominator
		denominator = sqrt(AsyVar * (1/Xrows) )
		return(denominator)
		
		real colvector t_stats
		t_stats = betahat / denominator
		return(t_stats)
		
		
		}
	
	olsbruteforce(X,y)
	
	st_matrix("betahat", betahat)
	st_matrix("denominator", denominator)
	st_matrix("tstat", tstat)
	
	}


end

matrix list betahat


forv i = 1(1)3 {
	
	g beta_`i' = betahat[1,`i']
	
}



/* PS 2 Q 5 */
/*

import delim "$data/LaLonde_1986", delim(",")

g educ2 = educ * educ
g black_earn74 = black * earn74

reg earn78 treat black age educ educ2 earn74 black_earn74 u74 u75
*/

