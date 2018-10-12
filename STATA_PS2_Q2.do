***Created by RM on 2018.09.30
***For ECON 675, PS 2, Q 2
*********************

/************
***Question 2: 5a and 5b
***********/

clear
set more off

global temp "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data/Temp"
global out "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Outputs"


clear all
/*
mata 
	void foobar(real scalar a) {
				clear
				X = st_data(.,("intercept", "x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10", "x_11", "x_12", "x_13", "x_14", "x_15", "x_16", "x_17", "x_18", "x_19", "x_20"))
				
				wmat = X' * luinv(X' * X) * X'

				st_matrix("wmat", wmat)
	}
end
*/

clear mata
mata:
real myfunction(c)
{
      X = st_data(.,("intercept", "x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10", "x_11", "x_12", "x_13", "x_14", "x_15", "x_16", "x_17", "x_18", "x_19", "x_20"))
				
		wmat = X' * luinv(X' * X) * X'
		
		wmatdiag = diagonal(wmat)

       return(wmatdiag)

}
end



*mata foobar(1)


global sims = 1

forv s = 1(1)$sims {

**Generate Data

	clear all

	global obs = 1000

	set obs $obs
	g x = runiform(-1,1)


	forv i = 1(1)5 {
		g rnorm`i'_sq = (rnormal())^2
	}

	g epsilon = x * x * (rnorm1_sq + rnorm2_sq + rnorm3_sq + rnorm4_sq + rnorm5_sq - 5)

	g y = exp(-.1 * (4*(x^2)-1)^2) * sin(5*x) + epsilon

	g intercept = 1
	
	forv k = 1(1)20 {
		g x_`k' = 0
	}
	
	*pause
	g obs = [_n]
	*X = st_data(.,("intercept", "x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "x_10", "x_11", "x_12", "x_13", "x_14", "x_15", "x_16", "x_17", "x_18", "x_19", "x_20"))
	 mata: matrix X = st_data(.,("intercept"))
	
	forv k = 1(1)20 {
		replace x_`k' = x^`k'
		
		di "creates x_`k'"
		
		mata: matrix add = st_data(.,("x_`k'"))
		matrix list add
		di "does first mata"
		mata: matrix X = (X,add)
		
		matrix list X
		
		reg y x_`k'
		predict resid_`k', resid
		
		capture eret clear
		
		di "right before mata"
		*mata: w = myfunction(1)
		
		pause
		
		g weight_i_`k' = wmat[obs,obs]
		
		*pause
		
		g cv_i_over_n_`k' = (1 / $obs) * (resid_`k'^2 /(1 - weight_i_`k')^2)
	
	}

}	



