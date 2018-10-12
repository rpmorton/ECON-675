***Created by RM on 2018.09.30
***For ECON 675, PS 2, Q 3
*********************

/************
***Question 3: 4
***********/

clear
set more off

**Generate Data

global n = 5000
global d = 5

set obs $n

g xnorm = 0

forv i = 1(1)$d {
	g x_`i' = runiform(-1,1)
	replace xnorm = xnorm + (x_`i' * x_`i')

}

replace xnorm = sqrt(xnorm)
g nu = rnormal()
g epsilon = .3637899*(1+xnorm^2) * nu
g g0 = exp(xnorm^2)
g u = rnormal()
g t = xnorm + u >= 0
