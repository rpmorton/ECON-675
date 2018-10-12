***Created by RM on 2018.09.30
***For ECON 675, PS 2, Q 2
*********************

/************
***Question 2: 5a
***********/

clear
set more off

**Generate Data

global obs = 1000

set obs $obs
g x = runiform(-1,1)


forv i = 1(1)5 {
	g rnorm`i'_sq = (rnormal())^2
}

g epsilon = x * x * (rnorm1_sq + rnorm2_sq + rnorm3_sq + rnorm4_sq + rnorm5_sq - 5)

g y = exp(-.1 * (4*(x^2)-1)^2) * sin(5*x) + epsilon

