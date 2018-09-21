*Created by RM on 2018.09.17
*For ECON 675, PS 1, Q 2

global data "/Users/russellmorton/Desktop/Coursework/Fall 2018/Econ 675/Problem Sets/Problem Set Data"

clear
set more off

import delim "$data/LaLonde_1986", delim(",")

g educ2 = educ * educ
g black_earn74 = black * earn74

reg earn78 treat black age educ educ2 earn74 black_earn74 u74 u75



