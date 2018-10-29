***Created by RM on 2018.10.15
***For ECON 675, PS 2, Q 3
*********************

/************
***Question 3: 4a
***********/

*Create a giant data set of all obs, all replications

clear
set more off

global obs = 500
global sims = 10

global theta0 = 1
global d = 5

set obs 1
g obs = $obs
g sims = $sims

g obs_sims = obs * sims

global obssims = obs_sims

clear

set obs $obssims

g x_norm_squared = 0
forv l = 1(1)$d {
	g x`l' = runiform()*2 - 1
	replace x_norm = x_norm + x`l'^2
}

g v = rnormal()

g epsilon = .3637899 * (1 + x_norm) * v

g g0 = exp(x_norm_squared)

g u = rnormal()
g t = sqrt(x_norm_squared) + u >= 0

g y = t * $theta0 + g0 + epsilon

g obs_n = [_n]
g group = ceil(obs_n / $obs)
g obs_in_group = obs_n - ( (group - 1) * $obs )


*Generate basis expansion without then with interactions

forv i = 1(1)$d {
	forv j = 1(1)10 {
		g p_`j'_x`i' = x`i' ^ `j'
	}
}

forv i = 1(1)$d {
	forv j = 1(1)$d {
		forv o1 = 1(1)4 {
			forv o2 = 1(1)4 {
				if `i' != `j' {
					g i_`o1'_`o2'_x`i'_x`j' = (x`i')^(`o1') * (x`j')^(`o2')
				}
			}
		}
	}
}

drop x*

/************
***Question 3: 4b
***********/

forv spec = 1(1)14 {
	g spec_`spec'_theta_hat_k = 0
	*g spec_`spec'_bias_theta_hat_k = 0
	*g spec_`spec'_var_theta_hat_k = 0
	g spec_`spec'_vhco_theta_hat_k = 0
	g spec_`spec'_ci_lb_theta_hat_k = 0
	g spec_`spec'_ci_ub_theta_hat_k = 0
}


/* for i through v compute in loop for each group */

forv g = 1(1)$sims {
	quietly reg y t p_1_* if group == `g'
	replace spec_1_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* if group == `g'
	replace spec_2_theta_hat_k = _b[t] if group == `g'

	quietly reg y t p_1_* p_2* i_1* if group == `g'
	replace spec_3_theta_hat_k = _b[t] if group == `g'

	quietly reg y t p_1_* p_2* i_1_1* p_3* if group == `g'
	replace spec_4_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* i_1* p_3* i_2_1* i_1_2* i_2_2* if group == `g'
	replace spec_5_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* i_1* p_3* i_2_1* i_1_2* i_2_2* p_4* if group == `g'
	replace spec_6_theta_hat_k = _b[t] if group == `g'

	quietly reg y t p_1_* p_2* i_1* p_3* i_2_1* i_1_2* i_2_2* p_4* i_3_1* i_3_2* i_1_3* i_2_3* i_3_3* if group == `g'
	replace spec_7_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* i_1* p_3* i_2_1* i_1_2* i_2_2* p_4* i_3_1* i_3_2* i_1_3* i_2_3* i_3_3* p_5* if group == `g'
	replace spec_8_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* p_3* p_4* p_5* i_* if group == `g'
	replace spec_9_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* p_3* p_4* p_5*  p_6* i_* if group == `g'
	replace spec_10_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* p_3* p_4* p_5*  p_6* p_7* i_* if group == `g'
	replace spec_11_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* p_3* p_4* p_5* p_6* p_7* p_8* i_* if group == `g'
	replace spec_12_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_1_* p_2* p_3* p_4* p_5*  p_6* p_7* p_8* p_9* i_* if group == `g'
	replace spec_13_theta_hat_k = _b[t] if group == `g'
	
	quietly reg y t p_* i_* if group == `g'
	replace spec_14_theta_hat_k = _b[t] if group == `g'
}

/* ii: compute bias */

forv spec = 1(1)14{
	g spec_`spec'_bias_theta_hat_k = ( $theta0 - spec_14_theta_hat_k)^2 
	
	
}

		

/** Now average over simulations */

egen tag_group = tag(group)
keep if tag_group == 1

forv spec = 1(1)14{
	egen avg_spec_`spec'_theta_hat_k = mean(spec_`spec'_theta_hat_k) 
	egen avg_spec_`spec'_bias_theta_hat_k = mean(spec_`spec'_bias_theta_hat_k )
	egen var_`spec'_theta_hat_k = sum((1 / ($sims - 1) ) * (spec_`spec'_theta_hat_k -avg_spec_`spec'_theta_hat_k  )^2 )
	
	
}


