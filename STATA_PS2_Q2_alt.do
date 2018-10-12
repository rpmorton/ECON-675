

clear

	global sims = 100
	global group = 1000

	set obs 1
	g sims = $sims
	g group = $group
	
	g total_obs = sims * group
	local totalobs = total_obs
	
	clear
	
	set obs `totalobs'
	g x = runiform(-1,1)


	forv i = 1(1)5 {
		g rnorm`i'_sq = (rnormal())^2
	}

	g epsilon = x * x * (rnorm1_sq + rnorm2_sq + rnorm3_sq + rnorm4_sq + rnorm5_sq - 5)

	g y = exp(-.1 * (4*(x^2)-1)^2) * sin(5*x) + epsilon

	*g intercept = 1
	drop rnorm*

	g obs_n = [_n]
	g obs_group = ceil(obs_n / $group)
	g obs_in_group = obs_n - ((obs_group-1) * $group)
	
	forv k = 1(1)20 {
		g x_`k' = x ^ `k'
		
		g cv_i_group_`k' = 0
		
		forv i = 1(1) $group {
			di "i is `i'"
			quietly bys obs_group: reg y x if obs_in_group != `i'
			predict resid_temp, resid
			replace cv_i_group_`k'= cv_i_group_`k' + resid_temp * resid_temp * (1/$obs)
			drop resid_temp
			
		}
		
		egen avg_cv_group_`k' = avg(cv_i_group_`k')
	}
	
			
			
		
		
		

	 
