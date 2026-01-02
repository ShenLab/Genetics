	
	
	h2 = 0.8
	z = 2
	s_max = 0.4
	
	k = 2
	
	
	### Results: 
	## 1.  mean genetic liability of siblings
	## 2.  mean genetic liability of parents
	## 3.  proportion of affected parents
	## 4.  mean genetic liability of unaffected siblings
	## 5.  mean genetic liability of unaffected parents
	## 6.  mean genetic liability of probands
	## 7.  mean genetic liability of probands with affected parents

	
	results = rep(0, 10)
	
	# population sample size for simulation
	N = 2000000
	  
	# parents from random population
	p1g = rnorm(N, 0, h2^0.5) 
	p2g = rnorm(N, 0, h2^0.5)
	p1e = rnorm(N, 0, (1-h2)^0.5)
	p2e = rnorm(N, 0, (1-h2)^0.5)
	
	p1 = p1g + p1e
	p2 = p2g + p2e
	
	
	
## fitness: 
## fitness impact of the condition, modeled as a sigmoid function:  f = 1 - s

## s  = s_max / (1 + exp( -1 * k( t - z) ))    #  t is the trait value  (liability)
 
 ##  t = 0 --> s = 0
 ## t -> inf, --> s = s_max
 ## z is the liability threshold for diagnosis. 
 ## k should be large to minimize 
  ## make it more modest for now. Power et al estimated s_max ~ 0.5 for female autism and 0.75 for males.
  ## but the autism cohort in the study might be much more severe than current cohorts
  
## default values: 
# s_max = 0.4
# k = 2
	
	p1fitness = 1 - s_max / (1 + exp(-1 * k * (p1 -z)))	
	p2fitness = 1 - s_max / (1 + exp(-1 * k * (p2 -z)))	
	
	# mean(p1)
	# var(p1)
	
## sig quantify the randomness of transmission of risk. estimate them to fit h2. 
	 ## because the variance of genetic load of the expected mid parent [(p1g+p2g)/2] is h2/2, the extra variance of genetic load from the randomness of transmission is ~ h2 - h2/2 = h2/2. 
	sig = (h2 / 2) ^ 0.5    
		
	## genetic load of child 1
	c1g = rnorm(N, p1g/2 + p2g/2, sig )    ## this model has problems with oligogenic component (ie. ultra-rare with large effect, which will much more discrete and cannot be modeled well by normal)
	
	## genetic load of child 2
	c2g = rnorm(N, p1g/2 + p2g/2, sig)

	## env load of children
	c1e = rnorm(N, 0, (1-h2)^0.5)
	c2e = rnorm(N, 0, (1-h2)^0.5)
	
	## total liability of child 1 and 2
	c1 = c1g + c1e
	c2 = c2g + c2e
	
	## take average of fitness of parents as the probablity of making the child in a Bernoulli: 
#	c1_chance = rbinom(N, 1, (p1fitness * p2fitness) ^ 0.5 )
#	c2_chance = rbinom(N, 1, (p1fitness * p2fitness) ^ 0.5 )
    c1_chance = rbinom(N, 1, (p1fitness + p2fitness) /2 )
	c2_chance = rbinom(N, 1, (p1fitness + p2fitness) /2 )
		
		
## reproduction events occured: 		
	c1_real = c1[as.logical(c1_chance)]   
	c2_real = c2[as.logical(c1_chance)]      ## this is iffy, as the second child should also be affected by fitness. But the additional impact is small, and model it will generate NA's in the data. For now, do not consider fitness impact of a second child in ascertained families. 
	p1_real = p1[as.logical(c1_chance)]
	p2_real = p2[as.logical(c1_chance)]
	
	c1g_real = c1g[as.logical(c1_chance)]
	c2g_real = c2g[as.logical(c1_chance)]
	p1g_real = p1g[as.logical(c1_chance)]
	p2g_real = p2g[as.logical(c1_chance)]

	
	## probands (cases): ascertain child 1's with z  > 2
	i <- c1_real > z 
	
	## total liability of ascertained families: 
	c1s = c1_real[i]
	p1s = p1_real[i]
	p2s = p2_real[i]
 	c2s = c2_real[i]  
	
	c1sg = c1g_real[i]
	p1sg = p1g_real[i]
	p2sg = p2g_real[i]
	c2sg = c2g_real[i]
	
	Nc1 = length(c1_real)
	Ncase = length(c1s)
	
		
	# mean genetic liability of siblings: 
#	results[1] = mean(c2sg)
	# cov(c1s, p1s + p2s)
	
	## prevalence in c1
	results[1] = Ncase / Nc1 
	
	## mean genetic liability of parents
	results[2]= mean(c(p1sg, p2sg))
	
	# var(p1s)
	
	## affected parents
	p1a <- p1s > z
	p2a <- p2s > z
	
	## unaffected parents
	p1u <- p1s < z
	p2u <- p2s < z
	
	## proportion of cases with >=1 affected parents
	Raffected = (length(c1s[p1s > z | p2s > z]) )/ Ncase 
	results[3] = Raffected
	
	## mean genetic liability of unaffected siblings: 
	results[4] = mean(c2sg[c2s < z])
	
	## mean genetic liability of unaffected parents:
	results[5] = mean(c(p1sg[p1s< z], p2sg[p2s < z])  )
	
	
	# mean genetic liability of probands: 
	results[6] = mean(c1sg)

	# mean genetic liability of probands with affected parents: 
	results[7] = mean(c1sg[p1s > z | p2s > z])

	# mean genetic liability of probands with both parents unaffected: 
	results[8] = mean(c1sg[p1s < z & p2s < z])

	# mean genetic liability of probands with affected sibs: 
	results[9] = mean(c1sg[c2s > z])
	
	# mean genetic liability of probands with no affected parent/sib (i.e. simplex): 
	results[10] = mean(c1sg[c2s < z & p1s < z & p2s < z])

simplex=c1sg[c2s < z & p1s < z & p2s < z]
multiplex = c1sg[c2s > z |  p1s > z | p2s > z]

par(family="sans")
hist(p1sg, freq=F, br=40, col='white', border='white', xlab="Genetic liability", main="", ylim=c(0,0.9))
grid()
abline(v=2, lty="dashed")
abline(v=0)

lines(density(p1g), lwd=2)
# lines(density(p1sg), col='salmon', lty="dashed")
unaffPar = c(p1sg[p1s< z], p2sg[p2s < z]) 
lines(density(unaffPar), col='blue', lwd=2)

# lines(density(c1sg), col="red", lty="dashed", lwd=1)
lines(density(simplex), col="orange", lwd=2)
lines(density(multiplex), col="red", lwd=2)
legend("topleft", c("General population", "Unaffected parents", "Simplex probands", "Multiplex probands"), col=c("black", "blue", "orange", "red"), lwd=c(2,2,2,2))
