# Genetic risks in families
#
#Questions:
#	1.	given heritability (h2) and fitness impact (S) of a binary condition, eg autism, what is the chance of parents being affected in randomly ascertained families with affected children?  i.e., Prob(affected parents | child is affected, h2, S)?
#	2.	Given the prevalence of the condition, and h2, S, what is the distribution of genetic load in unaffected parents of affected children or affected parents?

## how to simulate assortative mating? 




firstDegree  <- function(h2 = 0.8, z = 2, s_max_m = 0.99, s_max_f = 0.66, k = 2) {

	### Results: 
	## 1.  mean genetic liability of siblings
	## 2.  mean genetic liability of parents
	## 3.  proportion of affected parents
	## 4.  mean genetic liability of unaffected siblings
	## 5.  mean genetic liability of unaffected parents
	## 6.  mean genetic liability of probands
	## 7.  mean genetic liability of probands with affected parents

	
	results = rep(0, 16)
	
	
	# population sample size for simulation
	N = 5000000
	  
	# parents from random population, p1: fathers, p2: mothers
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
	
	p1fitness = 1 - s_max_m / (1 + exp(-1 * k * (p1 -z)))	
	p2fitness = 1 - s_max_f / (1 + exp(-1 * k * (p2 -z)))	
	
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

	## proportion of cases with >=1 affected parents
	Raffected = (length(c1s[p1s > z | p2s > z]) )/ Ncase 
	results[2] = Raffected

	# proportion of cases with affected fathers
	results[3] = (length(c1s[p1s > z]) )/ Ncase  	
	
	# proportion of cases with affected mothers	
	results[4] = (length(c1s[p2s > z]) )/ Ncase  


	## mean genetic liability of probands: 
	results[5] = mean(c1sg)
	
	## mean genetic liability of fathers
	results[6]= mean(p1sg)
	
	## mean genetic liability of mothers
	results[7]= mean(p2sg)	
	
	
	## mean genetic liability of affected fathers
	results[8]= mean(p1sg[p1s > z])
	
	## mean genetic liability of affected mothers
	results[9]= mean(p2sg[p2s > z])	


	# var(p1s)
	
	## affected parents
	p1a <- p1s > z
	p2a <- p2s > z
	
	## unaffected parents
	p1u <- p1s < z
	p2u <- p2s < z
	
		
	## mean genetic liability of unaffected siblings: 
	results[10] = mean(c2sg[c2s < z])
	
	## mean genetic liability of unaffected fathers:
	results[11] = mean(p1sg[p1s< z] )
	
	## mean genetic liability of unaffected mothers:	
	results[12] = mean(p2sg[p2s < z])
	
	# mean genetic liability of probands with affected parents: 
	results[13] = mean(c1sg[p1s > z | p2s > z])

	# mean genetic liability of probands with both parents unaffected: 
	results[14] = mean(c1sg[p1s < z & p2s < z])


simplex=c1sg[c2s < z & p1s < z & p2s < z]
multiplex = c1sg[c2s > z |  p1s > z | p2s > z]

	# mean genetic liability of probands with affected sibs or parents (multiplex): 
	results[15] = mean(multiplex)
	
	# mean genetic liability of probands with no affected parent/sib (i.e. simplex): 
	results[16] = mean(simplex)

	return(results)
} 



### heritability
h2v = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

## liability threshold (at z-score scale)
z = 2  ## --> 2.275% prevalence 

d = matrix(0, nrow = length(h2v), ncol = 16)

colnames(d) = c("prevalence", "frac_aff_parents", "frac_aff_father", "frac_aff_mother", "mu_proband", "mu_father", "mu_mother",  "mu_aff_father", "mu_aff_mother", "mu_unaff_sib", "mu_unaff_father",  "mu_unaff_mother",  "mu_case_w_aff_p", "mu_case_p_unaff", "mu_multiplex", "mu_simplex")


i = 1

sib_simplex_ratio = rep(0, length(h2v))
simplex_multi_ratio = rep(0, length(h2v))

for (h2 in h2v) {
	d[i,] = firstDegree(h2,z)
	
	sib_simplex_ratio[i] = d[i, 10]/ d[i, 16]
	simplex_multi_ratio[i] = d[i, 16]/ d[i, 15]
	i  = i + 1
	
}

par(family="sans")

plot(h2v, sib_simplex_ratio, xlab="h2", ylab="ratio of liability", xlim=range(c(0.3,1)), ylim=range(c(0,1)), type="l")

lines(h2v, simplex_multi_ratio, col='red')

grid()

# legend("topleft", c("General population", "Unaffected parents", "Simplex probands", "Multiplex probands"), col=c("black", "blue", "orange", "red"), lwd=c(2,2,2,2))

  
plot(h2v, d[,2], xlab="heritability (h^2)", ylab="fraction of probands with affected parents", xlim=range(c(0.3,0.9)), ylim=range(c(0,0.3)), type="b")
grid()



# xvec = c(5, 10, 20)   # number of rare risk sites in each individual
# varc1 = rep(0, length(xvec))

# children
# c1 =  (p1 + p2)  * (h2/2)^ 0.5 + rnorm(N, 0, (1-h2) ^ 0.5) 
# i = 1
# for (x in xvec) {
# c1g = p1g * rbinom(N, x, 0.5)/x +  p2g *rbinom(N, x, 0.5)/x

# c1e = rnorm(N, 0, (1-h2)^0.5)
# c2e = rnorm(N, 0, (1-h2)^0.5)

# c1 = c1e + c1g 

# c1 =  (p1g + p2g) / 2^0.5 + rnorm(N, 0, (1-h2) ^ 0.5) 
# c2 =  (p1g + p2g) / 2^0.5 + rnorm(N, 0, (1-h2) ^ 0.5) 

# c1a =  (p1 + p2) / 3 + rnorm(N, 0, (1-h2) ^ 0.5) 

