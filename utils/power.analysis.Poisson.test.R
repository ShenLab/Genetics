## analyze power of poisson test
# main question:  how power improves by including more cases that are related but not necessary the same phenotype
## assumption:  the effect size is 1/3 or 1/2 (or a range) in added cases (e.g. external cases) than ascertained cases


ptpower <- function(pcut  = 2e-6 ) {
	## pcut: type I error threshold
	
    # N <- c(1000, 2000, 3000)
    N = 2000
    N_ext <- c(0:10) * 1000
    
    effectDecay <- c(1/16, 1/8, 1/4, 1/3, 1/2)

    m = 3  ## 3 mutations in current data set
    rate = 9.45e-06   # median rate of LOF/Dmis mutation among HHE & RVIS genes (1428 of them)



    power = matrix(ncol = length(effectDecay), nrow = length(N_ext))
    for (i in 1:length(effectDecay)) {
    	decay = effectDecay[i]
	for (j in 1:length(N_ext)) {
	    n_external = N_ext[j]
	    cutoff = qpois(1 - pcut /2, (N + n_external) * 2 * rate ) 
	    power[j, i] = 1- ppois(cutoff + 1 , m + n_external * m / N * decay)	    
	}
    }
    return(power)
}
