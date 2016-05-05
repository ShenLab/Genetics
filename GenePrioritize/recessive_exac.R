## calculate likelihood ratio of haploinsufficiency, recessive, and neutral

library(dplyr)
recHISModel <- function(exp_lof, n_lof, i ) {
	ratioHIS = 0.09
	ratioRec = 0.5

	sizeN = 10
	sizeHIS = 5
	sizeRec = 5
	
	l0 = dnbinom(n_lof, mu=exp_lof, size = sizeN)
	l1 = dnbinom(n_lof, mu=exp_lof * ratioRec, size = sizeRec)
	l2 = dnbinom(n_lof, mu=exp_lof * ratioHIS, size = sizeHIS)
	LR1 = log10(l1 / l0)
	LR2 = log10(l1 / l2)
	if (i == 1) {
		return(LR1) 
	} else {
		return(LR2)
	}
	
}
# number of observed LGD variants is modeled as a gamma-poisson distribution, with \mu and \theta (negative binomial)

# for neutral genes, \mu = exp_lof (from ExAC calibration), \theta = 10

# for HIS genes, \mu = 0.09 * exp_lof, \theta = 10

# for recessive genes, \mu = 0.5 * exp_lof, \theta = 2  

# ideally, infer \theta jointly with \pi



exac = read.table("exac_r03_march16_z_data_pLI.txt", header=T, sep="\t")


y = mutate(exac, LR1 = recHISModel(exp_lof, n_lof, 1))

y = mutate(y,  LR2= recHISModel(exp_lof, n_lof, 2))

