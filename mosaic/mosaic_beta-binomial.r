library("bbmle")
library(emdbook) 


# MODEL 
# alt fraction is modeled as beta distribution
# rpararmeterize beta as p and theta:
# alpha = theta * p
# beta = theta * (1-p)

# Fit a beta-binomial model for each DP value and estimate p and theta. 
# then regress theta against DP 

# var = N * alpha * beta * (alpha + beta + N) / (alpha + beta) ^2 / (alpha + beta + 1) = N *  p (1-p) * (N + theta) / (theta + 1)


betaBinomEst <- function(counts, total) {

	mtmp <- function(prob,size,theta) { -sum(dbetabinom(counts,prob,size,theta,log=TRUE)) } 
 
	m0 <- mle2(mtmp,start=list(prob=0.5,theta=10),data=list(size=total))
	
	# MLE of theta
	t = coef(m0)["theta"]

}


fitReadCounts <- function(a) {
	## main function
	
	colnames(a)[1:5] = c("type", "chr", "pos", "ref", "alt")
# exploratory plot, showing the total number of reads (DP) versus the fraction of alt reads for each de novo variant
# plot(a$ref + a$alt, a$alt/ (a$alt+a$ref))

	N = a$ref + a$alt 
	x = cbind.data.frame(a[,1:5], N)

	results = estimateTheta(x)	
# look at theta vs DP:
	plot(results$N, results$theta)
## it seems theta is independent on DP
# just take the mean values: 
	

# to confirm fitting:

	results = results[results$m  > 0, ]
	
	# take the median value of estimation from the parts with large enough sample size
	
	mp = median(results[results$N > 12 & results$m > 20, ]$p)
	mtheta = median(results[results$N > 12 & results$m > 20, ]$theta)
	
	plot(results$N, results$var,  xlab = "DP", ylab = "var(Alt_counts)", log="y")

## if binomial
	lines(results$N, results$N * results$p * (1 - results$p), col='blue')
## if beta-binomial, with estimated p and theta
	lines(results$N, results$N * mp * (1 - mp) * (results$N + mtheta) / (mtheta + 1), col='red')

#	return(0)	
	pvalues <- testBetaBinomial(x, mp, mtheta)

#  visualize

	x = cbind.data.frame(x, pvalues)
	z  = x[x$pvalues < 0.005, ]
## plot 

	plot(a$ref + a$alt, a$alt/ (a$alt+a$ref),  xlab = "DP", ylab="alt fraction", pch=21)

# smoothScatter(a$ref + a$alt, a$alt/ (a$alt+a$ref),  xlab = "DP", ylab="alt fraction")

	
	mn = min(max(x$N), 500)


	ci = matrix(0, nrow = mn, ncol=3)
	for (i in 10:mn) {
		lci = 0
		hci = mn 
		d = 0
		for (j in 1:i) {
			d = sum(dbetabinom(0:j, prob = mp, size = i, theta = mtheta  ) )
			if (d > 0.025) {
				lci = j - 1
				break
			} 
		}
	
		for (j in i:1) {
			d = sum(dbetabinom(i:j, prob = mp, size = i, theta = mtheta  ) )
			if (d > 0.025) {
				hci = j + 1
				break
			} 
		}
		ci[i, 1] = lci
		ci[i, 2] = hci
		ci[i, 3] = i
	
	}

	lines(1:mn, ci[,1]/ci[,3], col='red')
	lines(1:mn, ci[,2]/ci[,3], col='red')
	abline(h=mp, col='blue')

# 	text(z$N, z$alt/z$N, label = z$gene, col='blue', cex=0.6, adj=-0.2)

	# write.table(z, "genes.p.lt.0.005.txt", quote=F)
	return(mtheta)

}




estimateTheta <- function(x) {

	mn = min(max(x$N), 500)
	results = matrix(0, nrow=length(levels(factor(x$N))) , ncol=5)
	results = as.data.frame(results)
	colnames(results) = c("N", "m", "p", "var", "theta")
	i = 1
	for (t in 1:mn) {
		v = x[x$N == t, ]$alt
		if (length(v) > 3 ) {
			theta = betaBinomEst(v, t)
	
			results[i,1] = t
			results[i,2] = length(v)
			results[i, 3] = mean(v)/t
			results[i,4] = var(v)
			results[i,5] = theta	
			i = i + 1
		}
	}
	return(results)
}


testBetaBinomial <- function(x, mp = 0.48, mtheta = 60) {
	### now given p and theta, test each variant against a null. The alternative hypothesis is that the variant is a mosaic

	pvalues = rep(0, nrow(x))
	for ( i in 1:nrow(x)) {
		p = sum(  dbetabinom(1:x[i,]$alt,   prob= mp, size = x[i,]$N,  theta = mtheta ) )
		pvalues[i] = p
	}
	return(pvalues)
	
}

### input:

#type	CHR	POS	REF_s1	ALT_s2	...

# only estimate parameters for the first subject
args<-commandArgs(TRUE)
ADfile = args[1]
a = read.table(ADfile, header=T)
a = a[!is.na(a[,4]), ]

tsnv = fitReadCounts(a[a[,1] == "SNV", ])
#tindel = fitReadCounts(a[a[,1] == "INDEL", ])
# print(t, row.names = FALSE)
cat(tsnv, sep="\n")

## todo:
# 1.  more complete candidate de novo list by removing allele fraction filter 
# 2. regress theta against DP and GC, p against DP and GC. 
# 3. confirmation


