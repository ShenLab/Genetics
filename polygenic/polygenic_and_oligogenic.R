## investigate polygenic + oligogenic  model of liability

library("truncnorm")

h2=0.8   ## total heritability
pi= 0.8  ## fraction of h2 by polygenic risk

N=20000

## firstGen <- function() {


no = 8000  ## choice of rare risk loci (binned by genes)
fo = 1e-4

# polygenic risk, modeled in aggregation by normal(0,sp)
sp=(h2*pi)^0.5  

so = 0.5  ## rare variants have bigger chance of large effect



## model of polygenic small effect
# polyg = rnorm(np, 0, sp)

## model of large effect  using  truncated normal
### !!! Poisson sampling of truncated normal is not quite stable.  Need  some  modification
### to ensure target variance 

oligo = rtruncnorm(no, a=-0.1, b=2, mean=0, sd=so)

#oligo = rnorm(no, 0, so)

## polygenic risk for each individual 
gp = rnorm(N, 0, sp)

## total oligogenic risk for each indivdual 
go = rep(0, N)

## individual risk rare variants for each individual 
oligoInd <- vector("list", N)

for (i in 1:N) {
#  co[i,] <- runif(no) < fo
  nr = rpois(1, no * fo)
  oligoInd[[i]] = sample(oligo, nr)  # Poisson sampling
  go[i] = sum(oligoInd[[i]]) 
  
}

g = gp + go

var(oligo)
var(go)


## }