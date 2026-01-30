## investigate polygenic + oligogenic  model of liability

library("truncnorm")

h2=0.8   ## total heritability
pi= 0.8  ## fraction of h2 by polygenic risk

N=20000

firstGen <- function() {


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

gp = rnorm(N, 0, sp)

go = rep(0, N)
co = matrix(0, nrow=N, ncol=no)
for (i in 1:N) {
#  co[i,] <- runif(no) < fo
  nr = rpois(1, no * fo)
 
  go[i] = sum(sample(oligo, nr))  # Poisson sampling
  
}

g = gp + go

#var(oligo)
# var(go)


}