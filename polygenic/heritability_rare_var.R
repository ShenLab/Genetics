## per SNP h2

# h2 / SNP ~ 2p(1-p)b^2   ### 

## rare variants, aggregated (s>0.01)
## ~4.8 D-mis / ind
##  ~1.22 LoF / ind


mean_from_z_p <- function(sigma, z, p) {
  q <- qnorm(p)
  mu <- z - q * sigma
  return(mu)
}


### 

z = 2 ## threshold for autism 
# N = 10000, then 
## total nvar ~ 6 * N, assume ultra rare, then freq = 1/N

N=30000
p = 1/N

# missense: 
## 5.11 / case, 4.81 / control 
m1 = 5.11
m0 = 4.81

# assuming 0.5 out of 4.8 are risk variants
# then 
x = c(0.3, 0.5, 0.7, 0.9, 1)
h2m = rep(0,length(x))
RRm = rep(0, length(x))
i = 1
for (mr in x) {
	RRm[i] = (m1 - m0) / mr + 1 
	b = mean_from_z_p(1, z, 1 - (1 - pnorm(z,0,1) ) * RRm[i])
	h2m[i] = 2 * b^2 
	i = i + 1

}


## LoF
n1 = 1.38
n0 = 1.22

y = c(0.1, 0.2, 0.3,  0.4, 0.5)
h2l = rep(0, length(y))
RRl = rep(0, length(y))

i = 1
for (mr in y) {
	RRl[i] = (n1 - n0) / mr + 1 
	b = mean_from_z_p(1, z, 1 - (1 - pnorm(z,0,1) ) * RRl[i])
	h2l[i] = 2 * b^2 
	i = i + 1

}


