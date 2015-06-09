library("VGAM")

# p-value threshold
p = 0.01

# beta distribtion parameters
a = 5
b = 5

total = c(10:200)

## binomial versus 
binomratio = rep(0, length(total))
bbratio = rep(0, length(total))

for (i in 1:length(total)) {
	b1 = qbinom(p, total[i], 0.5) 
	
	b2 = 0
	for (q in 0:total[i]) {
		if (pbetabinom.ab(q, total[i], a, b) >= p ) {
			b2 = q
			break
		}
	}
	binomratio[i] =  (total[i] - b1)  / b1
	bbratio[[i]] =   (total[i] - b2)  / b2
	
	
}


plot(total, bbratio, col='red', type="l", ylim=range(c(0,10)))
lines(total, binomratio )