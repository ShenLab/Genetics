 require(doParallel)
 cl<-makeCluster(2)  # use multiple cores
  registerDoParallel(cl)

 x = 1:10

 j <- {
      foreach(i=1:10, .combine = rbind) %dopar% {
      	 if( i %% 2 == 1) {
	     p = x[i]
 	  } else {
	     p = 0
         }
      }
    } [,1]


print(sum(j))