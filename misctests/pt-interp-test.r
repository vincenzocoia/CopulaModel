# pt cdf: table for monotone interpolation

library(CopulaModel)

#df=4
pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
nipol=length(pp) # number of points for interpolation
# nipol=111

# function to check accuracy of monotone interpolation
# df = degree of freedom parameter >0
# seed = seed for set.seed()
# nn = number of (random) checks of the interpolation
# Output within function: the interpolated and "true" values
chkacc=function(df,seed=125,nn=10)
{ cat("\ndf=", df,"\n")
  qq=qt(pp,df)
  pder=pcderiv(qq,pp)  
  #print(cbind(qq,pp,pder))
  set.seed(seed)
  qnew=runif(nn,qq[1],qq[nipol])
  pnew=pcinterpolate(qq,pp,pder,qnew)
  cat("compare interpolation and pt()\n")
  print(cbind(qnew,pnew[,1],pt(qnew,df)))
  invisible(0)
}

chkacc(4)
chkacc(7)
chkacc(10)
chkacc(20)
