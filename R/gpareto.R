# Generalized Pareto MLE

# xdata = data vector of exceedances above a threshold
# maxitn = max number of iterations
# This function tries to automatically come up with a good starting point;
# it has been extension tested for tail parameter in the range -0.5<xi<0.5
# Output: GEV log-likelihood, MLE and asymptotic covariance matrix:
#  (xi=shape/tail parameter, sigma=scale parameter)
gpmle=function(xdata,maxitn=20)
{ m=length(xdata)
  iconv=1
  av=rep(0,4)
  xi=0.1; s=1
  xdata=sort(xdata,decreasing=T)
  est= .C("gpmle",
  	as.double(xdata), as.integer(m), 
  	xi=as.double(xi), s=as.double(s),
  	iconv=as.integer(iconv), as.integer(maxitn), 
  	av=as.double(av))
  if(iconv==0) { warning("did not converge"); return (-1); }
  if(iconv==1) 
  { llk=est$av[1]
    acov=matrix(c(est$av[2], est$av[4], est$av[4], est$av[3]), 2,2)
    pm=matrix(c(est$xi, est$s, sqrt(acov[1,1]), sqrt(acov[2,2])), 2,2)
    dimnames(acov)= list(c("xi", "s"), c("xi", "s"))
    dimnames(pm)= list(c("xi", "s"), c("est", "s.e."))
    list(loglik=llk, params=pm, covar=acov)
  }
}

# For all functions below,
# xi = tail parameter, 
# sigma = scale parameter

# x = vector or scalar 
dgpareto=function(x,xi,sigma=1) 
{ tem=(1+(xi*x)/sigma)^((-1/xi) - 1)
  tem/sigma
}

# x = vector or scalar 
pgpareto=function(x,xi,sigma=1) 
{ tem=(1+(xi*x)/sigma)^(-1/xi)
  1-tem
}

# p = vector or scalar with values in the interval (0,1)
qgpareto=function(p,xi,sigma=1) 
{ (sigma/xi)*((1-p)^(-xi)-1) }


# checks
#xi=.5
#sigma=2
#pdf=dgpareto(1:5,xi,sigma)
#x=1:5; eps=1.e-5
#cdf1=pgpareto(x,xi,sigma)
#cdf2=pgpareto(x+eps,xi,sigma)
#print((cdf2-cdf1)/eps)
#qq=qgpareto(cdf1,xi,sigma)
#print(qq)

