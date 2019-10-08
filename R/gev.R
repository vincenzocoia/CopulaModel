# generalized extreme value MLE

# xdata = data vector of maxima
# maxitn = max number of iterations
# This function tries to automatically come up with a good starting point;
# it has been extension tested for tail parameter in the range -0.5<xi<0.5
# Output: GEV log-likelihood, MLE and asymptotic covariance matrix:
#  (xi=shape/tail parameter, sigma=scale parameter, mu=location parameter)
gevmle=function(xdata, maxitn=20)
{ m=length(xdata)
  iconv=1
  av=rep(0,7)
  xi=0.1
  mu=mean(xdata); s=sqrt(var(xdata))
  est= .C("gevmle",
      as.double(xdata), as.integer(m), 
      xi=as.double(xi), s=as.double(s), mu=as.double(mu),
      iconv=as.integer(iconv), as.integer(maxitn),
      av=as.double(av))
  if(iconv==0) stop("did not converge")
  if(iconv==1)
  { acov=diag(est$av[c(2,3,4)], 3,3)
    acov[1,2]= acov[2,1]= est$av[5]
    acov[1,3]= acov[3,1]= est$av[6]
    acov[2,3]= acov[3,2]= est$av[7]
    pm=c(est$xi,sqrt(acov[1,1]), est$s,sqrt(acov[2,2]), est$mu,sqrt(acov[3,3]))
    pm=matrix(pm, byrow=T, nrow=3)
    dimnames(pm)= list(c("xi", "s", "mu"), c("est", "s.e."))
    dimnames(acov)= list(c("xi", "s", "mu"), c("xi", "s", "mu"))
    list(loglik=est$av[1], params=pm, covar=acov)
  }
}

# revised from library(evir)
# For all functions below,
# xi = tail parameter, 
# mu = location parameter, 
# sigma = scale parameter

# x = vector or scalar 
dgev=function(x, xi=1, mu=0, sigma=1) 
{ tmp=1+(xi*(x-mu))/sigma
  tmpxi=tmp^(-1/xi)
  (as.numeric(tmp > 0) * (tmpxi/tmp) * exp(-tmpxi))/sigma
}

# x = vector or scalar
pgev=function(x, xi=1, mu=0, sigma=1) 
{ tem=pmax(0, 1+(xi*(x-mu))/sigma)
  exp(-tem^(-1/xi)) 
}

# p = vector or scalar with values in the interval (0,1)
qgev=function(p, xi=1, mu=0, sigma=1) 
{ mu+(sigma/xi)*((-log(p))^(-xi)-1) }

# This function assumes tmp>0 below
# x = vector or scalar 
logdgev=function(x, xi=1, mu=0, sigma=1)
{ tmp=1+(xi*(x-mu))/sigma
  lgtmp=log(tmp)
  tmpxi=exp(-lgtmp/xi)
  (-1/xi-1)*lgtmp -log(sigma) - tmpxi
}

# checks
#xi=.5
#mu=1
#sigma=2
#pdf=dgev(1:5,xi,mu,sigma)
#lgpdf=logdgev(1:5,xi,mu,sigma)
#print(cbind(pdf,lgpdf,exp(lgpdf)))
#x=1:5; eps=1.e-5
#cdf1=pgev(x,xi,mu,sigma)
#cdf2=pgev(x+eps,xi,mu,sigma)
#print((cdf2-cdf1)/eps)
#qq=qgev(cdf1,xi,mu,sigma)
#print(qq)
