
# bivariate Archimedean copula based on integrated positive stable LT
# initials are for integrated positive stable Archimedean

# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v,cpar are vectors of the same length,
#  or if only one of the three is a vector and the other two are scalars.
# cpar>0 is the copula parameter 

# copula cdf
# cpar = copula parameter >0
pipsA=function(u,v,cpar)
{ tem1=qgamma(1-u,cpar)
  tem2=qgamma(1-v,cpar)
  tem=(tem1^cpar+tem2^cpar)^(1/cpar)
  1-pgamma(tem,cpar)
}

# copula density
# cpar = copula parameter >0
dipsA=function(u,v,cpar)
{ tem1=qgamma(1-u,cpar)
  tem2=qgamma(1-v,cpar)
  con=gamma(1+cpar)/cpar
  sm=tem1^cpar+tem2^cpar
  tem=sm^(1/cpar)
  pdf=con*tem*exp(-tem+tem1+tem2)/sm
  pdf
}

# conditional copula cdf C_{2|1}(v|u;cpar) 
# cpar = copula parameter >0
# Output: conditional cdf
pcondipsA=function(v,u,cpar)
{ tem1=qgamma(1-u,cpar)
  tem2=qgamma(1-v,cpar)
  tem=(tem1^cpar+tem2^cpar)^(1/cpar)
  exp(-tem+tem1)
}

# conditional copula quantile C_{2|1}^{-1}(p|u;cpar), 
# 0<p<1, could be vector
# cpar = copula parameter >0
# Output: inverse of conditional cdf
qcondipsA=function(p,u,cpar)
{ tem1=qgamma(1-u,cpar)
  y=(tem1-log(p))^cpar - tem1^cpar
  y=y^(1/cpar)
  1-pgamma(y,cpar)
}

# log density
logdipsA=function(u,v,cpar)
{ tem1=qgamma(1-u,cpar)
  tem2=qgamma(1-v,cpar)
  #con=gamma(1+cpar)/cpar
  lcon=lgamma(cpar)
  sm=tem1^cpar+tem2^cpar
  lsm=log(sm)
  #tem=sm^(1/cpar)
  tem=exp(lsm/cpar)
  #pdf=con*tem*exp(-tem+tem1+tem2)/sm
  lpdf=lcon+lsm*(1/cpar-1)+(-tem+tem1+tem2)
  lpdf
}

# random pairs copula quantile C_{2|1}^{-1}(p|u;cpar), cpar>0
# n = sample size
# cpar = copula parameter >0
# Output: random sample of size n
ripsA=function(n,cpar)
{ u=runif(n)
  tem1=qgamma(1-u,cpar)
  p=runif(n)
  y=(tem1-log(p))^cpar - tem1^cpar
  y=y^(1/cpar)
  v=1-pgamma(y,cpar)
  cbind(u,v)
}

# Kendall's tau
# cpar = copula parameter >0
ipsA.cpar2tau=function(cpar)
{ ln2=log(2)
  tem=(2-2*cpar)*ln2+lgamma(2*cpar)-2*lgamma(1+cpar)
  1-cpar*exp(tem)
}

# Kendall's tau to copula parameter
# tau = vector or scalar with values in (-1,1)
# mxiter = maximum number of Newton-Raphson iterations 
# eps = tolerance for convergence of Newton-Raphson iterations 
# cparstart = starting point for Newton-Raphson iterations 
# iprint = print flag for Newton-Raphson iterations 
# Output: copula parameter(s)
ipsA.tau2cpar=function(tau, mxiter=20,eps=1.e-6,cparstart=0,iprint=F)
{ con=log((1-tau)*sqrt(pi)/2)
  cpar=cparstart
  if(cparstart<=0) cpar=tau+1
  iter=0
  diff=1
  while(iter<mxiter & max(abs(diff))>eps)
  { g=con+lgamma(1+cpar)-lgamma(cpar+.5)
    gp=digamma(1+cpar)-digamma(cpar+.5)
    iter=iter+1
    diff=g/gp
    cpar=cpar-diff
    while(min(cpar)<=0.) { diff=diff/2; cpar=cpar+diff }
    if(iprint) cat(iter," ",cpar," ",diff,"\n")
  }
  if(iter>=mxiter) cat("did not converge\n")
  cpar
}

#============================================================

# reflected/survival copula to have skewness to upper tail
# cpar = copula parameter >0

# copula density, 
dipsAr=function(u,v,cpar)
{ tem1=qgamma(u,cpar)
  tem2=qgamma(v,cpar)
  con=gamma(1+cpar)/cpar
  sm=tem1^cpar+tem2^cpar
  tem=sm^(1/cpar)
  pdf=con*tem*exp(-tem+tem1+tem2)/sm
  pdf
}

pipsAr=function(u,v,cpar)
{ u+v-1+pipsA(1-u,1-v,cpar) }

pcondipsAr=function(v,u,cpar)
{ 1-pcondipsA(1-v,1-u,cpar) }

qcondipsAr=function(p,u,cpar)
{ 1-qcondipsA(1-p,1-u,cpar) }

# log density of reflected copula
logdipsAr=function(u,v,cpar)
{ tem1=qgamma(u,cpar)
  tem2=qgamma(v,cpar)
  #con=gamma(1+cpar)/cpar
  lcon=lgamma(cpar)
  sm=tem1^cpar+tem2^cpar
  lsm=log(sm)
  #tem=sm^(1/cpar)
  tem=exp(lsm/cpar)
  #pdf=con*tem*exp(-tem+tem1+tem2)/sm
  lpdf=lcon+lsm*(1/cpar-1)+(-tem+tem1+tem2)
  lpdf
}

