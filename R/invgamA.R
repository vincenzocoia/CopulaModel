# density and conditional of Archimedean copula based on inverse gamma LT
# abbreviation is invgamA

# cpar = copula parameter = 1/alp to get increase in dependence,
#  where alp is shape parameter of inverse gamma distribution

# LT is psi(s)=2*s^{alp/2} K_alp(2*sqrt(s))/Gamma(alp) , K_alp is besselK(;alp)
# first derivative is  -2*s^{(alp-1)/2} K_(alp-1)(2*sqrt(s))/Gamma(alp)
# Inverse gamma is a special case of the generalized inverse Gaussian
# distribution which has besselK as a normalizing constant;
# see McNeil, Frey and Embrechts (2005), p 497, section A.2.5

# The first two derivatives are:
# psi'(s)=-2*s^{(alp-1)/2} K_(alp-1)(2*sqrt(s))/Gamma(alp)
# psi''(s)=2*s^{(alp-2)/2} K_(alp-2)(2*sqrt(s))/Gamma(alp)

# s = positive value
# alp = inverse gamma parameter,
# galp = gamma(alp), set  in calling routine

# inverse gamma LT 
iglt= function(s,alp,galp)
{ sroot=sqrt(s)
  tem=besselK(2*sroot,alp)
  2*s^(alp/2) * tem / galp
}

# first derivative of inverse gamma LT 
iglt1= function(s,alp,galp)
{ sroot=sqrt(s)
  tem=besselK(2*sroot,alp-1)
  -2*s^((alp-1)/2) * tem / galp
}

# second derivative of inverse gamma LT 
iglt2= function(s,alp,galp)
{ sroot=sqrt(s)
  tem=besselK(2*sroot,alp-2)
  2*s^((alp-2)/2) * tem / galp
}

# inverse LT via Newton-Raphson, 
# tt = value in (0,1)
# mxiter = max number of iterations
# iprint = print flag for iteration
# sstart = starting point for iterations
# tol = tolerance for convergence
igltinv= function(tt,alp,galp,mxiter=20,iprint=F,sstart=0,tol=1.e-8)
{ alp1=alp-1
  iter=0
  diff=1.
  s=sstart
  if(sstart<=0) s=1/tt^0.3-1
  while(max(abs(diff))>tol & iter<mxiter)
  { g=iglt(s,alp,galp)-tt
    sroot=sqrt(s)
    tem=besselK(2*sroot,alp1)
    gp= -2*s^(alp1/2) * tem / galp
    diff=g/gp
    s=s-diff
    while(min(s)<=0) { diff=diff/2; s=s+diff }
    iter=iter+1
    if(iprint) print(c(iter,s))
  }
  if(iter>=mxiter & iprint) cat("*** did not converge, tt=",tt," alp=",alp,"\n")
  s
}

# Archimedean copula based on inverse gamma LT
# 0<u<1, 0<v<1, one or both of these can be vectors
# cpar = copula parameter >0,  alp=1/cpar
# Output: copula cdf
pinvgamA=function(u,v,cpar)
{ alp=1/cpar
  galp=gamma(alp)
  tem1=igltinv(u,alp,galp)
  tem2=igltinv(v,alp,galp)
  iglt(tem1+tem2,alp,galp)
}

# conditional cdf C_{2|1}(v|u)
pcondinvgamA=function(v,u,cpar)
{ alp=1/cpar
  galp=gamma(alp)
  tem1=igltinv(u,alp,galp)
  tem2=igltinv(v,alp,galp)
  iglt1(tem1+tem2,alp,galp)/iglt1(tem1,alp,galp)
}

# copula density c_{12}(u,v)
dinvgamA=function(u,v,cpar)
{ alp=1/cpar
  galp=gamma(alp)
  tem1=igltinv(u,alp,galp)
  tem2=igltinv(v,alp,galp)
  iglt2(tem1+tem2,alp,galp)/(iglt1(tem1,alp,galp)*iglt1(tem2,alp,galp))
}


# copula parameter to Kendall's tau, 
# tau=cpar/(cpar+2) has been proved based on stochastic representation
# the function below provides a numerical check
invgamA.cpar2tau0= function(cpar)
{ alp=1/cpar
  galp=gamma(alp) 
  psider= function(s)
  { tem=iglt1(s,alp,galp)
    tem*tem*s
  }
  tem=integrate(psider,0,Inf, rel.tol=1e-06)
  tau=1-4*tem$value
  tau
}

invgamA.cpar2tau=function(cpar) { cpar/(cpar+2) }
invgamA.tau2cpar=function(tau) { 2*tau/(1-tau) }


# inverse of LT is obtained via monotone interpolation
# u = vector of values in (0,1) 
# v = vector of values in (0,1); same length as u
# cpar = copula parameter >0 
# Output: log copula density at each (u,v) pair
logdinvgamA=function(u,v,cpar,pgrid=0)
{ alp=1/cpar
  galp=gamma(alp)
  if(pgrid==0) 
  { p=seq(.01,.99,.01)
    p=c(.0001,.0002,.0005,.001,.002,.005,p)
  }
  else { p=pgrid }
  fn=igltinv(p,alp,galp)  # fast 
  p1=c(p,1)
  fn1=c(fn,0)
  der=pcderiv(p1,fn1)
  #tem1=igltinv(u,alp,galp)
  #tem2=igltinv(v,alp,galp)
  tem1=pcinterpolate(p1,fn1,der,u)[,1]
  tem2=pcinterpolate(p1,fn1,der,v)[,1]
  dens=iglt2(tem1+tem2,alp,galp)/(iglt1(tem1,alp,galp)*iglt1(tem2,alp,galp))
  log(dens)
}

# n = simulation sample size
# d = dimension
# cpar = copula parameter >0 
# Output: nxd matrix with random d-vectors from the d-dimensional Arch copula
rminvgamA=function(n,d,cpar) 
{ alp=1/cpar; galp=gamma(alp)
  r=1/rgamma(n,alp)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r
    uu[,j]=iglt(tem,alp,galp)
  }
  uu
}

