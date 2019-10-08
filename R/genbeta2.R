# Yang X, Frees EW, Zhang Z (2011). A generalized beta copula with
# applications in modeling multivariate long-tailed data.
# Insurance: Mathematics and Economics 49, 265-284.

#  [Y_j|Q=q)] ~ Gamma(eta_j,rate=q) conditionally independent for j=1,...,d, 
#  Q~ Gamma(ze,1).
#  f_{Y1...Yd}=
#   \int_0^oo \prod_{j=1}^d { y_j^{q-1}e^{-q y_j} q^{eta_j}\over Gamma(eta_j) } 
#          {e^{-q} q^{ze-1}\over Gamma(ze)} dq
#   = { Gamma(eta_+ +ze) \prod_{j=1}^d y_j^{eta_j-1} \over
#           Gamma(eta_1)...Gamma(eta_d) Gamma(ze) (1+y_+)^{eta_+ +ze} }. 

# y1>0, y2>0, param=c(eta1,eta2,zeta), all parameters>0 
# Output: bivariate generalized Beta2 density 
dbgb2=function(y1,y2,param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  lcon=lgamma(eta1+eta2+ze)-lgamma(eta1)-lgamma(eta2)-lgamma(ze)
  ltem=(eta1-1)*log(y1)+(eta2-1)*log(y2)-(eta1+eta2+ze)*log(1+y1+y2)
  exp(lcon+ltem)
}

# yvec = (y1,...,yd) >0 
# param = (eta1,...,etad, zeta), all >0 
# Output: multivariate generalized Beta2 density 
dmgb2=function(yvec,param)
{ d=length(param)-1
  eta=param[1:d]; ze=param[d+1]
  lcon=lgamma(sum(param))-sum(lgamma(param))
  ltem=sum((eta-1)*log(yvec))-sum(param)*log(1+sum(yvec))
  exp(lcon+ltem)
}

# y>0, eta>0, ze>0
# Output: univariate GB2 density
dugb2= function(y,eta,ze)
{ lcon=lgamma(eta+ze)-lgamma(eta)-lgamma(ze)
  ltem=(eta-1)*log(y)-(eta+ze)*log(1+y)
  exp(lcon+ltem)
}

# y>0, eta>0, ze>0
# Output: univariate GB2 cdf
ptbeta= function(y,eta,ze)
{ w=y/(1+y)
  pbeta(w,eta,ze)
}

# y>0, eta>0, ze>0
# Output: univariate GB2 quantile function
qtbeta= function(u,eta,ze)
{ w=qbeta(u,eta,ze)
  w/(1-w)
}

# 0<u,v<1 (could be vectorized)
# param = copula parameter (eta1,eta2,zeta)
# Output: bivariate copula density of generalized Beta2
dbgb2cop=function(u,v,param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  y1=qtbeta(u,eta1,ze)
  y2=qtbeta(v,eta2,ze)
  lcon=lgamma(eta1+eta2+ze)-lgamma(eta1+ze)-lgamma(eta2+ze)+lgamma(ze)
  ltem=(eta1+ze)*log(1+y1)+(eta2+ze)*log(1+y2)-(eta1+eta2+ze)*log(1+y1+y2)
  exp(lcon+ltem)
}

# uvec = (u1,...,ud) in (0,1) 
# Output: multivariate copula density of generalized Beta2
dmgb2cop=function(uvec,param)
{ d=length(param)-1
  eta=param[1:d]; ze=param[d+1]
  yvec=qtbeta(uvec,eta,ze)
  lcon=lgamma(sum(param))-sum(lgamma(eta+ze))+(d-1)*lgamma(ze)
  ltem=sum((eta+ze)*log(1+yvec))-(sum(eta)+ze)*log(1+sum(yvec))
  exp(lcon+ltem)
}

# n = simulation sample size
# param = (eta1,eta2,zeta), all parameters >0
# Output: random bivariate sample with U(0,1) margins
rbgb2cop= function(n,param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  th=rgamma(n,ze)
  y1=rgamma(n,eta1,rate=th)
  y2=rgamma(n,eta2,rate=th)
  u1=ptbeta(y1,eta1,ze)
  u2=ptbeta(y2,eta2,ze)
  cbind(u1,u2)
}

# n = simulation sample size
# param = (eta1,...,etad, zeta), all parameters >0 
# Output: random multivariate sample with U(0,1) margins
rmgb2cop=function(n,param)
{ d=length(param)-1
  eta=param[1:d]; ze=param[d+1]
  th=rgamma(n,ze)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { yy=rgamma(n,eta[j],rate=th) 
    uu[,j]=ptbeta(yy,eta[j],ze)
  }
  uu
}

# param = (eta1,eta2,zeta), all parameters>0 
# Output: upper tail dependence parameter of bivariate generalized Beta2
bgb2.cpar2lm=function(param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  tem1=(beta(eta1,ze))^(1/ze)
  tem2=(beta(eta2,ze))^(1/ze)
  w1=tem1/(tem1+tem2); w2=1-w1
  pbeta(w1,ze+eta2,eta1)+pbeta(w2,ze+eta1,eta2)
}

# uvec = (u1,...,ud) in (0,1) 
# Output: log copula density for one d-vector
logdmgb2cop=function(uvec,param)
{ d=length(param)-1
  eta=param[1:d]; ze=param[d+1]
  yvec=qtbeta(uvec,eta,ze)
  lcon=lgamma(sum(param))-sum(lgamma(eta+ze))+(d-1)*lgamma(ze)
  ltem=sum((eta+ze)*log(1+yvec))-(sum(eta)+ze)*log(1+sum(yvec))
  lcon+ltem
}

# param = (eta1,...,etad, zeta), all parameters >0 
# udat = nxd data matrix with U(0,1) margins
dmgb2copnllk= function(param,udat)
{ n=nrow(udat)
  if(any(param<=0)) return(1.e10)
  nllk=0
  for(i in 1:n) nllk=nllk-logdmgb2cop(udat[i,],param)
  nllk
}

# param = (eta1,...,etad, zeta), all parameters >0  and zeta>2
# Output: correlation matrix of multivariate GB2
mgb2.cpar2cor=function(param)
{ d=length(param)-1
  eta=param[1:d]; ze=param[d+1]
  if(ze<=2) return(NA)
  eq=1/(ze-1); eq2=eq/(ze-2); v=eq2=eq*eq
  r=eq2/v
  tem=eta/(r+eta)
  rmat=outer(tem,tem)
  diag(rmat)=1
  rmat
}

#============================================================

# bivariate cdf via integration, Gauss-Laguerre quadrature should be better
# y1>0, y2>0, param=c(eta1,eta2,zeta), all parameters>0 
# Output: bivariate generalized Beta2 cdf 
pbgb2=function(y1,y2,param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  intgfn= function(r)
  { pgamma(y1*r,eta1)*pgamma(y2*r,eta2)*dgamma(r,ze) }
  cdf=integrate(intgfn,0,Inf)
  cdf$value
}

# 0<u,v<1 
# param = copula parameter (eta1,eta2,zeta)
# Output: bivariate copula cdf of generalized Beta2
pbgb2cop=function(u,v,param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  y1=qtbeta(u,eta1,ze)
  y2=qtbeta(v,eta2,ze)
  pbgb2(y1,y2,param)
}

# 0<v,u<1 
# param = copula parameter (eta1,eta2,zeta)
# Output: copula conditional cdf of bivariate generalized Beta2
pcondbgb2cop=function(v,u,param)
{ eta1=param[1]; eta2=param[2]; ze=param[3]
  y1=qtbeta(u,eta1,ze)
  y2=qtbeta(v,eta2,ze)
  intgfnc= function(r)
  { r*dgamma(y1*r,eta1)*pgamma(y2*r,eta2)*dgamma(r,ze) }
  tem=integrate(intgfnc,0,Inf)
  ccdf=tem$value/dugb2(y1,eta1,ze)
  ccdf
}

