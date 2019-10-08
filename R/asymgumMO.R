# bivariate asymmetric Gumbel copulas with Marshall-Olkin in limit
# 0<u<1, 0<v<1 for all functions
# u,v are vectors of same length, or at least one is a scalar
# cpar = copula parameter = (de,pi1,pi2)
# with de>1,  pi1 in [0,1] and pi2 in [0,1] 

# 2-parameter bivariate MO copula cdf
# Cuadras-Auge copula is special case pi1=pi2
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (pi1,pi2), each in [0,1]
# Output: cdf
pbMO=function(u,v,cpar)
{ pi1=cpar[1]; pi2=cpar[2]
  temu=u^pi1
  temv=v^pi2
  cdf=pmin(temu,temv)*u*v/temu/temv
  cdf
}

# C_{2|1} for bivariate MO 
# v = scalar or vector of values in (0,1)
# u = scalar or vector of values in (0,1)
# cpar = copula parameter (pi1,pi2), each in [0,1]
# Output: conditional cdf
pcondbMO21=function(v,u,cpar)
{ if(length(u)!=length(v)) return(NA)
  pi1=cpar[1]; pi2=cpar[2]
  temu=u^pi1
  temv=v^pi2
  ii=(temu<=temv)
  nn=length(ii)
  ccdf=rep(0,nn)
  ccdf[ii]=v[ii]/temv[ii]
  ccdf[!ii]=(1-pi1)*v[!ii]/temu[!ii]
  ccdf
}

# asymmetric Gumbel with Marshall-Olkin at boundary
# exponent function A of the extreme value copula
# x = scalar or vector of positive values 
# y = scalar or vector of positive values 
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
# Output: A and its first/second order derivatives at x,y values
AasymgumMO=function(x,y,cpar)
{ de=cpar[1]; pi1=cpar[2]; pi2=cpar[3]
  temx=(pi1*x)^de; temy=(pi2*y)^de
  sm=temx+temy
  tem=sm^(1/de)
  Afn=tem+(1-pi1)*x+(1-pi2)*y
  Aderx=tem/sm*temx/x+1-pi1
  Adery=tem/sm*temy/y+1-pi2
  Aderxy=(tem/sm/sm)*(temx/x)*(temy/y)*(1-de)
  list(Afn=Afn,Aderx=Aderx,Adery=Adery,Aderxy=Aderxy)
}

# copula cdf for asymmetric Gumbel with Marshall-Olkin at boundary
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
# Output: cdf
pasymgumMO=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  Aval=AasymgumMO(x,y,cpar)$Afn
  cdf=exp(-Aval)
  cdf
}

# C_{2|1} for asymmetric Gumbel with Marshall-Olkin at boundary
# v = scalar or vector of values in (0,1)
# u = scalar or vector of values in (0,1)
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
# Output: conditional cdf
pcondasymgumMO21=function(v,u,cpar)
{ x=-log(u)
  y=-log(v)
  Aout=AasymgumMO(x,y,cpar)
  Aval=Aout$Afn
  cdf=exp(-Aval)
  ccdf=cdf*Aout$Aderx/u
  ccdf
}

# C_{1|2} for asymmetric Gumbel with Marshall-Olkin at boundary
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
# Output: conditional cdf
pcondasymgumMO12=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  Aout=AasymgumMO(x,y,cpar)
  Aval=Aout$Afn
  cdf=exp(-Aval)
  ccdf=cdf*Aout$Adery/v
  ccdf
}

# copula pdf for asymmetric Gumbel with Marshall-Olkin at boundary
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
# Output: copula density
dasymgumMO=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  Aout=AasymgumMO(x,y,cpar)
  Aval=Aout$Afn
  cdf=exp(-Aval)
  pdf=cdf*(Aout$Aderx*Aout$Adery-Aout$Aderxy)/u/v
  pdf
}

# Kendall's tau for asymmetric Gumbel
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
asymgumMO.cpar2tau=function(cpar)
{ de=cpar[1]; pi1=cpar[2]; pi2=cpar[3]
  tauasym= function(w)
  { w1=1-w
    wd=(pi1*w)^de
    w1d=(pi2*w1)^de
    sm=wd+w1d
    tem=sm^(1/de)
    B=tem+(1-pi1)*w+(1-pi2)*w1
    Bp=tem/sm*(wd/w-w1d/w1) + pi2-pi1
    ((2*w-1)*Bp*B + w*w1*Bp*Bp)/(B*B)
  }
  tem=integrate(tauasym,0,1, rel.tol=1e-06)
  tau=tem$value
  tau
}

# Spearman's rho for asymmetric Gumbel
# cpar = (de=delta, pi1,pi2) with de>1 and pi1,pi2 in (0,1]
asymgumMO.cpar2rhoS=function(cpar)
{ de=cpar[1]; pi1=cpar[2]; pi2=cpar[3]
  spasym= function(w)
  { w1=1-w
    wd=(pi1*w)^de
    w1d=(pi2*w1)^de
    sm=wd+w1d
    tem=sm^(1/de)
    B=tem+(1-pi1)*w+(1-pi2)*w1
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(spasym,0,1, rel.tol=1e-06)
  rho=12*tem$value-3
  rho
}
