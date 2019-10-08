# Functions for cdf of some multivariate copula families

# cdf of multivariate Frank copula
# uu = d-vector of values in (0,1) 
# cpar = copula parameter >0
pmfrk=function(uu,cpar)
{ d=length(uu)
  utem=1-exp(-cpar*uu)
  con=(1-exp(-cpar))^(d-1)
  tem=1-prod(utem)/con
  (-log(tem))/cpar
}

# cdf of multivariate Gumbel copula
# uu = d-vector of values in (0,1) 
# cpar = copula parameter >1
pmgum=function(uu,cpar)
{ d=length(uu)
  uu[uu<=0]=0.000001
  uu[uu>=1]=0.999999
  utem=-log(uu)
  sm=sum(utem^cpar)
  cdf=exp(-sm^(1/cpar))
  cdf
}

#print(pmgum(rep(.4,4),2))
#print(pmgum(rep(.6,4),2))
#print(pmgum(rep(.4,4),4))
#print(pmgum(rep(.4,5),2))
#print(pmgum(rep(.6,5),2))
#print(pmgum(rep(.6,6),2))
#pmgum(rep(0,4),3)
#pmgum(rep(1,4),3)
#pmgum(c(.5,.5,.5,1),3)
#pmgum(c(.5,.5,.5),3)
#pmgum(c(.5,.5,.5,0),3)

#print(pmfrk(rep(.4,4),4))
#print(pmfrk(rep(.6,4),4))
#print(pmfrk(rep(.4,4),5))
#print(pmfrk(rep(.4,5),4))
#print(pmfrk(rep(.6,5),4))
#print(pmfrk(rep(.6,6),4))
#pmfrk(rep(0,4),4)
#pmfrk(rep(1,4),4)
#pmfrk(c(.5,.5,.5,1),4)
#pmfrk(c(.5,.5,.5),4)
#pmfrk(c(.5,.5,.5,0),4)
# OK for boundaries

# cdf of multivariate Galambos copula
# uu = d-vector of values in (0,1) 
# cpar = copula parameter >0
pmgal=function(uu,cpar)
{ d=length(uu)
  uu[uu<=0]=0.000001
  uu[uu>=1]=0.999999
  xtem=-log(uu)
  xth=xtem^(-cpar)
  dd=2^d-1
  sm=sum(xtem)
  for(ii in 1:dd)
  { jj=d2b(d,ii)
    no1=sum(jj); sgn=2*(no1%%2)-1
    if(sum(jj)>1) sm=sm+sgn*((sum(xth[jj==1]))^(-1/cpar))
  }
  cdf=exp(-sm)
  cdf
}

