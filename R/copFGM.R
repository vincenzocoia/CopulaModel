# Farlie-Gumbel-Morgenstern, Morgenstern or FGM
# FGM copula: can used for some tests of copula functions
# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v are vectors of the same length
#   or if one of these is a vector and the other is a scalar

# FGM copula cdf
# cpar = copula parameter with -1<=cpar<=1
pfgm=function(u,v,cpar)
{ tem=1+cpar*(1-u)*(1-v)
  cdf=u*v*tem
  cdf
}

# FGM copula density
dfgm=function(u,v,cpar)
{ pdf=1+cpar*(1-2*u)*(1-2*v)
  pdf
}

# log copula density
logdfgm=function(u,v,cpar)
{ pdf=1+cpar*(1-2*u)*(1-2*v)
  log(pdf)
}

#============================================================

# cpar = copula parameter with -1<cpar<1
pcondfgm=function(v,u,cpar)
{ tem=1+cpar*(1-2*u)*(1-v)
  v*tem
}


# cpar = copula parameter with -1<cpar<1
#   cpar=0 input will not work for this function
# Output : conditional quantile
qcondfgm=function(p,u,cpar)
{ #if(u==.5) return(p)
  b=-cpar*(1-2*u)-1
  a=cpar*(1-2*u)
  discr=b*b-4*p*a
  #ifelse(u==0.5 | cpar==0, p, (-b-sqrt(discr))/(2*a))
  out=(-b-sqrt(discr))/(2*a)
  out[u==0.5]=p
  out
} 

