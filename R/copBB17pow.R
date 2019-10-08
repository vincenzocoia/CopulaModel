# Functions for BB7 and reflected BB1 with an extra power parameter

# BB1 reflected with power parameter eta
# BB7 with power parameter eta
# cpar is 3-dimensional
# 0<u<1, 0<v<1 for all functions
# these could be vectors of same length, or only one is a vector

# BB1rpow copula cdf
# cpar = copula parameter with theta>0, delta>1, eta>0
pbb1rpow=function(u,v,cpar)
{ if(is.matrix(cpar)) { par2=cpar[,1:2]; eta=cpar[,3] }
  else { par2=cpar[1:2]; eta=cpar[3] }
  ua=u^eta
  va=v^eta
  #tem=ua+va-1+pbb1(1-ua,1-va,par2)
  tem=pbb1r(ua,va,par2)
  cdf=tem^(1/eta)
  cdf
}

# conditional cdf
# cpar = copula parameter with theta>0, delta>1, eta>0
pcondbb1rpow=function(v,u,cpar)
{ if(is.matrix(cpar)) { par2=cpar[,1:2]; eta=cpar[,3] }
  else { par2=cpar[1:2]; eta=cpar[3] }
  ua=u^eta
  va=v^eta
  tem=pbb1r(ua,va,par2)
  cdf1=tem^(1/eta-1)
  tem1=pcondbb1r(va,ua,par2)
  ccdf=cdf1*tem1*ua/u
  ccdf
}

# conditional inverse cdf via pcinterpolate
# This seems quite good for p,u in midrange,
# worse when lml>>lmu
# p = scalar or vector in (0,1)
# u = scalar in (0,1)
# cpar = (theta,delta,eta)
# pvec = vector if quantile values in (0,1) for interpolation, 
#    use finer grid for more accuracy
# icheck = T to check pcond compose qcond
# Output: conditional cdf
qcondbb1rpow=function(p,u,cpar,pvec=c(0.01,seq(.02,.98,.02),.99),icheck=F)
{ tem=qcondbb1r(pvec,u,cpar[1:2])
  pp=pcondbb1rpow(tem,u,cpar)
  xx=c(0,pp,1)
  yy=c(0,tem,1)
  der=pcderiv(xx,yy)
  vv=pcinterpolate(xx,yy,der,p)
  if(is.matrix(vv)) { vv=vv[,1] } else { vv=vv[1] }
  if(icheck) { pnew=pcondbb1rpow(vv,u,cpar); print(cbind(p,vv,pnew)) }
  vv
}

# bisection method , can fail in uniroot in extreme cases
# above version is better
# p = scalar or vector in (0,1)
# u = scalar in (0,1)
# cpar = (theta,delta,eta)
#qcondbb1rpow1= function(p,u,cpar,low=0.00001,upp=.99999,icheck=F)
#{ gfn= function(v)
#  { pcondbb1rpow(v,u,cpar)-p }
#  v=uniroot(gfn,lower=low,upper=upp)
#  v=v$root
#  if(icheck) { pnew=pcondbb1rpow(v,u,cpar); cat(p,v,pnew,"\n") }
#  v
#}


# copula density
# cpar = copula parameter with theta>0, delta>1, eta>0
dbb1rpow=function(u,v,cpar)
{ if(is.matrix(cpar)) { par2=cpar[,1:2]; eta=cpar[,3] }
  else { par2=cpar[1:2]; eta=cpar[3] }
  ua=u^eta
  va=v^eta
  temcdf=pbb1r(ua,va,par2)
  cdf2=temcdf^(1/eta-2)
  tem1=pcondbb1r(va,ua,par2)
  tem2=pcondbb1r(ua,va,par2)
  tempdf=dbb1r(ua,va,par2)
  pdf=cdf2*((1-eta)*tem1*tem2+eta*temcdf*tempdf)
  pdf*ua/u*va/v
}

#============================================================

# BB7 with power parameter eta
# copula cdf
# cpar = copula parameter with theta>1, delta>0, eta>0
pbb7pow=function(u,v,cpar)
{ if(is.matrix(cpar)) { par2=cpar[,1:2]; eta=cpar[,3] }
  else { par2=cpar[1:2]; eta=cpar[3] }
  ua=u^eta
  va=v^eta
  tem=pbb7(ua,va,par2)
  cdf=tem^(1/eta)
  cdf
}

# conditional cdf
# cpar = copula parameter with theta>1, delta>0, eta>0
pcondbb7pow=function(v,u,cpar)
{ if(is.matrix(cpar)) { par2=cpar[,1:2]; eta=cpar[,3] }
  else { par2=cpar[1:2]; eta=cpar[3] }
  ua=u^eta
  va=v^eta
  tem=pbb7(ua,va,par2)
  cdf1=tem^(1/eta-1)
  tem1=pcondbb7(va,ua,par2)
  ccdf=cdf1*tem1*ua/u
  ccdf
}

# copula density
# cpar = copula parameter with theta>1, delta>0, eta>0
dbb7pow=function(u,v,cpar)
{ if(is.matrix(cpar)) { par2=cpar[,1:2]; eta=cpar[,3] }
  else { par2=cpar[1:2]; eta=cpar[3] }
  ua=u^eta
  va=v^eta
  temcdf=pbb7(ua,va,par2)
  cdf2=temcdf^(1/eta-2)
  tem1=pcondbb7(va,ua,par2)
  tem2=pcondbb7(ua,va,par2)
  tempdf=dbb7(ua,va,par2)
  pdf=cdf2*((1-eta)*tem1*tem2+eta*temcdf*tempdf)
  pdf*ua/u*va/v
}

#============================================================

