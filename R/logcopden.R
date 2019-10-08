# log of copula densities (for log-likelihoods)
# some 1-parameter and 2-parameter copula families

# cpar = copula parameter which could be a vector
# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v,cpar are vectors of the same length,
#  or if only one of the three is a vector and the other two are scalars.
# The boundary constraints on the functions are not checked on.


# log Plackett copula density, cpar>0
logdpla=function(u,v,cpar)
{ cpar1=cpar-1.;
  tem=1.+cpar1*(u+v); tem1=tem*tem-4.*cpar*cpar1*u*v;
  tem2=sqrt(tem1);
  pdf=cpar*(1.-u-v+2.*u*v+cpar*(u+v-2.*u*v))/(tem1*tem2);
  log(pdf)
}

# log bivariate Frank copula density, cpar>0 or cpar<0 
logdfrk=function(u,v,cpar)
{ t1=1.-exp(-cpar);
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  log(pdf)
}

# log MTCJ copula density, cpar>0
logdmtcj=function(u,v,cpar)
{ tem1=u^(-cpar); tem2=v^(-cpar);
  #pdf=(tem1+tem2-1)^(-1/cpar-2)*(1+cpar)*tem1*tem2/(u*v)
  # rewrite
  lpdf=(-1./cpar-2.)*log(tem1+tem2-1.)+log((1.+cpar)*tem1*tem2/(u*v))
  lpdf
}

# log reflected MTCJ copula density, cpar>0
logdmtcjr=function(u,v,cpar)
{ ubar=1-u; vbar=1-v
  tem1=ubar^(-cpar); tem2=vbar^(-cpar);
  lpdf=(-1./cpar-2.)*log(tem1+tem2-1.)+log((1.+cpar)*tem1*tem2/(ubar*vbar))
  lpdf
}

# log Joe/B5 copula density, cpar>1
logdjoe=function(u,v,cpar)
{ u1=1-u; v1=1-v;
  tem1=u1^cpar; tem2=v1^cpar;
  sm=tem1+tem2-tem1*tem2; 
  ltem=(-2+1/cpar)*log(sm)
  deriv12=tem1*tem2*((cpar-1.)+tem1+tem2-tem1*tem2)/(u1*v1)
  lpdf=ltem+log(deriv12)
  lpdf
}

# log Gumbel copula density, cpar>1
logdgum=function(u,v,cpar)
{ l1= -log(u); l2= -log(v);
  tem1=l1^cpar; tem2=l2^cpar; sm=tem1+tem2; tem=sm^(1./cpar);
  #cdf=exp(-tem);
  lcdf=-tem
  #pdf=cdf*tem*tem1*tem2*(tem+cpar-1.);
  #pdf=pdf/(sm*sm*l1*l2*u*v);
  deriv12=tem*tem1*tem2*(tem+cpar-1.)/(sm*sm*l1*l2)
  lpdf=lcdf+l1+l2+log(deriv12)
  lpdf
}

# log Galambos copula density, cpar>0
logdgal=function(u,v,cpar)
{ cpar1= -cpar;
  l1= -log(u); l2= -log(v);
  tem1=l1^(cpar1-1.); tem2=l2^(cpar1-1.); 
  sm=tem1*l1+tem2*l2; tem=sm^(1./cpar1-1.);
  lcdf=-(l1+l2-tem*sm);  # cdf=exp(lcdf);
  deriv12= 1.+ (tem/sm)*tem1*tem2*( 1.+cpar+tem*sm) -tem*(tem1+tem2);
  lpdf=lcdf+l1+l2+log(deriv12)
  lpdf
}

# log Huesler-Reiss copula density, cpar>0
logdhr=function(u,v,cpar)
{ z1=-log(u); z2=-log(v); 
  z=z1/z2; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;  
  p1=pnorm(tem1); p2=pnorm(tem2);
  #cdf=exp(-z1*p1-z2*p2);
  lcdf= -z1*p1-z2*p2;
  #pdf=cdf*(p1*p2+.5*cpar*dnorm(tem2)/z2)/(u*v); # error
  #deriv12=p1*p2+.5*cpar*dnorm(tem2)/z2 # error
  deriv12=p1*p2+.5*cpar*dnorm(tem1)/z2
  lpdf= lcdf+z1+z2+log(deriv12)
  lpdf
}


#============================================================

# log bivariate normal density, -1<rho<1
logdbvn=function(x1,x2,rho)
{ qf=x1^2+x2^2-2*rho*x1*x2
  qf=qf/(1-rho^2)
  con=sqrt(1-rho^2)*(2*pi)
  #pdf=exp(-.5*qf)/con
  lpdf=-.5*qf-log(con)
  lpdf
}

# log bivariate normal copula density, -1<cpar<1
logdbvncop=function(u,v,cpar)
{ x1=qnorm(u); x2=qnorm(v)
  qf=x1^2+x2^2-2*cpar*x1*x2
  qf=qf/(1-cpar^2)
  con=sqrt(1-cpar^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  #pdf=pdf/(dnorm(x1)*dnorm(x2))
  lpdf=-.5*qf-log(con)-dnorm(x1,log=T)-dnorm(x2,log=T)
  lpdf
}

# log Student t copula density 
# param = 2-vector (rho,df) or scalar rho
# rho = correlation parameter of length(param)=1
# df = degrees of freedom (if length(param)=1, set df=dfdefault before calling)
logdbvtcop=function(u,v,param,df=dfdefault)
{ if (length(param)==2) { rho=param[1]; df=param[2] }
  else rho=param
  lcon=lgamma((df+2)/2) -lgamma(df/2) -log(pi*df*sqrt(1-rho*rho))
  ex=-(df+2)/2
  xt1=qt(u,df)
  xt2=qt(v,df)
  r11=1/(1-rho*rho)
  r12=-rho*r11
  qf=xt1*xt1*r11 + xt2*xt2*r11 + 2*xt1*xt2*r12
  tem=1+qf/df
  lpdf=lcon +ex*log(tem) -dt(xt1,df,log=T) -dt(xt2,df,log=T)
  lpdf
}


#============================================================

# log BB1 density, cpar is 2-vector (th,de) with th>0, de>1
logdbb1=function(u,v,cpar)
{ th=cpar[1]; de=cpar[2] 
  de1=1/de; th1=1/th
  ut=(u^(-th)-1); vt=(v^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  lsm=log(sm); lut=log(ut); lvt=log(vt)
  #tem=(1+smd)^(-th1-2) * (th*(de-1)+(th*de+1)*smd)
  #pdf=tem*smd*x*y*(ut+1)*(vt+1)/sm/sm/ut/vt/u/v
  lpdf= (-th1-2)*log(1+smd)+ log(th*(de-1)+(th*de+1)*smd)
  lpdf= lpdf+(de1-2)*lsm+(de-1)*(lut+lvt) - (th+1)*log(u*v)
  lpdf
}

# log reflected BB1 density, cpar is 2-vector (th,de) with th>0, de>1
logdbb1r=function(u,v,cpar)
{ th=cpar[1]; de=cpar[2] 
  u1=1-u; v1=1-v;
  de1=1/de; th1=1/th
  ut=(u1^(-th)-1); vt=(v1^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  lsm=log(sm); lut=log(ut); lvt=log(vt)
  lpdf= (-th1-2)*log(1+smd)+ log(th*(de-1)+(th*de+1)*smd)
  lpdf= lpdf+(de1-2)*lsm+(de-1)*(lut+lvt) - (th+1)*log(u1*v1)
  lpdf
}

# log BB7 density, cpar is 2-vector (th,de) with th>1, de>0
logdbb7=function(u,v,cpar)
{ th=cpar[1]; de=cpar[2] 
  de1=1/de; th1=1/th
  ut=1.-(1.-u)^th; vt=1.-(1.-v)^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  #tem=(1.-smd)^(th1-2.);
  #tem=tem*(th*(de+1)-(th*de+1)*smd)
  #pdf=tem*smd*(x+1)*(y+1)*(1-ut)*(1-vt)/sm/sm/ut/vt/(1-u)/(1-v);
  lsm=log(sm); lut=log(ut); lvt=log(vt)
  lpdf=(th1-2.)*log(1.-smd)+ log(th*(de+1)-(th*de+1)*smd)
  lpdf= lpdf+ (-de1-2)*lsm+ (-de-1)*(lut+lvt) + (th-1)*log((1.-u)*(1.-v))
  lpdf
}

