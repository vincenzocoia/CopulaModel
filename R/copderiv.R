# copula derivatives for Gumbel, Frank, t(nu), BB1 for gradient vector
#  of vine loglik; cpar=copula parameter which could be a vector. 
# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v,cpar are vectors of the same length,
#  or if only one of the three is a vector and the other two are scalars.


# partial derivatives of logdfrk with respect to u, v and cpar
# cpar = copula parameter >0 or <0
logdfrk.deriv=function(u,v,cpar)
{ t1=1.-exp(-cpar); t0=1-t1
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  den=tem1+tem2-t0-tem1*tem2;
  den1=-u*tem1-v*tem2+t0+(u+v)*tem1*tem2;
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  lpdf=log(pdf)
  lpdfdpar = 1/cpar+t0/t1-(u+v)-2*den1/den;
  den1u=-cpar*tem1*(1-tem2);
  den1v=-cpar*tem2*(1-tem1);
  lpdfdu=-cpar-2*den1u/den
  lpdfdv=-cpar-2*den1v/den
  if(length(u)==1 & length(v)==1) { return(c(lpdf,lpdfdu,lpdfdv,lpdfdpar)) }
  else { return(cbind(lpdf,lpdfdu,lpdfdv,lpdfdpar)) }
}

# partial derivatives of pcondfrk  with respect to u, v and cpar
# deriv wrt v is c(u,v)
# deriv wrt u is d c(u,v)/du = pdf * d logc(u,v)/du
# cpar = copula parameter >0 or <0
pcondfrk.deriv=function(v,u,cpar)
{ t1=1.-exp(-cpar); t0=1-t1
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  dencond=t1/(1-tem2)-(1-tem1)
  ccdf=tem1/dencond;
  den=tem1+tem2-t0-tem1*tem2;
  den1=-u*tem1-v*tem2+t0+(u+v)*tem1*tem2;
  lcder1=-u+v*tem2/(1-tem2)-den1/den;
  ccdfdpar=ccdf*lcder1
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  ccdfdv=pdf
  ccdfdu=-cpar*ccdf+cpar*ccdf^2
  if(length(v)==1 & length(u)==1) { return(c(ccdf,ccdfdu,ccdfdv,ccdfdpar)) }
  else { return(cbind(ccdf,ccdfdu,ccdfdv,ccdfdpar)) }
}

# partial derivatives of logdgum with respect to u, v and cpar
# cpar = copula parameter >1 
logdgum.deriv=function(u,v,cpar)
{ u[u<=0]=1.e-7; v[v<=0]=1.e-7
  x=-log(u); y=-log(v);
  tx=log(x); ty=log(y)
  xd=x^cpar; yd=y^cpar; sm=xd+yd; tem=sm^(1./cpar);
  lcdf=-tem
  deriv12=tem*xd*yd*(tem+cpar-1.)/(sm*sm*x*y)
  lpdf=lcdf+x+y+log(deriv12)
  logs=log(sm)
  dlsq=cpar*cpar; 
  sder1=xd*tx+yd*ty;
  mder1=tem*sder1/(sm*cpar)-tem*logs/dlsq;
  den=tem+cpar-1.; 
  logm=log(tem); 
  lpdfdpar= -mder1+(mder1+1)/den-2*logm+(1-2*cpar)*mder1/tem+tx+ty;
  mu= -tem*xd/(u*sm*x)
  mv= -tem*yd/(v*sm*y)
  lpdfdu= -mu+mu/den+(1-2*cpar)*mu/tem-(cpar-1)/(u*x)-1/u;
  lpdfdv= -mv+mv/den+(1-2*cpar)*mv/tem-(cpar-1)/(v*y)-1/v;
  if(length(u)==1 & length(v)==1) { return(c(lpdf,lpdfdu,lpdfdv,lpdfdpar)) }
  else { return(cbind(lpdf,lpdfdu,lpdfdv,lpdfdpar)) }
}

# from gumtrvine.f90
# partial derivatives of pcondgum with respect to u, v and cpar
# deriv wrt v is c(u,v)
# deriv wrt u is d c(u,v)/du = pdf * d logc(u,v)/du
# cpar = copula parameter >1 
pcondgum.deriv=function(v,u,cpar)
{ u[u<=0]=1.e-7  # temporary fix (to be consistent with pgum)
  v[v<=0]=1.e-7
  x=-log(u); y=-log(v);
  tx=log(x); ty=log(y)
  xd=x^cpar; yd=y^cpar; sm=xd+yd; tem=sm^(1./cpar);
  deriv12=tem*xd*yd*(tem+cpar-1.)/(sm*sm*x*y)
  lcdf=-tem
  lpdf=lcdf+x+y+log(deriv12)
  logs=log(sm)
  dlsq=cpar*cpar; dlcu=dlsq*cpar
  sder1=xd*tx+yd*ty;
  mder1=tem*sder1/(sm*cpar)-tem*logs/dlsq;
  logm=log(tem); msq=tem*tem; usq=u*u;
  lccdf=x-tem+(1-cpar)*(logm-tx);      
  mu=-tem*xd/(u*sm*x);
  lcder1= -mder1+(1-cpar)*mder1/tem-logm+tx;  
  lcder1u= -(cpar-1)*mu/tem-mu-(cpar-1)/(x*u)-1/u
  ccdf=exp(lccdf)
  ccdfdu=ccdf*lcder1u
  ccdfdpar=ccdf*lcder1
  ccdfdv=exp(lpdf)
  if(length(v)==1 & length(u)==1) { return(c(ccdf,ccdfdu,ccdfdv,ccdfdpar)) }
  else { return(cbind(ccdf,ccdfdu,ccdfdv,ccdfdpar)) }
}

# partial derivatives of logdtcop with respect to u,v,rho, with df=dfdefault fixed
# param = 2-vector (rho,df) or scalar
# df = dfdefault defined externally if param is scalar
logdbvtcop.deriv=function(u,v,param,df=dfdefault)
{ if (length(param)==2) { rho=param[1]; df=param[2] }
  else rho=param
  con=0.5*log(pi*df);
  lgdif=lgamma(0.5*(df+1.))-lgamma(0.5*df);
  coef=1.-rho*rho; coef2=coef*coef; 
  t1=qt(u,df); t2=qt(v,df)
  t1sq=t1*t1; t2sq=t2*t2;
  den=1.+(t1sq-2.*rho*t1*t2+t2sq)/(df*coef);
  den1=2.*rho*(den-1.)/coef-2.*t1*t2/(df*coef);
  dt1=lgdif-con-0.5*(df+1.)*log(1.+t1sq/df);
  dt2=lgdif-con-0.5*(df+1.)*log(1.+t2sq/df);
  dd1=exp(dt1); dd2=exp(dt2) 
  reg12=t1-rho*t2; reg21=t2-rho*t1;
  coef3=den*coef*df
  lpdf= -log(2.*pi)-0.5*log(coef)-0.5*(df+2.)*log(den)-dt1-dt2;
  lpdfdpar=rho/coef-0.5*(df+2.)*den1/den;
  ltder1t1= -(df+2.)*reg12/coef3 + (df+1.)*t1/(df+t1sq);
  lpdfdu=ltder1t1/dd1;
  ltder1t2= -(df+2.)*reg21/coef3 + (df+1.)*t2/(df+t2sq);
  lpdfdv=ltder1t2/dd2;
  if(length(u)==1 & length(v)==1) { return(c(lpdf,lpdfdu,lpdfdv,lpdfdpar)) }
  else { return(cbind(lpdf,lpdfdu,lpdfdv,lpdfdpar)) }
}

# partial derivatives of pcondbvtcop with respect to u,v,rho, with df=dfdefault fixed
# deriv wrt v is c(u,v)
# deriv wrt u is d c(u,v)/du = pdf * d logc(u,v)/du
# param = 2-vector (rho,df) or scalar
# df = dfdefault defined externally if param is scalar
pcondbvtcop.deriv=function(v,u,param,df=dfdefault)
{ if (length(param)==2) { rho=param[1]; df=param[2] }
  else rho=param
  logpi=log(pi)
  rho2=rho*rho;
  t1=qt(u,df); t2=qt(v,df)
  r2=sqrt((df+t1*t1)/(df+1.));
  cr=sqrt(1.-rho2)
  xt21=(t2-rho*t1)/cr/r2;
  ccdf=pt(xt21,df+1.); 
  tem1=(t2*rho-t1)/cr^3/r2;
  tem2=xt21*xt21/(df+1.);
  const= exp(lgamma(0.5*df+1.)-lgamma(0.5*df+0.5)-0.5*logpi-0.5*log(df+1.));
  ccdfdpar= const*(1.+tem2)^(-0.5*df-1.) *tem1;
  t1sq=t1*t1; t2sq=t2*t2;
  coef=1.-rho*rho; 
  den=1.+(t1sq-2.*rho*t1*t2+t2sq)/(df*coef);
  lgdif=lgamma(0.5*(df+1.))-lgamma(0.5*df);
  con=0.5*log(pi*df);
  dt1= lgdif-con-0.5*(df+1.)*log(1.+t1sq/df);
  dt2= lgdif-con-0.5*(df+1.)*log(1.+t2sq/df);
  dd1=exp(dt1)
  lpdf= -log(2.*pi)-0.5*log(coef)-0.5*(df+2.)*log(den)-dt1-dt2;
  ccdfdu=dt(xt21,df+1.)/dd1*(-rho/r2/cr -xt21*t1/(df+t1*t1))
  ccdfdv=exp(lpdf)
  if(length(v)==1 & length(u)==1) { return(c(ccdf,ccdfdu,ccdfdv,ccdfdpar)) }
  else { return(cbind(ccdf,ccdfdu,ccdfdv,ccdfdpar)) }
}

# function used by logdbb1.deriv and pcondbb1.deriv
# u1 and u2 could be vectors of the same length
# theta>0, delta>1 are 2 parameters of BB1
bb1mderivs= function(u1,u2,theta,delta)
{ t1=u1^(-theta); t2=u2^(-theta)
  t1a=t1-1.; t2a=t2-1.
  tu1=-log(u1); tu2=-log(u2)
  ttu1=log(t1a); ttu2=log(t2a)
  td01=(t1a)^delta; td02=(t2a)^delta
  td11=td01/t1a; td12=td02/t2a;
  td21=td11/t1a; td22=td12/t2a;
  s=td01+td02
  sder1th= delta*(td11*t1*tu1+td12*t2*tu2)
  sder1dl= td01*ttu1+td02*ttu2
  m=s^(1./delta); m1=m/s
  ts=log(s)
  dlsq=delta*delta; dlcu=delta*dlsq
  mder1th=m1*sder1th/delta
  mder1dl=m1*sder1dl/delta - m*ts/dlsq
  m1der1dl= mder1dl/s - m*sder1dl/s^2
  list(sm=s,mexp=m,mder1th=mder1th,mder1dl=mder1dl)
}

# partial derivatives of logdbb1, with respect to u, v and cpar, 
# cpar = (theta,delta) with theta>0, delta>1
logdbb1.deriv=function(u,v,cpar)
{ theta=cpar[1]; delta=cpar[2]
  mder=bb1mderivs(u,v,theta,delta)
  m=mder$mexp; mder1th=mder$mder1th; mder1dl=mder$mder1dl
  sm=mder$sm
  mp1=1.+m; msq=m*m; thsq=theta*theta
  thtem=2.+1./theta
  logmp1=log(mp1)
  t10=-(thtem)*logmp1
  t1der1th=logmp1/thsq - thtem*mder1th/mp1
  t1der1dl=-thtem*mder1dl/mp1

  dltem=1.-2.*delta; dl1=delta-1.
  logm=log(m)
  t20=dltem*logm
  t2der1th=dltem*mder1th/m
  t2der1dl=-2.*logm+dltem*mder1dl/m

  coef=theta*delta+1.
  den=theta*dl1+coef*m
  den2=den*den
  t30=log(den)
  t3der1th=(dl1+delta*m+coef*mder1th)/den
  t3der1dl=(theta+theta*m+coef*mder1dl)/den

  t1=u^(-theta); t2=v^(-theta)
  t1a=t1-1.; t2a=t2-1.; smlog=log(t1a)+log(t2a)
  tu1=-log(u); tu2=-log(v)
  t40=dl1*smlog+(theta+1.)*(tu1+tu2)
  t4der1th=dl1*(t1*tu1/t1a+t2*tu2/t2a)+tu1+tu2  
  t4der1dl=smlog

  lpdf=t10+t20+t30+t40
  lpdfdth=t1der1th+t2der1th+t3der1th+t4der1th
  lpdfddl=t1der1dl+t2der1dl+t3der1dl+t4der1dl

  dmdx=m/sm/delta; x=t1a^delta
  dmdy=dmdx; y=t2a^delta
  dxdu=-delta*theta*x/t1a*t1/u
  dydv=-delta*theta*y/t2a*t2/v
  dlpdf= (-thtem/mp1+dltem/m+coef/den)
  lpdfdu= dlpdf*dmdx*dxdu + dl1/delta/x*dxdu - (theta+1)/u
  lpdfdv= dlpdf*dmdy*dydv + dl1/delta/y*dydv - (theta+1)/v

  if(length(u)==1 & length(v)==1) { return(c(lpdf,lpdfdu,lpdfdv,lpdfdth,lpdfddl)) }
  else { return(cbind(lpdf,lpdfdu,lpdfdv,lpdfdth,lpdfddl)) }
}

# partial derivatives of pcondbb1, with respect to u, v and cpar, 
# cpar = (theta,delta) with theta>0, delta>1
pcondbb1.deriv=function(v,u,cpar)
{ theta=cpar[1]; delta=cpar[2]
  mder=bb1mderivs(u,v,theta,delta)
  m=mder$mexp; mder1th=mder$mder1th; mder1dl=mder$mder1dl
  sm=mder$sm
  mp1=1.+m; msq=m*m; thsq=theta*theta
  lmp1=log(mp1); logm=log(m);
  t2=u^(-theta); tu2=-log(u); # t2 is same as  t1 in previous function
  t2a=t2-1.; lt2a=log(t2a)
  cf1=1.+1./theta;
  dl1n=1.-delta

  lcdf= -cf1*lmp1+(dl1n)*(logm-lt2a)+(theta+1.)*tu2;
  lcder1th= lmp1/thsq-cf1*mder1th/mp1;
  lcder1th= lcder1th+(dl1n)*(mder1th/m-t2*tu2/t2a)+tu2;
  lcder1dl= -cf1*mder1dl/mp1+(dl1n)*mder1dl/m-logm+lt2a;
  dmdx=m/sm/delta
  x=t2a^delta
  dxdu=-delta*theta*x/t2a*t2/u
  lcder1u= -cf1/mp1*dmdx*dxdu +dl1n/m*dmdx*dxdu - dl1n/delta/x*dxdu -(theta+1)/u
  ccdf=exp(lcdf);
  ccdfdth=ccdf*lcder1th;
  ccdfddl=ccdf*lcder1dl;
  ccdfdu=ccdf*lcder1u
  # pdf
  de1=1/delta; th1=1/theta
  ut=t2-1; vt=(v^(-theta) -1)
  x=ut^delta; y=vt^delta; sm=x+y
  smd=sm^(de1)
  tem= (1+smd)^(-th1-2) * (theta*(delta-1) + (theta*delta+1)*smd)
  pdf= tem*smd*x*y * (ut+1) * (vt+1)/sm/sm/ut/vt/u/v
  ccdfdv=pdf
  if(length(v)==1 & length(u)==1) { return(c(ccdf,ccdfdu,ccdfdv,ccdfdth,ccdfddl)) }
  else { return(cbind(ccdf,ccdfdu,ccdfdv,ccdfdth,ccdfddl)) }
}
