# Functions for copula cdfs and pdfs, 1-parameter bivariate families and t
# bvn = bivariate normal/Gaussian
# pla = Plackett
# frk = Frank
# mtcj = Mardia-Takahasi-Clayton-Cook-Johnson
# joe = Joe/B5
# gum = Gumbel
# gal = Galambos
# hr = Huelser-Reiss
# bvt = bivariate t

# copula cdfs, cpar=copula parameter which could be a vector
# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v,cpar are vectors of the same length,
#  or if only one of the three is a vector and the other two are scalars.
# The boundary constraints on the functions are not checked on.

# Plackett copula cdf
# cpar = copula parameter >0 
#   This function does not have the case of cpar=1 for independence
ppla=function(u,v,cpar)
{ cpar1=cpar-1.;
  tem=1.+cpar1*(u+v); tem1=tem*tem-4.*cpar*cpar1*u*v;
  tem2=sqrt(tem1);
  cdf=(tem-tem2)*.5/cpar1;
  #ifelse(cpar==1., u*v, cdf) # doesn't treat vectors as expected
  cdf
}

# Frank copula cdf
# cpar = copula parameter cpar>0 or cpar<0
#   This function does not have the case of cpar=0 for independence
pfrk=function(u,v,cpar)
{ cpar1=1.-exp(-cpar);
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  tem=cpar1-(1.-tem1)*(1.-tem2);
  cdf=(-log(tem/cpar1)/cpar);
  cdf
}

# MTCJ copula cdf
# cpar = copula parameter >0
pmtcj=function(u,v,cpar)
{ tem1=u^(-cpar); tem2=v^(-cpar);
  cdf=(tem1+tem2-1.)^(-1./cpar);
  cdf
}

# reflected/survival MTCJ copula cdf
# cpar = copula parameter >0
pmtcjr=function(u,v,cpar)
{ tem1=(1-u)^(-cpar); tem2=(1-v)^(-cpar);
  cdf=(tem1+tem2-1.)^(-1./cpar);
  u+v-1+cdf
}

# Joe/B5 copula cdf
# cpar = copula parameter >1
pjoe=function(u,v,cpar)
{ tem1=(1-u)^cpar; tem2=(1-v)^cpar;
  sm=tem1+tem2-tem1*tem2; tem=sm^(1./cpar);
  cdf=1-tem
  cdf
}

# Gumbel copula cdf
# cpar = copula parameter >1
pgum=function(u,v,cpar)
{ # below might not be best fix for boundary
  u[u<=0]=0.0000001; v[v<=0]=0.0000001
  # intermediate calculations can cause input to be 1+1.e-?? with roundoff
  u[u>=1]=1; v[v>=1]=1
  l1= -log(u); l2= -log(v);
  tem1=(l1^cpar); tem2=(l2^cpar); sm=tem1+tem2; tem=sm^(1./cpar);
  cdf=exp(-tem);
  cdf
}

# reflected/survival Gumbel copula cdf
# cpar = copula parameter >1
pgumr=function(u,v,cpar)
{ u[u>=1]=1-0.0000001; v[v>=1]=1-0.0000001
  u[u<=0]=0; v[v<=0]=0
  l1= -log(1-u); l2= -log(1-v);
  tem1=(l1^cpar); tem2=(l2^cpar); sm=tem1+tem2; tem=sm^(1./cpar);
  cdf=exp(-tem);
  cdf+u+v-1
}

# Galambos copula cdf
# cpar = copula parameter >0
pgal=function(u,v,cpar)
{ t1= -cpar;
  l1= -log(u); l2= -log(v);
  tem1=l1^t1; tem2=l2^t1; sm=tem1+tem2; tem=sm^(1./t1);
  cdf=exp(-(l1+l2-tem));
  cdf
}

# Huesler-Reiss copula cdf
# cpar = copula parameter >0
phr=function(u,v,cpar)
{ l1= -log(u); l2= -log(v);
  z=l1/l2; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;  # increasing in concordance
  p1=pnorm(tem1); p2=pnorm(tem2);
  lcdf=-l1*p1-l2*p2;  cdf=exp(lcdf);
  cdf
}

#============================================================
# copula density functions

# Plackett copula density
# cpar = copula parameter >0
dpla=function(u,v,cpar)
{ #if(cpar==1.) return(1.)
  cpar1=cpar-1.;
  tem=1.+cpar1*(u+v); tem1=tem*tem-4.*cpar*cpar1*u*v;
  tem2=sqrt(tem1);
  #pdf=cpar*(1.-u-v+2.*u*v+cpar*(u+v-2.*u*v))/(tem1*tem2);
  tem0=2.*cpar1*u*v; tem3=tem-tem0;
  pdf=cpar*tem3/tem1/tem2;
  #ifelse(cpar==1., 1,pdf)
  pdf
}

# bivariate Frank copula density
# cpar = copula parameter >0 or <0 (limit case of 0 is independence)
dfrk=function(u,v,cpar)
{ t1=1.-exp(-cpar);
  tem1=exp(-cpar*u); tem2=exp(-cpar*v);
  pdf=cpar*tem1*tem2*t1;
  tem=t1-(1.-tem1)*(1.-tem2);
  pdf=pdf/(tem*tem);
  pdf
}

# MTCJ copula density
# cpar = copula parameter >0 
dmtcj=function(u,v,cpar)
{ tem1=u^(-cpar); tem2=v^(-cpar);
  pdf=(tem1+tem2-1)^(-1/cpar-2)*(1+cpar)*tem1*tem2/(u*v)
  pdf
}

# reflected MTCJ copula density
# cpar = copula parameter >0 
dmtcjr=function(u,v,cpar)
{ tem1=(1-u)^(-cpar); tem2=(1-v)^(-cpar);
  pdf=(tem1+tem2-1)^(-1/cpar-2)*(1+cpar)*tem1*tem2/((1-u)*(1-v))
  pdf
}

# Joe/B5 copula density
# cpar = copula parameter >1 
djoe=function(u,v,cpar)
{ f1=1-u; f2=1-v;
  tem1=f1^cpar; tem2=f2^cpar;
  sm=tem1+tem2-tem1*tem2; tem=sm^(1./cpar);
  #cdf=1.-tem;
  pdf=tem*((cpar-1.)*tem1*tem2+tem1*tem1*tem2+tem1*tem2*tem2-tem1*tem1*tem2*tem2)
  pdf=pdf/(sm*sm);
  pdf=pdf/(f1*f2);
  pdf
}

# Gumbel copula density
# cpar = copula parameter >1 
dgum=function(u,v,cpar)
{ l1= -log(u); l2= -log(v);
  tem1=l1^cpar; tem2=l2^cpar; sm=tem1+tem2; tem=sm^(1./cpar);
  cdf=exp(-tem);
  pdf=cdf*tem*tem1*tem2*(tem+cpar-1.);
  pdf=pdf/(sm*sm*l1*l2*u*v);
  pdf
}

# reflected Gumbel copula density
# cpar = copula parameter >1 
dgumr=function(u,v,cpar)
{ u=1-u; v=1-v;
  l1= -log(u); l2= -log(v);
  tem1=l1^cpar; tem2=l2^cpar; sm=tem1+tem2; tem=sm^(1/cpar);
  cdf=exp(-tem);
  pdf=cdf*tem*tem1*tem2*(tem+cpar-1.);
  pdf=pdf/(sm*sm*l1*l2*u*v);
  pdf
}

# Galambos copula density
# cpar = copula parameter >0 
dgal=function(u,v,cpar)
{ cpar1= -cpar;
  l1= -log(u); l2= -log(v);
  tem1=l1^(cpar1-1.); tem2=l2^(cpar1-1.); 
  sm=tem1*l1+tem2*l2; tem=sm^(1./cpar1-1.);
  lcdf=-(l1+l2-tem*sm);  # cdf=exp(lcdf);
  deriv12= 1.+ (tem/sm)*tem1*tem2*( 1.+cpar+tem*sm) -tem*(tem1+tem2);
  pdf=lcdf+l1+l2;
  pdf=exp(pdf)*deriv12;
  pdf
}

# Huesler-Reiss copula density
# cpar = copula parameter >0 
dhr=function(u,v,cpar)
{ x=-log(u); y=-log(v); 
  z=x/y; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;  
  p1=pnorm(tem1); p2=pnorm(tem2);
  cdf=exp(-x*p1-y*p2);
  pdf=cdf*(p1*p2+.5*cpar*dnorm(tem1)/y)/(u*v);
  pdf
}


#============================================================

# bivariate normal density
# x = 2-vector or 2-column matrix
# rho = correlation parameter with -1<rho<1
dbvn=function(x,rho)
{ if(is.matrix(x)) { x1=x[,1]; x2=x[,2] }
  else { x1=x[1]; x2=x[2] }
  qf=x1^2+x2^2-2*rho*x1*x2
  qf=qf/(1-rho^2)
  con=sqrt(1-rho^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf
}

# bivariate normal density version 2, 
# x1,x2 = scalars or vectors
# rho = correlation parameter with -1<rho<1
dbvn2=function(x1,x2,rho)
{ qf=x1^2+x2^2-2*rho*x1*x2
  qf=qf/(1-rho^2)
  con=sqrt(1-rho^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf
}

# bivariate normal copula density
# cpar = copula parameter with -1<cpar<1
dbvncop=function(u,v,cpar)
{ x1=qnorm(u); x2=qnorm(v)
  qf=x1^2+x2^2-2*cpar*x1*x2
  qf=qf/(1-cpar^2)
  con=sqrt(1-cpar^2)*(2*pi)
  pdf=exp(-.5*qf)/con
  pdf=pdf/(dnorm(x1)*dnorm(x2))
  pdf
}

# bivariate t density
# x1,x2 = scalars or vectors
# param = (rho,nu) with -1<rho<1 and nu>0
# log = T to return log density, =F to return density
dbvt=function(x1,x2,param,log=FALSE)
{ rh=param[1]; nu=param[2]
  lcon=lgamma((nu+2.)/2.)-lgamma(nu/2.)-log(pi*nu*sqrt(1.-rh*rh));
  ex=-(nu+2.)/2.;
  r11=1./(1-rh*rh); r12=-rh*r11;
  qf=x1*x1*r11+x2*x2*r11+2.*x1*x2*r12;
  tem=1.+qf/nu;
  lpdf=lcon+ex*log(tem)
  if(log) { return(lpdf) } else { return(exp(lpdf)) }
}

# bivariate t copula density
# cpar = copula parameter (rho,nu) with -1<rho<1, nu>0
dbvtcop=function(u,v,cpar)
{ rho=cpar[1]; nu=cpar[2]
  con=exp(lgamma((nu+2.)/2.)-lgamma(nu/2.))/(pi*nu);
  con=con/sqrt(1.-rho*rho);
  ex=-(nu+2.)/2.;
  xt1=qt(u,nu); den1=dt(xt1,nu); 
  xt2=qt(v,nu); den2=dt(xt2,nu);
  r11=1./(1-rho*rho); r12=-rho*r11;
  qf=xt1*xt1*r11+xt2*xt2*r11+2.*xt1*xt2*r12;
  tem=1.+qf/nu;
  pdf=con*(tem^ex)/(den1*den2);
  #lpdf=log(con)+ex*log(tem)-log(den1)-log(den2);
  pdf
}

# bivariate t copula cdf
# cpar = copula parameter (rho,nu) with -1<rho<1, nu>0
#   nu is an integer here because input to pbvt is integer
pbvtcop=function(u,v,cpar)
{ rho=cpar[1]; nu=cpar[2]
  if(nu<1) nu=1
  nu=floor(nu) # input to pbvt is an integer >= 1.
  # need correction for boundary, check for better fix later (pbvn is OK)
  u[u>=1]=.9999999
  v[v>=1]=.9999999
  u[u<=0]=.0000001
  v[v<=0]=.0000001
  xt=qt(u,nu); yt=qt(v,nu); 
  # note that pbvt returns -1 for nu<1.
  out=pbvt(xt,yt,cpar) # function in this library
  out
}

# bivariate normal copula
# cpar = copula parameter with -1<cpar<1
pbvncop=function(u,v,cpar)
{ # endpoint corrections to prevent NaN
  u[1-u<1.e-9]=1-1.e-9
  v[1-v<1.e-9]=1-1.e-9
  u[u<1.e-9]=1.e-9
  v[v<1.e-9]=1.e-9
  x=qnorm(u)
  y=qnorm(v)
  cdf=pbnorm(x,y,cpar)  # function in this library
  cdf
}

#============================================================
