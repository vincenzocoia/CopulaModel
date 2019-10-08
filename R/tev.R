# Functions for bivariate t-EV (extreme value limit of bivariate t as in
#  Demarta-McNeil (2005) and Nikoloulopoulos-Joe-Li (2009).
# This code is modification of code written by David Lee.

# bivariate t-EV copula cdf, 
# u,v: scalars or vectors in (0,1) 
# cpar = (rho,nu) with |rho|<1, nu>0
# Output: cdf
ptev=function(u,v,cpar)
{ u[u<=0]=1e-10; v[v<=0]=1e-10;
  u[u>=1]=1-1e-10; v[v>=1]=1-1e-10;
  x=-log(u); y=-log(v); z=y/x;
  if(is.matrix(cpar)) { rho=cpar[,1]; nu=cpar[,2] }
  else { rho=cpar[1]; nu=cpar[2] }
  znu=z^(1/nu)
  tem1=sqrt((nu+1)/(1-rho^2))*(1/znu-rho);
  tem2=sqrt((nu+1)/(1-rho^2))*(znu-rho);
  cdf=exp(-x*pt(tem1,df=nu+1)-y*pt(tem2,df=nu+1));
  cdf
}

# t-EV copula conditional cdf C_{2|1}(v|u),
# u,v: scalars or vectors in (0,1) 
# cpar = (rho,nu) with |rho|<1, nu>0
# Output: conditional cdf
pcondtev=function(v,u,cpar)
{ u[u<=0]=1e-10; v[v<=0]=1e-10;
  u[u>=1]=1-1e-10; v[v>=1]=1-1e-10;
  x=-log(u); y=-log(v);
  if(is.matrix(cpar)) { rho=cpar[,1]; nu=cpar[,2] }
  else { rho=cpar[1]; nu=cpar[2] }
  ydx=y/x; ydxnu=ydx^(1/nu)
  xdy=x/y; xdynu=xdy^(1/nu); 
  ze=sqrt((nu+1)/(1-rho^2))
  tem1=ze*(xdynu-rho); tem2=ze*(ydxnu-rho);
  Taxy = pt(tem1,df=nu+1); Tayx = pt(tem2,df=nu+1);
  cdf=exp(-x*Taxy-y*Tayx);
  taxy = dt(tem1,df=nu+1); tayx = dt(tem2,df=nu+1);
  cond=cdf/u*(Taxy + taxy*ze*xdynu/nu - tayx*ze/nu*ydx*ydxnu)
  cond
}

#dtp= function(x,nu)
#{ return(-gamma((nu+1)/2)*x*(nu+1)*(1+x^2/nu)^(-(nu+3)/2)/gamma(nu/2)/sqrt(nu^3*pi));
#}

#t-EV copula density, 
# u,v: scalars or vectors in (0,1) 
# cpar = (rho,nu) with |rho|<1, nu>0
# Output: pdf
dtev=function(u,v,cpar)
{ u[u<=0]=1e-10; v[v<=0]=1e-10;
  u[u>=1]=1-1e-10; v[v>=1]=1-1e-10;
  x=-log(u); y=-log(v);
  if(is.matrix(cpar)) { rho=cpar[,1]; nu=cpar[,2] }
  else { rho=cpar[1]; nu=cpar[2] }
  ydx=y/x; ydxnu=ydx^(1/nu)
  xdy=x/y; xdynu=xdy^(1/nu); 
  ze=sqrt((nu+1)/(1-rho^2))
  tem1=ze*(xdynu-rho); tem2=ze*(ydxnu-rho);
  Taxy= pt(tem1,df=nu+1); Tayx= pt(tem2,df=nu+1);
  cdf=exp(-x*Taxy-y*Tayx);
  taxy= dt(tem1,df=nu+1); tayx= dt(tem2,df=nu+1);
  tpaxy= -taxy*(nu+2)*tem1/(nu+1)/(1+tem1^2/(nu+1)); 
  tpayx= -tayx*(nu+2)*tem2/(nu+1)/(1+tem2^2/(nu+1)); 
  brac1= (-Taxy - taxy*ze*xdynu/nu + tayx*ze/nu/xdynu/xdy) * 
          (-Tayx - tayx*ze/xdynu/nu + taxy*ze*xdynu*xdy/nu);
  brac2= taxy*ze*(nu+1)*xdynu/nu^2/y + tayx*ze*(nu+1)/nu^2/xdynu/x + 
      tpaxy*ze^2*xdynu^2/nu^2/y + tpayx*ze^2/xdynu^2/nu^2/x;
  pdf=cdf*(brac1+brac2)/u/v;
  pdf
}

# t-EV copula: copula parameter to Kendall's tau
# cpar = (rho,nu) with |rho|<1, nu>0
tev.cpar2tau=function(cpar)
{ tautev= function(w,rho,nu)
  { w1=1-w; a = sqrt((nu+1)/(1-rho^2)); w1w=w/w1;
    Tw= pt(a*(w1w^(1/nu)-rho),df=nu+1); Tw1= pt(a*(w1w^(-1/nu)-rho),df=nu+1);
    tw= dt(a*(w1w^(1/nu)-rho),df=nu+1); tw1= dt(a*(w1w^(-1/nu)-rho),df=nu+1);
    B=w*Tw+w1*Tw1;
    Bp= Tw + w*tw*a*w1w^(1/nu)/nu/w/w1 - Tw1 - w1*tw1*a*w1w^(-1/nu)/nu/w/w1;
    ((2*w-1)*Bp*B + w*w1*Bp*Bp )/(B*B)
  }
  rho=cpar[1]; nu=cpar[2]
  tem=integrate(tautev,0,1,rho=rho,nu=nu,rel.tol=1.e-6);
  tau=tem$value;
  tau
}

# t-EV copula: copula parameter to Spearman's rhoS
# cpar = (rho,nu) with |rho|<1, nu>0
tev.cpar2rhoS=function(cpar)
{ sptev= function(w,rho,nu)
  { w1=1-w; a=sqrt((nu+1)/(1-rho^2)); w1w=w/w1;
    Tw= pt(a*(w1w^(1/nu)-rho),df=nu+1); Tw1= pt(a*(w1w^(-1/nu)-rho),df=nu+1);
    B=w*Tw+w1*Tw1;
    tem=1/(B+1);
    tem^2;
  }
  rho=cpar[1]; nu=cpar[2]
  tem=integrate(sptev,0,1,rho=rho,nu=nu,rel.tol=1.e-6);
  rhos=12*tem$value-3;
  rhos;
}

# bivariate t-EV copula: copula parameter to Blomqvist beta
# cpar = (rho,nu) with |rho|<1, nu>0
tev.cpar2b=function(cpar)
{ if(is.matrix(cpar)) { rho=cpar[,1]; nu=cpar[,2] }
  else { rho=cpar[1]; nu=cpar[2] }
  ex=2*(log(2))*pt(sqrt((nu+1)*(1-rho)/(1+rho)),nu+1) 
  4*exp(-ex)-1
}

# tail dependence same as for bvt
tev.cpar2lm=bvt.cpar2lm

