# Functions for cdf and pdf for some 2-parameter copula families
# bb5 bb6 bb8 bb10; qcondbb to improve later for starting points
# The copula parameter cpar is a 2-vector (or a 2-column matrix in some cases)
# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v are vectors of the same length
#   or if one of these is a vector and the other is a scalar
# 0<p<1 for qcond functions

#============================================================

# functions for BB5 copula

# BB5 copula cdf
# cpar = copula parameter with th>=1, de>0
pbb5=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  x=-log(u); y=-log(v)
  xt=x^th; yt=y^th
  xdt=x^(-de*th); ydt=y^(-de*th)
  xydt=xdt+ydt
  xyp=xydt^(-1/de)
  w=xt+yt-xyp
  wth=w^(1/th)
  cdf=exp(-wth)
  cdf
}

# BB5 conditional cdf C_{2|1}(v|u;cpar)
# cpar = copula parameter with th>=1, de>0
pcondbb5=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  x=-log(u); y=-log(v)
  xt=x^th; yt=y^th
  xdt=x^(-de*th); ydt=y^(-de*th)
  xydt=xdt+ydt
  xyp=xydt^(-1/de)
  w=xt+yt-xyp
  wth=w^(1/th)
  cdf=exp(-wth)
  zx=xt/x-xyp/xydt*xdt/x
  ccdf=cdf*(wth/w)*zx/u
  ccdf
}

# BB5 density
# cpar = copula parameter with th>=1, de>0
dbb5=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  x=-log(u); y=-log(v)
  xt=x^th; yt=y^th
  xdt=x^(-de*th); ydt=y^(-de*th)
  xydt=xdt+ydt
  xyp=xydt^(-1/de)
  w=xt+yt-xyp
  wth=w^(1/th)
  cdf=exp(-wth)
  zx=xt/x-xyp/xydt*xdt/x
  zy=yt/y-xyp/xydt*ydt/y
  pdf=cdf*wth/w/w/u/v 
  pdf=pdf*((wth+th-1)*zx*zy + th*(1+de)*w*(xyp/xydt/xydt)*(xdt/x)*(ydt/y))
  pdf
}

# this code needs checking in C with 10^6 cases
# inverse of pcondbb5
# p = value in (0,1)
# u = value in (0,1)
# cpar = copula parameter with th>=1, de>0
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb5=function(p,u,cpar, eps=1.e-6,mxiter=30,iprint=F)
{ th=cpar[1]; de=cpar[2]
  de1=1./de; th1=1/th
  x=-log(u); 
  xt=x^th; 
  xdt=x^(-de*th); 
  con=x-log(p); 
  iter=0; diff=1.;
  y=.5*x; # what is good starting point?
  be=4*pbb5(.5,.5,cpar)-1
  if(be>=0.8) y=x
  while(iter<mxiter & max(abs(diff))>eps)
  { yt=y^th; ydt=y^(-de*th)
    xydt=xdt+ydt; xyp=xydt^(-1/de)
    w=xt+yt-xyp; wth=w^(1/th)
    zx=(xt-xyp/xydt*xdt)/x
    zy=(yt-xyp/xydt*ydt)/y
    h=-wth+(th1-1)*log(w)+log(zx)+con
    hp=(-wth+(1-th))*zy/w - th*(1+de)*(xyp/xydt/xydt)*(xdt/x)*(ydt/y)/zx
    diff=h/hp;
    y=y-diff;
    if(iprint) cat(iter, diff, y,"\n")
    while(min(y)<=0. | max(abs(diff))>5) { diff=diff/2.; y=y+diff;}
    iter=iter+1;
  }
  if(iprint & iter>=mxiter) cat("***did not converge\n");
  exp(-y)
}

#============================================================

# functions for BB6 copula

# BB6 copula cdf
# cpar = copula parameter with th>=1, de>=1
pbb6=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  ubar=1-u; vbar=1-v
  x=-log(1-ubar^th); y=-log(1-vbar^th)
  xd=x^de; yd=y^de; sm=xd+yd; tem=sm^(1/de)
  w=exp(-tem)
  cdf=(1-w)^(1/th)
  cdf=1-cdf
  cdf
}

# BB6 conditional cdf C_{2|1}(v|u;cpar) 
# cpar = copula parameter with th>=1, de>=1
pcondbb6=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  ubar=1-u; vbar=1-v
  zu=1-ubar^th
  x=-log(zu); y=-log(1-vbar^th)
  xd=x^de; yd=y^de; sm=xd+yd; tem=sm^(1/de)
  w=exp(-tem); #zu=exp(-x)
  ccdf=((1-w)/(1-zu))^(1/th-1) *(w/zu) *(tem/sm) * (xd/x)
  ccdf
}

# BB6 density
# cpar = copula parameter with th>=1, de>=1
dbb6=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  ubar=1-u; vbar=1-v
  zu=1-ubar^th; zv=1-vbar^th
  x=-log(zu); y=-log(zv)
  xd=x^de; yd=y^de; sm=xd+yd; tem=sm^(1/de)
  w=exp(-tem); 
  pdf=(1-w)^(1/th-2) *w*(tem/sm/sm)*(xd/x)*(yd/y) *((th-w)*tem+th*(de-1)*(1-w))
  pdf=pdf* (1-zu)*(1-zv)/zu/zv/ubar/vbar
  pdf
}

# this code needs checking in C with 10^6 cases
# this can get into NaN/Inf problems for th and de >=4.5 and u near 1
# inverse of pcondbb6
# p = value in (0,1)
# u = value in (0,1)
# cpar = copula parameter with th>=1, de>=1
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb6=function(p,u,cpar, eps=1.e-6,mxiter=30,iprint=F)
{ th=cpar[1]; de=cpar[2]
  de1=1./de; th1=1/th
  be=4*pbb6(.5,.5,cpar)-1
  ubar=1-u; 
  zu=1-ubar^th
  x=-log(zu); #zu=exp(-x)
  xd=x^de; 
  con=(de-1)*log(x)-(th1-1)*log(1-zu)+x-log(p); 
  iter=0; diff=1.;
  y=.5*x; # what is good starting point?
  if(be>=0.8) y=x
  while(iter<mxiter & max(abs(diff))>eps)
  { yd=y^de; sm=xd+yd; tem=sm^(1/de)
    w=exp(-tem); 
    h=(th1-1)*log(1-w)-tem+(de1-1)*log(sm)+con
    hp=(th1-1)*w*tem*yd/(1-w)/sm/y - tem*yd/sm/y +(1.-de)*yd/y/sm
    diff=h/hp;
    y=y-diff;
    if(iprint) cat(iter, diff, y,"\n")
    if(any(is.nan(y))) 
    { if(iprint) cat(p,u,cpar,x,"\n"); 
      return(u) 
    }
    while(min(y)<=0. | max(abs(diff))>5) { diff=diff/2.; y=y+diff;}
    iter=iter+1;
  }
  if(iprint & iter>=mxiter) cat("***did not converge\n");
  1-(1-exp(-y))^th1
}

#============================================================

# functions for BB8 copula

# BB8 copula cdf
# cpar = copula parameter with vth>=1, 0<de<=1
pbb8=function(u,v,cpar)
{ if(is.matrix(cpar)) { vth=cpar[,1]; de=cpar[,2] }
  else { vth=cpar[1]; de=cpar[2] }
  x=1-(1-de*u)^vth
  y=1-(1-de*v)^vth
  eta1=1/(1-(1-de)^vth)  # reciprocal of eta
  tem=(1-eta1*x*y)^(1/vth)
  cdf=(1-tem)/de
  cdf
}

# BB8 conditional cdf C_{2|1}(v|u;cpar)
# cpar = copula parameter with vth>=1, 0<de<=1
pcondbb8=function(v,u,cpar)
{ if(is.matrix(cpar)) { vth=cpar[,1]; de=cpar[,2] }
  else { vth=cpar[1]; de=cpar[2] }
  ut=(1-de*u)^vth
  x=1-ut
  y=1-(1-de*v)^vth
  eta1=1/(1-(1-de)^vth)
  tem=(1-eta1*x*y)^(1/vth-1)
  den=(1-de*u)/ut
  ccdf=eta1*y*tem/den
  ccdf
}

# BB8 density
# cpar = copula parameter with vth>=1, 0<de<=1
dbb8=function(u,v,cpar)
{ if(is.matrix(cpar)) { vth=cpar[,1]; de=cpar[,2] }
  else { vth=cpar[1]; de=cpar[2] }
  ut=(1-de*u)^vth; vt=(1-de*v)^vth
  x=1-ut; y=1-vt
  eta1=1/(1-(1-de)^vth)
  tem=(1-eta1*x*y)^(1/vth-2)
  pdf=eta1*de*tem*(vth-eta1*x*y)*ut*vt/(1-de*u)/(1-de*v)
  pdf
}


# this code needs checking in C with 10^6 cases
# inverse of pcondbb8
# p = value in (0,1)
# u = value in (0,1)
# cpar = copula parameter with vth>=1, 0<de<=1
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb8=function(p,u,cpar, eps=1.e-6,mxiter=30,iprint=F)
{ vth=cpar[1]; de=cpar[2]
  vth1=1/vth
  eta=(1-(1-de)^vth); eta1=1/eta
  be=4*pbb8(.5,.5,cpar)-1
  ut=(1-de*u)^vth
  x=1-ut
  con=log(eta1)+(1-vth1)*log(1-x)-log(p); 
  iter=0; diff=1.;
  y=.5*x; # what is good starting point?
  if(be>=0.8) y=x
  while(iter<mxiter & abs(diff)>eps)
  { tem=1-eta1*x*y
    h=log(y)+(vth1-1)*log(tem)+con
    hp=1/y+(1-vth1)*eta1*x/tem
    diff=h/hp;
    y=y-diff;
    if(iprint) cat(iter, diff, y,"\n")
    if(any(is.nan(y))) 
    { if(iprint) cat("***", p,u,cpar,x,"\n"); 
      return(u) 
    }
    while(min(y)<=0. | max(y)>=eta) { diff=diff/2.; y=y+diff;}
    iter=iter+1;
  }
  if(iprint & iter>=mxiter) cat("***did not converge\n");
  v=(1-y)^vth1
  v=(1-v)/de
  v
}

#============================================================

# functions for BB10 copula

# BB10 copula cdf based on NB LT
# cpar = copula parameter with th>0, 0<ppi<=1
pbb10=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; ppi=cpar[,2] }
  else { th=cpar[1]; ppi=cpar[2] }
  ut=u^th; vt=v^th
  tem=1-ppi*(1-ut)*(1-vt)
  cdf=u*v*tem^(-1/th)
  cdf
}

# BB10 conditional cdf 
# cpar = copula parameter with th>0, 0<ppi<=1
pcondbb10=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; ppi=cpar[,2] }
  else { th=cpar[1]; ppi=cpar[2] }
  ut=u^th; vt=v^th
  tem=1-ppi*(1-ut)*(1-vt)
  ttem=tem^(-1/th)
  ccdf=ttem/tem*v*(1-ppi+ppi*vt)
  ccdf
}

# BB10 density
# cpar = copula parameter with th>0, 0<ppi<=1
dbb10=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; ppi=cpar[,2] }
  else { th=cpar[1]; ppi=cpar[2] }
  ut=u^th; vt=v^th
  tem=1-ppi*(1-ut)*(1-vt)
  ttem=tem^(-1/th)
  pdf=ttem/tem/tem
  pdf=pdf*(1-ppi+ppi*(1+th)*ut*vt-ppi*(1-ppi)*(1-ut)*(1-vt))
  pdf
}


# inverse of pcondbb10
# p = value in (0,1)
# u = value in (0,1)
# cpar = copula parameter with th>0, 0<ppi<=1
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb10=function(p,u,cpar, eps=1.e-6,mxiter=30,iprint=F)
{ th=cpar[1]; ppi=cpar[2] 
  th1=1./th;
  be=4*pbb10(.5,.5,cpar)-1
  ut=u^th;
  mxdif=1.; iter=0; 
  diff=1.; 
  v=.7*u; # what is good starting point?
  if(be>=0.8) v=u
  while(mxdif>eps & iter<mxiter)
  { vt=v^th;
    tem=1.-ppi*(1.-ut)*(1.-vt);
    ttem=tem^(-th1);
    ccdf=ttem/tem*v*(1.-ppi+ppi*vt);
    pdf=ttem/tem/tem;
    pdf=pdf*(1.-ppi+ppi*(1.+th)*ut*vt-ppi*(1.-ppi)*(1.-ut)*(1.-vt));
    h=ccdf-p; hp=pdf;
    diff=h/hp;
    v=v-diff;
    iter=iter+1;
    while(min(v)<=0. | max(v)>=1.) { diff=diff/2.; v=v+diff; }
    if(iprint) cat(iter, diff, v, "\n")
    mxdif=abs(diff);
  }
  if(iprint & iter>=mxiter) 
  { cat("***did not converge\n");
    cat("p=", p, " u=", u, " theta=", th, " pi=", ppi, " lastv=", v, "\n")
  }
  v
}

