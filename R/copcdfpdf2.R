# Functions for cdf and pdf for some 2-parameter copula families
# bb1, bb2, bb3, bb4, bb7, bb9.

# The copula parameter cpar is a 2-vector (or a 2-column matrix in some cases)
# 0<u<1, 0<v<1 for all functions
# Most functions here should work if u,v are vectors of the same length
#   or if one of these is a vector and the other is a scalar

#  C(u,v)=G(x,y) for decreasing transforms x,y
#  C_{2|1}(v|u)= {\p G\over \p x} \cdot {\p x\over \p u} ,
#   {\p x\over \p u} =-\de\th(u^{-\th}-1)^{\de-1} u^{-\th-1} 
#  c(u,v) = {\p^2 G\over \p x \p x} \cdot {\p x\over \p u}
#  \cdot {\p y\over \p v}

# BB1 cdf
# x=(u^{-th}-1)^{de} and y=(v^{-\th}-1)^{\de}
# u=(x^{1/\de}+1)^{-1/\th} and v=(y^{1/\de}+1)^{-1/\th}
# C(u,v)=G(x,y) = [1+(x+y)^{1/de}]^{-1/th}
# cpar = copula parameter with th>0, de>=1
pbb1=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1/de; th1=1/th
  ut=(u^(-th)-1); vt=(v^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  cdf=(1+smd)^(-th1)
  cdf
}

# cdf of BB1 reflected 
pbb1r=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  u1=1-u; v1=1-v;
  de1=1/de; th1=1/th
  ut=(u1^(-th)-1); vt=(v1^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  cdf=(1+smd)^(-th1)
  cdf+u+v-1
}

# BB1 density : based on deriv of transform variables
# cpar = copula parameter with th>0, de>=1
dbb1=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1/de; th1=1/th
  ut=(u^(-th)-1); vt=(v^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  tem=(1+smd)^(-th1-2) * (th*(de-1)+(th*de+1)*smd)
  pdf=tem*smd*x*y*(ut+1)*(vt+1)/sm/sm/ut/vt/u/v
  pdf
}

# pdf of BB1 reflected
dbb1r=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  u1=1-u; v1=1-v;
  de1=1/de; th1=1/th
  ut=(u1^(-th)-1); vt=(v1^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  tem=(1+smd)^(-th1-2) * (th*(de-1)+(th*de+1)*smd)
  pdf=tem*smd*x*y*(ut+1)*(vt+1)/sm/sm/ut/vt/u1/v1
  pdf
}

#============================================================

# BB7 cdf
# x=(1-[1-u]^{\th})^{-\de}-1 and y=(1-[1-v]^{\th})^{-\de}-1
# u=1-(1-(1+x)^{-1/\de})^{1/\th} and v=1-(1-(1+y)^{-1/\de})^{1/\th}
# C(u,v)=G(x,y) = 1-[1-(x+y+1)^{-1/de}]^{1/th}
# cpar = copula parameter with th>1, de>0
pbb7=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1/de; th1=1/th
  ut=1.-(1.-u)^th; vt=1.-(1.-v)^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  cdf=1-(1.-smd)^th1;
  cdf
}

# cdf of BB7 reflected
pbb7r=function(u,v,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  u1=1-u; v1=1-v;
  de1=1/de; th1=1/th
  ut=1.-(1.-u1)^th; vt=1.-(1.-v1)^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  cdf=1-(1.-smd)^th1;
  cdf+u+v-1
}

# BB7 density : based on deriv of transform variables
# cpar = copula parameter with th>1, de>0
dbb7=function(u,v,cpar)
{ if(length(cpar)==2) { th=cpar[1]; de=cpar[2] }
  else { th=cpar[,1]; de=cpar[,2] }
  de1=1/de; th1=1/th
  ut=1.-(1.-u)^th; vt=1.-(1.-v)^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  tem=(1.-smd)^(th1-2.);
  tem=tem*(th*(de+1)-(th*de+1)*smd)
  pdf=tem*smd*(x+1)*(y+1)*(1-ut)*(1-vt)/sm/sm/ut/vt/(1-u)/(1-v);
  pdf
}

# pdf of BB7 reflected
dbb7r=function(u,v,cpar)
{ if(length(cpar)==2) { th=cpar[1]; de=cpar[2] }
  else { th=cpar[,1]; de=cpar[,2] }
  u1=1-u; v1=1-v;
  de1=1/de; th1=1/th
  ut=1.-u^th; vt=1.-v^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  tem=(1.-smd)^(th1-2.);
  tem=tem*(th*(de+1)-(th*de+1)*smd)
  pdf=tem*smd*(x+1)*(y+1)*(1-ut)*(1-vt)/sm/sm/ut/vt/u/v;
  pdf
}


#============================================================

# BB4 cdf
# cpar = copula parameter with th>1, de>=1
pbb4=function(u,v,cpar)
{ if(length(cpar)==2) { th=cpar[1]; de=cpar[2] }
  else { th=cpar[,1]; de=cpar[,2] }
  ut=u^(-th)-1; vt=v^(-th)-1
  x=ut^(-de); y=vt^(-de)
  tem=(x+y)^(-1/de)
  cdf=ut+vt+1-tem
  cdf=cdf^(-1/th)
  cdf
}

# cdf of BB4 reflected
# cpar = copula parameter with th>1, de>=1
pbb4r=function(u,v,cpar)
{ if(length(cpar)==2) { th=cpar[1]; de=cpar[2] }
  else { th=cpar[,1]; de=cpar[,2] }
  u1=1-u; v1=1-v;
  ut=u1^(-th)-1; vt=v1^(-th)-1
  x=ut^(-de); y=vt^(-de)
  tem=(x+y)^(-1/de)
  cdf=ut+vt+1-tem
  cdf=cdf^(-1/th)
  cdf+u+v-1
}

# BB4 pdf
# cpar = copula parameter with th>1, de>=1
dbb4=function(u,v,cpar)
{ if(length(cpar)==2) { th=cpar[1]; de=cpar[2] }
  else { th=cpar[,1]; de=cpar[,2] }
  ut=u^(-th)-1; vt=v^(-th)-1
  x=ut^(-de); y=vt^(-de); xy=x+y
  tem=xy^(-1/de)
  sm=ut+vt+1-tem
  xtem=ut/x-tem/xy
  ytem=vt/y-tem/xy
  pdf=sm^(-1/th-2)*x*y/ut/vt*(1+ut)*(1+vt)/u/v
  pdf=pdf*((th+1)*xtem*ytem+th*(1+de)*sm*tem/xy/xy)
  pdf
}

# pdf of BB4 reflected
# cpar = copula parameter with th>1, de>=1
dbb4r=function(u,v,cpar)
{ if(length(cpar)==2) { th=cpar[1]; de=cpar[2] }
  else { th=cpar[,1]; de=cpar[,2] }
  u1=1-u; v1=1-v;
  ut=u1^(-th)-1; vt=v1^(-th)-1
  x=ut^(-de); y=vt^(-de); xy=x+y
  tem=xy^(-1/de)
  sm=ut+vt+1-tem
  xtem=ut/x-tem/xy
  ytem=vt/y-tem/xy
  pdf=sm^(-1/th-2)*x*y/ut/vt*(1+ut)*(1+vt)/u1/v1
  pdf=pdf*((th+1)*xtem*ytem+th*(1+de)*sm*tem/xy/xy)
  pdf
}

#============================================================

# BB2 cdf
# x=\exp\{\de(u^{-\th}-1)\} and y=\exp\{\de(v^{-\th}-1)\}
# u=(1+\de^{-1}\log x)^{-1/\th} and v=(1+\de^{-1}\log y)^{-1/\th}
# C(u,v)=G(x,y) = [1+de^{-1}*log(x+y-1)]^{-1/th}
# cpar = copula parameter with th>0, de>0
pbb2=function(u,v,cpar,iprint=F)
{ th=cpar[1]; de=cpar[2] 
  de1=1/de; th1=1/th
  ut=(u^(-th)-1); vt=(v^(-th)-1)
  x=exp(ut*de)-1; y=exp(vt*de)-1
  if(iprint) cat(ut,vt,x,y,"\n")
  if(is.infinite(x) || is.infinite(y))
  { lr=de*(vt-ut); r=exp(lr);
    cdf=(1+ut+de1*log(1.+r))^(-th1)
  }
  else
  { sm=x+y+1; smd=de1*log(sm)
    cdf=(1+smd)^(-th1)
  }
  cdf
}

# BB2 density : based on deriv of transform variables
# cpar = copula parameter with th>0, de>0
dbb2=function(u,v,cpar,iprint=F)
{ th=cpar[1]; de=cpar[2] 
  de1=1/de; th1=1/th
  ut=(u^(-th)-1); vt=(v^(-th)-1)
  x=exp(ut*de)-1; y=exp(vt*de)-1
  if(iprint) cat(ut,vt,x,y,"\n")
  if(is.infinite(x) || is.infinite(y))
  { lr=de*(vt-ut); r=exp(lr); lr1=log(1.+r)
    tem=(1.+ut+de1*lr1)^(-th1-2)
    tem=tem*(1+th+th*lr1+th*de*(ut+1.))
    pdf=tem*(ut+1)*(vt+1)*r/u/v/(1.+r)^2   # y=rx => x/sm=1/(1+r) as x->oo
  }
  else
  { sm=x+y+1; smd=de1*log(sm)
    tem=(1+smd)^(-th1-2)
    tem=tem*(1+th+th*de*(1+smd))
    pdf=tem*(x+1)*(y+1)*(ut+1)*(vt+1)/sm/sm/u/v
  }
  pdf
}


#============================================================


# BB9 cdf  
# original alp -> 1/ga to get increasing in concordance
# x=ga1-\log u and y=ga1-\log v  , ga1=1/ga
# u=e^{-x+ga1} and v=e^{-y+\ga1}
# C(u,v)=G(x,y)=\exp\{-(x^\th+y^\th-ga1^\th)^{-1/\th}+ga1\}
# cpar = copula parameter with th>1, ga>=0
pbb9=function(u,v,cpar)
{ th=cpar[1]; ga=cpar[2]; ga1=1/ga 
  x= ga1-log(u); y= ga1-log(v);
  temx=x^th; temy=y^th; sm=temx+temy-ga1^th; smt=sm^(1./th);
  cdf=exp(-smt+ga1);
  cdf
}

# BB9 density 
# original alp -> 1/ga to get increasing in concordance
# cpar = copula parameter with th>1, ga>=0
dbb9=function(u,v,cpar)
{ th=cpar[1]; ga=cpar[2]; ga1=1/ga 
  x= ga1-log(u); y= ga1-log(v);
  temx=x^th; temy=y^th; sm=temx+temy-ga1^th; smt=sm^(1./th);
  cdf=exp(-smt+ga1);
  pdf=cdf*(smt+th-1)*smt*temx*temy/sm/sm/x/u/y/v;
  pdf 
}


#============================================================


# BB3 cdf  
# cpar = copula parameter with th>1, de>0
pbb3=function(u,v,cpar,iprint=F)
{ th=cpar[1]; de=cpar[2]
  de1=1./de; th1=1./th; dt=(de^th1);
  ul=-log(u); vl=-log(v);
  ut=(ul^th); vt=(vl^th);
  x=exp(ut*de)-1.; y=exp(vt*de)-1.;
  if(iprint) cat(ut,vt,x,y,"\n")
  if(is.infinite(x) || is.infinite(y) || x>1.e200 || y>1.e200) 
  { lr=de*(vt-ut); r=exp(lr); lx=de*ut;
    sm=r+1.; sml=log(sm)+lx; 
    tem=(sml^th1);
    cdf=exp(-tem/dt);
  }
  else if(x<=1.e-10 || y<=1.e-10) 
  { xx=de*ut; yy=vt*de; r=vt/ut;
    sm=xx+yy; 
    cdf=exp(-(sm^th1)/dt);
  }
  else
  { sm=x+y+1.; sml=log(sm);
    tem=(sml^th1);
    cdf=exp(-tem/dt);
  }
  cdf
}
  

# BB3 density 
# cpar = copula parameter with th>1, de>0
dbb3=function(u,v,cpar,iprint=F)
{ th=cpar[1]; de=cpar[2]
  de1=1./de; th1=1./th; dt=(de^th1);
  ul=-log(u); vl=-log(v);
  ut=(ul^th); vt=(vl^th);
  x=exp(ut*de)-1.; y=exp(vt*de)-1.;
  if(iprint) cat(ut,vt,x,y,"\n")
  if(is.infinite(x) | is.infinite(y) | x>1.e200 | y>1.e200) 
  { lr=de*(vt-ut); r=exp(lr); lx=de*ut;
    sm=r+1.; sml=log(sm)+lx; 
    tem=(sml^th1);
    cdf=exp(-tem/dt);
    pdf=cdf*((th-1)*dt +th*dt*sml + tem)
    pdf=pdf*tem*de*de*ut*vt*r/sml/sml/sm/sm/ul/vl/u/v/dt/dt;
  }
  else if(x<=1.e-10 | y<=1.e-10) 
  { xx=de*ut; yy=vt*de; r=vt/ut;
    sml=xx+yy; 
    tem=(sml^th1);
    cdf=exp(-tem/dt);
    pdf=cdf*((th-1)*dt +th*dt*sml +tem)
    pdf=pdf*tem*de*de*ut*vt*(xx+1)*(yy+1)/sml/sml/(1+sml)/(1+sml)/ul/vl/u/v/dt/dt;
  }
  else
  { sm=x+y+1.; sml=log(sm);
    tem=(sml^th1);
    cdf=exp(-tem/dt);
    pdf=cdf*((th-1)*dt +th*dt*sml +tem)
    pdf=pdf*tem*de*de*ut*vt*(x+1)*(y+1)/sml/sml/sm/sm/ul/vl/u/v/dt/dt;
  }
  pdf
}
