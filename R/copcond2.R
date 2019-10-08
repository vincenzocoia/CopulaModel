# pcond C_{2|1}, qcond C_{2|1}^{-1} for some 2-parameter copula families
# bb1, bb2, bb3, bb4, bb7, bb9

#  C_{2|1}(v|u) and C_{2|1}^{-1}(v|u)
# inputs v and u can be vectors of same length
# or u is vector, v is scalar 

# C_{2|1}(v|u) for BB1 : based on derivative of transform variables
# cpar = copula parameter with (th,de) : th>0, de>1
#    cpar can be a matrix with #row=length(u)
pcondbb1=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1/de; th1=1/th
  ut=(u^(-th)-1); vt=(v^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  tem=(1+smd)^(-th1-1)
  ccdf=tem*smd*x*(ut+1)/sm/ut/u
  ccdf
}

# reflected BB1
pcondbb1r=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1/de; th1=1/th
  u1=1-u; v1=1-v
  ut=(u1^(-th)-1); vt=(v1^(-th)-1)
  x=ut^de; y=vt^de
  sm=x+y; smd=sm^(de1)
  tem=(1+smd)^(-th1-1)
  ccdf=tem*smd*x*(ut+1)/sm/ut/u1
  1-ccdf
}

# qcondbb1 from C
# 0<p<1 
# cpar = copula parameter with (th,de) : th>0, de>1
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb1=function(p,u,cpar, eps=1.e-6, mxiter=30, iprint=F)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  if(length(p)>1 | length(u)>1 | length(th)>1 | length(de)>1)
  { # later add some checks to make p,u,th,de in range
    np=length(p)
    nu=length(u)
    nth=length(th)
    nde=length(de)
    nn=max(np,nu,nth,nde)
    if(np<nn & np>1) return(NA)
    if(nu<nn & nu>1) return(NA)
    if(nth<nn & nth>1) return(NA)
    if(nde<nn & nde>1) return(NA)
    if(nn>1)
    { if(np==1) p=rep(p,nn)
      if(nu==1) u=rep(u,nn)
      if(nth==1) th=rep(th,nn)
      if(nde==1) de=rep(de,nn)
    }
    out= .C("qcbb1", as.integer(nn), as.double(p), as.double(u),
         as.double(th), as.double(de), v=as.double(rep(0,nn)) )
    return(out$v)
  }
  else
  { de1=1./de; th1=1./th;
    ut=(u^(-th)-1); 
    x=ut^de; 
    den=(1.+ut)^(-th1-1)*ut/x; # density of univariate margin
    pden=p*den; # term 1/(th*de) cancels
    y=pden^(-1./(th1*de1+1))-x;
    if(x<1.e-5) # empirically based
    { #y=x*(p^(1.-de1));
      y=min( x* (p^(-de/(de-1.))-1.), 1.e-5)
      if(iprint) cat("\nsmall x case ",x,y,"\n");
    }
    if(x>1.e5) # empirically based
    { r=p^(-de*th/(1.+de*th))-1.;
      y=r*x;
      if(iprint) cat("\nlarge x case ",x,y,"\n");
      eps=eps*.0001; # good enough
    }
    #if th<0.1 or de<1.1 use boundary of Gumbel and MTCJ as starting point
    if(de<=1.1)  # MTCJ boundary
    { thr=-th/(1.+th); tem=(p^thr)-1.;
      tem=tem*(u^(-th))+1.; y=(tem-1.)^de;
    }
    else if (th<0.2) { v=qcondgum(p,u,de); y=((v^(-th))-1.)^de; }
    # v=(y^(1/de)+1)^(-1/th) so y=(v^(-th)-1)^de
    # modified Newton-Raphson
    diff=1; iter=0;
    while(abs(diff/y)> eps & iter<mxiter)
    { sm=x+y; smd=sm^de1;
      G21=(1.+smd)^(-th1-1.) * smd/sm; 
      gpdf=-G21; 
      gpdf=gpdf/(1.+smd)/sm/de/th;
      gpdf=gpdf*(th*(de-1.)+(th*de+1.)*smd);
      iter=iter+1;
      diff=(G21-pden)/gpdf;
      y=y-diff;
      if(iprint) cat(iter, y, diff,"\n");
      while(y<=0.) { diff=diff/2.; y=y+diff; }
    }
    v=(y^de1+1.)^(-th1);
    if(iprint) cat(v, "ended at iter. ", iter, "\n");
    return(v)
  }
}

# reflected BB1
qcondbb1r=function(p,u,cpar)
{ u1=1-u; p1=1-p
  v=qcondbb1(p1,u1,cpar) 
  1-v
}

#============================================================

# BB7 pcond
# cpar = copula parameter with (th,de) : th>1, de>0
pcondbb7=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1./de; th1=1./th;
  ut=1.-(1.-u)^th; vt=1.-(1.-v)^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  tem=(1.-smd)^(th1-1.);
  ccdf=tem*smd*(x+1)*(1.-ut)/sm/ut/(1.-u);
  ccdf;
}

# reflected BB7 
# cpar = copula parameter with (th,de) : th>1, de>0
pcondbb7r=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  de1=1./de; th1=1./th;
  u1=1-u; v1=1-v
  ut=1.-u^th; vt=1.-v^th; 
  x=ut^(-de)-1; y=vt^(-de)-1;
  sm=x+y+1; smd=sm^(-de1);
  tem=(1.-smd)^(th1-1.);
  ccdf=tem*smd*(x+1)*(1.-ut)/sm/ut/u;
  1-ccdf;
}

# BB7 qcond
# 0<p<1,  
# cpar = copula parameter with cpar=(th,de) : th>1, de>0
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb7=function(p,u,cpar, eps=1.e-6, mxiter=30, iprint=F)
{ th=cpar[1]; de=cpar[2] 
  de1=1./de; th1=1./th;
  ut=1.-(1.-u)^th; 
  x=ut^(-de)-1.; 
  if(x<=0.)  # might occur once in 1.e6
  { v=.999999; if(u>v) v=u;
    if(iprint) cat("\n**** x below 1 " , p,u,x,"\n");
    return(v);
  }
  den=(1.-ut)^(th1-1)*ut/(x+1.); # density of univariate margin
  pden=p*den; # term 1/(th*de) cancels
  # starting guess for y, 
  rhs=pden*(1.-(2*x+1.)^(-de1))^(1.-th1);
  y=rhs^(-de/(de+1))-1-x;
  if(y<=0.) y=0.1; # need better starting point if x<eps
  if(th>3 & (u>0.8 | p>0.9)) # switch to bisection method, scalar p,u
  { diff=1.
    v=u;
    v1=0.; v2=1.;
    vold=v;
    if(iprint) cat("\n (p,u)=", p,u,"\n")
    while(diff>eps)
    { ccdf=pcondbb7(v,u,cpar); di=ccdf-p;
      if(di<0) { v1=v; } else { v2=v; }
      diff=v2-v1; vold=v; v=(v1+v2)/2.;
      if(iprint) cat(vold,ccdf,di,"\n")
    }
    return(v)
  }

  if(x<1.e-5) # empirically based
  { epsx=de*(1.-ut);  # x
    tem=p*(1-(1+de1)*epsx);
    tem=tem^(-th/(th-1.))-1.;
    epsy=tem*epsx;
    if(epsy>1.e-5) epsy=1.e-5;
    y=epsy;
    if(iprint) cat("\nsmall x case " ,x,y,"\n");
  }
  # *** add boundary cases later for small th or de
  if(th<1.01) # MTCJ boundary
  { thr=-th/(1.+th); tem=(p^thr)-1.;
    tem=tem*(u^(-th))+1.; y=(tem-1.)^de;
  }
  else if(de<0.1)  # to check 
  { v=1.-qcondjoe(p,u,de)
    y=v^th; y=(1-y)^(-de)-1
  }

  # modified Newton-Raphson
  diff=1; iter=0;
  while(max(abs(diff/y))> eps & iter<mxiter)
  { sm=x+y+1; smd=sm^(-de1);
    G21=(1.-smd)^(th1-1.) * smd/sm; 
    gpdf=-G21; 
    gpdf=gpdf/(1.-smd)/sm/de/th;
    gpdf=gpdf*(th*(de+1.)-(th*de+1.)*smd);
    iter=iter+1;
    diff=(G21-pden)/gpdf;
    y=y-diff;
    while(min(y)<=0. | max(abs(diff))>5) { diff=diff/2.; y=y+diff; }
    if(iprint) cat(iter, y, diff, "\n");
  }
  v=1.-(1.-(y+1)^(-de1))^th1;
  if(iprint) cat(v, "ended at iter. ", iter, "\n");
  v
}

# reflected BB7
qcondbb7r=function(p,u,cpar)
{ u1=1-u; p1=1-p
  v=qcondbb7(p1,u1,cpar) 
  1-v
}

#============================================================

# BB2 pcond 
# cpar = copula parameter with th>0, de>0
pcondbb2=function(v,u,cpar, iprint=F)
{ th=cpar[1]; de=cpar[2] 
  de1=1./de; th1=1./th;
  ut=u^(-th)-1.; vt=v^(-th)-1.;
  x=exp(ut*de)-1.; y=exp(vt*de)-1.;
  if(iprint) cat(ut,vt,x,y,"\n")
  if(is.infinite(x) || is.infinite(y))
  { lr=de*(vt-ut); r=exp(lr);
    tem=(1.+ut+de1*log(1.+r))^(-th1-1);
    ccdf=tem*(ut+1.)/u/(1.+r);  # y=rx => x/sm=1/(1+r) as x->oo
  }
  else
  { sm=x+y+1.; smd=de1*log(sm);
    tem=(1.+smd)^(-th1-1.);
    ccdf=tem*(x+1.)*(ut+1.)/sm/u;
  }
  ccdf;
}

# BB2 qcond 
# 0<p<1
# cpar = copula parameter with th>0, de>0
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb2=function(p,u,cpar, eps=1.e-6, mxiter=30, iprint=F)
{ th=cpar[1]; de=cpar[2] 
  de1=1./de; th1=1./th;
  ut=u^(-th)-1.; 
  x=exp(ut*de)-1; 
  den=(1.+ut)^(-th1-1.)/(x+1.); # density of univariate margin
  if(is.infinite(x) | den<=1.e-100)
  { # v =pow(1+ ut+de1*log(1./p-1.),-th1); 
    #  solve a different equation based on y=rx as x->oo
    #  g(r)= 1+log(1+r) /lx -p(1+r); g'(r)=1/(1+r)/lx -p;
    # correction tha=-th/(1+th)
    #  g(r)= 1+log(1+r) /lx -p^tha (1+r)^tha; 
    #    g'(r)=1/(1+r)/lx -p^tha * tha *(1+r)^(tha-1);
    #  v=(1+de1*(log(r)+lx))^{-th1}
    if(iprint) cat("\n** infinite x ", p,u,th,de,ut,"\n");
    lx=ut*de; r=1.; tha=-th/(1+th)
    diff=1; iter=0;
    while(abs(diff)>eps  & iter<mxiter)
    { sm=1.+r; smd=de1*log(sm); rhs=(p*sm)^tha
      g=1.+log(sm)/lx-rhs;
      gp=1./sm/lx-tha*rhs/sm;
      iter=iter+1;
      diff=g/gp;
      r=r-diff;
      while(r<=0.) { diff=diff/2.; r=r+diff; }
      if(iprint) cat(iter, r, diff,"\n");
    }
    v=(1.+de1*(log(r)+lx))^(-th1);
    if(iprint) cat("** infinite x ", p,u,th,de,ut,v,"\n");
    return(v);
  }
  pden=p*den; # term 1/(th*de) cancels
  # starting guess for y, 
  rhs=(1.+log(2*x+1))^(-1.-th1)/pden;
  y=rhs-1+x;  
  if(y<=0.) y=0.1;
  if(y>1000.) y=1000.;
  # if x is large, y is large in need to replace by lower value
  # *** add boundary case later for small th or de
  # modified Newton-Raphson
  diff=1; iter=0;
  while((abs(diff/y)> eps) & iter<mxiter)
  { sm=x+y+1; smd=de1*log(sm);
    G21=((1.+smd)^(-th1-1.))/sm; 
    gpdf=-G21; 
    gpdf=gpdf/(1.+smd)/sm/de/th; 
    gpdf=gpdf*(1.+th+th*de*(1.+smd));
    iter=iter+1;
    diff=(G21-pden)/gpdf;
    y=y-diff;
    while(y<=0.) { diff=diff/2.; y=y+diff; }
    if(iprint) cat(iter, y, diff, "\n");
  }
  v=(1.+de1*log(y+1.))^(-th1);
  if(iprint) cat(v, "ended at iter. ", iter, "\n");
  v
}

#============================================================

# original alpha -> 1/ga to get increasing in concordance
# cpar = copula parameter with th>1, ga>=0
pcondbb9=function(v,u,cpar)
{ th=cpar[1]; ga=cpar[2]; ga1=1/ga 
  x= ga1-log(u); y= ga1-log(v);
  temx=x^th; temy=y^th; sm=temx+temy-ga1^th; smt=sm^(1./th);
  ccdf=exp(-smt+ga1);
  ccdf=ccdf*smt*temx/sm/x/u;
  ccdf
}

# 0<p<1 
# cpar = copula parameter with th>1, ga>=0
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb9=function(p,u,cpar, eps=1.e-6, mxiter=30, iprint=F)
{ th=cpar[1]; ga=cpar[2];  ga1=1/ga
  x=ga1-log(u); th1=th-1.;
  con=log(p)-x-th1*log(x); 
  z=(2.*x^th-ga1^th)^(1./th); 
  mxdif=1; iter=0; 
  diff=.1;  # needed in case first step leads to NaN??
  while(mxdif>eps & iter<mxiter)
  { h=z+th1*log(z)+con;
    hp=1.+th1/z;
    if(is.nan(h) | is.nan(hp) | is.nan(h/hp) ) { diff=diff/-2.; } # added for th>50
    else diff=h/hp;
    z=z-diff; iter=iter+1;
    while(z<=x) { diff=diff/2.; z=z+diff; }
    if(iprint) cat(iter, diff, z, "\n");
    mxdif=max(abs(diff));
  }
  if(iprint)
  { if(iter>=mxiter) 
    { cat("***did not converge ");
      cat("p, x, theta, ga, lastz :", p,x,th,ga,z, "\n");
    }
  }
  y=(z^th-x^th+ga1^th)^(1./th);
  vv=exp(-y+ga1);
  if(iprint) cat("p u v", p,u,vv,"\n");
  vv
}

#============================================================

# C_{2|1}(v|u) for BB3 : based on derivative of transform variables
# cpar = copula parameter with th>1, de>0
# iprint = print flag for intermediate calculations
pcondbb3=function(v,u,cpar,iprint=F)
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
    ccdf=cdf*tem*lx/sml/sm/ul/u/dt;   
  }
  else if(x<=1.e-10 || y<=1.e-10) 
  { xx=de*ut; yy=vt*de; r=vt/ut;
    sm=xx+yy; 
    tem=(1+r)^(th1-1);
    cdf=exp(-(sm^th1)/dt);
    ccdf=cdf*tem*(1+xx)/(1+xx+yy)/u;
  }
  else
  { sm=x+y+1.; sml=log(sm);
    tem=(sml^th1);
    cdf=exp(-tem/dt);
    ccdf=cdf*tem*de*ut*(x+1)/sml/sm/ul/u/dt;
  }
  ccdf
}

# 0<p<1 
# cpar = copula parameter with (th,de) : th>1, de>0
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb3=function(p,u,cpar, eps=1.e-6, mxiter=30, iprint=F)
{ th=cpar[1]; de=cpar[2]
  de1=1./de; th1=1./th; dt=(de^th1);
  ul=-log(u); ut=(ul^th); 
  x=exp(ut*de)-1.; 
  # *** add boundary case later for small th or de
  if(is.infinite(x) | x>1.e200)
  { if(iprint) cat("\n** infinite x " , p,u,th,de,ut,"\n");
    lx=ut*de; llx=log(lx); lxt=(lx^th1);
    con=log(p) -lx +(th1-1.)*llx -lxt/dt;
    r=1.; diff=1; iter=0;
    while(abs(diff)>eps  & iter<mxiter)
    { sm=r+1.; sml=log(sm)+lx; smlt=(sml^th1); 
      h=sml +(1.-th1)*log(sml) +smlt/dt +con;
      hp= 1./sm +(1.-th1)/sml/sm + th1*smlt/sml/sm/dt;
      iter=iter+1;
      diff=h/hp; r=r-diff;
      while(r<=0.) { diff=diff/2.; r=r+diff; }
      if(iprint) cat(iter, r, diff, "\n");
    }
    v=((log(r)+lx)^th1)/dt; v=exp(-v);
    if(iprint) cat("** infinite x ",   p,u,th,de,ut,v, "\n")
    return(v);
  }

  if(x<=1.e-10)
  { if(iprint) cat("\n** x near 0 " , p,u,th,de,ut,"\n");
    xx=ut*de; xxt=(xx^th1);
    con=log(p) -xxt/dt;
    r=1.; diff=1; iter=0;
    while(abs(diff)>eps  & iter<mxiter)
    { sm=r+1.; sml=log(sm); smt=(sm^th1); 
      h=xx*r +(1.-th1)*sml +smt*xxt/dt +con;
      hp= xx +(1.-th1)/sm + th1*smt*xxt/sm/dt;
      iter=iter+1;
      diff=h/hp; r=r-diff;
      while(r<=0.) { diff=diff/2.; r=r+diff; }
      if(iprint) cat(iter, r, diff, "\n");
    }
    v=(r*xx/de)^th1; v=exp(-v);
    if(iprint) cat("** x near 0 ",   p,u,th,de,ut,v, "\n")
    return(v);
  }

  lx=ut*de; 
  llx=log(lx);  lxt=(lx^th1);
  con=log(p) -lx +(th1-1.)*llx -lxt/dt;
  # starting guess for y 
  y=x;
  if(y<=0.) y=0.1;
  if(y>1000.) y=1000.;
  if(iprint) cat("p, u ", p,u, "\n");
  if(iprint) cat("x con and starting y: ", x,con,y, "\n");
  # *** add boundary case later for small th or de
  diff=1; iter=0;
  while((abs(diff/y)> eps) & iter<mxiter)
  { sm=x+y+1; sml=log(sm); smlt=(sml^th1);
    h=sml +(1.-th1)*log(sml) +smlt/dt + con;
    hp= 1./sm +(1.-th1)/sml/sm + th1*smlt/sml/sm/dt;
    iter=iter+1;
    diff=h/hp;
    y=y-diff;
    while(y<=0.) { diff=diff/2.; y=y+diff; }
    if(iprint) cat(iter, y, diff, "\n");
  }
  v=(de1*log(y+1.))^th1; v=exp(-v);
  if(iter>=mxiter) 
  { cat("** did not converge **\n");
    cat("p,u,th,de,lasty,v: ", p,u,th,de,y,v,"\n");
  }
  if(iprint) cat("v and #iter :", v,iter, "\n");
  v
}

#============================================================

# C_{2|1}(v|u) for BB4 : based on derivative of transform variables
# cpar = copula parameter with (th,de) : th>1, de>1
#   cpar can be a matrix with #row=length(u)
pcondbb4=function(v,u,cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  ut=u^(-th)-1; vt=v^(-th)-1
  x=ut^(-de); y=vt^(-de); xy=x+y
  tem=xy^(-1/de)
  ccdf=ut+vt+1-tem
  xtem=ut/x-tem/xy
  ccdf=ccdf^(-1/th-1)*xtem*x/ut*(1+ut)/u
  ccdf
}

# this code needs checking in C with 10^6 cases and better starting point
# 0<p<1 
# cpar = copula parameter with (th,de) : th>1, de>1
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondbb4=function(p,u,cpar, eps=1.e-6,mxiter=30,iprint=F)
{ th=cpar[1]; de=cpar[2]
  de1=1./de; th1=1/th
  be=4*pbb4(.5,.5,cpar)-1
  ut=u^(-th)-1; 
  x=ut^(-de); 
  con=(de1+1)*log(x)-(th+1)*log(u)-log(p); 
  iter=0; diff=1.;
  y=.5*x; # what is good starting point?
  if(de>1.7 | u>0.9 | p>0.9) # switch to bisection method, scalar p,u
  { diff=1.
    v=u;
    v1=0.; v2=1.;
    vold=v;
    if(iprint) cat("\n (p,u)=", p,u,"\n")
    while(diff>eps)
    { ccdf=pcondbb4(v,u,cpar); di=ccdf-p;
      if(di<0) { v1=v; } else { v2=v; }
      diff=v2-v1; vold=v; v=(v1+v2)/2.;
      if(iprint) cat(vold,ccdf,di,"\n")
    }
    return(v)
  }
  if(th<0.2) # one boundary 
  { v=qcondgal(p,u,de); vt=v^(-th)-1; y=vt^(-de) }
  if(be>=0.8) y=x
  while(iter<mxiter & max(abs(diff))>eps)
  { xy=x+y; vt=y^(-de1)
    tem=xy^(-de1)
    sm=ut+vt+1-tem
    xtem=ut/x-tem/xy
    ytem=vt/y-tem/xy
    h=-(th1+1)*log(sm)+log(xtem)+con
    hp=(th1+1)*ytem/sm/de + (1.+de1)*tem/xtem/xy/xy
    diff=h/hp;
    y=y-diff;
    if(iprint) cat(iter, diff, y,"\n")
    while(min(y)<=0. | max(abs(diff))>5) { diff=diff/2.; y=y+diff;}
    iter=iter+1;
  }
  if(iprint & iter>=mxiter) cat("***did not converge\n");
  (1+y^(-de1))^(-th1)
}

