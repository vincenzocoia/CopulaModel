# Functions for conditional cdf for 1-parameter bivariate copula families and t
# bvn = bivariate normal/Gaussian
# pla = Plackett
# frk = Frank
# mtcj = Mardia-Takahasi-Clayton-Cook-Johnson
# joe = Joe/B5
# gum = Gumbel
# gal = Galambos
# hr = Huelser-Reiss
# bvt = bivariate t

# conditional cdfs of bivariate copula C_{2|1}(v,u,cpar) 
# inverse of conditional cdfs C_{2|1}^{-1}(p,u,cpar)
# 0<v<1, 0<u<1, 0<p<1, cpar=copula parameter
# Most functions here should work if u,v,p,cpar are vectors of the same length,
#  or if only one of the three is a vector and the other two are scalars.
# The boundary constraints on the functions are not checked on.

# C_{2|1}(v|u} for bivariate  normal copula
# cpar = copula parameter with -1<cpar<1
pcondbvncop=function(v,u,cpar)
{ val=pnorm((qnorm(v)-cpar*qnorm(u))/sqrt(1-cpar^2))
  val[v <= 0 | u <= 0 | u >= 1]=0
  val[v == 1]=1
  val
}

# C_{2|1}^{-1}(p|v} for bivariate Gaussian copula
# cpar = copula parameter with -1<cpar<1, the correlation parameter
qcondbvncop=function(p,v,cpar)
{ x=qnorm(p)
  y=qnorm(v)
  tem=x*sqrt(1-cpar*cpar)+cpar*y
  pnorm(tem)
}

# alternative names of function (maybe remove later)
#pcondbvn=pcondbvncop
#qcondbvn=qcondbvncop


# C_{2|1}(u|v} for bivariate student t copula
# conditional t is Y2|Y1=y1 ~ t(nu+1,location=rho*y1, sigma(y1))
#   where sigma^2 = (1-rho^2)(nu+y1^2)/(nu+1)
# cpar = copula parameter: 2-vector with -1<rho<1, nu>0
#   or scalar with -1<cpar<1
# df = dfdefault is defined globally if cpar is scalar
pcondbvtcop=function(v,u,cpar,df=dfdefault)
{ #rho=cpar[1]; nu=cpar[2]
  if(length(cpar)==2) { rho=cpar[1]; df=cpar[2] }
  else rho=cpar  # not checking if length(cpar)>2
  u[u==0]=1.e-5  # temporary fix
  v[v==0]=1.e-6
  u[u==1]=1-1.e-5  # temporary fix
  v[v==1]=1-1.e-6
  y1=qt(u,df);
  y2=qt(v,df);
  mu=rho*y1; s2=(1.-rho*rho)*(df+y1*y1)/(df+1.);
  ccdf=pt((y2-mu)/sqrt(s2),df+1.);
  ccdf;
}

# cpar=(rho,df) or rho where -1<rho<1, df>0
#pcondt = function(u,v,cpar,df=dfdefault)
#{ if(length(cpar)==2) { rho=cpar[1]; df=cpar[2] }
#  else rho=cpar  # not checking if length(cpar)>2
#  u[u==0]=1.e-5  # temporary fix
#  v[v==0]=1.e-6
#  u[u==1]=1-1.e-5  # temporary fix
#  v[v==1]=1-1.e-6
#  x=qt(u,df); y=qt(v,df)
#  qq=(x-rho*y)/sqrt((df+y*y)*(1-rho*rho)/(df+1))
#  val=pt(qq,df+1)
#  val[u<=0]=0
#  val[u==1]=1
#  val
#}

#pcondbvt = pcondbvtcop
pcondt = pcondbvtcop

# C_{2|1}^{-1}(p|u} for bivariate student t copula
# 0<p<1, could be a vector
# cpar = copula parameter: 2-vector with -1<rho<1, nu>0
#   or scalar with -1<cpar<1
# df = dfdefault is defined globally if cpar is scalar
qcondbvtcop=function(p,u,cpar,df=dfdefault)
{ if(length(cpar)==2) { rho=cpar[1]; df=cpar[2] }
  else rho=cpar  # not checking if length(cpar)>2
  x2=qt(p,df+1)
  x1=qt(u,df)
  tem=x2*sqrt((df+x1*x1)*(1-rho*rho)/(df+1))+rho*x1
  pt(tem,df)
}

#qcondbvt = qcondbvtcop
qcondt = qcondbvtcop

#============================================================

# Plackett
# cpar = copula parameter >0
# cpar==1 will not work for this function?
pcondpla=function(v,u,cpar)
{ #if(cpar==1.) return(v)
  cpar[cpar==1]=1.+1.e-10
  eta=cpar-1.;
  tem=1.+eta*(u+v); tem1=tem*tem-4.*cpar*eta*u*v;
  tem2=sqrt(tem1);
  ccdf=(eta*u+1.-(eta+2.)*v)/tem2;
  ccdf=.5*(1.-ccdf);
  #ifelse(cpar==1., v, ccdf)
  ccdf
}

# Frank
# cpar = copula parameter: cpar>0 or cpar<0 (latter for negative dependence)
pcondfrk=function(v,u,cpar)
{ #if(cpar==0.) return(v)
  cpar[cpar==0]=1.e-10
  cpar1=1.-exp(-cpar);
  tem=1.-exp(-cpar*u);
  ccdf=(1.-tem)/(cpar1/(1.-exp(-cpar*v))-tem);
  ccdf
}

# MTCJ
# cpar = copula parameter >0
pcondmtcj=function(v,u,cpar)
{ tem=v^(-cpar)-1
  tem=tem*(u^cpar)+1
  ccdf=tem^(-1-1/cpar)
  ccdf
}

# reflected MTCJ
# cpar = copula parameter >0
pcondmtcjr=function(v,u,cpar)
{ eta=-cpar/(1+cpar)
  tem=(1-v)^(-cpar)-1.;
  tem=1.+((1-u)^cpar)*tem;
  1-tem^(1./eta)
}

# Joe/B5
# cpar = copula parameter >1
pcondjoe=function(v,u,cpar)
{ temv=(1.-v)^cpar;
  temu=(1.-u)^cpar;
  ccdf=1.+temv/temu-temv;
  ccdf=ccdf^(-1.+1./cpar);
  ccdf=ccdf*(1.-temv);
  ccdf
}

# Gumbel
# cpar = copula parameter >1
pcondgum=function(v,u,cpar)
{ u[u<=0]=1.e-7  # temporary fix (to be consistent with pgum)
  v[v<=0]=1.e-7
  x= -log(u); y= -log(v);
  tem1=x^cpar; tem2=y^cpar; sum=tem1+tem2; tem=sum^(1./cpar);
  ccdf=exp(-tem);
  ccdf=ccdf*(1+tem2/tem1)^(-1.+1./cpar);
  ccdf=ccdf/u;
  ccdf
}

# reflected Gumbel
# cpar = copula parameter >1
pcondgumr=function(v,u,cpar)
{ eta=-1+1/cpar
  u[u>=1]=1-1.e-7  # temporary fix (to be consistent with pgumr)
  v[v>=1]=1-1.e-7
  x=-log(1-u); y=-log(1-v)
  ud=x^cpar; vd=y^cpar
  tem1=(ud+vd)^(1/cpar)
  tem2=(1+vd/ud)^eta
  1-exp(-tem1)*tem2/(1-u)
}

# Galambos
# cpar = copula parameter >0
pcondgal=function(v,u,cpar)
{ x= -log(u); y= -log(v);
  tem1=x^(-cpar); tem2=y^(-cpar); sm=tem1+tem2; tem=sm^(-1./cpar);
  ccdf=exp(-(x+y-tem));
  ccdf=ccdf*(1.-(1+tem2/tem1)^(-1.-1./cpar));
  ccdf=ccdf/u;
  ccdf
}

# Huesler-Reiss
# cpar = copula parameter >0
pcondhr=function(v,u,cpar)
{ x= -log(u); y= -log(v);
  z=x/y; cpar1=1./cpar; lz=log(z);
  tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;  
  p1=pnorm(tem1); p2=pnorm(tem2);
  lcdf=-x*p1-y*p2;  cdf=exp(lcdf);
  ccdf=cdf*p1/u
  ccdf
}

# ============================================================
# qcondcop functions
# conditional quantile functions of parametric bivariate copula families
# C_{2|1}^{-1}(p|u;cpar)
# 0<p<1, 0<u<1, cpar=copula parameter
# p, u can be vectors of same length

# Plackett: C_{2|1}^{-1}(p|u;cpar) 
# cpar = copula parameter >0, must be scalar
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# Output : conditional quantile
qcondpla=function(p,u,cpar, eps=1.e-8,mxiter=30)
{ #if(cpar==1.) return(p);
  iter=0; diff=1.;
  v=u;
  #while(iter<mxiter & abs(diff)>eps)
  while(iter<mxiter & max(abs(diff))>eps)
  { num=pcondpla(v,u,cpar)-p; 
    den=dpla(u,v,cpar);
    diff=num/den;
    v=v-diff;
    #while(v<0. || v>1.) { diff=diff/2.; v=v+diff;}
    while(min(v)<0. | max(v)>1.) { diff=diff/2.; v=v+diff;}
    iter=iter+1;
  }
  v
}

# Frank
# cpar = copula parameter: cpar>0 or cpar<0; cpar=0 input will not work
# 1-exp(-cpar) becomes 1 in double precision for cpar>37.4
qcondfrk=function(p,u,cpar)
{ #if(cpar==0) return(p)
  cpar0=exp(-cpar)
  cpar1=1-cpar0
  #tem=1.-cpar1/((1./p-1.)*exp(-cpar*u)+1.);
  etem=exp(-cpar*u+log(1./p-1.))   
  tem=1.-cpar1/(etem+1.);
  v=(-log(tem))/cpar
  isinf=is.infinite(v)
  #print(cbind(v[isinf],tem[isinf],etem[isinf]))
  # v Inf, tem is 0 and etem < 1.e-16
  v[isinf]=(-log(cpar0+etem[isinf]))/cpar
  v
} 

# MTCJ
# cpar = copula parameter: cpar>0 
#   cpar=0 input will not work for this function
qcondmtcj=function(p,u,cpar)
{ #if(cpar==0) return(p)
  eta=-cpar/(1+cpar)
  tem=p^eta-1
  tem=tem*(u^(-cpar))+1
  #ifelse(cpar==0, p, tem^(-1/cpar))
  tem^(-1/cpar)
}

# reflected MTCJ
# cpar = copula parameter: cpar>0 
#   cpar=0 input will not work for this function
qcondmtcjr=function(p,u,cpar)
{ eta=-cpar/(1+cpar)
  tem=(1-p)^eta-1
  tem=tem*((1-u)^(-cpar))+1
  1-tem^(-1/cpar)
}

# Joe/B5
# solve C_{2|1}(v|u)=p for v given u,p 
# cpar = copula parameter >1
# eps = tolerance for convergence in Newton-Raphson iterations
# mxiter = maximum number of iterations
qcondjoe=function(p,u,cpar, eps=1.e-6,mxiter=30)
{ cpar1=cpar-1;  
  cpartem=-cpar1/(1.+cpar1); cpar1inv=-1./cpar1;
  ubar = 1.0-u; ud = ubar^cpar;
  cpari = 1./cpar;
  # Good starting point based on reflected B4 (June 2011)
  # A good starting point is crucial when cpar is large because
  # C_{2|1} will be steep
  tem=(1.-p)^cpartem-1.;
  tem=tem*(1.-u)^(-cpar1)+1.;
  v=tem^cpar1inv; v=1.-v;
  diff=1; iter=0;
  while(max(abs(diff))>eps & iter<mxiter)
  { vbar = 1.-v;
    vd = vbar^cpar;
    sm=ud+vd-ud*vd;
    smd = sm^cpari;
    c21=1.+vd/ud-vd;
    c21=c21^(-1.+cpari);
    c21=c21*(1.-vd);
    pdf=smd*ud*vd*(cpar1+sm)/sm/sm/ubar/vbar;
    iter=iter+1;
    if(any(is.nan(pdf)) | any(is.nan(c21)) ) { diff=diff/(-2.); } # for cpar>=30
    else diff=(c21-p)/pdf;
    v=v-diff;
    while(min(v)<=0 | max(v)>=1 | max(abs(diff))>0.25) 
    { diff=diff/2.; v=v+diff; }
  }
  v
}

# Gumbel
# vectorized version added on 130331
# G(x,y)=exp(-(x^d+y^cpar)^(1/cpar)), cpar>=1
# conditional 
#    P(Y>=y|x)= exp(x)*G(x,y)* (x^cpar+y^cpar)^(1/cpar-1) * x^(cpar-1)
#                  = p
# take logs, and let z=(x^cpar+y^cpar)^(1/cpar) >= x
# g(z) = log(p) - x + z + (cpar-1)*log(z) -(cpar-1)*log(x) =0
# g'(z) = 1 - (cpar-1)/z
# solve for z with NR, then solve for y=y(z,x,p,cpar)
# y = (z^cpar - x^cpar)^(1/cpar)
# x=-log(u1); y=-log(u2)
# cpar = copula parameter >1
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# Output : conditional quantile
qcondgum=function(p,u,cpar, eps=1.e-6, mxiter=30)
{ if(length(p)>1 | length(u)>1 | length(cpar)>1)
  { # later add some checks to make p,u,cpar in range?
    np=length(p)
    nu=length(u)
    ncpar=length(cpar)
    nn=max(np,nu,ncpar)
    if(np<nn & np>1) return(NA)
    if(nu<nn & nu>1) return(NA)
    if(ncpar<nn & ncpar>1) return(NA)
    if(nn>1)
    { if(np==1) p=rep(p,nn)
      if(nu==1) u=rep(u,nn)
      if(ncpar==1) cpar=rep(cpar,nn)
    }
    out= .C("qcgum", as.integer(nn), as.double(p), as.double(u),
         as.double(cpar), v=as.double(rep(0,nn)) )
    return(out$v)
  }
  else
  { x=-log(u); cpar1=cpar-1.;
    con=log(p)-x-cpar1*log(x);
    z=x*(2.^(1./cpar));
    mxdif=1; iter=0;
    while(mxdif>eps & iter<mxiter)
    { g=z+cpar1*log(z)+con;
      gp=1.+cpar1/z;
      diff=g/gp;
      z=z-diff; iter=iter+1;
      while(z<=x) { diff=diff/2.; z=z+diff; }
      mxdif=abs(diff);
    }
    y=((z^cpar)-(x^cpar))^(1./cpar);
    return(exp(-y))
  }
}

# reflected Gumbel 
# cpar = copula parameter >1
qcondgumr=function(p,u,cpar, eps=1.e-6, mxiter=30)
{ v1=qcondgum(1-p,1-u,cpar,eps,mxiter)
  1-v1
}

# Galambos 
# condgal Ok for de<.0001, minimal changes from de=.04 and smaller 
# G(x,y)=exp(-x-y+(x^dn+y^dn)^(1/dn)), de>=0, dn=-de
#   conditional P(Y>=y|x)= p
#     = exp(x)*G(x,y)*[1-x^(dn-1)*(x^dn+y^dn)^(1/dn-1)]
#     = exp(-y+z)*[1-x^(dn-1)*(x^dn+y^dn)^(1/dn-1)]
#     = exp(-y+z)* x^(dn-1)* [x^(de+1)-(x^dn+y^dn)^(1/dn-1)]
#   take logs, and let z=(x^dn+y^dn)^(1/dn)
#   g(y) = log(p) + y -(x^dn+y^dn)^(1/dn)  +
#           - log[1-x^(dn-1)*(x^dn+y^dn)^(1/dn-1)]
# cpar = copula parameter >0
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondgal=function(p,u,cpar, eps=1.e-6, mxiter=30,iprint=F)
{ if(length(p)>1 | length(u)>1 | length(cpar)>1)
  { # later add some checks to make p,u,cpar in range?
    np=length(p)
    nu=length(u)
    ncpar=length(cpar)
    nn=max(np,nu,ncpar)
    if(np<nn & np>1) return(NA)
    if(nu<nn & nu>1) return(NA)
    if(ncpar<nn & ncpar>1) return(NA)
    if(nn>1)
    { if(np==1) p=rep(p,nn)
      if(nu==1) u=rep(u,nn)
      if(ncpar==1) cpar=rep(cpar,nn)
    }
    out= .C("qcgal", as.integer(nn), as.double(p), as.double(u),
         as.double(cpar), v=as.double(rep(0,nn)) )
    return(out$v)
  }
  else
  { x=-log(u); 
    con=log(p); dn=-cpar; dn1=1./dn;
    xd=x^dn;
    y=(x*cpar-con)/(1+cpar);
    if(cpar>2.7) y=x
    # different equation r=[(-log(u))^de]/[(-log(v))^de]; 
    # or y=x * r^(1/de) when de is large and u is close to 1
    if(xd>1.e6) # empirical
    { con=log(p)
      r=1
      #if(iprint) cat("entering loop solving for r=(y/x)^de\n")
      mxdif=1; iter=0;
      while(mxdif>eps & iter<mxiter)
      { rd=r^(1/cpar); r1=1+1/r; r1d=r1^dn1
        h=con+rd*x - x*r1d- log(1-r1d/r1)
        hp=rd*x/r/cpar - x*r1d/r1/cpar/r/r + (1+1/cpar)*r1d/r1/r1/r/r/(1-r1d/r1)
        dif=h/hp;
        if(iprint) cat(iter, dif, r,"\n")
        r=r-dif; iter=iter+1;
        while(r<=0.) { dif=dif/2.; r=r+dif; }
        mxdif=abs(dif);
      }
      y=x* r^(1/cpar)
      return(exp(-y))
    }
    mxdif=1; iter=0;
    while(mxdif>eps & iter<mxiter)
    { yd=y^dn; z=(xd+yd)^dn1;
      tem=z^(cpar+1.);
      den=1-tem*xd/x;
      g=con-z+y-log(den);
      gp=1.-tem*yd/y+(1.+cpar)*tem*xd*yd/(den*x*y*(xd+yd));
      dif=g/gp;
      y=y-dif; iter=iter+1;
      while(y<=0.) { dif=dif/2.; y=y+dif; }
      mxdif=abs(dif);
    }
    #if(iter>=mxiter) cat("***did not converge\n");
    return(exp(-y))
  }
}

# Huesler-Reiss
# cpar = copula parameter >0
# eps = tolerance for convergence
# mxiter = maximum number of Newton-Raphson iterations
# iprint = print flag for intermediate calculations
# Output : conditional quantile
qcondhr=function(p,u,cpar, eps=1.e-5,mxiter=30,iprint=F)
{ if(length(p)>1 | length(u)>1 | length(cpar)>1)
  { # later add some checks to make p,u,cpar in range?
    np=length(p)
    nu=length(u)
    ncpar=length(cpar)
    nn=max(np,nu,ncpar)
    if(np<nn & np>1) return(NA)
    if(nu<nn & nu>1) return(NA)
    if(ncpar<nn & ncpar>1) return(NA)
    if(nn>1)
    { if(np==1) p=rep(p,nn)
      if(nu==1) u=rep(u,nn)
      if(ncpar==1) cpar=rep(cpar,nn)
    }
    out= .C("qchr", as.integer(nn), as.double(p), as.double(u),
         as.double(cpar), v=as.double(rep(0,nn)) )
    return(out$v)
  }
  else # scalar u and p
  { diff=1.;
    # bisection method for larger cpar
    if(cpar>1.8) # assumes p and u are scalar
    { v=u;
      v1=0.; v2=1.;
      vold=v;
      if(iprint) cat("\n (p,u)=", p,u,"\n")
      while(diff>eps)
      { ccdf=pcondhr(v,u,cpar); di=ccdf-p;
        if(di<0) { v1=v; } else { v2=v; }
        diff=v2-v1; vold=v; v=(v1+v2)/2.;
        if(iprint) cat(vold,ccdf,di,"\n")
      }
      return(v)
    }
    # cpar<=1.8
    con=-log(p); cpar1=1./cpar
    iter=0; 
    x= -log(u); 
    y=.5*x; 
    while(iter<mxiter & abs(diff)>eps)
    #while(iter<mxiter & max(abs(diff))>eps)
    { z=x/y; lz=log(z)
      tem1=cpar1+.5*cpar*lz; tem2=cpar1-.5*cpar*lz;
      p1=pnorm(tem1); p2=pnorm(tem2);
      g=x*(1-p1)-y*p2+log(p1)+con
      gp= -p2-0.5*dnorm(tem1)/p1/y
      diff=g/gp;
      y=y-diff;
      if(iprint) cat(iter, diff, y,"\n")
      while(y<=0.) { diff=diff/2.; y=y+diff;}
      #while(min(y)<=0.) { diff=diff/2.; y=y+diff;}
      iter=iter+1;
    }
    if(iprint & iter>=mxiter) cat("***did not converge\n");
    return(exp(-y))
  }
}



