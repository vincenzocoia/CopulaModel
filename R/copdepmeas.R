# Functions for copula parameter (cpar) to/from rhoS, tau, beta=b
# and lm = lambda = tail dependence parameter(s) 
# for bivariate copula families.

# cpar is abbreviation for copula parameter

# Spearman rho for Plackett copula
# In book of Hutchinson and Lai (originally due to Mardia) 
# rho= (cpar+1)/(cpar-1) - (2*cpar*log(cpar))/(cpar-1)^2
# rho'(cpar)= 1/(cpar-1) - [(cpar+1)+2*cpar*log(cpar)+2]/(cpar-1)^2
#        +4*cpar*log(cpar)/(cpar-1)^3
# cpar = copula parameter >0 except cpar=1
pla.cpar2rhoS=function(cpar)
{ ifelse(cpar==1, 0, (cpar+1)/(cpar-1) - (2*cpar*log(cpar))/(cpar-1)^2 ) }

# copula parameter given rhoS for Plackett
# vectorized input rhoS is OK
# rhoS = value of Spearman's rho
# cpar0 = starting point for Newton-Raphson iterations
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for iterations
# mxstep = bound on step size for Newton-Raphson iterations
pla.rhoS2cpar=function(rhoS, cpar0=1.5,mxiter=25,eps=1.e-6,iprint=F,mxstep=10)
{ cpar=cpar0
  iter=0
  diff=1
  while(iter<mxiter & max(abs(diff))>eps)
  { cpar1=cpar-1
    cpar2=cpar1*cpar1; lgcpar=log(cpar)
    g=(cpar+1)/cpar1 - (2*cpar*lgcpar)/cpar2 -rhoS
    gp = 1/cpar1- ((cpar+1)+2*lgcpar+2)/cpar2 + 4*cpar*lgcpar/(cpar1*cpar2)
    diff=g/gp; cpar=cpar-diff
    iter=iter+1
    while(max(abs(diff))>mxstep | min(cpar)<=0) { diff=diff/2; cpar=cpar+diff }
    if(iprint) cat("iter=", iter, " cpar=", cpar, " ",g, " ",gp, "\n")
  }
  if(iter>=mxiter) { cat("did not converge\n") }
  cpar
}

# integrans for Kendall's tau for Frank copula
debye= function(t) { t/(exp(t)-1) }
debye2= function(t) { t^2/(exp(t)-1) }
frk.cpar2tau=function(cpar)
{ tem=integrate(debye,0,cpar,rel.tol=1.e-6)$value
  tem=tem/cpar
  tau=1+4*(tem-1)/cpar
  tau
}

# reference is Nelsen (1986). Comm Stat, 15(11), 3277-3285
# cpar = copula parameter >0 or <0
frk.cpar2rhoS=function(cpar)
{ tem=integrate(debye,0,cpar,rel.tol=1.e-6)$value
  tem2=integrate(debye2,0,cpar,rel.tol=1.e-6)$value
  tem=tem/cpar  # debye1
  tem2=tem2*2/cpar^2  # debye2
  rho=1+12*(tem2-tem)/cpar
  rho
}

# Kendall's tau for Gumbel is tau=(cpar-1)/cpar, for MTCJ is tau=cpar/(cpar+2)
# (also inverse functions)
# cpar = copula parameter >1, 0<tau<1
gum.cpar2tau=function(cpar) { (cpar-1)/cpar }
gum.tau2cpar=function(tau) { 1/(1-tau) }

# cpar = copula parameter >0, 0<tau<1
mtcj.cpar2tau=function(cpar) { cpar/(cpar+2) }
mtcj.tau2cpar=function(tau) { 2*tau/(1-tau) }

# function in integrand for Kendall's tau for Archimedean copula with Sibuya LT
# s* [psi'(s)]^2 for Sibuya LT
# s = positive value
# cpar = copula parameter >1
sibuyader= function(s,cpar)
{ es=exp(-s)
  psid=(1-es)^(1/cpar-1)
  psid=psid*es/cpar
  s*psid^2
}

# old version with integration has been replaced
#joe.cpar2tau= function(cpar)
#{ tem=integrate(sibuyader,0,Inf,cpar=cpar,rel.tol=1.e-6)
#  tem=tem$value
#  tau=1-4*tem
#  tau
#}

# simpler version of Kendall's tau for Joe/B5 copula
# also simpler than formula in Schepsmeier's Master's thesis
# cpar = copula parameter >1
joe.cpar2tau=function(cpar)
{ cpar1=2/cpar+1
  tem=digamma(2)-digamma(cpar1)
  tau=1+tem*2/(2-cpar)
  tau[cpar==2]=1-trigamma(2)
  tau
}

# cpar = copula parameter >1
gum.cpar2rhoS=function(cpar)
{ spgum= function(w,cpar)
  { w1=1-w
    wd=w^cpar; w1d=w1^cpar
    sm=wd+w1d
    B=sm^(1/cpar)
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(spgum,0,1,cpar=cpar,rel.tol=1.e-6)
  rho=12*tem$value-3
  rho
}

# cpar = copula parameter >0
gal.cpar2tau=function(cpar)
{ # B(w) = 1- (w^cpar1+w1^cpar1)^(1/cpar1),   cpar1=-cpar
  # B'(w) =  (w^cpar1+w1^cpar1)^(1/cpar1-1) (-w^(cpar1-1) + w1^(cpar1-1) )
  taugal= function(w,cpar)
  { w1=1-w; cpar1=-cpar
    wd=w^cpar1; w1d=w1^cpar1
    sm=wd+w1d
    B=1-sm^(1/cpar1)
    Bp=(1+(w1/w)^cpar)^(1/cpar1-1) - (1+(w/w1)^cpar)^(1/cpar1-1)
    ((2*w-1)*Bp*B + w*w1*Bp*Bp )/(B*B)
  }
  tem=integrate(taugal,0,1,cpar=cpar,rel.tol=1.e-6)
  tau=tem$value
  tau
}

# cpar = copula parameter >0
gal.cpar2rhoS=function(cpar)
{ spgal= function(w,cpar)
  { w1=1-w; cpar1=-cpar
    wd=w^cpar1; w1d=w1^cpar1
    sm=wd+w1d
    B=1-sm^(1/cpar1)
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(spgal,0,1,cpar=cpar,rel.tol=1.e-6)
  rho=12*tem$value-3
  rho
}

# cpar = copula parameter >0
hr.cpar2tau=function(cpar)
{ # B(w) =w\Phi\Bigl(\cpar^{-1}+\half\cpar \log {w\over 1-w}\Bigr)
  # +(1-w)\Phi\Bigl(\cpar^{-1}+\half\cpar \log {1-w\over w}\Bigr)
  # B'(w)=\Phi\Bigl(\cpar^{-1}+\half\cpar \log {w\over 1-w}\Bigr)
  #   -\Phi\Bigl(\cpar^{-1}+\half\cpar \log {1-w\over w}\Bigr)
  tauhr= function(w,cpar)
  { w1=1-w
    cpar1=1/cpar
    rat=w/w1; lograt=log(rat)
    tem1=pnorm(cpar1+.5*cpar*lograt); tem2=pnorm(cpar1-.5*cpar*lograt)
    B=w*tem1+w1*tem2
    Bp=tem1-tem2
    ((2*w-1)*Bp*B + w*w1*Bp*Bp )/(B*B)
  }
  tem=integrate(tauhr,0,1,cpar=cpar,rel.tol=1.e-6)
  tau=tem$value
  tau
}

# cpar = copula parameter >0
hr.cpar2rhoS=function(cpar)
{ sphr= function(w,cpar)
  { w1=1-w; cpar1=1/cpar
    rat=w/w1; lograt=log(rat)
    tem1=pnorm(cpar1+.5*cpar*lograt); tem2=pnorm(cpar1-.5*cpar*lograt)
    B=w*tem1+w1*tem2
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(sphr,0,1,cpar=cpar,rel.tol=1.e-6)
  rho=12*tem$value-3
  rho
}

# lm = tail dependence parameter lambda in (0,1), possibly a vector
joe.cpar2lm=function(cpar) { 2-2^(1/cpar) }
joe.lm2cpar=function(lm) { log(2)/log(2-lm) }
gum.cpar2lm=function(cpar) { 2-2^(1/cpar) }
gum.lm2cpar=function(lm) { log(2)/log(2-lm) }
mtcj.cpar2lm=function(cpar) { 2^(-1/cpar) }
mtcj.lm2cpar=function(lm) { log(2)/(-log(lm)) }
gal.cpar2lm=function(cpar) { 2^(-1/cpar) }
gal.lm2cpar=function(lm) { log(2)/(-log(lm)) }
hr.cpar2lm=function(cpar) { 2*(1-pnorm(1/cpar)) }
hr.lm2cpar=function(lm) { 1/qnorm(1-lm/2) }

#============================================================
# Blomqvist beta for bivariate 1-parameter copulas, and  bivariate t, BB1
# given beta, find copula parameter

# BVN and bivariate t(nu)
# rho=sin(.5*pi*b), -1<b<1, -1<rho<1, -1<rhoS<1
bvn.b2cpar=function(b) { sin(.5*pi*b) }
bvn.tau2cpar=bvn.b2cpar
bvn.cpar2b=function(rho) { (2/pi)*asin(rho) }
bvn.cpar2tau=bvn.cpar2b
# Spearman is different for biv t
bvn.cpar2rhoS=function(rho) { (6/pi)*asin(rho/2) }
bvn.rhoS2cpar=function(rhoS) { 2*sin(pi*rhoS/6) }

# cpar = copula parameter (rho,nu)
bvt.cpar2b=function(cpar) 
{ if(is.matrix(cpar)) { rho=cpar[,1] } else { rho=cpar[1] }
  (2/pi)*asin(rho) 
}  
bvt.cpar2tau=bvt.cpar2b
bvt.b2cpar=function(b,nu) { sin(.5*pi*b)  }
bvt.tau2cpar=bvt.b2cpar

bvt.cpar2lm=function(cpar)
{ if(is.matrix(cpar)) { rho=cpar[,1]; nu=cpar[,2] }
  else { rho=cpar[1]; nu=cpar[2] }
  2*pt(-sqrt((nu+1)*(1-rho)/(1+rho)),nu+1) 
}

# lm = tail dependent lamdba  0<lm<1
# nu = degree of freedom parameter >0
bvt.lm2cpar=function(lm,nu)
{ tem=-qt(lm/2,1+nu)
  tem=tem^2/(1+nu)
  rho=(1-tem)/(1+tem)
  rho
}

# b = value of Blomqvist's beta, -1<b<1
pla.b2cpar=function(b) { ((1+b)/(1-b))^2 }

# Frank copula, solve equation to get cpar given Blomqvist's beta b
# vectorized input b is OK, b=0 fails
# b = value of Blomqvist's beta, -1<b<1
# cpar0 = starting point for Newton-Raphson iterations
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for iterations
# Output: copula parameter with given b value
frk.b2cpar=function(b, cpar0=0,mxiter=20,eps=1.e-8,iprint=F)
{ iter=0
  #b[b==0]=.00000001
  b[abs(b)<1.e-8]=.00000001
  diff=1.
  cpar=cpar0
  if(cpar0<=0) cpar=10*b
  ln2=log(2)
  while(iter<mxiter & max(abs(diff))>eps)
  { tem=exp(cpar/2)
    g=log(1+tem)-ln2-cpar*(1+b)/4
    gp=.5/(1+1/tem)-(1+b)/4
    iter=iter+1
    diff=g/gp
    cpar=cpar-diff
    if(iprint) cat(iter," ",cpar," ",diff,"\n")
  }
  if(iter>=mxiter) cat("did not converge\n")
  cpar
}

# MTCJ copula, solve equation to get cpar given Blomqvist's beta b
# vectorized input b is OK
# b = value of Blomqvist's beta, 0<b<1
# cpar0 = starting point for Newton-Raphson iterations
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for iterations
# Output: copula parameter with given b value
mtcj.b2cpar=function(b, cpar0=0,mxiter=20,eps=1.e-8,iprint=F)
{ iter=0
  diff=1.
  #b[b==0]=.00000001
  b[abs(b)<1.e-8]=.00000001
  cpar=cpar0
  if(cpar0<=0) cpar=2*b/(1-b)  # value of cpar with given tau
  ln2=log(2)
  btem=log((1+b)/4)
  while(iter<mxiter & max(abs(diff))>eps)
  { tem=2^(cpar+1)
    g=log(tem-1)+cpar*btem
    gp=ln2*tem/(tem-1)+btem
    iter=iter+1
    diff=g/gp
    cpar=cpar-diff
    if(iprint) cat(iter," ",cpar," ",diff,"\n")
  }
  if(iter>=mxiter) cat("did not converge\n")
  cpar
}

# Joe/B5, solve equation to get cpar given Blomqvist's beta b
# vectorized input b is OK
# b = value of Blomqvist's beta, 0<b<1
# cpar0 = starting point for Newton-Raphson iterations
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for iterations
# Output: copula parameter with given b value
joe.b2cpar=function(b, cpar0=0,mxiter=20,eps=1.e-8,iprint=F)
{ iter=0
  diff=1.
  #b[b==0]=.00000001
  b[abs(b)<1.e-8]=.00000001
  cpar=cpar0
  if(cpar0<=0) cpar=2*b/(1-b)+1  # value of cpar with given tau for b4
  ln2=log(2)
  btem=log((3-b)/2)
  while(iter<mxiter & max(abs(diff))>eps)
  { tem=2^(-cpar)
    g=log(2-tem)-cpar*btem
    gp=ln2*tem/(2-tem)-btem
    iter=iter+1
    diff=g/gp
    cpar=cpar-diff
    while(min(cpar)<=1) { diff=diff/2; cpar=cpar+diff }
    if(iprint) cat(iter," ",cpar," ",diff,"\n")
  }
  if(iter>=mxiter) cat("did not converge\n")
  cpar
}

# Gumbel
# b = value of Blomqvist's beta, 0<b<1
gum.b2cpar=function(b)
{ ln2=log(2)
  ln2/log(2-(log(1+b))/ln2)
}

# Galambos
# b = value of Blomqvist's beta, 0<b<1
gal.b2cpar=function(b)
{ ln2=log(2)
  ln2/log(ln2/log(1+b))
}

# Huesler-Reiss
# b = value of Blomqvist's beta, 0<b<1
hr.b2cpar=function(b)
{ ln2=log(2)
  1/qnorm(1-(log(1+b))/(2*ln2))
}


# BB1 cpar=(th,de)
# b = value of Blomqvist's beta, 0<b<1
# de = delta parameter of copula
# thstart = starting point for theta for Newton-Raphson iterations
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# Output: theta parameter given Blomqvist's beta and copula delta parameter
#   also tail dependence parameters  lml, lmu
bb1.b2cpar=function(b,de,thstart=1,mxiter=30,eps=1.e-6)
{ iter=0; 
  diff=1.;
  ln2=log(2.)
  #cat("beta=", b, " de=", de,"\n")
  lmu=2.-2^(1./de);
  bp=(1.+b)/4.; lbp=log(bp); 
  th=thstart;
  while(iter<mxiter & abs(diff)>eps)
  { tem1=2^th; tem2=bp^(-th); 
    g=2.*tem1-lmu*(tem1-1.)-tem2;
    gp=2.*tem1*ln2-lmu*tem1*ln2+tem2*lbp;
    iter=iter+1;
    diff=g/gp;
    th=th-diff;
    while(th<=0. | abs(diff)>5.) { diff=diff/2.; th=th+diff; }
    #print(c(iter,th,diff))
  }
  if(iter>=mxiter | th<1.e-6) cat("did not converge, use larger thstart\n");
  lml=2^(-1./(th*de))
  list(th=th,lmu=lmu,lml=lml)
}

# cpar = copula parameter with th>0, de>1
bb1.cpar2tau=function(cpar) 
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  1-2/(de*(th+2.)) 
}

# lmpar = tail dependence parameters (lml,lmu) in (0,1)^2 
# Output copula parameter (th,de)
bb1.lm2cpar=function(lmpar)  # order is lml,lmu
{ de=log(2)/log(2-lmpar[2]); th=log(2)/(-de*log(lmpar[1]))
  c(th,de)
}

# cpar = copula parameter with th>0, de>1
bb1.cpar2lm=function(cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  lml=2^(-1/(th*de)); lmu=2-2^(1/de)
  cbind(lml,lmu)
}

# BB1, given 0<tau<1, find theta and delta with lml=lmu
# tau = Kendall tau value
# destart = starting point for delta
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for iterations
# Output: copula parameter (th,de) with lml=lmu anf given tau
bb1.tau2eqlm=function(tau,destart=1.5,mxiter=30,eps=1.e-6,iprint=F)
{ iter=0; 
  diff=1.;
  ln2=log(2.)
  de=destart
  rhs=2/(1-tau)
  while(iter<mxiter & max(abs(diff))>eps)
  { tem=2-2^(1/de); 
    ltem= -log(tem)
    g= ln2/ltem +2*de -rhs
    gp= (ln2/ltem/de)^2 * (2-tem)/tem +2
    iter=iter+1;
    diff=g/gp;
    de=de-diff;
    while(min(de)<=1. | max(abs(diff))>5.) { diff=diff/2.; de=de+diff; }
    if(iprint) print(c(iter,de,diff))
  }
  if(iter>=mxiter) cat("did not converge\n");
  th=ln2/(de*ltem)
  lml=2^(-1./(th*de))
  lmu=2.-2^(1./de)
  if(length(tau)==1)
  { out=c(th,de,lml,lmu)
    names(out)=c("th","de","lml","lmu")
  }
  else
  { out=cbind(th,de,lml,lmu) }
  out
}

# lmpar = tail dependence parameters (lml,lmu) in (0,1)^2 
# Output copula parameter (th,de)
bb7.lm2cpar=function(lmpar)  # order is lml,lmu
{ de=log(2)/(-log(lmpar[1])); th=log(2)/log(2-lmpar[2])
  c(th,de)
}

# cpar = copula parameter with th>1, de>0
bb7.cpar2lm=function(cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  lml=2^(-1/de); lmu=2-2^(1/th)
  cbind(lml,lmu)
}

# lmpar = tail dependence parameters (lml,lmu) in (0,1)^2 
# Output copula parameter (th,de)
bb4.lm2cpar=function(lmpar)  # order is lml,lmu
{ de=log(2)/(-log(lmpar[2])); th=log(2-lmpar[2])/(-log(lmpar[1]))
  c(th,de)
}

# cpar = copula parameter with th>0, de>0
bb4.cpar2lm=function(cpar)
{ if(is.matrix(cpar)) { th=cpar[,1]; de=cpar[,2] }
  else { th=cpar[1]; de=cpar[2] }
  lml=(2-2^(-1/de))^(-1/th); lmu=2^(-1/de)
  cbind(lml,lmu)
}

#============================================================
