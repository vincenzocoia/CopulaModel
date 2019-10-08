# Functions for 
# (a) bivariate dependence measures for copulas including semi-correlations 
# (b) empirical dependence measures for bivariate data  
# (c) mappings from dependence measures to copula parameters 

# rho = correlation in (-1,1)
# output: Cor(Z1,Z2| Z1>0,Z2>0) when (Z1,Z2)~BVN(rho)
bvnsemic=function(rho)
{ bp=.25+asin(rho)/(2*pi)
  r1=sqrt(1-rho^2)
  denom=2*bp*sqrt(2*pi)
  v10=(1+rho)/denom
  v20=1+rho*r1/(2*bp*pi)
  v11=rho+r1/(2*bp*pi)
  sc=(v11-v10^2)/(v20-v10^2)
  sc
}

# semi-correlations (lower and upper) applied to normal scores
# bivdat = nx2 data set
# inscore = T if bivdat has already been converted to normal scores
# Output: a 3-vector with rhoN and lower/ upper semi-correlations
semicor=function(bivdat,inscore=T)
{ if(!inscore) bivdat=nscore(bivdat)
  ii1=(bivdat[,1]<0 & bivdat[,2]<0)
  ii2=(bivdat[,1]>0 & bivdat[,2]>0)
  lcorr=cor(bivdat[ii1,1],bivdat[ii1,2])
  ucorr=cor(bivdat[ii2,1],bivdat[ii2,2])
  ncorr=cor(bivdat[,1],bivdat[,2])
  out=c(ncorr,lcorr,ucorr)
  names(out)=c("ncorr","lcorr","ucorr")
  out
}

# semi-correlation table of normal scores
# mdat = nxd data set with d>=2 columns
# inscore = T if mdat has already been converted to normal scores
# Output: d*(d-1)/2 by 6 matrix with columns j1,j2,ncorr,lcorr,ucorr,bvnsemic
semicortable=function(mdat,inscore=F)
{ if(!inscore) mdat=nscore(mdat)
  d=ncol(mdat)
  d2=(d*(d-1))/2
  out=matrix(0,d2,6)
  ii=0
  for(j2 in 2:d)
  { for(j1 in 1:(j2-1))
    { scs=semicor(mdat[,c(j1,j2)])
      bvnsc=bvnsemic(scs[1])
      ii=ii+1
      out[ii,]=c(j1,j2,scs, bvnsc)
    }
  }
  out=as.data.frame(out)
  names(out)=c("j1","j2","ncorr","lcorr","ucorr","bvnsemic")
  out
}

# Blomqvist beta, efficient computation
# bivdat = nx2 dat set
# iunif = T if bivdat has already been converted to uniform scores
# Output: empirical Blomqvist beta value
blomq=function(bivdat,iunif=F)
{ n=nrow(bivdat)
  udat=uscore(bivdat)
  ii=(((udat[,1]-0.5)*(udat[,2]-0.5))>=0)
  tem=mean(ii)
  2*tem-1
}

# Spearman's rho for bivariate copula with parameter cpar
# cpar = copula parameter
# cop = function name of joint cdf or conditional cdf C_{2|1}
# zero = 0 or something like 1.e-6 (to avoid endpoint problems)
# if icond = T, cop = conditional cdf pcondcop 
# if icond = F, cop = joint cdf pcop 
#   default is to integrate on [zero,1-zero]^2 for icond=T
# tol = accuracy for 2-dimensional integration
# Output: Spearman rho value for copula
rhoS=function(cpar,cop,zero=0,icond=F,tol=.0001)
{ # use of pdf has been eliminated
  if(icond)
  { # cop=pcond 
    lb=zero; ub=1-zero
    spearfn= function(u) 
    { u1=u[1]; u2=u[2];
      u1*cop(u2,u1,cpar)
    }
    spear=adaptIntegrate(spearfn,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
    spear=3-12*spear$integ
  }
  else
  { # cop=pcop
    lb=zero; ub=1-zero
    spearfn= function(u) 
    { u1=u[1]; u2=u[2]; cop(u1,u2,cpar) }
    spear=adaptIntegrate(spearfn,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
    spear=12*spear$integ-3
  }
  spear
}

# Kendall's tau for bivariate copula with parameter cpar
# cpar = copula parameter
# if icond = T,
#  pcond12 = C_{1|2}, pcond21=C_{2|1}, they are same if C is permutation symmetric
#  zero = 0 or something like 1.e-6
# if icond = F,
#  dcop = copula density;
#  pcop = copula cdf 
#  B = upper limit of numerical integration
#  default is to integrate on [-B,B]^2 for icond=F (with split to 4 regions)
# tol = accuracy for 2-dimensional integration
# Output: Kendall tau value for copula
ktau=function(cpar,icond=T,pcond12,pcond21,zero=0,dcop,pcop,B=6,tol=.0001)
{ # option that uses pcond12, pcond21
  if(icond)
  { lb=zero; ub=1-zero
    taufn= function(u) 
    { u1=u[1]; u2=u[2];
      pcond12(u1,u2,cpar)*pcond21(u2,u1,cpar)
    }
    tau=adaptIntegrate(taufn,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
    tau=1-4*tau$integ
  }
  else
  { lb=0; ub=B
    taufn2= function(z) 
    { z1=z[1]; z2=z[2]; u1=pnorm(z1); u2=pnorm(z2)
      tem=dcop(u1,u2,cpar)*pcop(u1,u2,cpar) +
        dcop(1-u1,1-u2,cpar)*pcop(1-u1,1-u2,cpar) +
        dcop(1-u1,u2,cpar)*pcop(1-u1,u2,cpar)  +
        dcop(u1,1-u2,cpar)*pcop(u1,1-u2,cpar) 
      dnorm(z1)*dnorm(z2)*tem
    }
    tau=adaptIntegrate(taufn2,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
    tau=4*tau$integ-1
  }
  tau
}

# rhoN (correlation of normal scores) for bivariate copula 
# compute for pcond, or cdf or density (3 choices)
# cpar = copula parameter
# icond = T if numerical integration via conditional cdf,
# pcond = conditonal cdf C_{2|1}
# icdf = T and icond = F if numerical integration via copula cdf,
# icdf = F and icond = F if numerical integration via copula density
# dcop = copula density
# pcop = copula cdf
# B = upper limit of integration
#   default is to integrate on [0,B]^2, with split into sum of 4 integrands 
# tol = accuracy for 2-dimensional integration
rhoN=function(cpar,icond=T,pcond,icdf=F,pcop,dcop,B=6,tol=.0001)
{ lb=0; ub=B
  # using conditional cdf
  if(icond)
  { nsfn= function(z) 
    { z1=z[1]; z2=z[2]; u1=pnorm(z1); u2=pnorm(z2)
      tem= 1-pcond(u2,u1,cpar)+ pcond(1-u2,1-u1,cpar) -
           pcond(1-u2,u1,cpar) - (1-pcond(u2,1-u1,cpar))
      z1*dnorm(z1)*tem
    }
  }
  # using cdf via Hoeffding's identity
  else if(icdf)
  { nsfn= function(z) 
    { z1=z[1]; z2=z[2]; u1=pnorm(z1); u2=pnorm(z2)
      tem= (pcop(u1,u2,cpar)-u1*u2) + (pcop(1-u1,1-u2,cpar)-(1-u1)*(1-u2)) +
           (pcop(1-u1,u2,cpar)-(1-u1)*u2) + (pcop(u1,1-u2,cpar)-u1*(1-u2))
      tem
    }
  }
  # using pdf 
  else
  { nsfn= function(z) 
    { z1=z[1]; z2=z[2]; u1=pnorm(z1); u2=pnorm(z2)
      tem= dcop(u1,u2,cpar)+ dcop(1-u1,1-u2,cpar) -
           dcop(1-u1,u2,cpar) - dcop(u1,1-u2,cpar)
      z1*z2*dnorm(z1)*dnorm(z2)*tem
    }
  }
  ns=adaptIntegrate(nsfn,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
  ns$integral
}

# semi-correlations of normal scores for bivariate copula 
# cpar = copula parameter
# dcop = copula density
# pcop = copula cdf
# pcond = conditional C_{2|1}(v|u)
# B = upper bound for 2-dimensional numerical integration
# isym = T means reflection symmetric and lower semi-corr is same as upper
# iinfbd = T means use upper bound Inf for 1-dimensional integral. otherwise B
# iprint = print flag for intermediate calculations
# Output: 2-vector with lcorr = rhoN^-, ucorr = rhoN^+
rhoNsemic=function(cpar,dcop,pcop,pcond,B=6,isym=F,iinfbd=T,iprint=F)
{ 
  denom=pcop(.5,.5,cpar)
  # cpar is parameter value of copula with conditional pcond()
  bcopu1= function(z) { z*dnorm(z)*(1.-pcond(.5,pnorm(z),cpar)) }
  bcopu2= function(z) { z*z*dnorm(z)*(1.-pcond(.5,pnorm(z),cpar)) }
  bcopl1= function(z) { z*dnorm(z)*pcond(.5,pnorm(-z),cpar) }
  bcopl2= function(z) { z*z*dnorm(z)*pcond(.5,pnorm(-z),cpar) }
  if(iprint) cat("\n", " cpar=",cpar, " denom=", denom," bound=",B,"\n")
  # next integrand can be used for range -B,B to get ncorr or 0,B for ucorr
  g12u= function(z) 
  { z1=z[1]; z2=z[2]; 
    z1*z2*dnorm(z1)*dnorm(z2)*dcop(pnorm(z1),pnorm(z2),cpar)
  }
  # B=6 might be too large
  B1=B-1
  ex12=adaptIntegrate(g12u,lowerLimit=c(-B1,-B1),upperLimit=c(B1,B1),absError=.0001)
  if(iprint)
  { cat("corr of nscores= ", ex12$integ, "\n")
    cat("upper\n")
  }
  Binf=B; if(iinfbd) { Binf=Inf }
  ex1=integrate(bcopu1,0,Binf)
  #print(ex1)
  ex2=integrate(bcopu2,0,Binf)
  #print(ex2)
  ex12=adaptIntegrate(g12u,lowerLimit=c(0,0),upperLimit=c(B,B),absError=.0001)
  #print(ex12)
  ex1=ex1$value/denom; ex2=ex2$value/denom; ex12=ex12$integ/denom
  varr=ex2-ex1*ex1; ucorr=(ex12-ex1*ex1)/varr
  if(iprint) print(c(ex1,ex2,varr,ex12,ucorr))
  if(!isym)
  { if(iprint) cat("lower\n")
    # this might fail for bivariate t if Binf=Inf
    ex1=integrate(bcopl1,0,Binf)
    ex2=integrate(bcopl2,0,Binf)
    g12l= function(z) 
    { z1=z[1]; z2=z[2]; 
      z1*z2*dnorm(z1)*dnorm(z2)*dcop(pnorm(-z1),pnorm(-z2),cpar)
    }
    ex12=adaptIntegrate(g12l,lowerLimit=c(0,0),upperLimit=c(B,B),absError=.0001)
    ex1=ex1$value/denom; ex2=ex2$value/denom; ex12=ex12$integ/denom
    varr=ex2-ex1*ex1; lcorr=(ex12-ex1*ex1)/varr
    if(iprint) print(c(ex1,ex2,varr,ex12,lcorr))
  }
  else { lcorr=ucorr }
  c(lcorr,ucorr)
}


# dependence measure(s) to copula parameter(s)
# values =  vector, each element is between 0 and 1
# type = "beta", "tau", "rhoS" or "rhoN"
# copname = name of copula family such as 
#   "plackett", "frank", "mtcj", "joe", "gumbel", "galambos", "huesler-reiss"
#   string variables are converted to lower case letters
# Output: (vector of) copula parameter value(s)
depmeas2cpar=function(values,type="beta",copname="gumbel")
{ data(deptabder) 
  type=tolower(type)
  copname=tolower(copname)
  if(copname=="gumbel" | copname=="gum") { deptab=gum.deptab; bfn=gum.b2cpar }
  else if(copname=="galambos" | copname=="gal") { deptab=gal.deptab; bfn=gal.b2cpar }
  else if(copname=="plackett" | copname=="pla") { deptab=pla.deptab; bfn=pla.b2cpar }
  else if(copname=="frank" | copname=="frk") { deptab=frk.deptab; bfn=frk.b2cpar }
  else if(copname=="mtcj" | copname=="mtcjr") { deptab=mtcj.deptab; bfn=mtcj.b2cpar }
  else if(copname=="joe") { deptab=joe.deptab; bfn=joe.b2cpar }
  else if(copname=="hr" | copname=="huesler-reiss" | copname=="husler-reiss") { deptab=hr.deptab; bfn=hr.b2cpar }
  if(type=="beta") { cpar=bfn(values) }
  else
  { 
    if(type=="tau") { bvec=pcinterpolate(deptab$tau,deptab$beta,deptab$tauder,values) }
    else if(type=="rhos") { bvec=pcinterpolate(deptab$rhoS,deptab$beta,deptab$rhoSder,values) }
    else if(type=="rhon") { bvec=pcinterpolate(deptab$rhoN,deptab$beta,deptab$rhoNder,values) }
    bvec=bvec[,1]  # second column has derivatives
    cpar=bfn(bvec)  
  }
  cpar
}

# 2-parameter family BB1 
# (this function can serve as a template for other 2-parameter copula families)
# This code can fail for type="rhoN" and large dependence 
# value = value of dependence measure in (0,1)
# lmU = value of upper tail dependence lambda
# type = "tau", "beta", "rhoS" or "rhoN".
# iprint = print flag for intermediate results
# Output: value of (th,de) given lmU=upper tail dependence, and type=value
bb1.dep2cpar=function(value,lmU,type="tau",iprint=F)
{ type=tolower(type)
  ln2=log(2)
  de=ln2/log(2-lmU)
  if(lmU>=1 | lmU <=0) return(NA)
  if(value>=1 | value <=0) return(NA)
  if(lmU<0.9) ublim=6 else ublim=3
  # min possible value given lmU
  if(type=="tau") { valmin=1-1/de }
  else if(type=="beta") { valmin=be=4*pgum(.5,.5,de)-1 }
  else if(type=="rhos") { valmin=gum.cpar2rhoS(de) }
  else if(type=="rhon") { valmin=rhoN(de,pcond=pcondgum,B=ublim) }
  if(iprint) cat("\n",type, value, "lmU=", lmU,"\n")
  if(value<valmin) return(NA)
  if(type=="tau") { th=2/(de*(1-value))-2; return(c(th,de)) }
  ugrid=seq(.05,.95,.1)
  nn=length(ugrid)
  zevec=rep(0,nn)
  tauvec=rep(0,nn)
  for(i1 in 1:nn)
  { tau=(1-1/de)+ugrid[i1]/de
    tauvec[i1]=tau
    th=2/((1-tau)*de)-2
    cpar=c(th,de)
    if(tau>=0.9 & type=="rhon") break
    if(type=="beta") { ze=4*pbb1(.5,.5,cpar)-1 }
    else if(type=="rhos") { ze=rhoS(cpar,cop=pbb1) }
    else if(type=="rhon") { ze=rhoN(cpar,pcond=pcondbb1,B=ublim) }
    zevec[i1]=ze
    if(iprint) print(c(ugrid[i1],de,th,tau,ze))
  }
  # interpolate step 1
  nn1=sum(zevec>0)
  tauvec=c(1-1/de,tauvec[1:nn1],1)
  zevec=c(valmin,zevec[1:nn1],1)
  nn1=nn1+2
  der=pcderiv(zevec,tauvec)
  tau0=pcinterpolate(zevec,tauvec,der,value)
  tau0=tau0[1]
  if(iprint) cat("tau0=", tau0, " (th0,de)=", 2/(de*(1-tau0))-2,de,"\n")
  if(value>=0.9) return(c(2/(de*(1-tau0))-2,de))
  #c(2/(de*(1-tau0))-2,de)
  # interpolate step 2 : refining 
  taulb=max(1-1/de,tau0-0.02)
  tauub=min(.99,tau0+0.02)
  tauvec=seq(taulb,tauub,length=10)
  nn=length(tauvec)
  zevec=rep(0,nn)
  for(i2 in 1:nn)
  { tau=tauvec[i2]
    th=2/((1-tau)*de)-2
    cpar=c(th,de)
    if(tau>=0.9 & type=="rhon") break
    if(type=="beta") { ze=4*pbb1(.5,.5,cpar)-1 }
    else if(type=="rhos") { ze=rhoS(cpar,cop=pbb1) }
    else if(type=="rhon") { ze=rhoN(cpar,pcond=pcondbb1,B=ublim) }
    zevec[i2]=ze
    #if(iprint) print(c(de,th,tau,ze))
  }
  nn1=sum(zevec>0)
  tauvec=c(1-1/de,tauvec[1:nn1],1)
  zevec=c(valmin,zevec[1:nn1],1)
  nn1=nn1+2
  der=pcderiv(zevec,tauvec)
  tau0=pcinterpolate(zevec,tauvec,der,value)
  tau0=tau0[1]
  if(iprint) cat("tau0=", tau0, " (th0,de)=", 2/(de*(1-tau0))-2,de,"\n")
  c(2/(de*(1-tau0))-2,de)
}

# bivariate t copula family  
# value = value of dependence measure in (-1,1)
# nu = degree of freedom parameter >0
# type = "tau", "beta", "rhoS" or "rhoN".
# ngrid = number of grid points for monotone interpolation
#   if ngrid is larger, could be non-monotone because of integration roundoff
# iprint = print flag for intermediate results
# Output: value of rho given nu and type=value
bvt.dep2cpar=function(value,nu,type="tau",ngrid=15,iprint=F)
{ if(value>=1 | value <=-1) return(NA)
  if(nu<=0) return(NA)
  type=tolower(type)
  if(type=="beta" | type=="tau" ) { ze=bvn.b2cpar(value); return(ze) }
  else if(type=="lambda") { ze=bvt.lm2cpar(value,nu); return(ze) }
  if(abs(value)<0.9) ublim=6 else ublim=4
  if(iprint) cat("\n",type, value, "nu=", nu,"\n")
  # rhoS or rhoN; these are larger as nu decreases
  if(value>=0) { rholb=max(-.99,value-0.02); rhoub=min(.99,value+0.2) }
  else { rholb=max(-.99,value-0.2); rhoub=min(.99,value+0.02) }
  rhovec=seq(rholb,rhoub,length=ngrid)
  nn=ngrid
  zevec=rep(0,nn)
  for(ii in 1:nn)
  { rho=rhovec[ii]
    cpar=c(rho,nu)
    if(type=="rhos") { ze=rhoS(cpar,cop=pcondbvtcop,icond=T) }
    else if(type=="rhon") { ze=rhoN(cpar,pcond=pcondbvtcop,B=ublim) }
    zevec[ii]=ze
    if(iprint) print(c(rhovec[ii],ze))
  }
  der=pcderiv(zevec,rhovec)
  rho0=pcinterpolate(zevec,rhovec,der,value)
  rho0=rho0[1]
  if(iprint) cat("rho=", rho0, "\n")
  rho0
}

