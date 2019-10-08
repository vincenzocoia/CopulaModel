# Markov order 1 and 2 time series models with negative binomial (NB), 
#  generalized Poission (GP) and Poisson (Po) margins

# functions are
#(gp1mcnllk,gp2mcnllk)
#(nb1mcnllk,nb2mcnllk,pomcnllk)
#(rmspe.mccop1,mltscop1,rmspe.mccop2,mltscop2)
#(nb1mc2nllk,nb2mc2nllk,pomc2nllk)
#(gp1mc2nllk,gp2mc2nllk)
#(nb1ar2nllk,nb2ar2nllk,poar2nllk,gp1ar2nllk,gp2ar2nllk,rmspe.tvn)

# see coptrivmxid.R for some trivariate mixture of max-id copulas for
#  Markov order 2

#============================================================

# NB count Markov ts with copula model for consecutive observations

# form is log P(Y1=y1) + sum log P(Y_t=y_t| Y_{t-1}=y_{t-1})
# P(Y_t=y) is dnbinom(y,th_t,p)   mu_t = exp(b0+b*x_t), p constant for NB1
#   also have the alternative version of NB2 regression.  
# P(Y_{t-1}=y1, Y_t=y2) = C(pnb(y1,tht,p),pnb(y2,tht1,p))
# -C(pnb(y1-1,tht,p),pnb(y2,tht1,p)) -C(pnb(y1,tht,p),pnb(y2-1,tht1,p))
#  +C(pnb(y1-1,tht,p),pnb(y2-1,tht1,p))
# for copula C, use Gumbel, reflectedGumbel , BVN , Frank for 4 tail patterns

# NB2 negative log-likelihood function 
# NB2 parametrization with convolution parameter th fixed over covariates
# param = b0, bvec[], th, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop = bivariate copula cdf, default is BVN copula
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
nb2mcnllk=function(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  #nc=ncol(xdat)  # should be same length as bvec
  b0=param[1]
  np=length(param)
  th=param[np-1]; #pp=1/(1+xi)
  cpar=param[np]
  if(th<=0 | cpar<=cparlb | cpar>=cparub) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  xiv=muvec/th
  ppv=1/(1+xiv)
  cdf1=pnbinom(yy,size=th,prob=ppv)
  pmf=dnbinom(yy,size=th,prob=ppv)
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  nllk=-log(pmf[1])
  for(i in 2:n)
  { tem=pcop(cdf1[i-1],cdf1[i],cpar) - pcop(cdf0[i-1],cdf1[i],cpar) -
      pcop(cdf1[i-1],cdf0[i],cpar) + pcop(cdf0[i-1],cdf0[i],cpar)
    condpr=tem/pmf[i-1]
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
  }
  nllk
}

# NB1 negative log-likelihood function 
# NB1 parametrization with p and xi fixed over covariates
# param = b0, bvec[], xi=1/p-1, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop = bivariate copula cdf, default is BVN copula
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
nb1mcnllk=function(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  #nc=ncol(xdat)  # should be same length as bvec
  b0=param[1]
  np=length(param)
  xi=param[np-1]; pp=1/(1+xi)
  cpar=param[np]
  if(xi<=0 | cpar<=cparlb | cpar>=cparub) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  thv=muvec*pp/(1-pp)
  cdf1=pnbinom(yy,size=thv,prob=pp)
  pmf=dnbinom(yy,size=thv,prob=pp)
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  nllk=-log(pmf[1])
  for(i in 2:n)
  { tem=pcop(cdf1[i-1],cdf1[i],cpar) - pcop(cdf0[i-1],cdf1[i],cpar) -
      pcop(cdf1[i-1],cdf0[i],cpar) + pcop(cdf0[i-1],cdf0[i],cpar)
    condpr=tem/pmf[i-1]
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
  }
  nllk
}

# Poisson count Markov ts with copula model for consecutive observations
# form is log P(Y1=y1) + sum log P(Y_t=y_t| Y_{t-1}=y_{t-1})
# P(Y_t=y) is dpois(y,th_t)   mu_t = exp(b0+b*x_t), 
# P(Y_{t-1}=y1, Y_t=y2) = C(ppois(y1,tht),ppois(y2,tht1))
# -C(ppois(y1-1,tht),ppois(y2,tht1)) -C(ppois(y1,tht),ppois(y2-1,tht1))
#  +C(ppois(y1-1,tht),ppois(y2-1,tht1))
# for copula C, use Gumbel, reflectedGumbel , BVN , Frank for 4 tail patterns

# Poisson negative log-likelihood function 
# param = b0, bvec[], cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop = bivariate copula cdf, default is BVN copula
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
pomcnllk=function(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  #nc=ncol(xdat)  # should be same length as bvec
  b0=param[1]
  np=length(param)
  cpar=param[np]
  if(cpar<=cparlb | cpar>=cparub) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  #thv=muvec
  cdf1=ppois(yy,muvec)
  pmf=dpois(yy,muvec)
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  nllk=-log(pmf[1])
  for(i in 2:n)
  { tem=pcop(cdf1[i-1],cdf1[i],cpar) - pcop(cdf0[i-1],cdf1[i],cpar) -
      pcop(cdf1[i-1],cdf0[i],cpar) + pcop(cdf0[i-1],cdf0[i],cpar)
    condpr=tem/pmf[i-1]
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
  }
  nllk
}

# reparametrization for GP regression
# theta=convolution parameter, vrh is second parameter linked to overdispersion
# mu=theta/(1-vrh), sigma2=theta/(1-vrh)^3,
# sigma2/mu=1/(1-vrh)^2 = 1+xi
# xi=1/(1-vrh)^2-1 >0,  vrh=1-sqrt(1/(1+xi))
# GP1: vrh,xi is fixed and theta is function of covariates through mu
# GP2: theta is fixed and vrh is function of covariates through mu

# GP2 negative log-likelihood function 
# GP2 parametrization with convolution parameter th fixed over covariates
# param = b0,bvec[],th, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop = bivariate copula cdf, default is BVN copula
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
gp2mcnllk=function(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  th=param[np-1]; 
  cpar=param[np]
  if(th<=0 | cpar<=cparlb | cpar>=cparub) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  vrh=1-th/muvec  # should be in 0,1 # is it necessary to check?
  if(any(vrh<=0) | any(vrh>=1)) return(1.e10)
  cdf1=rep(0,n); pmf=rep(0,n);
  #cdf1=pnbinom(yy,size=th,prob=ppv)
  #pmf=dnbinom(yy,size=th,prob=ppv)
  # pgpois currently does not have vectorized interface
  for(i in 1:n) 
  { cdf1[i]=pgpois(yy[i],c(th,vrh[i]));
    pmf[i]=dgpois(yy[i],c(th,vrh[i]));
  }
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  nllk=-log(pmf[1])
  for(i in 2:n)
  { tem=pcop(cdf1[i-1],cdf1[i],cpar) - pcop(cdf0[i-1],cdf1[i],cpar) -
      pcop(cdf1[i-1],cdf0[i],cpar) + pcop(cdf0[i-1],cdf0[i],cpar)
    condpr=tem/pmf[i-1]
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
  }
  nllk
}

# GP1 negative log-likelihood function 
# GP1 parametrization with xi and vrh fixed over covariates
# param = b0,bvec[],xi=1/(1-vrh)^2-1>0, cpar for copula parameter   
#    where vrh=1-sqrt(1/(1+xi)), 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop = bivariate copula cdf, default is BVN copula
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
gp1mcnllk=function(param,yy,xdat=0,pcop=pbvncop,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  xi=param[np-1]; vrh=1-sqrt(1/(1+xi))
  cpar=param[np]
  if(xi<=0 | cpar<=cparlb | cpar>=cparub) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  thv=muvec*(1-vrh)
  cdf1=rep(0,n); pmf=rep(0,n);
  for(i in 1:n) 
  { cdf1[i]=pgpois(yy[i],c(thv[i],vrh));
    pmf[i]=dgpois(yy[i],c(thv[i],vrh));
  }
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  nllk=-log(pmf[1])
  for(i in 2:n)
  { tem=pcop(cdf1[i-1],cdf1[i],cpar) - pcop(cdf0[i-1],cdf1[i],cpar) -
      pcop(cdf1[i-1],cdf0[i],cpar) + pcop(cdf0[i-1],cdf0[i],cpar)
    condpr=tem/pmf[i-1]
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
  }
  nllk
}


# root mean square prediction error
# Markov order 1 model with copula model for 2 consecutive observations
# param = b0,bvec[], constant parameter xi or th, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop = bivariate copula cdf, default is BVN copula
# family = NB1, NB2, GP1, GP2 or Po
# iprint = print flag for intermediate results
# Output: root mean square prediction error plus predictions
rmspe.mccop1=function(param,yy,xdat=0,pcop=pbvncop,family="Po",iprint=F)
{ xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  n=length(yy)
  b0=param[1]
  np=length(param)
  family=tolower(family)
  if(family=="nb1") { xi=param[np-1]; pp=1/(1+xi) }
  else if(family=="nb2") { th=param[np-1] }  
  else if(family=="gp1") { xi=param[np-1]; vrh=1-sqrt(1/(1+xi)) }
  else if(family=="gp2") { th=param[np-1] }  
  cpar=param[np]
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  if(family=="nb1") { th=muvec*pp/(1-pp) }
  else if(family=="nb2") { xi=muvec/th; pp=1/(1+xi) }
  else if(family=="gp1") { th=muvec*(1-vrh); }
  else if(family=="gp2") { vrh=1-th/muvec; xi=1/(1-vrh)^2-1 }
  if(family=="nb1" | family=="nb2")
  { cdf1=pnbinom(yy,size=th,prob=pp)
    pmf=dnbinom(yy,size=th,prob=pp)
  }
  else if(family=="gp1")
  { cdf1=rep(0,n); pmf=rep(0,n);
    for(i in 1:n) 
    { cdf1[i]=pgpois(yy[i],c(th[i],vrh));
      pmf[i]=dgpois(yy[i],c(th[i],vrh));
    }
  }
  else if(family=="gp2")
  { cdf1=rep(0,n); pmf=rep(0,n);
    for(i in 1:n) 
    { cdf1[i]=pgpois(yy[i],c(th,vrh[i]));
      pmf[i]=dgpois(yy[i],c(th,vrh[i]));
    }
  }
  else # Poisson
  { cdf1=ppois(yy,muvec)
    pmf=dpois(yy,muvec)
    xi=0
  }
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  predss=0
  ub=floor(muvec+6*sqrt(muvec*(1+xi)))
  predv=rep(0,n)
  for(i in 2:n)
  { condexp=0
    if(family=="nb1") { cdf=pnbinom(0:ub[i],size=th[i],prob=pp) }
    else if(family=="nb2") { cdf=pnbinom(0:ub[i],size=th,prob=pp[i]) }
    else if(family=="gp1") 
    { tbl=gpoispmfcdf(ub[i],th[i],vrh); cdf=tbl[,3] }
    else if(family=="gp2") 
    { tbl=gpoispmfcdf(ub[i],th,vrh[i]); cdf=tbl[,3] }
    else # family=="po") 
    { cdf=ppois(0:ub[i],muvec[i]) }
    for(j in 1:ub[i])
    { tem=pcop(cdf1[i-1],cdf[j+1],cpar) - pcop(cdf0[i-1],cdf[j+1],cpar) -
        pcop(cdf1[i-1],cdf[j],cpar) + pcop(cdf0[i-1],cdf[j],cpar)
      condexp=condexp+j*tem/pmf[i-1]
    }
    predv[i]=condexp
    predss=predss+(condexp-yy[i])^2
  }
  if(iprint) print(predv)
  list(rmse=sqrt(predss/(n-1)), pred=predv)
}

# wrapper function for Markov order 1 time series with copula
# y = count times series vector
# x = matrix of covariates, could be 0 for NULL (no covariates)
# pcop = function for bivariate copula cdf
# start = starting vector for param for maximum likelihood 
# family = univariate count regression model: "Po", "NB1", "NB2", "GP1", "GP2"
# iprint = print flag for extra summary information including SEs
# cparlb = lower bound for copula parameter
# cparub = upper bound for copula parameter
# prlevel = print.level for nlm() for numerical minimization
#  max number of iterations is the default for nlm (add later??)
# Output: MLE as nlm object plus acov, rmspe
mltscop1=function(y,x,pcop=pbvncop,start,family="Po",iprint=F,
  cparlb=0,cparub=30, prlevel=1)
{ family=tolower(family)
  if(family=="nb1")
  { out=nlm(nb1mcnllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop=pcop,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="nb2")
  { out=nlm(nb2mcnllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop=pcop,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="gp1")
  { out=nlm(gp1mcnllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop=pcop,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="gp2")
  { out=nlm(gp2mcnllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop=pcop,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="po") # poisson
  { out=nlm(pomcnllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop=pcop,cparlb=cparlb,cparub=cparub) 
  }
  else { cat("family not supported\n"); return(NA); }
  acov=solve(out$hessian)
  #print(acov)
  if(iprint)
  { cat("famtype=", family,"\n")
    cat("nllk=", out$min, "\n")
    cat("MLE=", out$estimate,"\n")
    cat("SEs", sqrt(diag(acov)),"\n\n")
  }
  outpred=rmspe.mccop1(out$estimate,y,x,pcop=pcop,family=family,iprint=iprint)
  if(iprint) cat("RMSPE=",outpred$rmse,"\n\n")
  list(nllk=out$min,mle=out$estimate,acov=acov,rmspe=outpred$rmse)
}



#============================================================

# NB count Markov order2 ts with copula model for 3 consecutive observations
# trivariate Gaussian copula is handled in a separate function

# NB2 parametrization with conv param th fixed over covariates
# param = b0,bvec[],th, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop3 = trivarate copula cdf, 
# pcop2 = its (1,2) bivariate margin
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
nb2mc2nllk=function(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,
   cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  th=param[np-2]; 
  de=param[np-1]; de13=param[np]
  cpar=c(de,de13)
  if(th<=0 | any(cpar<=cparlb) | any(cpar>=cparub)) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  xiv=muvec/th
  ppv=1/(1+xiv)
  cdf1=pnbinom(yy,size=th,prob=ppv)
  pmf=dnbinom(yy,size=th,prob=ppv)
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  # first probability 
  #pmf12=pcop2(cdf1[1:2],cpar)-pcop2(c(cdf1[1],cdf0[2]),cpar)-
  #      pcop2(c(cdf0[1],cdf1[2]),cpar)+ pcop2(cdf0[1:2],cpar)
  pmf12=pcop2(cdf1[1],cdf1[2],cpar)-pcop2(cdf1[1],cdf0[2],cpar)-
      pcop2(cdf0[1],cdf1[2],cpar)+ pcop2(cdf0[1],cdf0[2],cpar)
  nllk=-log(pmf12)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pcop3(cdf1[ii],cpar) - pcop3(c(cdf0[i-2],cdf1[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf0[i-1],cdf1[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf0[i]),cpar) - pcop3(cdf0[ii],cpar)
    condpr=tem/pmf12
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    #ii=(i-1):i
    #pmf12=pcop2(cdf1[ii],cpar)-pcop2(c(cdf1[i-1],cdf0[i]),cpar)-
    #  pcop2(c(cdf0[i-1],cdf1[i]),cpar)+ pcop2(cdf0[ii],cpar)
    pmf12=pcop2(cdf1[i-1],cdf1[i],cpar)-pcop2(cdf1[i-1],cdf0[i],cpar)-
      pcop2(cdf0[i-1],cdf1[i],cpar)+ pcop2(cdf0[i-1],cdf0[i],cpar)
  }
  nllk
}

# NB1 parametrization with p fixed over covariates
# param = b0,b[],xi=1/p-1, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop3 = trivarate copula cdf, 
# pcop2 = its (1,2) bivariate margin
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
nb1mc2nllk=function(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  xi=param[np-2]; pp=1/(1+xi)
  de=param[np-1]; de13=param[np]
  cpar=c(de,de13)
  if(xi<=0 | any(cpar<=cparlb) | any(cpar>=cparub)) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  thv=muvec*pp/(1-pp)
  cdf1=pnbinom(yy,size=thv,prob=pp)
  pmf=dnbinom(yy,size=thv,prob=pp)
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  # first probability 
  #pmf12=pcop2(cdf1[1:2],cpar)-pcop2(c(cdf1[1],cdf0[2]),cpar)-
  #      pcop2(c(cdf0[1],cdf1[2]),cpar)+ pcop2(cdf0[1:2],cpar)
  pmf12=pcop2(cdf1[1],cdf1[2],cpar)-pcop2(cdf1[1],cdf0[2],cpar)-
      pcop2(cdf0[1],cdf1[2],cpar)+ pcop2(cdf0[1],cdf0[2],cpar)
  nllk=-log(pmf12)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pcop3(cdf1[ii],cpar) - pcop3(c(cdf0[i-2],cdf1[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf0[i-1],cdf1[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf0[i]),cpar) - pcop3(cdf0[ii],cpar)
    condpr=tem/pmf12
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    #ii=(i-1):i
    #pmf12=pcop2(cdf1[ii],cpar)-pcop2(c(cdf1[i-1],cdf0[i]),cpar)-
    #  pcop2(c(cdf0[i-1],cdf1[i]),cpar)+ pcop2(cdf0[ii],cpar)
    pmf12=pcop2(cdf1[i-1],cdf1[i],cpar)-pcop2(cdf1[i-1],cdf0[i],cpar)-
      pcop2(cdf0[i-1],cdf1[i],cpar)+ pcop2(cdf0[i-1],cdf0[i],cpar)
  }
  nllk
}

# Poisson count Markov order2 ts with copula model for 3 consecutive observations
# param = b0, b[], cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop3 = trivarate copula cdf, 
# pcop2 = its (1,2) bivariate margin
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
pomc2nllk=function(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  de=param[np-1]; de13=param[np]
  cpar=c(de,de13)
  if(any(cpar<=cparlb) | any(cpar>=cparub)) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  cdf1=ppois(yy,muvec)
  pmf=dpois(yy,muvec)
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  # first probability 
  pmf12=pcop2(cdf1[1],cdf1[2],cpar)-pcop2(cdf1[1],cdf0[2],cpar)-
      pcop2(cdf0[1],cdf1[2],cpar)+ pcop2(cdf0[1],cdf0[2],cpar)
  nllk=-log(pmf12)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pcop3(cdf1[ii],cpar) - pcop3(c(cdf0[i-2],cdf1[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf0[i-1],cdf1[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf0[i]),cpar) - pcop3(cdf0[ii],cpar)
    condpr=tem/pmf12
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    pmf12=pcop2(cdf1[i-1],cdf1[i],cpar)-pcop2(cdf1[i-1],cdf0[i],cpar)-
      pcop2(cdf0[i-1],cdf1[i],cpar)+ pcop2(cdf0[i-1],cdf0[i],cpar)
  }
  nllk
}

# GP2 parametrization with conv param th fixed over covariates
# param = b0, b[], th, cpar for copula parameter, 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop3 = trivarate copula cdf, 
# pcop2 = its (1,2) bivariate margin
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
gp2mc2nllk=function(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  th=param[np-2]; 
  de=param[np-1]; de13=param[np]
  cpar=c(de,de13)
  if(th<=0 | any(cpar<=cparlb) | any(cpar>=cparub)) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  vrh=1-th/muvec  
  cdf1=rep(0,n); pmf=rep(0,n);
  for(i in 1:n) 
  { cdf1[i]=pgpois(yy[i],c(th,vrh[i]));
    pmf[i]=dgpois(yy[i],c(th,vrh[i]));
  }
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  # first probability 
  pmf12=pcop2(cdf1[1],cdf1[2],cpar)-pcop2(cdf1[1],cdf0[2],cpar)-
      pcop2(cdf0[1],cdf1[2],cpar)+ pcop2(cdf0[1],cdf0[2],cpar)
  nllk=-log(pmf12)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pcop3(cdf1[ii],cpar) - pcop3(c(cdf0[i-2],cdf1[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf0[i-1],cdf1[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf0[i]),cpar) - pcop3(cdf0[ii],cpar)
    condpr=tem/pmf12
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    pmf12=pcop2(cdf1[i-1],cdf1[i],cpar)-pcop2(cdf1[i-1],cdf0[i],cpar)-
      pcop2(cdf0[i-1],cdf1[i],cpar)+ pcop2(cdf0[i-1],cdf0[i],cpar)
  }
  nllk
}

# GP1 parametrization with xi and vrh fixed over covariates
# param = b0, b[], xi=1/(1-vrh)^2-1, cpar for copula parameter 
#    where vrh=1-sqrt(1/(1+xi))
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop3 = trivarate copula cdf, 
# pcop2 = its (1,2) bivariate margin
# cparlb, cparub = lower / upper bound on cpar
# Output: negative log-likelihood
gp1mc2nllk=function(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,cparlb=0,cparub=30)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  xi=param[np-2]; vrh=1-sqrt(1/(1+xi))
  de=param[np-1]; de13=param[np]
  cpar=c(de,de13)
  if(xi<=0 | any(cpar<=cparlb) | any(cpar>=cparub)) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  thv=muvec*(1-vrh)
  cdf1=rep(0,n); pmf=rep(0,n);
  for(i in 1:n) 
  { cdf1[i]=pgpois(yy[i],c(thv[i],vrh));
    pmf[i]=dgpois(yy[i],c(thv[i],vrh));
  }
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  # first probability 
  pmf12=pcop2(cdf1[1],cdf1[2],cpar)-pcop2(cdf1[1],cdf0[2],cpar)-
      pcop2(cdf0[1],cdf1[2],cpar)+ pcop2(cdf0[1],cdf0[2],cpar)
  nllk=-log(pmf12)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pcop3(cdf1[ii],cpar) - pcop3(c(cdf0[i-2],cdf1[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf1[i]),cpar) -
      pcop3(c(cdf1[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf0[i-1],cdf1[i]),cpar) +
      pcop3(c(cdf0[i-2],cdf1[i-1],cdf0[i]),cpar) +
      pcop3(c(cdf1[i-2],cdf0[i-1],cdf0[i]),cpar) - pcop3(cdf0[ii],cpar)
    condpr=tem/pmf12
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    pmf12=pcop2(cdf1[i-1],cdf1[i],cpar)-pcop2(cdf1[i-1],cdf0[i],cpar)-
      pcop2(cdf0[i-1],cdf1[i],cpar)+ pcop2(cdf0[i-1],cdf0[i],cpar)
  }
  nllk
}


# root mean square prediction error
# Markov order 2 model with copula model for 3 consecutive observations
# param = b0, bvec[], constant parameter xi or th, cpar for copula parameter 
# yy = vector of counts, 
# xdat = matrix of covariates
# pcop3 = trivarate copula cdf, 
# pcop2 = its (1,2) bivariate margin
# family = NB1, NB2, GP1, GP2 or Po
# iprint = print flag for intermediate results
# Output: root mean square prediction error plus predictions
rmspe.mccop2=function(param,yy,xdat=0,pcop3=pmxid3ls,pcop2=pmxid2ls,family="Po",iprint=F)
{ xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  n=length(yy)
  b0=param[1]
  np=length(param)
  family=tolower(family)
  if(family=="nb1") { xi=param[np-2]; pp=1/(1+xi) }
  else if(family=="nb2") { th=param[np-2] }  
  else if(family=="gp1") { xi=param[np-2]; vrh=1-sqrt(1/(1+xi)) }
  else if(family=="gp2") { th=param[np-2] }
  # to adapt for Poisson similar to rmspe.mccop
  de=param[np-1]; de13=param[np]
  cpar=c(de,de13)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  if(family=="nb1") { th=muvec*pp/(1-pp) }
  else if(family=="nb2") { xi=muvec/th; pp=1/(1+xi) }
  else if(family=="gp1") { th=muvec*(1-vrh); }
  else if(family=="gp2") { vrh=1-th/muvec; xi=1/(1-vrh)^2-1 }
  if(family=="nb1" | family=="nb2")
  { cdf1=pnbinom(yy,size=th,prob=pp)
    pmf=dnbinom(yy,size=th,prob=pp)
  }
  else if(family=="gp1")
  { cdf1=rep(0,n); pmf=rep(0,n);
    for(i in 1:n) 
    { cdf1[i]=pgpois(yy[i],c(th[i],vrh));
      pmf[i]=dgpois(yy[i],c(th[i],vrh));
    }
  }
  else if(family=="gp2")
  { cdf1=rep(0,n); pmf=rep(0,n);
    for(i in 1:n) 
    { cdf1[i]=pgpois(yy[i],c(th,vrh[i]));
      pmf[i]=dgpois(yy[i],c(th,vrh[i]));
    }
  }
  else # Poisson
  { cdf1=ppois(yy,muvec)
    pmf=dpois(yy,muvec)
    xi=0
  }
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=0
  predss=0
  ub=floor(muvec+6*sqrt(muvec*(1+xi)))
  predv=rep(0,n)
  for(i in 3:n)
  { condexp=0
    # probability for the Y[i]=j given Y[i-2], Y[i-1]
    if(family=="nb1") { cdf=pnbinom(0:ub[i],size=th[i],prob=pp) }
    else if(family=="nb2") { cdf=pnbinom(0:ub[i],size=th,prob=pp[i]) }
    else if(family=="gp1") 
    { tbl=gpoispmfcdf(ub[i],th[i],vrh); cdf=tbl[,3] }
    else if(family=="gp2") 
    { tbl=gpoispmfcdf(ub[i],th,vrh[i]); cdf=tbl[,3] }
    else # family=="po" 
    { cdf=ppois(0:ub[i],muvec[i]) }
    pmf12=pcop2(cdf1[i-2],cdf1[i-1],cpar)-pcop2(cdf1[i-2],cdf0[i-1],cpar)-
      pcop2(cdf0[i-2],cdf1[i-1],cpar)+ pcop2(cdf0[i-2],cdf0[i-1],cpar)
    for(j in 1:ub[i])
    { tem=pcop3(c(cdf1[i-2],cdf1[i-1],cdf[j+1]),cpar) - 
        pcop3(c(cdf0[i-2],cdf1[i-1],cdf[j+1]),cpar) -
        pcop3(c(cdf1[i-2],cdf0[i-1],cdf[j+1]),cpar) -
        pcop3(c(cdf1[i-2],cdf1[i-1],cdf[j]),cpar) +
        pcop3(c(cdf0[i-2],cdf0[i-1],cdf[j+1]),cpar) +
        pcop3(c(cdf0[i-2],cdf1[i-1],cdf[j]),cpar) +
        pcop3(c(cdf1[i-2],cdf0[i-1],cdf[j]),cpar) - 
        pcop3(c(cdf0[i-2],cdf0[i-1],cdf[j]),cpar)
      condpr=tem/pmf12
      condexp=condexp+j*condpr
    }
    predv[i]=condexp
    predss=predss+(condexp-yy[i])^2
  }
  if(iprint) print(predv)
  list(rmse=sqrt(predss/(n-2)), pred=predv)
}

# wrapper function for Markov order 2 time series based on trivariate copula
# y = count times series vector
# x = matrix of covariates, could be 0 for NULL (no covariates)
# pcop3 = function for trivariate copula cdf
# pcop2 = function for marginal bivariate copula cdf (for first two)
# start = starting vector for param for maximum likelihood 
# family = univariate count regression model: "Po", "NB1", "NB2", "GP1", "GP2"
# iprint = print flag for extra summary information including SEs
# cparlb = lower bound for copula parameter
# cparub = upper bound for copula parameter
# prlevel = print.level for nlm() for numerical minimization
#  max number of iterations is the default for nlm (add later??)
# Output: MLE as nlm object plus acov, rmspe
mltscop2=function(y,x,pcop3=pmxid3ls,pcop2=pmxid2ls,start,family="Po",iprint=F,cparlb=0,cparub=30,prlevel=1)
{ family=tolower(family) 
  if(family=="nb1")
  { out=nlm(nb1mc2nllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop3=pcop3,pcop2=pcop2,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="nb2")
  { out=nlm(nb2mc2nllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop3=pcop3,pcop2=pcop2,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="gp1")
  { out=nlm(gp1mc2nllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop3=pcop3,pcop2=pcop2,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="gp2")
  { out=nlm(gp2mc2nllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop3=pcop3,pcop2=pcop2,cparlb=cparlb,cparub=cparub) 
  }
  else if(family=="po")
  { out=nlm(pomc2nllk,p=start,hessian=T,print.level=prlevel,yy=y,xdat=x,
      pcop3=pcop3,pcop2=pcop2,cparlb=cparlb,cparub=cparub) 
  }
  else { cat("family not supported\n"); return(NA); }
  acov=solve(out$hessian)
  #print(acov)
  if(iprint)
  { cat("famtype=", family,"\n")
    cat("nllk=", out$min, "\n")
    cat("MLE=", out$estimate,"\n")
    cat("SEs", sqrt(diag(acov)),"\n\n")
  }
  outpred=rmspe.mccop2(out$estimate,y,x,pcop3,pcop2,family=family,iprint=iprint)
  if(iprint) cat("RMSPE=",outpred$rmse,"\n\n")
  list(nllk=out$min,mle=out$estimate,acov=acov,rmspe=outpred$rmse)
}


#============================================================

# Markov order 2 based on trivariate Gaussian

# AR(2) with NB1/NB2 margins
# pmnorm tkane from library(mprobit) for this package 

# trivariate Gaussian for Markov order 2 transition
# NB2 parametrization with convolution param th fixed over covariates
# This uses pmnorm for trivariate Gaussian rectangle probabilities 
# yy = vector of counts, 
# xdat = matrix of covariates
# param = b0, bvec[], th, cpar=(rh1,rh2) 
# Output: negative log-likelihood
nb2ar2nllk=function(param,yy,xdat)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  nc=ncol(xdat)  # should be same length as bvec
  b0=param[1]
  np=length(param)
  th=param[np-2]; #pp=1/(1+xi)
  rh1=param[np-1]; rh2=param[np]  # acf lag1 and lag2
  cpar=c(rh1,rh2)
  if(th<=0 | any(cpar<=-1) | any(cpar>=1)) return(1.e10)
  dt=1+2*rh1*rh1*rh2-2*rh1*rh1-rh2*rh2
  if(dt<=0) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  xiv=muvec/th
  ppv=1/(1+xiv)
  cdf1=pnbinom(yy,size=th,prob=ppv)
  pmf=dnbinom(yy,size=th,prob=ppv)
  cdf1[cdf1>=1]=1-1.e-9
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=1.e-9
  z1=qnorm(cdf1); z0=qnorm(cdf0)  # convert to z values 
  # first probability 
  rmat12=matrix(c(1,rh1,rh1,1),2,2)
  pmf12=pmnorm(lb=z0[1:2],ub=z1[1:2],rep(0,2),rmat12,eps=1.e-05)
  nllk=-log(pmf12$pr)
  n=length(yy)
  rmat123=matrix(c(1,rh1,rh2,rh1,1,rh1,rh2,rh1,1),3,3)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pmnorm(lb=z0[ii], ub=z1[ii], rep(0,3),rmat123)
    tem=tem$pr
    condpr=tem/pmf12$pr
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    ii=(i-1):i
    pmf12=pmnorm(lb=z0[ii],ub=z1[ii],rep(0,2),rmat12,eps=1.e-05)
  }
  nllk
}


# trivariate Gaussian for Markov order 2 transition
# NB1 parametrization with p or xi fixed over covariates
# This uses pmnorm for trivariate Gaussian rectangle probabilities 
# yy = vector of counts, 
# xdat = matrix of covariates
# param = b0, bvec[], xi, cpar=(rh1,rh2) 
# Output: negative log-likelihood
nb1ar2nllk=function(param,yy,xdat)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  xi=param[np-2]; pp=1/(1+xi)
  rh1=param[np-1]; rh2=param[np]  # acf lag1 and lag2
  cpar=c(rh1,rh2)
  if(xi<=0 | any(cpar<=-1) | any(cpar>=1)) return(1.e10)
  dt=1+2*rh1*rh1*rh2-2*rh1*rh1-rh2*rh2
  if(dt<=0) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  thv=muvec*pp/(1-pp)
  cdf1=pnbinom(yy,size=thv,prob=pp)
  pmf=dnbinom(yy,size=thv,prob=pp)
  cdf1[cdf1>=1]=1-1.e-9
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=1.e-9
  z1=qnorm(cdf1); z0=qnorm(cdf0)  # convert to z values 
  # first probability 
  rmat12=matrix(c(1,rh1,rh1,1),2,2)
  pmf12=pmnorm(lb=z0[1:2],ub=z1[1:2],rep(0,2),rmat12,eps=1.e-05)
  nllk=-log(pmf12$pr)
  n=length(yy)
  rmat123=matrix(c(1,rh1,rh2,rh1,1,rh1,rh2,rh1,1),3,3)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pmnorm(lb=z0[ii], ub=z1[ii], rep(0,3),rmat123)
    tem=tem$pr
    condpr=tem/pmf12$pr
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    ii=(i-1):i
    pmf12=pmnorm(lb=z0[ii],ub=z1[ii],rep(0,2),rmat12,eps=1.e-05)
  }
  nllk
}

# AR(2) with Poisson margin
# trivariate Gaussian for Markov order 2 transition
# This uses pmnorm for trivariate Gaussian rectangle probabilities 
# yy = vector of counts, 
# xdat = matrix of covariates
# param= b0, bvec[], cpar=(rh1,rh2) 
# Output: negative log-likelihood
poar2nllk=function(param,yy,xdat)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  rh1=param[np-1]; rh2=param[np]  # acf lag1 and lag2
  cpar=c(rh1,rh2)
  if(any(cpar<=-1) | any(cpar>=1)) return(1.e10)
  dt=1+2*rh1*rh1*rh2-2*rh1*rh1-rh2*rh2
  if(dt<=0) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  #thv=muvec*pp/(1-pp)
  cdf1=ppois(yy,muvec)
  pmf=dpois(yy,muvec)
  cdf1[cdf1>=1]=1-1.e-9
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=1.e-9
  z1=qnorm(cdf1); z0=qnorm(cdf0)  # convert to z values 
  # first probability 
  rmat12=matrix(c(1,rh1,rh1,1),2,2)
  pmf12=pmnorm(lb=z0[1:2],ub=z1[1:2],rep(0,2),rmat12,eps=1.e-05)
  nllk=-log(pmf12$pr)
  n=length(yy)
  rmat123=matrix(c(1,rh1,rh2,rh1,1,rh1,rh2,rh1,1),3,3)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pmnorm(lb=z0[ii], ub=z1[ii], rep(0,3),rmat123)
    tem=tem$pr
    condpr=tem/pmf12$pr
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    ii=(i-1):i
    pmf12=pmnorm(lb=z0[ii],ub=z1[ii],rep(0,2),rmat12,eps=1.e-05)
  }
  nllk
}

# trivariate Gaussian for Markov order 2 transition
# GP2 parametrization with convolution param th fixed over covariates
# This uses pmnorm for trivariate Gaussian rectangle probabilities 
# yy = vector of counts, 
# xdat = matrix of covariates
# param = b0, bvec[], th, cpar=(rh1,rh2),
# Output: negative log-likelihood
gp2ar2nllk=function(param,yy,xdat)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  th=param[np-2]; 
  rh1=param[np-1]; rh2=param[np]  # acf lag1 and lag2
  cpar=c(rh1,rh2)
  if(th<=0 | any(cpar<=-1) | any(cpar>=1)) return(1.e10)
  dt=1+2*rh1*rh1*rh2-2*rh1*rh1-rh2*rh2
  if(dt<=0) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  vrh=1-th/muvec  # 
  if(any(vrh<=0) | any(vrh>=1)) return(1.e10)
  cdf1=rep(0,n); pmf=rep(0,n);
  for(i in 1:n) 
  { cdf1[i]=pgpois(yy[i],c(th,vrh[i]));
    pmf[i]=dgpois(yy[i],c(th,vrh[i]));
  }
  cdf1[cdf1>=1]=1-1.e-9
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=1.e-9
  z1=qnorm(cdf1); z0=qnorm(cdf0)  # convert to z values 
  # first probability 
  rmat12=matrix(c(1,rh1,rh1,1),2,2)
  pmf12=pmnorm(lb=z0[1:2],ub=z1[1:2],rep(0,2),rmat12,eps=1.e-05)
  nllk=-log(pmf12$pr)
  n=length(yy)
  rmat123=matrix(c(1,rh1,rh2,rh1,1,rh1,rh2,rh1,1),3,3)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pmnorm(lb=z0[ii], ub=z1[ii], rep(0,3),rmat123)
    tem=tem$pr
    condpr=tem/pmf12$pr
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    ii=(i-1):i
    pmf12=pmnorm(lb=z0[ii],ub=z1[ii],rep(0,2),rmat12,eps=1.e-05)
  }
  nllk
}


# trivariate Gaussian for Markov order 2 transition
# GP1 parametrization with vrh or xi fixed over covariates
# This uses pmnorm for trivariate Gaussian rectangle probabilities 
# yy = vector of counts, 
# xdat = matrix of covariates
# param = b0, bvec[], xi, cpar=(rh1,rh2)
# Output: negative log-likelihood
gp1ar2nllk=function(param,yy,xdat)
{ n=length(yy)
  xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  b0=param[1]
  np=length(param)
  xi=param[np-2]; vrh=1-sqrt(1/(1+xi))
  rh1=param[np-1]; rh2=param[np]  # acf lag1 and lag2
  cpar=c(rh1,rh2)
  if(xi<=0 | any(cpar<=-1) | any(cpar>=1)) return(1.e10)
  dt=1+2*rh1*rh1*rh2-2*rh1*rh1-rh2*rh2
  if(dt<=0) return(1.e10)
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  thv=muvec*(1-vrh)
  cdf1=rep(0,n); pmf=rep(0,n);
  for(i in 1:n) 
  { cdf1[i]=pgpois(yy[i],c(thv[i],vrh));
    pmf[i]=dgpois(yy[i],c(thv[i],vrh));
  }
  cdf1[cdf1>=1]=1-1.e-9
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=1.e-9
  z1=qnorm(cdf1); z0=qnorm(cdf0)  # convert to z values 
  # first probability 
  rmat12=matrix(c(1,rh1,rh1,1),2,2)
  pmf12=pmnorm(lb=z0[1:2],ub=z1[1:2],rep(0,2),rmat12,eps=1.e-05)
  nllk=-log(pmf12$pr)
  n=length(yy)
  rmat123=matrix(c(1,rh1,rh2,rh1,1,rh1,rh2,rh1,1),3,3)
  for(i in 3:n)
  { ii=(i-2):i
    tem=pmnorm(lb=z0[ii], ub=z1[ii], rep(0,3),rmat123)
    tem=tem$pr
    condpr=tem/pmf12$pr
    if(condpr<=0. | is.na(condpr)) condpr=1.e-15
    nllk=nllk-log(condpr)
    ii=(i-1):i
    pmf12=pmnorm(lb=z0[ii],ub=z1[ii],rep(0,2),rmat12,eps=1.e-05)
  }
  nllk
}


# root mean square prediction error for trivariate Gaussian/normal transition
# param = b0,bvec[], constant parameter xi or th, cpar=(rh1,rh2) 
# yy = vector of counts, 
# xdat = matrix of covariates
# family = NB1, NB2, GP1, GP2 or Po
# iprint = print flag for intermediate results
# Output: root mean square prediction error plus predictions
rmspe.tvn=function(param,yy,xdat=0,family="Po",iprint=F)
{ xdat=as.matrix(xdat)
  if(nrow(xdat)==1) nc=0 else nc=ncol(xdat)
  n=length(yy)
  b0=param[1]
  np=length(param)
  family=tolower(family)
  if(family=="nb1") { xi=param[np-2]; pp=1/(1+xi) }
  else if(family=="nb2") { th=param[np-2] }
  else if(family=="gp1") { xi=param[np-2]; vrh=1-sqrt(1/(1+xi)) }
  else if(family=="gp2") { th=param[np-2] }
  rh1=param[np-1]; rh2=param[np]  # acf lag1 and lag2
  if(nc>0) 
  { bvec=param[2:(nc+1)];  
    muvec=b0+xdat%*%bvec; muvec=exp(muvec) 
  }
  else muvec=rep(exp(b0),n)
  if(family=="nb1") { th=muvec*pp/(1-pp) }
  else if(family=="nb2") { xi=muvec/th; pp=1/(1+xi) }
  else if(family=="gp1") { th=muvec*(1-vrh) }
  else if(family=="gp2") { vrh=1-th/muvec; xi=1/(1-vrh)^2-1 }
  if(family=="nb1" | family=="nb2")
  { cdf1=pnbinom(yy,size=th,prob=pp)
    pmf=dnbinom(yy,size=th,prob=pp)
  }
  else if(family=="gp1")
  { cdf1=rep(0,n); pmf=rep(0,n);
    for(i in 1:n) 
    { cdf1[i]=pgpois(yy[i],c(th[i],vrh));
      pmf[i]=dgpois(yy[i],c(th[i],vrh));
    }
  }
  else if(family=="gp2")
  { cdf1=rep(0,n); pmf=rep(0,n);
    for(i in 1:n) 
    { cdf1[i]=pgpois(yy[i],c(th,vrh[i]));
      pmf[i]=dgpois(yy[i],c(th,vrh[i]));
    }
  }
  else # Poisson
  { cdf1=ppois(yy,muvec)
    pmf=dpois(yy,muvec)
    xi=0
  }
  cdf1[cdf1>=1]=1-1.e-9
  cdf0=cdf1-pmf
  cdf0[cdf0<=0]=1.e-9
  z1=qnorm(cdf1); z0=qnorm(cdf0)  # convert to z values 
  rmat12=matrix(c(1,rh1,rh1,1),2,2)
  rmat123=matrix(c(1,rh1,rh2,rh1,1,rh1,rh2,rh1,1),3,3)
  predss=0
  ub=floor(muvec+6*sqrt(muvec*(1+xi)))
  predv=rep(0,n)
  for(i in 3:n)
  { condexp=0
    # probability for the Y[i]=j given Y[i-2], Y[i-1]
    if(family=="nb1") { cdf=pnbinom(0:ub[i],size=th[i],prob=pp) }
    else if(family=="nb2") { cdf=pnbinom(0:ub[i],size=th,prob=pp[i]) }
    else if(family=="gp1") 
    { tbl=gpoispmfcdf(ub[i],th[i],vrh); cdf=tbl[,3] }
    else if(family=="gp2") 
    { tbl=gpoispmfcdf(ub[i],th,vrh[i]); cdf=tbl[,3] }
    else  # family=="po"
    { cdf=ppois(0:ub[i],muvec[i]) }
    z=qnorm(cdf)
    ii=(i-2):(i-1)
    pmf12=pmnorm(lb=z0[ii],ub=z1[ii],rep(0,2),rmat12,eps=1.e-05)
    for(j in 1:ub[i])
    { tem=pmnorm(lb=c(z0[ii],z[j]), ub=c(z1[ii],z[j+1]), rep(0,3),rmat123)
      tem=tem$pr
      condpr=tem/pmf12$pr
      condexp=condexp+j*condpr
    }
    #print(c(i,condexp))
    predv[i]=condexp
    predss=predss+(condexp-yy[i])^2
  }
  
  if(iprint) print(predv)
  list(rmse=sqrt(predss/(n-2)), pred=predv)
}

#============================================================

