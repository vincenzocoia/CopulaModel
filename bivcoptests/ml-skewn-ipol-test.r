# compare 1-parameter, 2-parameter copula families with skewness to upper tail
# bivariate skew-normal copula using monotone interpolation for qsn

library(CopulaModel)
library(sn)

data(alae600)
alae=alae600$alae 
loss=alae600$loss
n=length(alae) # 600
jj=rank(alae)
ualae.emp=(jj-0.5)/n
jj=rank(loss)
uloss.emp=(jj-0.5)/n

#============================================================

# Pareto cdf with F(x)=1-[1+(x/s)]^{-alp}
# x = positive value
# param = (alp, s), both positive
ppareto=function(x,param)
{ alp=param[1]; s=param[2]
  u=1-(1+(x/s))^{-alp}
  u
}

# density= [1+(x/s)]^{-alp-1} * (alp/s)
# log of Pareto density
logdpareto=function(x,param)
{ alp=param[1]; s=param[2]
  lpdf=log(alp/s)-(alp+1)*log(1+(x/s))
  lpdf
}

# Burr cdf with F(x)= 1- [1+(x/s)^zeta]^{-alp}
# x = positive value
# param = (alp, zeta, s), all positive
pburr=function(x,param)
{ alp=param[1]; zeta=param[2]; s=param[3]
  u=1-(1+(x/s)^zeta)^(-alp)
  u
}

# density= [1+(x/s)^zeta]^{-alp-1} *(x/s)^(zeta-1) * (zeta*alp/s)
# log of Burr density
logdburr=function(x,param)
{ alp=param[1]; zeta=param[2]; s=param[3]
  xs=x/s
  lpdf=log(zeta*alp/s)-(alp+1)*log(1+xs^zeta)+(zeta-1)*log(xs)
  lpdf
}

# negative log-likelihood for Pareto density
# Fbar(x)= [1+(x/s)]^{-alp}
# density= [1+(x/s)]^{-alp-1} * (alp/s)
paretonllk=function(param,xdat)
{ alp=param[1]; s=param[2]
  if(alp<=0 | s<=0) return(1.e10)
  xs=xdat/s
  nllk=(alp+1)*log(1+xs) - log(alp/s)
  sum(nllk)
}

# negative log-likelihood for Burr density
# Fbar(x)= [1+(x/s)^zeta]^{-alp}
# density= [1+(x/s)^zeta]^{-alp-1} *(x/s)^(zeta-1) * (zeta*alp/s)
burrnllk=function(param,xdat)
{ alp=param[1]; zeta=param[2]; s=param[3]
  if(alp<=0 | zeta<=0 | s<=0) return(1.e10)
  xs=xdat/s; tem=xs^zeta
  nllk=(alp+1)*log(1+tem) - (zeta-1)*log(xs) -log(zeta*alp/s)
  sum(nllk)
}

#============================================================

# univariate analysis
# variable 1 is alae/10000, variable 2 is loss/10000
# alae/10000 scaled so that parameter estimates have same order of magnitude
#    (this makes nlm perform better)

cat("\nunivariate for alae, Burr\n")
alae.burr=nlm(burrnllk,p=c(2.39,1,1.6),hessian=T,print.level=1,xdat=alae/10000)
#acov=solve(alae.burr$hess) # inverse Hessian, inverse Fisher information
#SE=sqrt(diag(acov))  # standard errors

cat("\nunivariate for loss, Pareto\n")
loss.pareto=nlm(paretonllk,p=c(1.5,1.5),hessian=T,print.level=1,xdat=loss/10000)
#acov=solve(loss.pareto$hess) # inverse Hessian, inverse Fisher information
#SE=sqrt(diag(acov))  # standard errors
#print(cbind(loss.pareto$estimate,SE))
cat("\n============================================================\n")

#============================================================

# Functions for bivariate skew-normal copula

# log of univariate skew-normal density
# y = real value
# par1 = parameter xi in (-1.1), then ze=xi/sqrt(1-xi^2)
logdsn=function(y,par1)
{ ze=par1/sqrt(1-par1^2)
  dsn(y,0,1,shape=ze,log=TRUE) 
}

# quantile function of univariate skew-normal density
# p = value in (0,1)
# par1 = parameter xi in (-1.1), then ze=xi/sqrt(1-xi^2)
qsnshape=function(p,par1)
{ ze=par1/sqrt(1-par1^2)
  qsn(p,0,1,shape=ze)
}

# log of bivariate skew-normal density
# y1, y2 = real values, could be vectors
# param = c(rho,xi), both in (-1,1), xi appears in the univariate parameter 
#  there is a constraint shown within the function
logdbsn=function(y1,y2,param)
{ rho=param[1]; xi=param[2]
  if(1+2*xi^2*(rho-1)-rho^2<=0) return(-1.e10)
  #if(any(param<=-1) | any(param>=1.)) return(-1.e10)
  ze=xi/sqrt(1-xi^2)
  s=sqrt(1-2*xi^2/(1+rho))
  be=xi/(1+rho)
  bpdf=2*dbvn2(y1,y2,rho)*pnorm(be*(y1+y2)/s)
  log(bpdf)
}

#============================================================

# Burr for alae, Pareto for loss
# 2-stage ifm to estimate dependence parameter 
mle1b=alae.burr$estimate
mle2p=loss.pareto$estimate
ualae.burr=pburr(alae/10000,mle1b)
uloss.pareto=ppareto(loss/10000,mle2p)
udat=cbind(ualae.burr,uloss.pareto)

umin=min(udat)
umax=max(udat)
#ngrid=50 # not reliable enough
ngrid=100
cat("min/max=", umin,umax,"\n")
ugrid=seq(umin,umax,length=ngrid)
param=c(.5,.1)
out=bivcopnllk.ipol(param,udat,logbpdf=logdbsn,logupdf=logdsn,uquant=qsnshape,
  iunivar=2,ppvec=ugrid, LB=-1,UB=1)
print(out)
# -69.87156

copname="bivariate skew-normal copula"
#param=c(0.4602755,0.05)
cat("start1\n")
param=c(0.5,0.2)
cat("\nIFM, Burr for alae, Pareto for loss", copname, "\n")
ifm.bp=nlm(bivcopnllk.ipol,p=param,hessian=T,print.level=1,udat=udat,
  logbpdf=logdbsn,logupdf=logdsn,uquant=qsnshape,
  iunivar=2,ppvec=ugrid, LB=-1,UB=1)

cat("\nstart2\n")
param=c(0.8,0.7)
cat("\nIFM, Burr for alae, Pareto for loss", copname, "\n")
ifm.bp=nlm(bivcopnllk.ipol,p=param,hessian=T,print.level=1,udat=udat,
  logbpdf=logdbsn,logupdf=logdsn,uquant=qsnshape,
  iunivar=2,ppvec=ugrid, LB=-1,UB=1)

# check contour plot with above parameters

#============================================================

# simulation of bivariate skew-normal
# n = simulation sample size
# rho = value in (-1,1)
# xi = value in (-1,1)
#  xi1=xi2=xi is a univariate skew parameter 
#  there is a constraint on rho,xi shown within the function
rbskncop=function(n,rho,xi)
{ if(xi^2>(1+rho)/2) 
  { cat("not positive definite\n"); return(matrix(c(0,1,0,1),2,2)) }
  rmat=matrix(c(1,xi,xi,xi,1,rho,xi,rho,1),3,3)
  amat=chol(rmat)
  zmat=matrix(rnorm(n*3),n,3)
  zmat=zmat%*%amat
  z0=zmat[,1]
  ipos=(z0>0)
  zmat=zmat[ipos,2:3]
  ze=xi/sqrt(1-xi^2) # parameter for univariate skew-normal
  # convert to uniform(0,1)
  umat=zmat
  umat[,1]=psn(zmat[,1],shape=ze)
  umat[,2]=psn(zmat[,2],shape=ze)
  umat
}

# simulation for bivariate skew-normal copula and checking semi-correlations
cat("============================================================\n\n")
cat("semi-correlations of simulated data with MLE of skew-normal copula\n")
rho=ifm.bp$estimate[1]; xi=ifm.bp$estimate[2]
cat("trial 1\n")
set.seed(123)
udat=rbskncop(1350,rho,xi)  
zdat=nscore(udat) 
# plot(zdat)
print(semicor(zdat))
cat("trial 2\n")
set.seed(1234)
udat=rbskncop(1350,rho,xi)  
zdat=nscore(udat) 
print(semicor(zdat))

cat("\nsemi-correlations of alae600\n")
zzz=nscore(cbind(alae600$alae,alae600$loss))
print(semicor(zzz))


