# Dependence tables (beta, tau, rhoS, rhoN, lm) for bivariate copula familes

# 1-parameter copula families: independence to comonotonicity
# Table of dependence measures beta, tau, rhoS, rhoN, lm, cpar,
#  start with beta on a grid, invert to get cpar, then evaluate others
# bvec = vector of beta values (increasing order, first is 0, last is 1)
# bfn = function to get copula parameter cpar given Blomqvist beta
# dcop = function for copula density 
# pcop = function for copula cdf 
# pcond12 = function for conditional cdf C_{1|2}
# pcond21 = function for conditional cdf C_{2|1}, 
#     above conditional cdfs are the same if C is permutation symmetric
# LBcpar = parameter for independence
# UBcpar = parameter for comonotonicity
# itaildep = T to compute lm=lambda
# lmfn = function to get lambda given cpar if itaildep=T
# zero = 0 or something like 1.e-6 (used with integration)
# zbd = integration bound with respect to N(0,1) margins for rhoN
# Output: table with "cpar","beta","tau","rhoS","rhoN","lambda"
makedeptable=function(bvec,bfn,dcop,pcop,pcond12,pcond21,LBcpar=0,UBcpar=Inf,
   itaildep=F,lmfn,zero=0,zbd=6,iprint=F)
{ np=length(bvec)
  cpar=rep(0,np); cpar[np]=UBcpar; cpar[1]=LBcpar
  tau=rep(0,np); tau[np]=1
  rhs=rep(0,np); rhs[np]=1
  rhn=rep(0,np); rhn[np]=1
  #rhn2=rep(0,np); rhn2[np]=1
  lm=rep(0,np); 
  np1=np-1
  for(i in 2:np1)
  { be=bvec[i]; de=bfn(be)
    cpar[i]=de
    rhs[i]=rhoS(de,cop=pcop) # integrate on [0,1]^2
    rhn[i]=rhoN(de,icond=T,pcond=pcond21,B=zbd,tol=1.e-5) 
    #rhn2[i]=rhoN(de,icond=F,icdf=T,pcop=pcop,B=zbd,tol=1.e-5) 
    tau[i]=ktau(de,icond=T,pcond12=pcond12,pcond21=pcond21,zero=zero,tol=1.e-5)
    if(iprint) cat(i,de,be,tau[i],rhs[i],rhn[i],"\n")  
    #if(iprint) cat(i,de,be,tau[i],rhs[i],rhn[i],rhn2[i],"\n")  
  }
  if(itaildep) lm=lmfn(cpar)
  lm[np]=1
  out=cbind(cpar,bvec,tau,rhs,rhn,lm)
  out=data.frame(out)
  names(out)=c("cpar","beta","tau","rhoS","rhoN","lambda")
  out
}

