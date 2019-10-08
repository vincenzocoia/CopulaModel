# compare 1-parameter, 2-parameter copula families with skewness to upper tail
# Archimedean copula based on inverse gamma LT

library(CopulaModel)

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
acov=solve(alae.burr$hess) # inverse Hessian, inverse Fisher information
SE=sqrt(diag(acov))  # standard errors
#print(cbind(alae.burr$estimate,SE))
# alphat  2.124780 0.48157019
# zetahat 1.095413 0.06377033
# sighat  1.423125 0.43074198

cat("\nunivariate for loss, Pareto\n")
loss.pareto=nlm(paretonllk,p=c(1.5,1.5),hessian=T,print.level=1,xdat=loss/10000)
acov=solve(loss.pareto$hess) # inverse Hessian, inverse Fisher information
SE=sqrt(diag(acov))  # standard errors
print(cbind(loss.pareto$estimate,SE))
#    alphat  1.320371 0.1289922
#    sighat  1.666202 0.2598457
cat("\n============================================================\n")

#============================================================

# wrapper for estimation with IFM and full log-likelihood
#     Burr for alae, Pareto for loss
# lgdcop = function for log copula density
# copname = string for copula name
# cpar = starting point for copula parameter
# cparlb = lower bound for copula parameter
# cparub = upper bound for copula parameter
# prlev = print.level for nlm()
# Output: nlm objects from IFM, pseudo likelihood and full likelihood 
ifmandfull=function(lgdcop,copname,cpar,cparlb=0,cparub=100,prlev=1)
{ cat("\nIFM, Burr for alae, Pareto for loss", copname, "\n")
  ifm.bp=nlm(bivcopnllk,p=cpar,hessian=T,print.level=prlev,udat=cbind(ualae.burr,uloss.pareto),logdcop=lgdcop,LB=cparlb,UB=cparub)
  cat("cpar=", ifm.bp$estimate,"\n") 
  if(length(cpar)>1) ifmcov=solve(ifm.bp$hess) else ifmcov=1/ifm.bp$hess
  ifmSE= sqrt(diag(ifmcov))
  cat("ifmSE=", ifmSE,"\n")

  cat("\npseudo, Burr for alae, Pareto for loss\n")
  emp.bp=nlm(bivcopnllk,p=cpar,hessian=T,print.level=prlev,udat=cbind(ualae.emp,uloss.emp),logdcop=lgdcop,LB=cparlb,UB=cparub)
  cat("cpar=", emp.bp$estimate,"\n") 
  if(length(cpar)>1) empcov=solve(emp.bp$hess) else empcov=1/emp.bp$hess
  empSE= sqrt(diag(empcov))
  cat("empSE=", empSE,"\n")

  cat("\nfull log-likelihood, Burr for alae, Pareto for loss\n")
  parambp=c(mle1b,mle2p,ifm.bp$estimate)
  full.bp=nlm(bivmodnllk,p=parambp,hessian=T,print.level=prlev,
    xdat=cbind(alae/10000,loss/10000),logdcop=lgdcop,
    logpdf1=logdburr,cdf1=pburr,np1=3,logpdf2=logdpareto,cdf2=ppareto,np2=2,
    LB=c(rep(0,5),cparlb),UB=c(rep(100,5),cparub) )
  acov=solve(full.bp$hess) 
  SE=sqrt(diag(acov))  
  print(cbind(full.bp$estimate,SE))

  cat("\nFor Comparisons ", copname,"\n")
  cat("bivmodnllk(Burr,Pareto,full)=", full.bp$min,"\n")
  cat("bivcopnllk(Burr,Pareto,IFM) =", ifm.bp$min,"\n")
  cat("empnllk(Burr,Pareto)=        ", emp.bp$min,"\n")
  cat("\n============================================================\n\n")
  list(ifm=ifm.bp, emp=emp.bp, full=full.bp)
}


# 2-stage ifm to estimate dependence parameter 
mle1b=alae.burr$estimate
mle2p=loss.pareto$estimate
ualae.burr=pburr(alae/10000,mle1b)
uloss.pareto=ppareto(loss/10000,mle2p)

fit.invgamA=ifmandfull(logdinvgamA,"invgamA",.8,cparlb=0,cparub=10,prlev=1)
fit.gum=ifmandfull(logdgum,"gumbel",1.8,cparlb=1,cparub=20,prlev=1)

