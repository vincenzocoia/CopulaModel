# numerical check of Kendall's tau for (C+Chat)/2, where 
#   Chat is the reflected copula of C, and C is bivariate 

library(CopulaModel)

# u,v = values in (0,1)
# param = copula parameter

# cdf for symmetrized Gumbel
psymgum=function(u,v,param)
{ .5*(pgum(u,v,param)+u+v-1+pgum(1-u,1-v,param)) }

# conditional cdf for symmetrized Gumbel
pcondsymgum=function(v,u,param)
{ .5*(pcondgum(v,u,param)+1-pcondgum(1-v,1-u,param)) }

# cdf for symmetrized MTCJ
pcondsymmtcj=function(v,u,param)
{ .5*(pcondmtcj(v,u,param)+1-pcondmtcj(1-v,1-u,param)) }


cat("Gumbel and symmetrized version\n")
for(cpar in 2:4)
{ cat("cpar=",cpar,"\n")
  tauori=gum.cpar2tau(cpar)
  tau=ktau(cpar,icond=T,pcondsymgum,pcondsymgum)
  cat(tauori,tau,"\n")
}

cat("\nMTCJ and symmetrized version\n")
for(cpar in 2:4)
{ cat("cpar=",cpar,"\n")
  tauori=mtcj.cpar2tau(cpar)
  tau=ktau(cpar,icond=T,pcondsymmtcj,pcondsymmtcj)
  cat(tauori,tau,"\n")
}

# close but not the same
