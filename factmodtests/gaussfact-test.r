# testing of f90gaussfactor.R (p-factor and bi-factor)
#   with gamma factor model
#source("../R/f90gaussfactor.R")

# random gamma convolution bi-factor model
# n = sample size
# th0 = scalar for shape parameter of the shared/common component
# thgrp = vector of shape parameters for group components (length mgrp)
# thvec = vector of shape parameters of individual components, length d
# grsize = vector of grouop sizes (length mgrp)
# Output: nxd random sample; d=length(thvec)=sum(grsize)
rgammaconv.bi=function(n,th0,thgrp,thvec,grsize)
{ d=length(thvec) # also sum(grsize)
  y=matrix(0,n,d)
  mgrp=length(grsize)
  thtmp=NULL
  for(jg in 1:mgrp) thtmp=c(thtmp,rep(thgrp[jg],grsize[jg]))
  for(i in 1:n)
  { yindiv=rgamma(d,thvec,1)
    ycom=rgamma(1,th0,1)
    ygrp=rgamma(mgrp,thgrp,1)
    ytmp=NULL
    for(jg in 1:mgrp) ytmp=c(ytmp,rep(ygrp[jg],grsize[jg]))
    y[i,]=yindiv+ytmp+ycom
  }
  y
}

#============================================================

library(CopulaModel)
set.seed(123)
n=100
th0=3.5
thgrp=c(2.2,2.4)
d=6
thvec=seq(1,2,length=d)
grsz=c(floor(d/2),ceiling(d/2))
ydat=rgammaconv.bi(n,th0,thgrp,thvec,grsz)
Robs=cor(ydat)
cat("correlation matrix\n")
print(Robs)

st=c(rep(.5,d),rep(.4,d))
stp=c(rep(.5,d),rep(.4,grsz[1]),rep(.05,grsz[2]),rep(.05,grsz[1]),rep(.4,grsz[2]))
cat("\nbi-factor nllk\n")
mlb=factanal.bi(grsz,st,cormat=Robs,n=n,prlevel=1)
stb=mlb$rhmat
stp=c(stb[,1],stb[1:grsz[1],2],rep(.0,grsz[2]),rep(.0,grsz[1]),
  stb[(grsz[1]+1):d,2])
cat("============================================================\n")
cat("\n3-factor nllk\n")
ml3=factanal.co(factors=3,stp,cormat=Robs,n=n,prlevel=1)
cat("============================================================\n")
cat("\n2-factor nllk\n")
ml2=factanal.co(factors=2,stp[1:(2*d)],cormat=Robs,n=n,prlevel=1)

# can sometimes fail to convergence for 3-factor

