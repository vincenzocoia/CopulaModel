# Multivariate Huesler-Reiss (based on NJL 2009, Extremes)

# Conditional correlations in multivariate Huesler-Reiss model
# thmat = d*(d-1)/2 parameter vector as a symmetric matrix
#  (the diagonal can be anything on input as it will be changed)
# j = index in 1:d
# Output: matrix for the jth term of exponent A, conditioned on variable j
hrcondcor=function(thmat,j)
{ d=nrow(thmat)
  if(j<1) j=1
  if(j>d) j=d
  j=floor(j)
  d1=d-1
  thinv=1/thmat
  diag(thinv)=0
  pc=matrix(1,d,d)
  for(i2 in 2:d)
  { for(i1 in 1:(i2-1))
    { pc[i1,i2]=(thinv[i1,j]^2+thinv[i2,j]^2-thinv[i1,i2]^2)/2/thinv[i1,j]/thinv[i2,j]
      pc[i2,i1]=pc[i1,i2]
    }
  }
  pc=pc[-j,]; pc=pc[,-j]
  pc
}


# Exponent function of multivariate Huesler-Reiss for d>=3
# This function depends on library mvtnorm
# ww = d-vector of positive values 
# param = d*(d-1)/2 -vector with order 12, 13, 23, 14, ... for parameters
# Output: exponent A function 
mhrA=function(ww,param)
{ d=length(ww)
  thmat=corvec2mat(param)
  asum=0
  lb=rep(-Inf,d-1)
  for(j in 1:d)
  { pcmat=hrcondcor(thmat,j)
    wwj=ww[j]/ww; wwj=wwj[-j]
    thj=thmat[j,]; thj=thj[-j]
    ub= 1/thj+.5*thj*log(wwj)
    tem=pmvnorm(lb,ub,mean=rep(0,d-1),corr=pcmat,algorithm=GenzBretz())
    #print(tem) # accuracy of order 5.e-5
    asum=asum+ww[j]*c(tem)
  }
  asum
}

# Multivariate Huesler-Reiss copula cdf
# uu = d-vector of values in (0,1)
# cpar = d*(d-1)/2 -vector with order 12, 13, 23, 14, ... for parameters
# Output: multivariate cdf
pmhr=function(uu,cpar)
{ ww=-log(uu)
  aa=mhrA(ww,cpar)
  exp(-aa)
}

