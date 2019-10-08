# utilities for bivariate discrete

# bpmf = matrix of probability masses
# Output: matrix with bivariate cdf
bivpmf2cdf=function(bpmf)
{ d1=nrow(bpmf)
  d2=ncol(bpmf)
  cdf=apply(bpmf,2,cumsum)
  for(j in 2:d2) cdf[,j]=cdf[,j-1]+cdf[,j]
  cdf
}

# bpmf = matrix of probability masses
# Output: means, variances, covariance, correlation of bivariate count distribution
corbivpmf=function(bpmf)
{ d1=nrow(bpmf)
  d2=ncol(bpmf)
  pmf1=apply(bpmf,1,sum)  # sum over cols
  pmf2=apply(bpmf,2,sum)  # sum over rows
  i1=0:(d1-1)
  i2=0:(d2-1)
  mu1=sum(i1*pmf1)
  ss1=sum(i1*i1*pmf1)
  v1=ss1-mu1^2
  mu2=sum(i2*pmf2)
  ss2=sum(i2*i2*pmf2)
  v2=ss2-mu2^2
  ii=outer(i1,i2)
  ss12=sum(ii*bpmf)
  cv=ss12-mu1*mu2
  r=cv/sqrt(v1*v2)
  out=c(mu1,v1,mu2,v2,cv,r)
  names(out)=c("mu1","var1","mu2","var2","cov","cor")
  out
}

