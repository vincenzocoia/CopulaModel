
# Inputs:
#  rmat = sample correlation matrix
#  mxfactor = max number of factors to use
#  n = sample size
#  iprint = print flag for intermediate results
# Output : -log-likelihoods for 1 to mxfactor factors
#  note that upper bound of mxfactor depends on dim(rmat)
factanal2nllk=function(rmat,mxfactor,n,iprint=F)
{ nllk=rep(0,mxfactor)
  d=nrow(rmat)
  for(k in 1:mxfactor) 
  { fa=factanal(covmat=rmat,factors=k)
    load=fa$loadings
    rk=load%*%t(load)
    diag(rk)=1
    if(iprint) cat("max abs diff ", max(abs(rmat-rk))," for ",k, "factors\n")
    tr=sum(diag(solve(rk,rmat)))
    nllk[k]=.5*n*(d*log(2*pi)+log(det(rk))+tr)
  }
  nllk
}

