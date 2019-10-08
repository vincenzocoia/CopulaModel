# functions for correlation matrix

# rvec = correlations in vector in the order r12,r13,r23,r14,...
# Output: dxd correlation matrix
corvec2mat=function(rvec)
{ dd=length(rvec)
  d=(1+sqrt(1+8*dd))/2
  rmat=matrix(1,d,d)
  k=0
  for(i in 2:d)
  { for(j in 1:(i-1))
    { k=k+1; rmat[i,j]=rvec[k]; rmat[j,i]=rvec[k] }
  }
  rmat
}

# rmat = dxd correlation matrix
# Output: vector from rmat  in the order r12,r13,r23,r14,...
cormat2vec=function(rmat)
{ d=nrow(rmat)
  dd=d*(d-1)/2
  rvec=NULL
  for(i in 2:d)
  rvec=c(rvec,rmat[i,1:(i-1)])
  rvec
}

# discrepancy of model-based correlation matrix and observed one based on
#  log-likelihood
# Rmod = model-based correlation matrix
# Robs = empirical correlation matrix (could be observed or polychoric)
# n = sample size (if positive integer)
# npar = #parameters in the correlation structure
# Output: discrepancy Dfit
#  also nllk, BIC, AIC if n and npar are inputted
corDis=function(Rmod,Robs,n=0,npar=0)
{ lgdetmod=log(det(Rmod))
  d=nrow(Robs)
  tem=sum(diag(solve(Rmod,Robs)))
  Dfit=lgdetmod-log(det(Robs))+tem-d
  if(n==0) return(Dfit)
  nllk2=n*d*(log(2*pi))+n*lgdetmod + n*tem
  bic=nllk2+log(n)*npar
  aic=nllk2+2*npar
  out=c(Dfit,nllk2,aic,bic)
  names(out)=c("Dfit","nllk2","AIC","BIC")
  out
}

