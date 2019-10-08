# R interface to pfactor/bifactor nllk+grad 

# rhvec = vector of length d*p with partial corr representation of loadings
# Robs = dxd empirical correlation matrix
# nsize = sample size
# Output: negative log-likelihood for Gaussian p-factor model
pfactnllk=function(rhvec,Robs,nsize)
{ d=nrow(Robs)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  p=length(rhvec)/d
  out= .Fortran("pfactnllk",
      as.integer(d), as.integer(p), as.double(rhvec), as.double(Robs),
      as.integer(nsize),
      nllk=as.double(0), grad=as.double(rep(0,d*p))  )
  nllk=out$nllk; 
  attr(nllk,"gradient") = out$grad;
  nllk
}

# rhvec = vector of length d*2 for partial corr representation of loadings,
#   first d correlations with common factor, then
#   partial correlations with group factor given common factor
# grsize = vector of group sizes for bi-factor model
# Robs = dxd empirical correlation matrix
# nsize = sample size
# Output: negative log-likelihood for Gaussian bi-factor model
bifactnllk=function(rhvec,grsize,Robs,nsize)
{ d=nrow(Robs)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  mgrp=length(grsize)
  out= .Fortran("bifactnllk",
      as.integer(d), as.integer(mgrp), as.double(rhvec), as.integer(grsize),
      as.double(Robs), as.integer(nsize),
      nllk=as.double(0), grad=as.double(rep(0,d*2))  )
  nllk=out$nllk; 
  attr(nllk,"gradient") = out$grad;
  nllk
}

# front end like factanal() 
# grsize = vector of group sizes for bi-factor model
# start = starting point should have dimension 2*d
# data = nsize x d data set to compute the correlation matrix if
#        cormat not given
# cormat = dxd empirical correlation matrix
# n = sample size
# prlevel = print.level for nlm()
# mxiter = maximum number of iterations for nlm()
# Output: a list with
# $nllk, $rhmat= dx2 matrix of correlations and partial correlations
factanal.bi=function(grsize,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
{ if(is.null(cormat))
  { n=nrow(data); d=ncol(data);
    cormat=cor(data)
  }
  else { d=sum(grsize) }
  if(length(start)!=2*d) { cat("start should have length 2*d\n"); return(0) }
  mle=nlm(bifactnllk,p=start,grsize=grsize,Robs=cormat,nsize=n,hessian=F,
    iterlim=mxiter,print.level=prlevel,check.analyticals=F)
  list(nllk=mle$minimum,rhmat=matrix(mle$estimate,d,2))
}

# front end like factanal() 
# factors = p = #factors
# start = starting point should have dimension 2*d
# data = nsize x d data set to compute the correlation matrix if
#        cormat not given
# cormat = dxd empirical correlation matrix
# n = sample size
# prlevel = print.level for nlm()
# mxiter = maximum number of iterations for nlm()
# Output: a list with
# $nllk, $rhmat = dxp matrix of partial correlations,
#  $loading = dxp loading matrix after varimax,
#  $rotmat = pxp rotation matrix used by varimax
factanal.co=function(factors,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
{ if(is.null(cormat))
  { n=nrow(data); d=ncol(data);
    cormat=cor(data)
  }
  else { d=nrow(cormat) }
  p=factors
  if(length(start)!=p*d) 
  { cat("start should have length factors*d\n"); return(0) }
  mle=nlm(pfactnllk,p=start,Robs=cormat,nsize=n,hessian=F,
    iterlim=mxiter,print.level=prlevel,check.analyticals=F)
  rhmat=matrix(mle$estimate,d,p)
  if(p>=2)
  { amat=pcor2load(rhmat)
    load=varimax(amat)
    loading=as.matrix(load$loadings[,1:p])
    rotmat=load$rotmat
  }
  else # count #negative loadings and reverse sign if needed
  { nneg=sum(rhmat<0)
    if(nneg>d/2) { rotmat=-1; loading=-rhmat } else { rotmat=1; loading=rhmat }
  }
  list(nllk=mle$minimum,rhmat=rhmat,loading=loading,rotmat=rotmat)
}

#============================================================

# tri-factor

# grsize = vector of group sizes for tri-factor model
# sbgrsize = vector of subgroup sizes for tri-factor model
# Output: TRUE if grsize,sbgrsize are consistent, FALSE otherwise
subgr.consistent=function(grsize,sbgrsize)
{ if(any(grsize<0) | any(sbgrsize<0)) 
  { cat("negative values not allowed\n"); return(FALSE) }
  if(any(grsize%%1!=0) | any(sbgrsize%%1!=0)) 
  { cat("non-integer values not allowed\n"); return(FALSE) }
  cmgr=cumsum(grsize)
  cmsbgr=cumsum(sbgrsize)
  ii=(cmgr %in% cmsbgr)
  (sum(ii)==length(grsize))
}

# R interface to trifactor nllk+grad 

# rhvec = vector of length d*2 for partial corr representation of loadings,
#   first d correlations with common factor, then
#   partial correlations with group factor given common factor
# grsize = vector of group sizes for tri-factor model
# sbgrsize = vector of subgroup sizes for tri-factor model
# Robs = dxd empirical correlation matrix
# nsize = sample size
# Output: negative log-likelihood for Gaussian tri-factor model
trifactnllk=function(rhvec,grsize,sbgrsize,Robs,nsize)
{ d=nrow(Robs)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  mgrp=length(grsize)
  msbgrp=length(sbgrsize)
  out= .Fortran("trifactnllk",
      as.integer(d), as.integer(mgrp),  as.integer(msbgrp),
      as.double(rhvec), as.integer(grsize), as.integer(sbgrsize),
      as.double(Robs), as.integer(nsize),
      nllk=as.double(0), grad=as.double(rep(0,d*3))  )
  nllk=out$nllk; 
  attr(nllk,"gradient") = out$grad;
  nllk
}

# front end like factanal() 
# grsize = vector of group sizes for tri-factor model
# sbgrsize = vector of subgroup sizes for tri-factor model
# start = starting point should have dimension 3*d
# data = nsize x d data set to compute the correlation matrix if
#        cormat not given
# cormat = dxd empirical correlation matrix
# n = sample size
# prlevel = print.level for nlm()
# mxiter = maximum number of iterations for nlm()
# Output: a list with
# $nllk, $rhmat= dx3 matrix of correlations and partial correlations
factanal.tri=function(grsize,sbgrsize,start,data=1,cormat=NULL,n=100,
  prlevel=0,mxiter=150)
{ 
  if(length(start)!=3*d) { cat("start should have length 3*d\n"); return(0) }
  if(!subgr.consistent(grsize,sbgrsize))
  { cat("grsize and sbgrsize not consistent\n"); return(-1.e10) }
  if(is.null(cormat))
  { n=nrow(data); d=ncol(data);
    cormat=cor(data)
  }
  else { d=sum(grsize) }
  mle=nlm(trifactnllk,p=start,grsize=grsize,sbgrsize=sbgrsize,
    Robs=cormat,nsize=n,hessian=F,
    iterlim=mxiter,print.level=prlevel,check.analyticals=F)
  list(nllk=mle$minimum,rhmat=matrix(mle$estimate,d,3))
}


