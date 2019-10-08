# discrete R-vine with bivariate normal/Gaussian for pair-copulas

# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ncl = #clusters, each with d repeated measures
# zzdat = matrix of dimension nclx(2d) with corners of rectangle, N(0,1) scale
# A = dxd vine array with 1:d on diagonal 
# Output: negative log-likelihood
rvinediscbvnnllk=function(parvec,zzdat,A)
{ d=ncol(zzdat)/2
  ncl=nrow(zzdat)
  out=varray2M(A) # maximal array
  M=out$mxarray
  if(any(parvec<=-1) | any(parvec>=1) ) return(1.e10)
  out= .Fortran("rvinediscbvnnllk",
    as.integer(ncl), as.integer(d), as.double(parvec), 
    as.double(zzdat),
    as.integer(A), as.integer(M), nllk=as.double(0) )
  nllk=out$nllk 
  nllk
}

# Full log-likelihood, ordinal response with covariates
#    univariate regression parameters estimated
# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ncl = #clusters, each with d=nrep repeated ordinal measures
# A = dxd vine array with 1:d on diagonal 
# xmat = nn*npred covariate matrix, nn = nrep*ncl
# yvec = integer-valued vector (nnx1) : with values in 0:(ncateg-1) or 1:ncateg
# nrep = d = cluster size or #repeated ordinal measures
# ncateg = number of ordinal categories
# Output: negative log-likelihood
rvinediscbvnfullnllk=function(parvec,A,xmat,yvec,nrep,ncateg)
{ if(is.matrix(xmat)) { nn=nrow(xmat); npred=ncol(xmat) }
  else { nn=length(xmat); npred=1 }
  d=nrep
  #ncl=nn/d  # number of clusters
  #np=length(parvec)
  out=varray2M(A) # maximal array
  M=out$mxarray
  # f90 code assumes 1:ncateg
  if(min(yvec)==0) { yvec=yvec+1 }
  out= .Fortran("rvineOrdwPrednllk",
    as.integer(nn), as.integer(d), as.integer(npred), as.integer(ncateg),
    as.double(parvec), as.integer(A), as.integer(M), 
    as.double(xmat), as.integer(yvec), nllk=as.double(0) )
  nllk=out$nllk 
  nllk
}

