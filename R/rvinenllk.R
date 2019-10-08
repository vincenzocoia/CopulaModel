# Functions for negative log-likelihood and loglik vector at MLE for R-vine,
#  the loglik vector (for each observation) is used for Vuong's procedure

# Note that if pair-copula families are permutation asymmetric,
# the code below needs to replace pcond (pcondnames)
# with pcond21 (pcond21names) and pcond12 (pcond12names),
# with small changes to the code, replacing pcond by either pcond21 or pcond12

# parvec = vector of parameters to be optimized in nllk
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# logdcopnames = vector of names of log copula densities for trees 1,...,ntrunc
# pcondnames = vector of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
# np = dxd matrix where np[ell,j] is #parameters for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# ifixed = vector with length equal to length(parvec)+length(parfixed),
#  position is F if parameter is free to vary and T if fixed in the optimization
# parfixed = vector with dimension equal to sum(ifixed) for the positions where
#   the ifixed vector is T, parfixed[1] goes to the first position
#   where ifixed is T, etc.
# LB,UB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: negative log-likelihood for R-vine model
rvinenllk.trunc=function(parvec,udat,A,logdcopnames,pcondnames,np,ifixed,
  parfixed, LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) return(1.e10)
  ntrunc=length(logdcopnames) # assume >=1
  if(sum(ifixed)>0)
  { npar=length(ifixed); param0=rep(0,npar)
    ii=1; jj=1
    for(ip in 1:npar)
    { if(!ifixed[ip]) { param0[ip]=parvec[ii]; ii=ii+1 } 
      else { param0[ip]=parfixed[jj]; jj=jj+1 }
    }
    #print(param0)
  }
  else { param0=parvec } 
  d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0;
  # create (ragged) array th[1:np[ell,j],ell,j]
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=param0[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  logdcop=match.fun(logdcopnames[1])
  for(j in 2:d) 
  { ipp=1:np[1,j]
    nllk=nllk-sum(logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j]))
  }
  #if(iprint) print(c(1,llk))
  # tree 2
  if(ntrunc>=2)
  { pcond=match.fun(pcondnames[1])
    for(j in 2:d) 
    { ipp=1:np[1,j]
      if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[ipp,1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,1,j])
    }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    logdcop=match.fun(logdcopnames[2])
    for(j in 3:d) 
    { ipp=1:np[2,j]
      nllk=nllk-sum(logdcop(s[,j],v[,j],th[ipp,2,j]))
    }
    #if(iprint) print(c(2,llk))
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcond=match.fun(pcondnames[ell-1]) 
      logdcop=match.fun(logdcopnames[ell]) 
      #for(j in ell:d) 
      #{ if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ell-1,j]) }
      #for(j in ell:d) v[,j]=pcond(w[,j],s[,j],th[ell-1,j])
      for(j in ell:d) 
      { ipp=1:np[ell-1,j]
        if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) 
      { ipp=1:np[ell,j]
        nllk=nllk-sum(logdcop(s[,j],v[,j],th[ipp,ell,j]))
      }
      #if(iprint) print(c(ell,llk))
      w=v; wp=vp
    }
  }
  nllk
}

# This is a simpler version of the R-vine negative log-likelihood where 
#  there is a scalar parameter for each pair-copula.
# parvec = vector of parameters of pair-copulas
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# logdcopnames = vector of names of log copula densities for trees 1,...,ntrunc
# pcondnames = vector of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
# LB,UB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: negative log-likelihood for R-vine model
rvinenllk1.trunc=function(parvec,udat,A,logdcopnames,pcondnames,LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) return(1.e10)
  ntrunc=length(logdcopnames) # assume >=1
  d=ncol(A)  # or ncol(udat)
  ii=0; th=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { th[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  logdcop=match.fun(logdcopnames[1])
  for(j in 2:d) nllk=nllk-sum(logdcop(udat[,A[1,j]],udat[,j],th[1,j]))
  # tree 2
  if(ntrunc>=2)
  { pcond=match.fun(pcondnames[1])
    for(j in 2:d) 
    { if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[1,j])
    }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    logdcop=match.fun(logdcopnames[2])
    for(j in 3:d) nllk=nllk-sum(logdcop(s[,j],v[,j],th[2,j]))
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcond=match.fun(pcondnames[ell-1]) 
      logdcop=match.fun(logdcopnames[ell]) 
      for(j in ell:d) 
      { if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ell-1,j])
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) nllk=nllk-sum(logdcop(s[,j],v[,j],th[ell,j]))
      w=v; wp=vp
    }
  }
  nllk
}

# This is version of rvine log-likelihood vector where there is a scalar
#   parameter for each pair-copula.
# parvec = vector of parameters of pair-copulas
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# logdcopnames = vector of names of log copula densities for trees 1,...,ntrunc
# pcondnames = vector of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
# Output: log-likelihood vector for R-vine model (for Vuong's procedure)
rvinellkv1.trunc= function(parvec,udat,A,logdcopnames,pcondnames)
{ ntrunc=length(logdcopnames) # assume >=1
  d=ncol(A)  # or ncol(udat)
  ii=0; th=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { th[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  llkv=rep(0,n)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  logdcop=match.fun(logdcopnames[1])
  for(j in 2:d) llkv=llkv+logdcop(udat[,A[1,j]],udat[,j],th[1,j])
  # tree 2
  if(ntrunc>=2)
  { pcond=match.fun(pcondnames[1])
    for(j in 2:d) 
    { if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[1,j])
    }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    logdcop=match.fun(logdcopnames[2])
    for(j in 3:d) llkv=llkv+logdcop(s[,j],v[,j],th[2,j])
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcond=match.fun(pcondnames[ell-1]) 
      logdcop=match.fun(logdcopnames[ell]) 
      for(j in ell:d) 
      { if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ell-1,j])
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) llkv=llkv+logdcop(s[,j],v[,j],th[ell,j])
      w=v; wp=vp
    }
  }
  llkv
}

# parvec = vector of parameters for pair-copulas
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# logdcopnames = vector of names of log copula densities for trees 1,...,ntrunc
# pcondnames = vector of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
# np = dxd matrix where np[ell,j] is size for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# Output: log-likelihood vector for R-vine model (for Vuong's procedure)
rvinellkv.trunc=function(parvec,udat,A,logdcopnames,pcondnames,np)
{ ntrunc=length(logdcopnames) # assume >=1
  d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0;
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=parvec[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  llkv=rep(0,n)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  logdcop=match.fun(logdcopnames[1])
  for(j in 2:d) 
  { ipp=1:np[1,j]
    llkv=llkv+logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j])
  }
  # tree 2
  if(ntrunc>=2)
  { pcond=match.fun(pcondnames[1])
    for(j in 2:d) 
    { ipp=1:np[1,j]
      if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[ipp,1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,1,j])
    }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    logdcop=match.fun(logdcopnames[2])
    for(j in 3:d) 
    { ipp=1:np[2,j]
      llkv=llkv+logdcop(s[,j],v[,j],th[ipp,2,j])
    }
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcond=match.fun(pcondnames[ell-1]) 
      logdcop=match.fun(logdcopnames[ell]) 
      for(j in ell:d) 
      { ipp=1:np[ell-1,j] 
        if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) 
      { ipp=1:np[ell,j] 
        llkv=llkv+logdcop(s[,j],v[,j],th[ipp,ell,j])
      }
      w=v; wp=vp
    }
  }
  llkv
}

#============================================================

# truncated R vine where simplifying assumption need not hold,
# This is an illustration/template of how above R-vine likelihood could be 
#  modified for non-simplfiying assumption.
# Regression on sum of values of conditioning variables
# for tree 2 and above let log(cpar)=b0+b1*ugiven 
# 2-trunc: number of parameters is (d-1)+(d-2)*2
# 3-trunc: number of parameters is (d-1)+(d-2)*2+(d-3)*2 etc
# parvec = vector of parameters to be optimized in nllk 
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# logdcopnames = vector of names of log copula densities for trees 1,...,ntrunc
# pcondnames = vector of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
# ib0fixed = flag for whether intercept b0 is fixed
# bOfixed = fixed intercept when copula parameter has form 
#  cpar=exp(b0+b1*sum(ucond)) for trees 2,3.. 
# iprint = print flag for intermediate calculations
# LB,UB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: negative log-likelihood for R-vine model
rvinenllk1.nonsimpl=function(parvec,udat,A,logdcopnames,pcondnames,ib0fixed=F,
  b0fixed=0,iprint=F,LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) return(1.e10)
  ntrunc=length(logdcopnames) # assume >=1
  d=ncol(A)  # or ncol(udat)
  ii=0; th=matrix(0,d,d); b0=matrix(0,d,d); b1=matrix(0,d,d);
  th[1,2:d]=parvec[1:(d-1)]
  ii=d-1
  if(ntrunc>=2)
  { if(!ib0fixed)
    { for(ell in 2:ntrunc)
      { b0[ell,(ell+1):d]=parvec[ii+seq(1,by=2,length=d-ell)] 
        b1[ell,(ell+1):d]=parvec[ii+seq(2,by=2,length=d-ell)] 
        ii=ii+2*(d-ell)
      }
    }
    else # fixed b0s
    { for(ell in 2:ntrunc)
      { b0[ell,(ell+1):d]=b0fixed
        b1[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
        ii=ii+(d-ell)
      }
    }
    if(iprint) { print(b0); print(b1) }
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  d1=d-1
  n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  logdcop=match.fun(logdcopnames[1])
  for(j in 2:d) nllk=nllk-sum(logdcop(udat[,A[1,j]],udat[,j],th[1,j]))
  # tree 2
  if(ntrunc>=2)
  { pcond=match.fun(pcondnames[1])
    for(j in 2:d) 
    { if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[1,j]) }
    for(j in 2:d) v[,j]=pcond(udat[,j],udat[,A[1,j]],th[1,j])
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    logdcop=match.fun(logdcopnames[2])
    for(j in 3:d) 
    { cpar=exp(b0[2,j]+b1[2,j]*udat[,A[1,j]])
      nllk=nllk-sum(logdcop(s[,j],v[,j],cpar))
    }
    w=v; wp=vp
  }
  # remaining trees 
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcond=match.fun(pcondnames[ell-1]) 
      logdcop=match.fun(logdcopnames[ell]) 
      for(j in ell:d) 
      { cpar=exp(b0[ell-1,j]+b1[ell-1,j]*sum(udat[,A[1:(ell-2),j]]))
        if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],cpar) 
        v[,j]=pcond(w[,j],s[,j],cpar)
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) 
      { cpar=exp(b0[ell,j]+b1[ell,j]*sum(udat[,A[1:(ell-1),j]]))
        nllk=nllk-sum(logdcop(s[,j],v[,j],cpar))
      }
      w=v; wp=vp
    }
  }
  nllk
}


# ============================================================

# Versions of -loglik and loglik vector for R-vine
#   with different pair-copula family for each edge of the vine

# parvec = vector of parameters of pair-copulas
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# ntrunc = truncated level, assume >=1
# logdcopmat =  matrix of names of log copula densities for trees 1,...,ntrunc
# pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
#   logdcopmat and pcondmat are empty for diagonal and lower triangle,
#    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
# np = dxd where np[ell,j] is #parameters for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# ifixed = vector with length equal to length(parvec)+length(parfixed),
#   F is parameter is free to vary and T if fixed in the optimization
# parfixed = vector with dimension equal to sum(ifixed) for the positions where
#   the ifixed vector is T, parfixed[1] goes to the first position
#   where ifixed is T, etc.
# UB,LB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: negative log-likelihood for R-vine model
rvinenllk.trunc2=function(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np,
  ifixed,parfixed, LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) return(1.e10)
  if(sum(ifixed)>0)
  { npar=length(ifixed); param0=rep(0,npar)
    ii=1; jj=1
    for(ip in 1:npar)
    { if(!ifixed[ip]) { param0[ip]=parvec[ii]; ii=ii+1 } 
      else { param0[ip]=parfixed[jj]; jj=jj+1 }
    }
    #print(param0)
  }
  else { param0=parvec } 
  d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0;
  # create (ragged) array th[1:np[ell,j],ell,j]
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=param0[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  for(j in 2:d) 
  { ipp=1:np[1,j]
    logdcop=match.fun(logdcopmat[1,j])
    nllk=nllk-sum(logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j]))
  }
  #if(iprint) print(c(1,llk))
  # tree 2
  if(ntrunc>=2)
  { for(j in 2:d) 
    { ipp=1:np[1,j]
      pcond=match.fun(pcondmat[1,j])
      if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[ipp,1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,1,j])
    }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    for(j in 3:d) 
    { ipp=1:np[2,j]
      logdcop=match.fun(logdcopmat[2,j])
      nllk=nllk-sum(logdcop(s[,j],v[,j],th[ipp,2,j]))
    }
    #if(iprint) print(c(2,llk))
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { for(j in ell:d) 
      { ipp=1:np[ell-1,j]
        pcond=match.fun(pcondmat[ell-1,j]) 
        if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) 
      { ipp=1:np[ell,j]
        logdcop=match.fun(logdcopmat[ell,j]) 
        nllk=nllk-sum(logdcop(s[,j],v[,j],th[ipp,ell,j]))
      }
      #if(iprint) print(c(ell,llk))
      w=v; wp=vp
    }
  }
  nllk
}

# parvec = vector of parameters of pair-copulas
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# ntrunc = truncated level, assume >=1
# logdcopmat =  matrix of names of log copula densities for trees 1,...,ntrunc
# pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
#   logdcopmat and pcondmat are empty for diagonal and lower triangle,
#    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
# np = dxd where np[ell,j] is size for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# Output: log-likelihood vector for R-vine model (for Vuong's procedure)
rvinellkv.trunc2=function(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np)
{ d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0;
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=parvec[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  n=nrow(udat)
  llkv=rep(0,n)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  for(j in 2:d) 
  { ipp=1:np[1,j]
    logdcop=match.fun(logdcopmat[1,j])
    llkv=llkv+logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j])
  }
  # tree 2
  if(ntrunc>=2)
  { for(j in 2:d) 
    { ipp=1:np[1,j]
      pcond=match.fun(pcondmat[1,j])
      if(icomp[1,j]==1) vp[,j]=pcond(udat[,A[1,j]],udat[,j],th[ipp,1,j]) 
      v[,j]=pcond(udat[,j],udat[,A[1,j]],th[ipp,1,j])
    }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] } 
    for(j in 3:d) 
    { ipp=1:np[2,j]
      logdcop=match.fun(logdcopmat[2,j])
      llkv=llkv+logdcop(s[,j],v[,j],th[ipp,2,j])
    }
    w=v; wp=vp
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { for(j in ell:d) 
      { ipp=1:np[ell-1,j] 
        pcond=match.fun(pcondmat[ell-1,j]) 
        if(icomp[ell-1,j]==1) vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
        v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] } 
      for(j in (ell+1):d) 
      { ipp=1:np[ell,j] 
        logdcop=match.fun(logdcopmat[ell,j]) 
        llkv=llkv+logdcop(s[,j],v[,j],th[ipp,ell,j])
      }
      w=v; wp=vp
    }
  }
  llkv
}

