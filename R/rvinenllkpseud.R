# Functions to output pseudo-observations after tree ntrunc of a vine
#  for diagnostic checks for pair-copulas of next tree.

# Version that outputs pseudo-observations after tree ntrunc
#   can be used for sequential tree estimation
# parvec = vector of parameters to be optimized in nllk 
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# logdcopnames =  vector of names of logcopula densities for trees 1,...,ntrunc
# pcondnames = vector of names of cond cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric bivariate copula)
# np = dxd where np[ell,j] is #parameters for parameter th[ell,j]
#   for bivariate copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# Output:  a list with $nllk
#   condforw = matrix with C_{j|a_{ntrunc,j};S} in the forward direction
#   condbackw= matrix with C_{a_{ntrunc,j}|j;S} in the backward direction
#   The dimension of the matrices is nx(d-ntrunc).
rvinenllkpseud=function(parvec,udat,A,logdcopnames,pcondnames,np)
{ 
  ntrunc=length(logdcopnames) # assume >=1
  d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0; 
  # create (ragged) array th[1:np[ell,j],ell,j]
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=parvec[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  #icomp=out$icomp
  icomp=matrix(1,d,d)  # output all conditionals as pseudo-observations
  n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  logdcop=match.fun(logdcopnames[1])
  for(j in 2:d) 
  { ipp=1:np[1,j]
    nllk=nllk-sum(logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j]))
  }
  # define for tree 1 (this line not in rvinenllk1)
  w=v
  for(j in 2:d) { s[,j]=udat[,A[1,j]]; w[,j]=udat[,j] }
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
  #if(ipseud) 
  pcond=match.fun(pcondnames[ntrunc]) 
  ell=ntrunc+1
  for(j in ell:d) 
  { ipp=1:np[ell-1,j]
    vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
    v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
  }
  list(nllk=nllk, condforw=v[,ell:d], condbackw=vp[,ell:d])  
}

#============================================================

# This version allows different pair-copula family for each edge of the vine

# Pseudo-observations are output after tree ntrunc and
#    can be used for sequential tree estimation
# parvec = vector of parameters to be optimized in nllk 
# udat = nxd matrix with uniform scores
# A = dxd vine array with 1:d on diagonal
# ntrunc = truncated level, assume >=1
# logdcopmat =  matrix of names of logcopula densities for trees 1,...,ntrunc
# pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
#   logdcopmat and pcondmat are empty for diagonal and lower triangle,
#    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
# np = dxd where np[ell,j] is #parameters for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# UB,LB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output:  a list with $nllk
#   condforw = matrix with C_{j|a_{ntrunc,j};S} in the forward direction
#   condbackw= matrix with C_{a_{ntrunc,j}|j;S} in the backward direction
#   The dimension of the matrices is nx(d-ntrunc).
rvinenllkpseud2=function(parvec,udat,A,ntrunc,logdcopmat,pcondmat,np)
{ 
  d=ncol(A)  # or ncol(udat)
  npmax=max(np); th=array(0,c(npmax,d,d)) # format of parameter for vine
  ii=0; 
  # create (ragged) array th[1:np[ell,j],ell,j]
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ipp=1:np[ell,j]
      th[ipp,ell,j]=parvec[ii+ipp]
      ii=ii+np[ell,j]
    }
  }
  out=varray2M(A)
  M=out$mxarray
  #icomp=out$icomp
  icomp=matrix(1,d,d)  # output all conditionals as pseudo-observations
  n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d);
  nllk=0
  # tree 1
  for(j in 2:d) 
  { ipp=1:np[1,j]
    logdcop=match.fun(logdcopmat[1,j])
    nllk=nllk-sum(logdcop(udat[,A[1,j]],udat[,j],th[ipp,1,j]))
  }
  # define for tree 1 (this line not in rvinenllk1)
  w=v
  for(j in 2:d) { s[,j]=udat[,A[1,j]]; w[,j]=udat[,j] }
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
      w=v; wp=vp
    }
  }
  ell=ntrunc+1
  for(j in ell:d) 
  { ipp=1:np[ell-1,j]
    pcond=match.fun(pcondmat[ell-1,j]) 
    vp[,j]=pcond(s[,j],w[,j],th[ipp,ell-1,j]) 
    v[,j]=pcond(w[,j],s[,j],th[ipp,ell-1,j])
  }
  list(nllk=nllk, condforw=v[,ell:d], condbackw=vp[,ell:d])  
}

