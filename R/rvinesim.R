# Functions for simulation from C-vine, D-vine and R-vine copulas

# The simpler versions of the code can be used to simulate one d-vector
# with U(0,1) margins; they are used for illustration of the algorithms.

# C-vine simulation as in Joe (2011, Chapter 7, Vine Copula Handbook),
# common copula family with scalar parameter for all edges of the vine.
# p = d-vector with values in (0,1)
# parmat = dxd parameter matrix, parmat[l,j] for tree l, variable j
# qcond = function for conditional quantile C_{2|1}^{-1}
# Output: vector of length d in (0,1)
cvinesim=function(p,parmat,qcond,pcond)
{ d=length(p)
  u=rep(0,d)
  u[1]=p[1]; 
  u[2]=qcond(p[2],p[1],parmat[1,2]); 
  # the main loop 
  for(j in 3:d)  # variable index
  { ttem=p[j]
    for(ell in seq(j-1,1)) # tree index
    { ttem=qcond(ttem, p[ell], parmat[ell,j]); }
    u[j]=ttem
  }
  u
}

# D-vine simulation,
# common copula family with scalar parameter for all edges of the vine.
# assumes pair-copula is permutation symmetric.
# p = d-vector with values in (0,1)
# parmat = dxd parameter matrix, parmat[l,j] for tree l, variable j
# pcond = function for conditional cdf C_{2|1}
# qcond = function for conditional quantile C_{2|1}^{-1}
# iprint = T to print intermediate calculations
# Output: vector of length d in (0,1)
dvinesim=function(p,parmat,qcond,pcond,iprint=F)
{ d=length(p)
  qq=matrix(0,d,d); v=matrix(0,d,d)
  u=rep(0,d)
  u[1]=p[1]; qq[1,1]=p[1]; qq[2,2]=p[2];
  u[2]=qcond(p[2],p[1],parmat[1,2]); qq[1,2]=u[2]
  v[1,2]=pcond(u[1],u[2],parmat[1,2])
  # the main loop 
  for(j in 3:d)  # variable index
  { qq[j,j]=p[j];
    for(ell in seq(j-1,2)) # tree index
    { qq[ell,j]=qcond(qq[ell+1,j], v[ell-1,j-1], parmat[ell,j]); }
    qq[1,j]=qcond(qq[2,j],u[j-1], parmat[1,j])
    u[j]=qq[1,j] 
    # set up for next iteration (not needed for last j=d)
    v[1,j]=pcond(u[j-1],u[j],parmat[1,j])
    for(ell in 2:(j-1))
    { v[ell,j]=pcond(v[ell-1,j-1], qq[ell,j], parmat[ell,j]); }
  }
  if(iprint) { print(qq); print(v) }
  u
}

# R-vine simulation with a common pcop (scalar parameter) for all edges
# version 0 with notation from recursion and use of icomp
# assumes pair-copula is permutation symmetric
# p = d-vector with values in (0,1)
# A = dxd vine array with 1:d on diagonal
# parmat = dxd parameter matrix, parmat[l,j] for tree l, variable j
# pcond = function for conditional cdf C_{2|1}
# qcond = function for conditional quantile C_{2|1}^{-1}
# iprint = print flag for intermediate calculations
# Output: vector of length d in (0,1)
rvinesim0=function(p,A,parmat,qcond,pcond,iprint=F)
{ out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  d=length(p)
  qq=matrix(0,d,d); v=matrix(0,d,d)
  u=rep(0,d)
  u[1]=p[1]; qq[1,1]=p[1]; qq[2,2]=p[2];
  u[2]=qcond(p[2],p[1],parmat[1,2]); qq[1,2]=u[2]
  if(icomp[1,2]==1) v[1,2]=pcond(u[1],u[2],parmat[1,2])
  # the main loop 
  for(j in 3:d)  # variable
  { qq[j,j]=p[j];
    for(ell in seq(j-1,2)) # tree
    { s=ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
      qq[ell,j]=qcond(qq[ell+1,j], s, parmat[ell,j]); 
    }
    qq[1,j]=qcond(qq[2,j],u[A[1,j]], parmat[1,j])
    u[j]=qq[1,j] 
    # set up for next iteration (not needed for last j=d)
    v[1,j]=pcond(u[A[1,j]],u[j],parmat[1,j])
    for(ell in 2:(j-1))
    { s=ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
      if(icomp[ell,j]==1) v[ell,j]=pcond(s, qq[ell,j], parmat[ell,j]); 
    }
  }
  if(iprint) { print(qq); print(v) }
  u
}

# R-vine simulation, different copula family for each tree.
# version 1 with notation from recursion and use of icomp,
# assumes pair-copulas are permutation symmetric
# p = vector of length d (independent random U(0,1))
# A = dxd vine array with 1:d on diagonal 
#      if truncated vine, only rows 1 to ntrunc are used.
# parmat = dxd parameter matrix, parmat[l,j] for tree l, variable j
# qcondnames = string vector of length ntrunc, 1<=ntrunc<=d-1, of inverse
#                conditional cdfs
# pcondnames = string vector of length ntrunc of conditional cdfs
#  In position ell, qcond and pcond should be functional inverses of each other
# iprint = print flag for intermediate calculations in the qq, v matrices
# Output: vector of length d in (0,1)
rvinesim1=function(p,A,parmat,qcondnames,pcondnames,iprint=F)
{ ntrunc=length(qcondnames)
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  d=length(p)
  qq=matrix(0,d,d); v=matrix(0,d,d)
  u=rep(0,d)
  u[1]=p[1]; qq[1,1]=p[1]; qq[2,2]=p[2];
  qcond=match.fun(qcondnames[1])
  u[2]=qcond(p[2],p[1],parmat[1,2]); qq[1,2]=u[2]
  if(icomp[1,2]==1) 
  { pcond=match.fun(pcondnames[1])
    v[1,2]=pcond(u[1],u[2],parmat[1,2])
  }
  # the main loop 
  for(j in 3:d)  # variable index
  { #qq[j,j]=p[j];
    #for(ell in seq(j-1,2)) # tree
    tt=min(ntrunc,j-1)
    qq[tt+1,j]=p[j]
    if(tt>1)
    { for(ell in seq(tt,2)) # tree index
      { s=ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
        qcond=match.fun(qcondnames[ell])
        qq[ell,j]=qcond(qq[ell+1,j], s, parmat[ell,j]); 
      }
    }
    qcond=match.fun(qcondnames[1])
    qq[1,j]=qcond(qq[2,j],u[A[1,j]], parmat[1,j])
    u[j]=qq[1,j] 
    # set up for next iteration (not needed for last j=d)
    pcond=match.fun(pcondnames[1])
    v[1,j]=pcond(u[A[1,j]],u[j],parmat[1,j])
    #for(ell in 2:(j-1))
    if(tt>1)
    { for(ell in 2:tt)
      { s=ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
        if(icomp[ell,j]==1) 
        { pcond=match.fun(pcondnames[ell])
          v[ell,j]=pcond(s, qq[ell,j], parmat[ell,j]); 
        }
      }
    }
  }
  if(iprint) { print(qq); print(v) }
  u
}

# checks
#pcond=pcondfrk; qcond=qcondfrk
#set.seed(123)
#d=6
#A=genVineArray(d)
# random theta parameters 
#th=matrix(0,d,d)
#for(j in 2:d)
#{ for(ell in 1:(j-1)) th[ell,j]=runif(1,1,4) }
#pp=runif(d)
#out=rvinesim0(pp,A,th,qcond,pcond,iprint=F)

#============================================================

# More flexible versions below

# Note that if pair-copula families are not permutation asymmetric,
# the code below needs to replace pcond (pcondnames) and qcond (qcondnames)
# with pcond21 (pcond21names) and pcond12 (pcond12names),
# and qcond21 (qcond21names) and qcond12 (qcond12names).
# Also small changes to the code, replacing pcond by either pcond21 or pcond12
# and replacing qcond by either qcond21 or qcond12.

# Version with dxd matrix np indicate parameter size for each edge of vine.
# Assumes pair-copulas are permutation symmetric.
# ntrunc = truncation level
# p = vector of length d (with independent random U(0,1))
# A = dxd vine array with 1:d on diagonal
#      if truncated vine, only rows 1 to ntrunc are used.
# parvec = vector of copula parameters; length is sum(np)
# np = dxd matrix with positive integers in the upper triangle
# np[ell,j] has number of parameters for tree ell, variable j
#    for j>ell, ell=1,...,ntrunc ; np=0 otherwise
# qcondnames = string vector of length ntrunc, 1<=ntrunc<=d-1, of inverse
#                conditional cdfs
# pcondnames = string vector of length ntrunc of conditional cdfs
#  In position ell, qcond and pcond should be functional inverses of each other
# iprint = T to print intermediate calculations in the qq, v matrices
# Output: vector of length d in (0,1)
rvinesim2=function(p,A,parvec,np,qcondnames,pcondnames,iprint=F)
{ ntrunc=length(qcondnames)
  d=ncol(A)
  # get matrix ip1,ip2 of indices
  ii=0
  ip1=matrix(0,d,d); ip2=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ip1[ell,j]=ii+1; ip2[ell,j]=ii+np[ell,j]
      ii=ii+np[ell,j]
    }
  }
  if(iprint) { print(ip1); print(ip2) }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  qq=matrix(0,d,d); v=matrix(0,d,d)
  u=rep(0,d)
  u[1]=p[1]; qq[1,1]=p[1]; qq[2,2]=p[2];
  qcond=match.fun(qcondnames[1])
  #u[2]=qcond(p[2],p[1],th[1,2]); 
  u[2]=qcond(p[2],p[1],parvec[ip1[1,2]:ip2[1,2]])
  qq[1,2]=u[2]
  if(icomp[1,2]==1) 
  { pcond=match.fun(pcondnames[1])
    #v[1,2]=pcond(u[1],u[2],th[1,2])
    v[1,2]=pcond(u[1],u[2],parvec[ip1[1,2]:ip2[1,2]])
  }
  # the main loop 
  for(j in 3:d)  # variable
  { #qq[j,j]=p[j];
    #for(ell in seq(j-1,2)) # tree
    tt=min(ntrunc,j-1)
    qq[tt+1,j]=p[j]
    if(tt>1)
    { for(ell in seq(tt,2)) # tree
      { s=ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
        qcond=match.fun(qcondnames[ell])
        #qq[ell,j]=qcond(qq[ell+1,j], s, th[ell,j]); 
        qq[ell,j]=qcond(qq[ell+1,j], s, parvec[ip1[ell,j]:ip2[ell,j]]); 
      }
    }
    qcond=match.fun(qcondnames[1])
    #qq[1,j]=qcond(qq[2,j],u[A[1,j]], th[1,j])
    qq[1,j]=qcond(qq[2,j],u[A[1,j]], parvec[ip1[1,j]:ip2[1,j]])
    u[j]=qq[1,j] 
    # set up for next iteration (not needed for last j=d)
    pcond=match.fun(pcondnames[1])
    #v[1,j]=pcond(u[A[1,j]],u[j],th[1,j])
    v[1,j]=pcond(u[A[1,j]],u[j],parvec[ip1[1,j]:ip2[1,j]])
    #for(ell in 2:(j-1))
    if(tt>1)
    { for(ell in 2:tt)
      { s=ifelse(A[ell,j]==M[ell,j], qq[ell,A[ell,j]], v[ell-1,M[ell,j]] )
        if(icomp[ell,j]==1) 
        { pcond=match.fun(pcondnames[ell])
          #v[ell,j]=pcond(s, qq[ell,j], th[ell,j]); 
          v[ell,j]=pcond(s, qq[ell,j], parvec[ip1[ell,j]:ip2[ell,j]]); 
        }
      }
    }
  }
  if(iprint) { print(qq); print(v) }
  u
}

#============================================================

# vectorized version of rvinesim2
# rvinesim2= function(p,A,parvec,np,qcondnames,pcondnames,iprint=F)
# vectorized version of above:
# Difference in input
#   nsim = #replications or sample size (replaces p)
# Difference in output 
#   matrix nsim x d : nsim replications of the given R-vine distribution
# ntrunc = truncation level
# nsim = simulation sample size
# A = dxd vine array with 1:d on diagonal
#      if truncated vine, only rows 1 to ntrunc are used.
# parvec = vector of copula parameters; length is sum(np)
# np = dxd matrix with positive integers in the upper triangle
# np[ell,j] has number of parameters for tree ell, variable j
#    for j>ell, ell=1,...,ntrunc ; np=0 otherwise
# qcondnames = string vector of length ntrunc, 1<=ntrunc<=d-1, of inverse
#                conditional cdfs
# pcondnames = string vector of length ntrunc of conditional cdfs
#  In position ell, qcond and pcond should be functional inverses of each other
# iprint = T to print intermediate calculations in the qq, v matrices
#  seed should be set before using this function.
# Output: nsim x d matrix with U(0,1) margins
rvinesimvec=function(nsim,A,parvec,np,qcondnames,pcondnames,iprint=F)
{ ntrunc=length(qcondnames)
  d=ncol(A)
  # get matrix ip1,ip2 of indices
  ii=0
  ip1=matrix(0,d,d); ip2=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ip1[ell,j]=ii+1; ip2[ell,j]=ii+np[ell,j]
      ii=ii+np[ell,j]
    }
  }
  if(iprint) { print(ip1); print(ip2) }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  p=matrix(runif(nsim*d),nsim,d)
  qq=array(0,c(nsim,d,d)); v=array(0,c(nsim,d,d))
  u=matrix(0,nsim,d)
  u[,1]=p[,1]; qq[,1,1]=p[,1]; qq[,2,2]=p[,2];
  qcond=match.fun(qcondnames[1])
  #u[2]=qcond(p[2],p[1],th[1,2]); 
  u[,2]=qcond(p[,2],p[,1],parvec[ip1[1,2]:ip2[1,2]])
  qq[,1,2]=u[,2]
  if(icomp[1,2]==1) 
  { pcond=match.fun(pcondnames[1])
    #v[1,2]=pcond(u[1],u[2],th[1,2])
    v[,1,2]=pcond(u[,1],u[,2],parvec[ip1[1,2]:ip2[1,2]])
  }
  # the main loop 
  for(j in 3:d)  # variable
  { tt=min(ntrunc,j-1)
    qq[,tt+1,j]=p[,j]
    if(tt>1)
    { for(ell in seq(tt,2)) # tree
      { if(A[ell,j]==M[ell,j]) { s= qq[,ell,A[ell,j]] } 
        else { s=v[,ell-1,M[ell,j]] }
        qcond=match.fun(qcondnames[ell])
        #qq[ell,j]=qcond(qq[ell+1,j], s, th[ell,j]); 
        qq[,ell,j]=qcond(qq[,ell+1,j], s, parvec[ip1[ell,j]:ip2[ell,j]]); 
      }
    }
    qcond=match.fun(qcondnames[1])
    #qq[1,j]=qcond(qq[2,j],u[A[1,j]], th[1,j])
    qq[,1,j]=qcond(qq[,2,j],u[,A[1,j]], parvec[ip1[1,j]:ip2[1,j]])
    u[,j]=qq[,1,j] 
    # set up for next iteration (not needed for last j=d)
    pcond=match.fun(pcondnames[1])
    #v[1,j]=pcond(u[A[1,j]],u[j],th[1,j])
    v[,1,j]=pcond(u[,A[1,j]],u[,j],parvec[ip1[1,j]:ip2[1,j]])
    #for(ell in 2:(j-1))
    if(tt>1)
    { for(ell in 2:tt)
      { if(A[ell,j]==M[ell,j]) { s=qq[,ell,A[ell,j]] }
        else { s=v[,ell-1,M[ell,j]] }
        if(icomp[ell,j]==1) 
        { pcond=match.fun(pcondnames[ell])
          #v[ell,j]=pcond(s, qq[ell,j], th[ell,j]); 
          v[,ell,j]=pcond(s, qq[,ell,j], parvec[ip1[ell,j]:ip2[ell,j]]); 
        }
      }
    }
  }
  u
}

#============================================================

# Simulation with different copula family for each edge
# nsim = #replications or sample size
# parvec = vector of parameters to be optimized in nllk 
# A = dxd vine array with 1:d on diagonal
# ntrunc = truncated level, assume >=1
# pcondmat = matrix of names of conditional cdfs for trees 1,...,ntrunc
#  (assuming only one needed for permutation symmetric pair-copulas)
#   logdcopmat and pcondmat are empty for diagonal and lower triangle,
#    and could have ntrunc rows or be empty for rows ntrunc+1 to d-1
# qcondmat = matrix of names of conditional quantile functions for 
#        trees 1,...,ntrunc
# np = dxd where np[ell,j] is #parameters for parameter th[ell,j]
#   for pair-copula in tree ell, variables j and A[ell,j] 
#   np=0 on and below diagonal
# Output: nsim x d matrix with U(0,1) margins
rvinesimvec2=function(nsim,A,ntrunc,parvec,np,qcondmat,pcondmat,iprint=F)
{ d=ncol(A)
  # get matrix ip1,ip2 of indices
  ii=0
  ip1=matrix(0,d,d); ip2=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { for(j in (ell+1):d)
    { ip1[ell,j]=ii+1; ip2[ell,j]=ii+np[ell,j]
      ii=ii+np[ell,j]
    }
  }
  if(iprint) { print(ip1); print(ip2) }
  out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  p=matrix(runif(nsim*d),nsim,d)
  qq=array(0,c(nsim,d,d)); v=array(0,c(nsim,d,d))
  u=matrix(0,nsim,d)
  u[,1]=p[,1]; qq[,1,1]=p[,1]; qq[,2,2]=p[,2];
  #qcond=match.fun(qcondnames[1])
  qcond=match.fun(qcondmat[1,2])
  u[,2]=qcond(p[,2],p[,1],parvec[ip1[1,2]:ip2[1,2]])
  qq[,1,2]=u[,2]
  if(icomp[1,2]==1) 
  { #pcond=match.fun(pcondnames[1])
    pcond=match.fun(pcondmat[1,2])
    v[,1,2]=pcond(u[,1],u[,2],parvec[ip1[1,2]:ip2[1,2]])
  }
  # the main loop 
  for(j in 3:d)  # variable
  { tt=min(ntrunc,j-1)
    qq[,tt+1,j]=p[,j]
    if(tt>1)
    { for(ell in seq(tt,2)) # tree
      { if(A[ell,j]==M[ell,j]) { s= qq[,ell,A[ell,j]] } 
        else { s=v[,ell-1,M[ell,j]] }
        #qcond=match.fun(qcondnames[ell])
        qcond=match.fun(qcondmat[ell,j])
        qq[,ell,j]=qcond(qq[,ell+1,j], s, parvec[ip1[ell,j]:ip2[ell,j]]); 
      }
    }
    #qcond=match.fun(qcondnames[1])
    qcond=match.fun(qcondmat[1,j])
    qq[,1,j]=qcond(qq[,2,j],u[,A[1,j]], parvec[ip1[1,j]:ip2[1,j]])
    u[,j]=qq[,1,j] 
    # set up for next iteration (not needed for last j=d)
    #pcond=match.fun(pcondnames[1])
    pcond=match.fun(pcondmat[1,j])
    v[,1,j]=pcond(u[,A[1,j]],u[,j],parvec[ip1[1,j]:ip2[1,j]])
    if(tt>1)
    { for(ell in 2:tt)
      { if(A[ell,j]==M[ell,j]) { s=qq[,ell,A[ell,j]] }
        else { s=v[,ell-1,M[ell,j]] }
        if(icomp[ell,j]==1) 
        { #pcond=match.fun(pcondnames[ell])
          pcond=match.fun(pcondmat[ell,j])
          v[,ell,j]=pcond(s, qq[,ell,j], parvec[ip1[ell,j]:ip2[ell,j]]); 
        }
      }
    }
  }
  u
}

#============================================================
