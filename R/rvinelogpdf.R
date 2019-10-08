# functions for vines : log pdf for continuous variables

# varray2M is in file varray.R

# log pdf for an R-vine density with
# same 1-parameter copula family for all edges; 
# this is version 1 for implementing and checking algorithms.
# uvec = d-vector with values in (0,1)
# A = dxd vine array with 1:d on diagonal
# parmat = dxd parameter array (1-parameter common copula for all edges)
# logdcop =log(dcop) = function for log copula density
# pcond = function for conditional cdf 
#    (only one needed for permutation symmetric pair-copula)
# iprint = print flag for intermediate calculations
# Output: log density of R-vine
rvinelogpdf=function(uvec,A,parmat,logdcop,pcond,iprint=F)
{ out=varray2M(A)
  M=out$mxarray
  icomp=out$icomp
  d=ncol(A)  # or length(uvec)
  d1=d-1
  v=rep(0,d); vp=rep(0,d)    
  s=rep(0,d);
  llk=0
  # tree 1
  for(j in 2:d) llk=llk+logdcop(uvec[A[1,j]],uvec[j],parmat[1,j])
  if(iprint) print(c(1,llk))
  # tree 2
  for(j in 2:d) 
  { if(icomp[1,j]==1) vp[j]=pcond(uvec[A[1,j]],uvec[j],parmat[1,j]) }
  for(j in 2:d) v[j]=pcond(uvec[j],uvec[A[1,j]],parmat[1,j])
  for(j in 3:d) s[j]=ifelse(A[2,j]<M[2,j],vp[M[2,j]],v[A[2,j]])  
  for(j in 3:d) llk=llk+logdcop(s[j],v[j],parmat[2,j])
  if(iprint) print(c(2,llk))
  w=v; wp=vp
  # remaining trees
  for(ell in 3:d1)
  { for(j in ell:d) 
    { if(icomp[ell-1,j]==1) vp[j]=pcond(s[j],w[j],parmat[ell-1,j]) }
    for(j in ell:d) v[j]=pcond(w[j],s[j],parmat[ell-1,j])
    for(j in (ell+1):d) s[j]=ifelse(A[ell,j]<M[ell,j],vp[M[ell,j]],v[A[ell,j]]) 
    for(j in (ell+1):d) llk=llk+logdcop(s[j],v[j],parmat[ell,j])
    if(iprint) print(c(ell,llk))
    w=v; wp=vp
  }
  llk
}

# C-vine (as check of R-vine code)
# this is like rvinelogpdf without array A as an argument
# uvec = d-vector with values in (0,1)
# parmat = dxd parameter array (1-parameter common copula for all edges)
# logdcop =log(dcop) = function for log copula density
# pcond = function for conditional cdf 
#    (only one needed for permutation symmetric pair-copula)
# iprint = print flag for intermediate calculations
# Output: log density of C-vine
cvinelogpdf=function(uvec,parmat,logdcop,pcond,iprint=F)
{ d=length(uvec)
  d1=d-1
  v=rep(0,d); 
  llk=0
  # tree 1
  for(j in 2:d) llk=llk+logdcop(uvec[1],uvec[j],parmat[1,j])
  if(iprint) print(c(1,llk))
  # tree 2
  for(j in 2:d) v[j]=pcond(uvec[j],uvec[1],parmat[1,j])
  for(j in 3:d) llk=llk+logdcop(v[2],v[j],parmat[2,j])
  if(iprint) print(c(2,llk))
  # remaining trees
  for(ell in 3:d1)
  { for(j in ell:d) v[j]=pcond(v[j],v[ell-1],parmat[ell-1,j])
    for(j in (ell+1):d) llk=llk+logdcop(v[ell],v[j],parmat[ell,j])
    if(iprint) print(c(ell,llk))
  }
  llk
}

# D-vine (as check of R-vine code)
# this is like rvinelogpdf without array A as an argument
# uvec = d-vector with values in (0,1)
# parmat = dxd parameter array (1-parameter common copula for all edges)
# logdcop =log(dcop) = function for log copula density
# pcond = function for conditional cdf 
#    (only one needed for permutation symmetric pair-copula)
# iprint = print flag for intermediate calculations
# Output: log density of D-vine
dvinelogpdf=function(uvec,parmat,logdcop,pcond,iprint=F)
{ d=length(uvec)
  d1=d-1
  v=rep(0,d); vp=rep(0,d)
  llk=0
  # tree 1
  for(j in 2:d) llk=llk+logdcop(uvec[j-1],uvec[j],parmat[1,j])
  if(iprint) print(c(1,llk))
  # tree 2
  for(j in 2:d1) vp[j]=pcond(uvec[j-1],uvec[j],parmat[1,j])
  for(j in 3:d) v[j]=pcond(uvec[j],uvec[j-1],parmat[1,j])
  for(j in 3:d) llk=llk+logdcop(vp[j-1],v[j],parmat[2,j])
  if(iprint) print(c(2,llk))
  w=v; wp=vp
  # remaining trees
  for(ell in 3:d1)
  { for(j in ell:d1) vp[j]=pcond(wp[j-1],w[j],parmat[ell-1,j])
    for(j in (ell+1):d) v[j]=pcond(w[j],wp[j-1],parmat[ell-1,j])
    for(j in (ell+1):d) llk=llk+logdcop(vp[j-1],v[j],parmat[ell,j])
    if(iprint) print(c(ell,llk))
    w=v; wp=vp
  }
  llk
}

#uvec=c(.1,.3,.4,.5,.7)
#parmat= matrix(c(0,0,0,0,0, 1.5,0,0,0,0, 1.5,1.2,0,0,0, 1.5,1.2,1.3,0,0, 1.5,1.2,1.3,1.4,0), 5,5)
# A=vbin2array(5,6)

# R-vine (not C or D)
#cat("\nR-vine\n")
#out2=rvinelogpdf(uvec,A,parmat,logdcop=logdgum,pcond=pcondgum,iprint=F)
#print(out2)

#cat("\nC-vine\n")
#out2=rvinelogpdf(uvec,C,parmat,logdcop=logdgum,pcond=pcondgum,iprint=F)
#out=cvinelogpdf(uvec,parmat,logdcop=logdgum,pcond=pcondgum,iprint=T)

#cat("\nD-vine\n")
#out2=rvinelogpdf(uvec,D,parmat,logdcop=logdgum,pcond=pcondgum,iprint=F)
#print(out2)
#out=dvinelogpdf(uvec,parmat,logdcop=logdgum,pcond=pcondgum,iprint=T)

