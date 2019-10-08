
# Mappings from correlation matrix to partial correlation vine
#  to regression coefficients

#============================================================

# partial correlation function for D-vine
# rr = positive definite dxd correlation matrix 
# Output: matrix with partial correlations in upper triangle 
cor2pcor.dvine=function(rr)
{ d=nrow(rr)
  if(d<=2) return(rr)
  pc=rr
  # tree 2
  for(j in 3:d) 
  { a1=j-1; a2=j-2 
    pc[j-2,j]=(rr[j,a2]-rr[j,a1]*rr[a1,a2])/sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
    pc[j,j-2]=pc[j-2,j]
  }
  # remaining trees
  for(ell in 3:(d-1))
  { for(j in (ell+1):d)
    { given=(j-1):(j-ell+1)
      pc[j-ell,j]=partcor(rr,given,j-ell,j)  
      pc[j,j-ell]=pc[j-ell,j]
    }
  }
  pc
}

# pc = dxd array with rho_{jk;(j+1):(k-1)} in position (j,k)
# Output: dxd correlation matrix
pcor2cor.dvine=function(pc)
{ d=nrow(pc)
  if(d<=2) { return(pc) }
  rr=matrix(0,d,d)
  diag(rr)=1
  for(j in 2:d) { a1=j-1; rr[a1,j]=pc[j-1,j]; rr[j,a1]=pc[j-1,j]; }
  # tree 2
  for(j in 3:d) 
  { a1=j-1; a2=j-2 
    rr[j,a2]=rr[j,a1]*rr[a1,a2]+pc[j-2,j]*sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
    rr[a2,j]=rr[j,a2]
  }
  # remaining trees
  for(ell in 3:(d-1))
  { for(j in (ell+1):d)
    { #given=A[1:(ell-1),j]
      given=(j-1):(j-ell+1)
      S11=rr[given,given]
      anew=j-ell
      jk=c(anew,j)
      S12=rr[given,jk]; S21=rr[jk,given]; S22=rr[jk,jk]
      tem=solve(S11,S12)
      Om212=S21%*%tem
      om11=1-Om212[1,1]; om22=1-Om212[2,2]
      tem12=pc[j-ell,j]*sqrt(om11*om22)
      rr[anew,j]=tem12+Om212[1,2]
      rr[j,anew]=rr[anew,j]
    }
  }
  rr
}


# test cases
#rr=toeplitz(c(1,.5,.25,.125))
#rr=toeplitz(c(1,.5,.25,.125,.1))
#pc=cor2pcor.dvine(rr)
#print(rr)
#print(pc)
#print(pcor2cor.dvine(pc))

#============================================================

# partial correlation function for C-vine
# rr = positive definite dxd correlation matrix 
# Output: matrix with partial correlations by rows in upper triangle
cor2pcor.cvine=function(rr)
{ d=nrow(rr)
  if(d<=2) return(rr)
  pp=matrix(0,d,d)
  for(i in 1:(d-1)) pp[i,(i+1):d]=rr[i,(i+1):d]
  for(m in 1:(d-2))
  { for(j in (m+1):(d-1))
    { for(k in (j+1):d)
      { pp[j,k]=(pp[j,k]-pp[m,j]*pp[m,k])/sqrt((1-pp[m,j]^2)*(1-pp[m,k]^2)) }
    }
  }
  pp
}

# generate correlation matrix based on partial correlations of C-vine
# pc = dxd array with partial cor rho_{jk;1:(j-1)} in position (j,k) for j<k
# Output: dxd correlation matrix
pcor2cor.cvine=function(pc)
{ d=nrow(pc)
  rr=matrix(0,d,d)
  diag(rr)=1
  rr[1,2:d]=pc[1,2:d]
  rr[2:d,1]=pc[1,2:d]
  for(j in 2:(d-1))
  { for(k in (j+1):d)
    { tem=pc[j,k]
      for(m in (j-1):1)
      { tem=pc[m,j]*pc[m,k]+tem*sqrt((1-pc[m,j]^2)*(1-pc[m,k]^2)) }
      rr[j,k]=tem; rr[k,j]=rr[j,k]
    }
  }
  rr
}

# example 
#rmat=matrix(c(1,.4,.5,.6,.7, .4,1,.2,.24,.28, .5,.2,1,.3,.35, .6,.24,.3,1,.42,
#      .7,.28,.35,.42,1),5,5)
#pmat=matrix(c(0,.4,.5,.6,.7, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0), 5,5,byrow=T)

#pcor2cor.cvine(pmat)
#cor2pcor.cvine(rmat)

#============================================================

# correlations to partial correlations or regression coefficients and vice
# versa for general R-vine with vine array A

# rr = positive definite dxd correlation matrix 
# A = dxd vine array with 1:d on diagonal (only upper triangle is used)
# Output: List with $pctree: tree ell in row ell
#         and $pcmat: rho_{jk;S} in position (j,k)
cor2pcor.rvine=function(rr,A)
{ d=nrow(rr)
  if(d<=2) return(rr)
  pp=matrix(0,d,d); pcmat=matrix(0,d,d); diag(pcmat)=1
  for(j in 2:d) pp[1,j]=rr[A[1,j],j]
  # tree 2
  for(j in 3:d) 
  { a1=A[1,j]; a2=A[2,j] 
    pp[2,j]=(rr[j,a2]-rr[j,a1]*rr[a1,a2])/sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
  }
  # remaining trees
  if(d>3)
  { for(ell in 3:(d-1))
    { for(j in (ell+1):d)
      { given=A[1:(ell-1),j]
        pp[ell,j]=partcor(rr,given,A[ell,j],j)  # assuming A[j,j]=j
      }
    }
  }
  for(ell in 1:(d-1))
  { for(j in (ell+1):d) 
    { alj=A[ell,j]; pcmat[alj,j]=pp[ell,j]; pcmat[j,alj]=pcmat[alj,j] }
  }
  list(pctree=pp, pcmat=pcmat)
}


# generate correlation matrix based on partial correlations of R-vine
# pc = dxd array with partial cor rho_{jk;S} for j<k
# A = dxd vine array with 1:d on diagonal (only upper triangle is used)
# if byrow=T, pc is matrix of d*(d-1)/2 partial correlations in upper triangle
#       with row matching tree
# if byrow=F, pc has rho_{jk;S} in position (j,k)
# Output: dxd correlation matrix
pcor2cor.rvine=function(pc,A,byrow=T)
{ d=nrow(A)
  if(d<=2) { rr=matrix(c(1,pc[1,2],pc[1,2],1)); return(rr) }
  if(!byrow) # convert to pc by tree (by row)
  { pp=matrix(0,d,d)
    for(ell in 1:(d-1))
    { for(j in (ell+1):d) 
      { alj=A[ell,j]; pp[ell,j]=pc[alj,j] }
    }
  }
  else { pp=pc }
  rr=matrix(0,d,d)
  diag(rr)=1
  for(j in 2:d) { a1=A[1,j]; rr[a1,j]=pp[1,j]; rr[j,a1]=pp[1,j]; }
  # tree 2
  for(j in 3:d) 
  { a1=A[1,j]; a2=A[2,j] 
    rr[j,a2]=rr[j,a1]*rr[a1,a2]+pp[2,j]*sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
    rr[a2,j]=rr[j,a2]
  }
  #print(rr)
  # remaining trees
  for(ell in 3:(d-1))
  { for(j in (ell+1):d)
    { given=A[1:(ell-1),j]
      S11=rr[given,given]
      anew=A[ell,j]
      jk=c(anew,j)
      S12=rr[given,jk]; S21=rr[jk,given]; #S22=rr[jk,jk]
      tem=solve(S11,S12)
      Om212=S21%*%tem
      om11=1-Om212[1,1]; om22=1-Om212[2,2]
      tem12=pp[ell,j]*sqrt(om11*om22)
      rr[anew,j]=tem12+Om212[1,2]
      rr[j,anew]=rr[anew,j]
      #print(rr)
    }
  }
  rr
}

#============================================================

# correlation matrix > matrix of regression coefficients for vine array A

# rr = dxd correlation matrix, 
# A = dxd vine array with 1:d on diagonal
# iprint = print flag for intermediate results
# Output: matrix of regression coefficients
cor2reg=function(rr,A,iprint=F)
{ d=nrow(A)
  phm=matrix(0,d,d)
  for(j in 2:d)
  { lh=matrix(0,j-1,j-1)
    for(k in  1:(j-1)) 
    { for(ell in 1:(j-1) ) lh[k,ell]=rr[A[ell,j],k] }
    if(iprint) print(lh)
    ph=solve(lh,rr[1:(j-1),j])
    phm[j,1:(j-1)]=ph
  }
  phm
}

# phm = lower triangular dxd matrix of regression coefficients
# A = dxd vine array with 1:d on diagonal
# Output: dxd correlation matrix
reg2cor=function(phm,A)
{ d=nrow(A)
  rr=matrix(0,d,d)
  diag(rr)=1
  for(j in 2:d)
  { for(k in  1:(j-1)) 
    { tem=0
      for(ell in 1:(j-1)) tem=tem+rr[A[ell,j],k]*phm[j,ell]
      rr[k,j]=tem; rr[j,k]=tem
    }
  }
  rr
}

#============================================================

# truncated vine: partial correlations to correlation matrix

# pp = dx1 or dxd with correlations for tree 1,
#     the correlations are in pp[2],...,pp[d];
# A = dxd vine array with 1:d on diagonal, but only the first row is used
# Output: dxd correlation matrix
pcor2cor.1tr=function(pp,A)
{ d=ncol(A)
  if(is.matrix(pp)) pp=pp[1,]
  if(length(pp)<d-1) pp=c(pp,rep(0,d-1-length(pp)))
  if(length(pp)==d-1) pp=c(0,pp)
  ph=pp
  rr=matrix(NA,d,d)
  diag(rr)=1
  for(j in 2:d) 
  { a1=A[1,j]; rr[a1,j]=pp[j]; rr[j,a1]=pp[j] } 
  # first row/column
  for(j in 2:d)
  { if(is.na(rr[1,j])) { a1=A[1,j]; rr[1,j]=ph[j]*rr[1,a1]; rr[j,1]=rr[1,j] }}
  # remaining correlations
  for(i in 2:(d-1))
  { for(j in (i+1):d)
    { if(is.na(rr[i,j])) 
      { a1=A[1,j]; rr[i,j]=ph[j]*rr[i,a1]; rr[j,i]=rr[i,j] }
    }
  }
  return(rr)
}

# pp = dx2 or dxd with correlations for tree 1, and partials for tree 2
#   the correlations of tree 1 are in pp[1,2],...,pp[1,d]
#   partial correlations of tree 2 are in pp[2,3],...,pp[2,d]
# A = dxd vine array with 1:d on diagonal but only the first two rows are used, d>=3
# Output: rmat=correlation matrix and phmat=matrix of regression coefficients
pcor2cor.2tr=function(pp,A)
{ d=ncol(A)
  if(!is.matrix(pp)) return(0)
  if(nrow(pp)<2 | ncol(pp)!=d) return(0)
  rr= matrix(NA,d,d)
  diag(rr)=1
  # correlations from tree 1
  for(j in 2:d) 
  { a1=A[1,j]; rr[a1,j]=pp[1,j]; rr[j,a1]=rr[a1,j] } 
  # correlations from tree 2
  for(j in 3:d) 
  { a1=A[1,j]; a2=A[2,j]; 
    rr[j,a2]=rr[j,a1]*rr[a1,a2] + pp[2,j]*sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
    rr[a2,j]=rr[j,a2]
  }
  # ph1, ph2 dx1 vectors of regression coefficients
  ph1=rep(0,d); ph2=rep(0,d);
  ph1[2]=rr[1,2]
  for(j in 3:d)
  { a1=A[1,j]; a2=A[2,j];
    ph2[j]=pp[2,j]*sqrt((1-rr[j,a1]^2)/(1-rr[a1,a2]^2))
    ph1[j]=rr[j,a1]-ph2[j]*rr[a1,a2]
  }
  # remaining correlations
  for(i in 1:(d-1))
  { for(j in (i+1):d)
    { if(is.na(rr[i,j])) 
      { a1=A[1,j]; a2=A[2,j];
        rr[i,j]=ph1[j]*rr[i,a1]+ph2[j]*rr[i,a2]; 
        rr[j,i]=rr[i,j] 
      }
    }
  }
  list(rmat=rr,phmat=cbind(ph1,ph2))
}


# pp = dx3 or dxd, with correlations for tree 1, 
#    and partials for trees 2,3
#   the correlations of tree 1 are in pp[1,2],...,pp[1,d]
#   partial correlations of tree 2 are in pp[2,3],...,pp[2,d]
#   partial correlations of tree 3 are in pp[3,4],...,pp[3,d]
# A = dxd vine array with 1:d on diagonal, but only the first 3 rows are used, d>=4
# iprint = print flag for intermediate results
# Output: rmat=correlation matrix and phmat=matrix of regression coefficients
pcor2cor.3tr=function(pp,A,iprint=F)
{ d=ncol(A)
  if(!is.matrix(pp)) return(0)
  if(nrow(pp)<3 | ncol(pp)!=d) return(0)
  rr= matrix(NA,d,d)
  diag(rr)=1
  # correlations from tree 1
  for(j in 2:d) 
  { a1=A[1,j]; rr[a1,j]=pp[1,j]; rr[j,a1]=rr[a1,j] } 
  # correlations from tree 2
  for(j in 3:d) 
  { a1=A[1,j]; a2=A[2,j]; 
    rr[a2,j]=rr[j,a1]*rr[a1,a2] + pp[2,j]*sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
    rr[j,a2]=rr[a2,j]
  }
  # correlations from tree 3; r[a1,a2],r[a1,a3],r[a2,a3] known after tree2
  #          r[a1,j],r[a2,j] known after tree2, need r[a3,j] 
  # loop ell=3 to m for m-truncation
  for(j in 4:d) 
  { a1=A[1,j]; a2=A[2,j]; a3=A[3,j]
    # r[a3,j] similar to code in pcor2cor.rvine
    given=c(a1,a2); S11=rr[given,given]; jk=c(a3,j)
    S12=rr[given,jk]; S21=rr[jk,given]; 
    tem=solve(S11,S12)
    Om212=S21%*%tem; om11=1-Om212[1,1]; om22=1-Om212[2,2]
    tem12=pp[3,j]*sqrt(om11*om22)
    rr[a3,j]=tem12+Om212[1,2]
    rr[j,a3]=rr[a3,j]
  }
  # ph1, ph2, ph3 dx1 vectors of regression coefficients
  ph1=rep(0,d); ph2=rep(0,d); ph3=rep(0,d);
  ph1[2]=rr[1,2]
  j=3; a1=A[1,j]; a2=A[2,j];
  ph2[j]=pp[2,j]*sqrt((1-rr[j,a1]^2)/(1-rr[a1,a2]^2))
  ph1[j]=rr[j,a1]-ph2[j]*rr[a1,a2]
  for(j in 4:d)
  { a1=A[1,j]; a2=A[2,j]; a3=A[3,j]
    vnum=1-rr[a1,a2]^2-rr[a1,j]^2-rr[a2,j]^2+2*rr[a1,a2]*rr[a1,j]*rr[a2,j]
    vden=1-rr[a1,a2]^2-rr[a1,a3]^2-rr[a2,a3]^2+2*rr[a1,a2]*rr[a1,a3]*rr[a2,a3]
    ph3[j]=pp[3,j]*sqrt(vnum/vden)
    pcnum=(rr[a2,a3]-rr[a1,a2]*rr[a1,a3])
    ph2[j]=pp[2,j]*sqrt((1-rr[j,a1]^2)/(1-rr[a1,a2]^2)) -
           ph3[j]*pcnum/(1-rr[a1,a2]^2)
    ph1[j]=rr[j,a1]-ph2[j]*rr[a1,a2]-ph3[j]*rr[a1,a3]
  }
  if(iprint) print(cbind(ph1,ph2,ph3))
  # remaining correlations
  for(i in 1:(d-1))
  { for(j in (i+1):d)
    { if(is.na(rr[i,j])) 
      { a1=A[1,j]; a2=A[2,j]; a3=A[3,j]
        rr[i,j]=ph1[j]*rr[i,a1]+ph2[j]*rr[i,a2]+ph3[j]*rr[i,a3]; 
        rr[j,i]=rr[i,j] 
      }
    }
  }
  list(rmat=rr,phmat=cbind(ph1,ph2,ph3))
}

# pp = dxntrunc or dxd, with correlations for tree 1, 
#    and partials for trees 2,...,ntrunc,  ntrunc>=1;
#   the correlations of tree 1 are in pp[1,2],...,pp[1,d],
#   partial correlations of tree ell are in pp[ell,ell+1],...,pp[ell,d].
#   for ell=2,...,ntrunc.
# A = dxd vine array with 1:d on diagonal, but only the first ntrunc rows 
#       are used for d>=ntrunc+1
# ntrunc = the truncation level
# iprint = print flag for intermediate results
# Output: rmat=correlation matrix 
#   phmat=matrix of regression coefficients 
# see vineresidvar() below  for output of psi2 = vector of residual variances: 
pcor2cor.truncvine=function(pp,A,ntrunc,iprint=F)
{ if(ntrunc<1) return(0) 
  if(ntrunc==1) return(pcor2cor.1tr(pp,A))
  d=ncol(A)
  rr=matrix(NA,d,d)
  diag(rr)=1
  for(j in 2:d) { a1=A[1,j]; rr[a1,j]=pp[1,j]; rr[j,a1]=pp[1,j]; }
  # tree 2
  for(j in 3:d) 
  { a1=A[1,j]; a2=A[2,j] 
    rr[j,a2]=rr[j,a1]*rr[a1,a2]+pp[2,j]*sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
    rr[a2,j]=rr[j,a2]
  }
  # step (i) trees 3 to ntrunc
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { for(j in (ell+1):d)
      { given=A[1:(ell-1),j]
        S11=rr[given,given]; anew=A[ell,j]; jk=c(anew,j)
        S12=rr[given,jk]; S21=rr[jk,given]; 
        tem=solve(S11,S12)
        Om212=S21%*%tem; om11=1-Om212[1,1]; om22=1-Om212[2,2]
        tem12=pp[ell,j]*sqrt(om11*om22)
        rr[anew,j]=tem12+Om212[1,2]
        rr[j,anew]=rr[anew,j]
      }
    }
  }
  # step (ii) solve for regression coefficient phi's  -- from cor2reg()
  phm=matrix(0,d,ntrunc)
  phm[1,1]=1
  phm[2,1]=rr[1,2]
  for(j in 3:d)
  { mn=min(ntrunc,j-1)
    apast=A[1:mn,j]
    lh=rr[apast,apast] 
    #if(iprint) print(lh)
    ph=solve(lh,rr[apast,j])
    if(iprint) print(ph)
    phm[j,1:mn]=ph
  }
  # step (iii) remaining correlations
  for(i in 1:(d-1))
  { for(j in (i+1):d)
    { if(is.na(rr[i,j])) 
      { apast=A[1:ntrunc,j]; 
        rr[i,j]=sum(phm[j,]*rr[i,apast]); 
        rr[j,i]=rr[i,j] 
      }
    }
  }
  list(rmat=rr,phmat=phm)
}

#============================================================

# rmatobj = object from  pcor2cor.truncvine()
#   with $rmat=correlation matrix and 
#      $phmat=matrix of regression coefficients 
# A = dxd vine array with 1:d on diagonal
# ntrunc = the truncation level
# Output: psi2 = vector of residual variances
vineResidVar=function(rmatobj,A,ntrunc) 
{ phmat=rmatobj$phmat
  rmat=rmatobj$rmat
  d=ncol(A)
  psi2=rep(1,d)
  for(j in 2:d)
  { j1=min(j-1,ntrunc)
    tem=0
    for(ell in 1:j1)
    { for(k in 1:j1) 
      { a1=A[ell,j]; a2=A[k,j]
        tem=tem+phmat[j,ell]*phmat[j,k]*rmat[a1,a2]
      }
    }
    psi2[j]=1-tem
  }
  psi2
}

