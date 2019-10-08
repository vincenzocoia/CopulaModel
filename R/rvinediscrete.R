# nllk and joint pmf for R-vine with discrete/count variables

#============================================================

# D-vine (as check of R-vine code and implementation of algorithm in DMwC)
# D-vine probability for multivariate discrete 
# This assumes 1-dimensional parameter for each edge of vine.
# parmat = dxd matrix with parameters in standard D-vine
# u1vec = lower vector of hyperrectangle
# u2vec = upper vector of hyperrectangle
# pcopnames = string vector of names of pair-copula cdfs of length ntrunc
# iprint = print flag for intermediate steps in the probability calculations.
# Output: joint pmf of discrete D-vine
dvinepmf.discrete=function(parmat,u1vec,u2vec,pcopnames,iprint=F)
{ d=length(u1vec)
  d1=d-1
  upr=rep(0,d); Fp=rep(0,d); Fm=rep(0,d);
  Vp=rep(0,d); Vbp=rep(0,d); Vm=rep(0,d); Vbm=rep(0,d);
  vbpr=rep(0,d); vpr=rep(0,d)
  ovbpr=rep(0,d); ovpr=rep(0,d)
  bpp=rep(0,d); bpm=rep(0,d); bmp=rep(0,d); bmm=rep(0,d);
  mult=matrix(0,d,d)  # for f_{j-l..j}, j=1,...,d, l=0,...,d-1
  lk=0
  # univariate probabilities
  Fp=u2vec; Fm=u1vec; upr=Fp-Fm
  if(iprint) { print(u2vec); print(u1vec); print(upr) }
  mult[1,]=upr
  # tree 1
  pcop=match.fun(pcopnames[1])
  for(j in 2:d)
  { bpp[j]=pcop(Fp[j-1],Fp[j],parmat[1,j])
    bpm[j]=pcop(Fp[j-1],Fm[j],parmat[1,j])
    bmp[j]=pcop(Fm[j-1],Fp[j],parmat[1,j])
    bmm[j]=pcop(Fm[j-1],Fm[j],parmat[1,j])
  }
  for(j in 2:d) mult[2,j]=bpp[j]-bpm[j]-bmp[j]+bmm[j] # f_{j-1,j}
  lk=mult[2,2]   # f_{1,2}
  if(iprint) print(c(1,lk))
  # j=d not used for vb,  j=2 not used for v?
  for(j in 2:d) 
  { Vbp[j]=(bpp[j]-bpm[j])/upr[j] # F_{j-1|j}
    Vbm[j]=(bmp[j]-bmm[j])/upr[j]
    vbpr[j]=Vbp[j]-Vbm[j]  # f_{j-1|j}
    Vp[j]=(bpp[j]-bmp[j])/upr[j-1] # F_{j|j-1}
    Vm[j]=(bpm[j]-bmm[j])/upr[j-1]
    vpr[j]=Vp[j]-Vm[j]    # f_{j|j-1}
  }
  # tree 2 and on
  # main loop
  for(ell in 2:d1) # tree level
  { ell1=ell+1
    if(ell<=length(pcopnames)) { pcop=match.fun(pcopnames[ell]) }
    else { pcop=pindepcop }
    for(j in ell1:d) # Vb in  first arg, V in second arg
    { bpp[j]=pcop(Vbp[j-1],Vp[j],parmat[ell,j])  
      # C_{j-ell,j;}(F_{j-ell|j-ell+1,j-1},F_{j|.})
      bpm[j]=pcop(Vbp[j-1],Vm[j],parmat[ell,j])
      bmp[j]=pcop(Vbm[j-1],Vp[j],parmat[ell,j])
      bmm[j]=pcop(Vbm[j-1],Vm[j],parmat[ell,j])
      #if(iprint) cat("***", ell,j,Vbp[j-1],Vbm[j-1],Vp[j],Vm[j],"\n")
    }
    ovbpr=vbpr; ovpr=vpr
    for(j in ell1:d) 
    { Vbp[j]=(bpp[j]-bpm[j])/ovpr[j] # F_{j-ell|j,j-1,...j-ell+1} den prev vpr
      Vbm[j]=(bmp[j]-bmm[j])/ovpr[j]
      vbpr[j]=Vbp[j]-Vbm[j]  # f_{j-ell|j,j-1...j-ell+1}
      Vp[j]=(bpp[j]-bmp[j])/ovbpr[j-1] # F_{j|j-ell...j-1} den prev vbpr
      Vm[j]=(bpm[j]-bmm[j])/ovbpr[j-1]
      vpr[j]=Vp[j]-Vm[j]    # f_{j|j-ell,..j-1}
    }
    for(j in ell1:d) mult[ell1,j]=mult[ell,j-1]*vpr[j]
    lk=mult[ell1,ell1] # vpr[ell1]*mult[ell,ell]  # ell1=d at the end
    if(iprint) print(c(ell,lk))
  }
  # There needs to be a better way to handle this
  # Sometimes sum of enumeration is not 1 if there are boundary problems
  #  with pcop
  if(is.na(lk))
  { lk=0.;
    if(iprint) { print(parmat); print(mult); }
  }
  lk
}

# R-vine probability for multivariate discrete 
# This assumes 1-dimensional parameter for each edge of vine.
# parmat = dxd matrix with parameters in same position as vine array A
# u1vec = lower vector of hyperrectangle
# u2vec = upper vector of hyperrectangle
# A = dxd vine array with 1:d on diagonal
# M = dxd max array corresponding to A
# pcopnames = string vector of names of pair-copula cdfs of length ntrunc
# iprint = print flag for intermediate steps in the probability calculations.
# Output: joint pmf of discrete R-vine
rvinepmf.discrete=function(parmat,u1vec,u2vec,A,M,pcopnames,iprint=F)
{ d=length(u1vec)
  d1=d-1
  upr=rep(0,d); Fp=rep(0,d); Fm=rep(0,d);
  Vp=rep(0,d); Vbp=rep(0,d); Vm=rep(0,d); Vbm=rep(0,d);
  vbpr=rep(0,d); vpr=rep(0,d)
  ovbpr=rep(0,d); ovpr=rep(0,d)
  bpp=rep(0,d); bpm=rep(0,d); bmp=rep(0,d); bmm=rep(0,d);
  mult=matrix(0,d,d)  # for f_{j-l..j}, j=1,...,d, l=0,...,d-1
  lk=0
  # univariate probabilities
  #for(j in 1:d) { Fp[j]=u2vec[j]; Fm[j]=u1vec[j]; upr[j]=Fp[j]-Fm[j] } 
  Fp=u2vec; Fm=u1vec; upr=Fp-Fm
  if(iprint) { print(u2vec); print(u1vec); print(upr) }
  mult[1,]=upr
  # tree 1
  pcop=match.fun(pcopnames[1])
  for(j in 2:d)
  { bpp[j]=pcop(Fp[A[1,j]],Fp[j],parmat[1,j])
    bpm[j]=pcop(Fp[A[1,j]],Fm[j],parmat[1,j])
    bmp[j]=pcop(Fm[A[1,j]],Fp[j],parmat[1,j])
    bmm[j]=pcop(Fm[A[1,j]],Fm[j],parmat[1,j])
  }
  for(j in 2:d) mult[2,j]=bpp[j]-bpm[j]-bmp[j]+bmm[j] # f_{j-1,j}
  lk=mult[2,2]   # f_{1,2}
  if(iprint) print(c(1,lk))
  for(j in 2:d) 
  { Vbp[j]=(bpp[j]-bpm[j])/upr[j] # F_{a_{1j}|j}
    Vbm[j]=(bmp[j]-bmm[j])/upr[j]
    vbpr[j]=Vbp[j]-Vbm[j]  # f_{a_{1j}|j}
    Vp[j]=(bpp[j]-bmp[j])/upr[A[1,j]] # F_{j|a_{1j}}
    Vm[j]=(bpm[j]-bmm[j])/upr[A[1,j]]
    vpr[j]=Vp[j]-Vm[j]    # f_{j|a_{1j}}
  }
  # tree 2 and on
  # main loop
  for(ell in 2:d1) # tree level
  { ell1=ell+1
    if(ell<=length(pcopnames)) { pcop=match.fun(pcopnames[ell]) }
    else { pcop=pindepcop }
    for(j in ell1:d) # Vb in  first arg, V in second arg
    { sp=ifelse(A[ell,j]<M[ell,j],Vbp[M[ell,j]],Vp[M[ell,j]])
      sm=ifelse(A[ell,j]<M[ell,j],Vbm[M[ell,j]],Vm[M[ell,j]])
      bpp[j]=pcop(sp,Vp[j],parmat[ell,j])  
      # C_{a_{ell,j},j;}
      bpm[j]=pcop(sp,Vm[j],parmat[ell,j])
      bmp[j]=pcop(sm,Vp[j],parmat[ell,j])
      bmm[j]=pcop(sm,Vm[j],parmat[ell,j])
    }
    ovbpr=vbpr; ovpr=vpr
    for(j in ell1:d) 
    { 
      { Vbp[j]=(bpp[j]-bpm[j])/ovpr[j] # F_{a_{ell,j}|j,S} den prev vpr
        Vbm[j]=(bmp[j]-bmm[j])/ovpr[j]
        vbpr[j]=Vbp[j]-Vbm[j]  # f_{a_{ell,j}|j,S}
      }
      tt=ifelse(A[ell,j]<M[ell,j],ovbpr[M[ell,j]],ovpr[M[ell,j]])
      Vp[j]=(bpp[j]-bmp[j])/tt # F_{j|a_{ell,j},S}
      Vm[j]=(bpm[j]-bmm[j])/tt
      vpr[j]=Vp[j]-Vm[j]    # f_{j|a_{ell,j},S}
    }
    for(j in ell1:d) mult[ell1,j]=mult[ell,M[ell,j]]*vpr[j]
    lk=mult[ell1,ell1] # ell1=d at the end
    if(iprint) print(c(ell,lk))
  }
  # There needs to be a better way to handle this
  # Sometimes sum of enumeration is not 1 if there are boundary problems
  #  with pcop
  if(is.na(lk)) 
  { lk=0.; 
    if(iprint) { print(parmat); print(mult); }
  }
  lk
}

#============================================================

# R-vine nllk for multivariate discrete 
# This assumes 1-dimensional parameter for each edge of vine
# parvec = vector of parameters of pair-copula cdfs pcop,
#   current version: scalar parameter for each pair-copula
# uudat = matrix dimension nx(2d) with corners of rectangle 
# A = dxd vine array with 1:d on diagonal
# pcopnames = string vector of names of pair-copula cdfs of length ntrunc
# LB, UB = lower/upper bounds on parvec, could be scalars or 
#        same dimension as parvec
# iprint = print flag for intermediate steps in the probability calculations.
# Output: nllk of discrete R-vine
rvinediscretenllk=function(parvec,uudat,A,pcopnames,LB=0,UB=10,iprint=F)
{ ntrunc=length(pcopnames) # assume >=1
  d=ncol(uudat)/2
  n=nrow(uudat)
  if(any(parvec<=LB) | any(parvec>=UB) ) return(1.e10)
  parmat=matrix(0,d,d)
  ii=0
  for(ell in 1:ntrunc)
  { parmat[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  out=varray2M(A) # maximal array
  M=out$mxarray
  nllk=0
  jj1=1:d; jj2=(d+1):(2*d)
  for(i in 1:n)
  { pr=rvinepmf.discrete(parmat,uudat[i,jj1],uudat[i,jj2],A,M,pcopnames,iprint)
    if(pr<=0) pr=1.e-15
    nllk=nllk-log(pr)
  }
  nllk
}

#============================================================

# full log-likelihood : univariate parameters and dependence parameters of vine
# common parameters for univariate margins
# This assumes 1-dimensional parameter for each edge of vine
# parvec = vector of parameters of pair-copula cdfs pcop,
# ydat = (n*d)x1 vector of counts
# xdat = (n*d)xq matrix, n=#subjects, d=#repeated measurements, q=#covariates
#   xdat has d=nrep consecutive rows for each subject/unit
# nrep = d = number of repeated measurements per subject
# upmfcdf = function that computes the pmf,cdf up to the value mx
# npar1 = number of parameters in upmfcdf
# A = dxd vine array with 1:d on diagonal
# pcopnames = string vector of names of pair-copula cdfs of length ntrunc
# LB, UB = lower/upper bounds on vec, could be scalars or 
#        same dimension as parvec
# zero = either 0 or epsilon depending on the copula
# one = either 1 or 1-epsilon depending on the copula
# iprint = print flag for intermediate steps in the probability calculations.
# Output: full nllk of discrete R-vine model
rvinediscfullnllk=function(parvec,ydat,xdat,nrep,upmfcdf,npar1,
   A,pcopnames,LB=0,UB=10,zero=0.00001,one=0.99999,iprint=F)
{ if(any(parvec<=LB) | any(parvec>=UB)) return(1.e10)
  np=length(parvec)
  uparam=parvec[1:npar1]
  nrec=nrow(xdat)
  d=nrep
  n=nrec/d
  # step 1: get univariate cdfs
  u1mat=matrix(zero,n,d)
  u2mat=matrix(zero,n,d)
  ij=0
  for(i in 1:n)
  { for(j in 1:d)
    { ij=ij+1 
      x=xdat[ij,]
      y=ydat[ij] 
      out=upmfcdf(y,uparam,x)
      # cdf in column 3, first row is 0
      u2mat[i,j]=out[y+1,3]
      if(y>0) u1mat[i,j]=out[y,3]
    }
  }
  u2mat[u2mat>=1]=one
  # step 2: get vine probability 
  cparam=parvec[(npar1+1):np]
  ntrunc=length(pcopnames) # assume >=1
  parmat=matrix(0,d,d)
  ii=0
  for(ell in 1:ntrunc)
  { parmat[ell,(ell+1):d]=cparam[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  out=varray2M(A) # maximal array
  M=out$mxarray
  nllk=0
  for(i in 1:n)
  { pr=rvinepmf.discrete(parmat,u1mat[i,],u2mat[i,],A,M,pcopnames,iprint)
    if(pr<=0) pr=1.e-15
    nllk=nllk-log(pr)
  }
  nllk
}

# R-vine probability for multivariate ordinal, no covariates 
# parmat = dxd matrix with parameters in same position as vine array A
# yvec = integer-valued d-vector, values in 0:(ncateg-1) 
# A = dxd vine array with 1:d on diagonal
# M = dxd max array corresponding to A
# pcopnames = string vector of names of pair-copula cdfs of length ntrunc
#   copula for levels ntrunc+1,...,d-1 are assumed to be pindepcop
# ucuts =  (ncateg+1)xd matrix of cut points for ordinal, computed from unifcuts
#  via ucuts=unifcuts(y), ucuts=rbind(rep(0,d),ucuts,rep(1,d))
# iprint = print flag for intermediate steps in the probability calculations.
# Output: joint pmf of discrete R-vine with ordinal response
rvinepmf.ordinal=function(parmat,yvec,A,M,pcopnames,ucuts,iprint=F)
{ d=length(yvec)
  d1=d-1
  upr=rep(0,d); Fp=rep(0,d); Fm=rep(0,d);
  Vp=rep(0,d); Vbp=rep(0,d); Vm=rep(0,d); Vbm=rep(0,d);
  vbpr=rep(0,d); vpr=rep(0,d)
  ovbpr=rep(0,d); ovpr=rep(0,d)
  bpp=rep(0,d); bpm=rep(0,d); bmp=rep(0,d); bmm=rep(0,d);
  mult=matrix(0,d,d)  # for f_{j-l..j}, j=1,...,d, l=0,...,d-1
  lk=0
  # univariate probabilities
  for(j in 1:d) 
  { icateg=yvec[j]; Fp[j]=ucuts[icateg+2,j]; Fm[j]=ucuts[icateg+1,j]; 
    upr[j]=Fp[j]-Fm[j] 
  } 
  if(iprint) { print(yvec); print(upr) }
  mult[1,]=upr
  # tree 1
  pcop=match.fun(pcopnames[1])
  for(j in 2:d)
  { bpp[j]=pcop(Fp[A[1,j]],Fp[j],parmat[1,j])
    bpm[j]=pcop(Fp[A[1,j]],Fm[j],parmat[1,j])
    bmp[j]=pcop(Fm[A[1,j]],Fp[j],parmat[1,j])
    bmm[j]=pcop(Fm[A[1,j]],Fm[j],parmat[1,j])
  }
  for(j in 2:d) mult[2,j]=bpp[j]-bpm[j]-bmp[j]+bmm[j] # f_{j-1,j}
  lk=mult[2,2]   # f_{1,2}
  if(iprint) print(c(1,lk))
  # j=d not used for vb,  j=2 not used for v?
  for(j in 2:d) 
  { Vbp[j]=(bpp[j]-bpm[j])/upr[j] # F_{a_{1j}|j}
    Vbm[j]=(bmp[j]-bmm[j])/upr[j]
    vbpr[j]=Vbp[j]-Vbm[j]  # f_{a_{1j}|j}
    Vp[j]=(bpp[j]-bmp[j])/upr[A[1,j]] # F_{j|a_{1j}}
    Vm[j]=(bpm[j]-bmm[j])/upr[A[1,j]]
    vpr[j]=Vp[j]-Vm[j]    # f_{j|a_{1j}}
  }
  # tree 2 and on
  # main loop
  for(ell in 2:d1) # tree level
  { ell1=ell+1
    if(ell<=length(pcopnames)) { pcop=match.fun(pcopnames[ell]) }
    else { pcop=pindepcop }
    for(j in ell1:d) # Vb in  first arg, V in second arg
    { sp=ifelse(A[ell,j]<M[ell,j],Vbp[M[ell,j]],Vp[M[ell,j]])
      sm=ifelse(A[ell,j]<M[ell,j],Vbm[M[ell,j]],Vm[M[ell,j]])
      bpp[j]=pcop(sp,Vp[j],parmat[ell,j])  
      # C_{a_{ell,j},j;}
      bpm[j]=pcop(sp,Vm[j],parmat[ell,j])
      bmp[j]=pcop(sm,Vp[j],parmat[ell,j])
      bmm[j]=pcop(sm,Vm[j],parmat[ell,j])
    }
    ovbpr=vbpr; ovpr=vpr
    for(j in ell1:d) 
    { 
      { Vbp[j]=(bpp[j]-bpm[j])/ovpr[j] # F_{a_{ell,j}|j,S} den prev vpr
        Vbm[j]=(bmp[j]-bmm[j])/ovpr[j]
        vbpr[j]=Vbp[j]-Vbm[j]  # f_{a_{ell,j}|j,S}
      }
      tt=ifelse(A[ell,j]<M[ell,j],ovbpr[M[ell,j]],ovpr[M[ell,j]])
      Vp[j]=(bpp[j]-bmp[j])/tt # F_{j|a_{ell,j},S}
      Vm[j]=(bpm[j]-bmm[j])/tt
      vpr[j]=Vp[j]-Vm[j]    # f_{j|a_{ell,j},S}
    }
    for(j in ell1:d) mult[ell1,j]=mult[ell,M[ell,j]]*vpr[j]
    lk=mult[ell1,ell1] # vpr[ell1]*mult[ell,ell]  # ell1=d at the end
    if(iprint) print(c(ell,lk))
  }
  if(is.na(lk)) { print(parmat); print(mult); lk=0.; }
  #if(is.na(lk)) { lk=0.; }
  lk
}

# R-vine nllk for ordinal response and no covariates
# parvec = vector of parameter of pair-copula cdfs pcop,
#   current version: scalar parameter for each pair-copula
# ydat = integer-valued 0:(ncateg-1) or 1:(ncateg)
# ncateg = number of categories
# A = dxd vine array with 1:d on diagonal
# pcopnames = string vector of names of pair-copula cdfs of length ntrunc,
#   copula for levels ntrunc+1,...,d-1 are assumed to be pindepcop
# LB, UB = lower/upper bounds on vec, could be scalars or 
#        same dimension as parvec
# iprint = print flag for intermediate steps in the probability calculations.
# Output: nllk of ordinal R-vine model
rvineordinalnllk=function(parvec,ydat,ncateg,A,pcopnames,LB=0,UB=10,iprint=F)
{ ntrunc=length(pcopnames) # assume >=1
  d=ncol(ydat)
  if(min(ydat)==1) ydat=ydat-1  # convert to 0:(ncateg-1)
  if(any(parvec<=LB) | any(parvec>=UB) ) return(1.e10)
  th=matrix(0,d,d)
  ii=0
  for(ell in 1:ntrunc)
  { th[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  # may be better to do next few lines outside this function?
  odatfr=ordinal2fr(ydat,ncateg)
  ucuts=unifcuts(ydat)
  ucuts=rbind(rep(0,d),ucuts,rep(1,d))
  datfr=odatfr[,1:d]; fr=odatfr[,d+1]
  if(iprint) print(sum(fr))
  nrec=nrow(datfr)
  out=varray2M(A) # maximal array
  M=out$mxarray
  nllk=0
  for(i in 1:nrec)
  { if(i>3) iprint=F
    pr=rvinepmf.ordinal(th,datfr[i,],A,M,pcopnames,ucuts,iprint)
    if(pr<=0) pr=1.e-15
    nllk=nllk-fr[i]*log(pr)
  }
  nllk
}

