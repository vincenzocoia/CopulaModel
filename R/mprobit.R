# Rectangle probabilities for multivariate ordinal probit.
# KL divergence and sample size for multivariate probit
#   versus discrete D-vine and R-vine 

# This function depends on library(mvtnorm).
# It can be used for dimensions d=3,4,5.
# MVN rectangle probabilities for multivariate ordinal.
# zcuts = (ncateg+1)xd matrix of cutpts in N(0,1) scale for ordinal responses
# rmat = latent correlation matrix
# iprint = print flag for intermediate calculations and some diagnostics
# ifixseed = flag for fixed seed, 
#   if T: seed is 12345 before each call to pmvnorm with algorithm=GenzBretz
# Output : pmf vector in ordinal categories in lexicographic order.
pmfmordprobit=function(zcuts,rmat,iprint=F,ifixseed=F)
{ d=ncol(zcuts)
  ncateg=nrow(zcuts)-1
  nn=ncateg^d
  pr=rep(0,nn)
  lb=rep(0,d)
  ub=rep(0,d)
  for(i in 1:nn)
  { jj=d2v(d,ncateg,i-1)
    # categories are indexed as 0, 1, ... n^categ-1 in d2v()
    for(j in 1:d) { lb[j]=zcuts[jj[j],j]; ub[j]=zcuts[jj[j]+1,j] }
    # option to fix seed
    if(ifixseed) set.seed(12345)
    pr[i]=pmvnorm(lb,ub,mean=rep(0,d),corr=rmat,algorithm=GenzBretz())
    if(pr[i]<=0) pr[i]=1.e-10
    if(iprint) print(c(i,pr[i]))
  }
  if(iprint) print(sum(pr))
  pr
}

# KL divergence for discrete D-vine model
# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ucuts = (ncateg+1)xd matrix of cutpts in U(0,1) scale for ordinal responses,
#   includes the boundary cutpts.
# pr = vector from pmfmordprobit()
# Output: KL divergence
dvineKLfn=function(parvec,ucuts,pr)
{ if(any(parvec>=1) | any(parvec<=-1)) return(1.e10)
  d=ncol(ucuts)
  ncateg=nrow(ucuts)-1
  D=Dvinearray(d)
  out=varray2M(D)
  M=out$mxarray
  pcopnames=rep("pbvncop",d-1)
  nn=ncateg^d
  prv=rep(0,nn)
  u1vec=rep(0,d)
  u2vec=rep(0,d)
  parmat=matrix(0,d,d)
  d2=(d*(d-1))/2
  ii=0
  for(ell in 1:(d-1))
  { for(j in (ell+1):d)
    { ii=ii+1; parmat[ell,j]=parvec[ii] }
  }
  kl=0
  for(i in 1:nn)
  { jj=d2v(d,ncateg,i-1)
    for(j in 1:d) 
    { u1vec[j]=ucuts[jj[j],j]; u2vec[j]=ucuts[jj[j]+1,j] }
    prv[i]=rvinepmf.discrete(parmat,u1vec,u2vec,D,M,pcopnames,iprint=F)
    if(prv[i]<=0) prv[i]=1.e-10
    if(is.na(prv[i]) | is.infinite(prv[i])) return(1.e10)
    kl=kl+pr[i]*log(pr[i]/prv[i])
  }
  kl
}

# KL sample size for discrete D-vine model assume true model is
#   discretized MVN for ordinal response.
# min KL vector from dvineKLfn via nlm
# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ucuts = (ncateg+1)xd matrix of cutpts in U(0,1) scale for ordinal responses
# pr = vector from pmfmordprobit()
# iprint = print flag for intermediate results
# Output:
#  KLdiv = KL divergence of mult probit and D-vine
#  KLss = KL sample size
#  vinepr = probability vector from D-vine approximation
dvineKLss=function(parvec,ucuts,pr,iprint=F)
{ d=ncol(ucuts)
  ncateg=nrow(ucuts)-1
  D=Dvinearray(d)
  out=varray2M(D)
  M=out$mxarray
  pcopnames=rep("pbvncop",d-1)
  nn=ncateg^d
  prv=rep(0,nn)
  u1vec=rep(0,d)
  u2vec=rep(0,d)
  parmat=matrix(0,d,d)
  d2=(d*(d-1))/2
  ii=0
  for(ell in 1:(d-1))
  { for(j in (ell+1):d)
    { ii=ii+1; parmat[ell,j]=parvec[ii] }
  }
  kl=0; kl2=0
  for(i in 1:nn)
  { jj=d2v(d,ncateg,i-1)
    for(j in 1:d) 
    { u1vec[j]=ucuts[jj[j],j]; u2vec[j]=ucuts[jj[j]+1,j] }
    prv[i]=rvinepmf.discrete(parmat,u1vec,u2vec,D,M,pcopnames,iprint=F)
    if(prv[i]<=0) prv[i]=1.e-10
    tem=log(pr[i]/prv[i])
    if(iprint) print(c(i,prv[i]))
    kl=kl+pr[i]*tem
    kl2=kl2+pr[i]*tem*tem
  }
  if(iprint) print(sum(prv))
  sigma2=kl2-kl*kl
  zcv=qnorm(.95)
  klss=(zcv/kl)^2*sigma2
  list(KLdiv=kl,KLss=klss,vinepr=prv)
}

# wrapper function for dvineKLfn and dvineKLss
# ucuts = (ncateg-1)xd matrix of cutpoints on U(0,1) scale, without boundaries
#    d = #ordinal variables (3<=d<=5)
#    ncateg = #ordinal categories
# rmat = latent correlation matrix (AR-Toeplitz structure preferable)
# iprint = print flag for intermediate results and diagnostic info
# prlevel = printlevel for nlm() optimization
# mxiter = max iterations for nlm()
# ifixseed = flag for fixed seed, 
#   if T: seed is 12345 before each call to pmvnorm with algorithm=GenzBretz
# Output:
#  mordprobitpr = probability vector from multivariate ordinal probit
#  dvineparam = parameter of D-vine approximation
#  KLdiv = KL divergence of mult probit and D-vine
#  KLss = KL sample size
#  vinepr = probability vector from D-vine approximation
ARprobitvsDvine=function(ucuts,rmat,iprint=F,prlevel=1,mxiter=50,ifixseed=F)
{ d=ncol(ucuts)
  if(d>=6) return(0)
  ncateg=nrow(ucuts)-1
  zcuts=qnorm(ucuts)
  zcuts[1,]=-6; zcuts[ncateg+1]=6
  if(iprint)
  { cat("\ncut points in U(0,1) scale\n")
    print(ucuts)
    cat("rmat\n"); print(rmat)
  }
  pr=pmfmordprobit(zcuts,rmat,iprint=iprint,ifixseed=ifixseed)
  if(d>3) 
  { D=Dvinearray(d)
    pcmat=cor2pcor.rvine(rmat,D)$pctree 
    stpar=NULL
    for(ell in 1:(d-1)) stpar=c(stpar,pcmat[ell,(ell+1):d])
  }
  else 
  { stpar=rep(0,3); stpar[1:2]=(rmat[1,2]+rmat[1,3]+rmat[2,3])/3 }
  if(iprint) { cat("start for nlm:\n"); print(stpar) }
  out=nlm(dvineKLfn,stpar,hessian=T,iterlim=mxiter,print.level=prlevel,
    ucuts=ucuts,pr=pr)
  klobj=dvineKLss(out$estimate,ucuts,pr,iprint=iprint)
  if(iprint) cat("KLss=",klobj$KLss,"\n")
  list(mordprobitpr=pr, dvineparam=out$estimate,
    KLdiv=klobj$KLdiv, KLss=klobj$KLss, vinepr=klobj$vinepr)
}

#============================================================

# KL divergence for discrete R-vine model
# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ucuts = (ncateg+1)xd matrix of cutpts in U(0,1) scale for ordinal responses
# A = dxd vine array with 1:d on diagonal
# pr = vector from pmfmordprobit()
# Output: KL divergence
rvineKLfn=function(parvec,ucuts,A,pr)
{ if(any(parvec>=1) | any(parvec<=-1)) return(1.e10)
  #if(any(parvec>=0.99999) | any(parvec<=-0.99999)) return(1.e10)
  d=ncol(ucuts)
  ncateg=nrow(ucuts)-1
  out=varray2M(A)
  M=out$mxarray
  pcopnames=rep("pbvncop",d-1)
  nn=ncateg^d
  prv=rep(0,nn)
  u1vec=rep(0,d)
  u2vec=rep(0,d)
  parmat=matrix(0,d,d)
  d2=(d*(d-1))/2
  ii=0
  for(ell in 1:(d-1))
  { for(j in (ell+1):d)
    { ii=ii+1; parmat[ell,j]=parvec[ii] }
  }
  kl=0
  for(i in 1:nn)
  { jj=d2v(d,ncateg,i-1)
    for(j in 1:d) 
    { u1vec[j]=ucuts[jj[j],j]; u2vec[j]=ucuts[jj[j]+1,j] }
    prv[i]=rvinepmf.discrete(parmat,u1vec,u2vec,A,M,pcopnames,iprint=F)
    if(prv[i]<=0) prv[i]=1.e-10
    if(is.na(prv[i]) | is.infinite(prv[i])) return(1.e10)
    kl=kl+pr[i]*log(pr[i]/prv[i])
  }
  kl
}

# KL sample size for discrete R-vine model assume true model is
#   discretized MVN for ordinal response.
# min KL vector from dvineKLfn via nlm
# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ucuts = (ncateg+1)xd matrix of cutpts in U(0,1) scale for ordinal responses
# A = dxd vine array with 1:d on diagonal
# pr = vector from pmfmordprobit()
# iprint = print flag for intermediate results
# Output:
#  KLdiv = KL divergence of mult probit and R-vine
#  KLss = KL sample size
#  vinepr = probability vector from R-vine approximation
rvineKLss=function(parvec,ucuts,A,pr,iprint=F)
{ d=ncol(ucuts)
  ncateg=nrow(ucuts)-1
  out=varray2M(A)
  M=out$mxarray
  pcopnames=rep("pbvncop",d-1)
  nn=ncateg^d
  prv=rep(0,nn)
  u1vec=rep(0,d)
  u2vec=rep(0,d)
  parmat=matrix(0,d,d)
  d2=(d*(d-1))/2
  ii=0
  for(ell in 1:(d-1))
  { for(j in (ell+1):d)
    { ii=ii+1; parmat[ell,j]=parvec[ii] }
  }
  kl=0; kl2=0
  for(i in 1:nn)
  { jj=d2v(d,ncateg,i-1)
    for(j in 1:d) 
    { u1vec[j]=ucuts[jj[j],j]; u2vec[j]=ucuts[jj[j]+1,j] }
    prv[i]=rvinepmf.discrete(parmat,u1vec,u2vec,A,M,pcopnames,iprint=F)
    if(prv[i]<=0) prv[i]=1.e-10
    tem=log(pr[i]/prv[i])
    if(iprint) print(c(i,prv[i]))
    kl=kl+pr[i]*tem
    kl2=kl2+pr[i]*tem*tem
  }
  if(iprint) print(sum(prv))
  sigma2=kl2-kl*kl
  zcv=qnorm(.95)
  klss=(zcv/kl)^2*sigma2
  list(KLdiv=kl,KLss=klss,vinepr=prv)
}

# wrapper function for rvineKLfn and rvineKLss
# This function also requires library(combinat)
# ucuts = (ncateg+1)xd matrix of cutpoints on U(0,1) scale, with boundaries
#    d = #ordinal variables (3<=d<=5)
#    ncateg = #ordinal categories
# rmat = latent correlation matrix 
# A = dxd vine array with 1:d on diagonal
# iprint = print flag for intermediate results and diagnostic info
# prlevel = printlevel for nlm() optimization
# mxiter = max iterations for nlm()
# ifixseed = flag for fixed seed, 
#   if T: seed is 12345 before each call to pmvnorm with algorithm=GenzBretz
# Output: 
#  mordprobitprmat = prob vectors from mult ordinal probit (ncateg^d x d!/2)
#  parmat = parameter of R-vine approximation for each perm (C(d,2) x d!/2)
#  vKLdiv = vector of KL divergences of mult probit and R-vine
#  vKLss = vector of KL sample size
#  rvineprmat = probability vectors from R-vine approximation (ncateg^d x d!/2) 
mprobitvsRvine=function(ucuts,rmat,A,iprint=F,prlevel=1,mxiter=50,ifixseed=F)
{ d=ncol(ucuts)
  ncateg=nrow(ucuts)-1
  zcuts=qnorm(ucuts)
  zcuts[1,]=-6; zcuts[ncateg+1]=6
  if(iprint)
  { cat("\ncut points in U(0,1) scale\n")
    print(ucuts[-c(1,ncateg+1),])
    cat("rmat\n"); print(rmat)
  }
  pr=pmfmordprobit(zcuts,rmat,iprint=iprint,ifixseed=ifixseed)
  nn=length(pr)
  # set up d!/2 permutations
  permd=permn(d)
  permd=t(array(unlist(permd), dim = c(d,prod(1:d))))
  ii=(permd[,d-1]<permd[,d])
  permd=permd[ii,]
  nperm=nrow(permd)
  if(iprint) cat("nperm=",nperm,"\n")
  parmat=matrix(0,d*(d-1)/2,nperm)
  vKLdiv=rep(0,nperm)
  vKLss=rep(0,nperm)
  mprobmat=matrix(0,ncateg^d ,nperm) 
  vineprmat=matrix(0,ncateg^d ,nperm) 
  # permute all
  for(ii in 1:nperm)
  { perm=permd[ii,]
    if(iprint) cat("\nperm=", perm,"\n")
    ucutsp=ucuts[,perm]
    rmatp=rmat[perm,perm]
    # need to permute pr also
    prp=pr
    for(i in 1:nn)
    { jj=d2v(d,ncateg,i-1)
      jj2=jj[perm]
      i2=v2d(jj2-1,ncateg)
      prp[i2+1]=pr[i]
    }
    pcmat=cor2pcor.rvine(rmatp,A)$pctree 
    stpar=NULL
    for(ell in 1:(d-1)) stpar=c(stpar,pcmat[ell,(ell+1):d])
    if(iprint) { cat("start for nlm:\n"); print(stpar) }
    out=nlm(rvineKLfn,stpar,hessian=T,iterlim=mxiter,print.level=prlevel,
       ucuts=ucutsp,pr=prp,A=A)
    klobj=rvineKLss(out$estimate,ucutsp,A,prp,iprint=iprint)
    if(iprint) cat("KLss=",klobj$KLss,"\n")
    parmat[,ii]=out$estimate
    mprobmat[,ii]=prp
    vineprmat[,ii]=klobj$vinepr
    vKLdiv[ii]=klobj$KLdiv
    vKLss[ii]=klobj$KLss
  }
  list(mordprobitprmat=mprobmat, parmat=parmat,
    vKLdiv=vKLdiv, vKLss=vKLss, rvineprmat=vineprmat)
}

#============================================================

# f90 interface for KL divergence 
# parvec = parameter vector of partial correlations with length d*(d-1)/2
# ucuts = (ncateg+1)xd matrix of cutpts in U(0,1) scale for ordinal responses
# A = dxd vine array with 1:d on diagonal
# M = maximal array obtained from A
#  out=varray2M(A); M=out$mxarray (before call)
# pr = vector from pmfmordprobit()
# Output: KL divergence
f90rvineKL=function(parvec,ucuts,A,M,pr)
{ if(any(parvec>=1) | any(parvec<=-1)) return(1.e10)
  d=ncol(ucuts); ncateg=nrow(ucuts)-1
  out= .Fortran("rvineklfn",
    as.integer(d), as.integer(ncateg), as.double(parvec), 
    as.double(ucuts), as.integer(A), as.integer(M), 
    as.double(pr),kl=as.double(0) )
  kl=out$kl 
  kl
}


# wrapper function for rvineKLfn and rvineKLss with call to f90
# This function also requires library(combinat)
# ucuts = (ncateg+1)xd matrix of cutpoints on U(0,1) scale, with boundaries
#    d = #ordinal variables (3<=d<=5)
#    ncateg = #ordinal categories
# rmat = latent correlation matrix 
# A = dxd vine array with 1:d on diagonal
# iprint = print flag for intermediate results and diagnostic info
# prlevel = printlevel for nlm() optimization
# mxiter = max iterations for nlm()
# ifixseed = flag for fixed seed, 
#   if T: seed is 12345 before each call to pmvnorm with algorithm=GenzBretz
# Output: 
#  mordprobitprmat = prob vectors from mult ordinal probit (ncateg^d x d!/2)
#  parmat = parameter of R-vine approximation for each perm (C(d,2) x d!/2)
#  vKLdiv = vector of KL divergences of mult probit and R-vine
#  vKLss = vector of KL sample size
#  rvineprmat = probability vectors from R-vine approximation (ncateg^d x d!/2) 
f90mprobitvsRvine=function(ucuts,rmat,A,iprint=F,prlevel=1,mxiter=50,ifixseed=F)
{ d=ncol(ucuts)
  ncateg=nrow(ucuts)-1
  zcuts=qnorm(ucuts)
  zcuts[1,]=-6; zcuts[ncateg+1]=6
  if(iprint)
  { cat("\ncut points in U(0,1) scale\n")
    print(ucuts[-c(1,ncateg+1),])
    cat("rmat\n"); print(rmat)
  }
  pr=pmfmordprobit(zcuts,rmat,iprint=iprint,ifixseed=ifixseed)
  nn=length(pr)
  # set up d!/2 permutations
  permd=permn(d)
  permd=t(array(unlist(permd), dim = c(d,prod(1:d))))
  ii=(permd[,d-1]<permd[,d])
  permd=permd[ii,]
  nperm=nrow(permd)
  if(iprint) cat("nperm=",nperm,"\n")
  parmat=matrix(0,d*(d-1)/2,nperm)
  vKLdiv=rep(0,nperm)
  vKLss=rep(0,nperm)
  mprobmat=matrix(0,ncateg^d ,nperm) 
  vineprmat=matrix(0,ncateg^d ,nperm) 
  vout=varray2M(A); M=vout$mxarray
  # permute all
  for(ii in 1:nperm)
  { perm=permd[ii,]
    if(iprint) cat("\nperm=", perm,"\n")
    ucutsp=ucuts[,perm]
    rmatp=rmat[perm,perm]
    # need to permute pr also
    prp=pr
    for(i in 1:nn)
    { jj=d2v(d,ncateg,i-1)
      jj2=jj[perm]
      i2=v2d(jj2-1,ncateg)
      prp[i2+1]=pr[i]
    }
    pcmat=cor2pcor.rvine(rmatp,A)$pctree 
    stpar=NULL
    for(ell in 1:(d-1)) stpar=c(stpar,pcmat[ell,(ell+1):d])
    if(iprint) { cat("start for nlm:\n"); print(stpar) }
    out=nlm(f90rvineKL,stpar,hessian=T,iterlim=mxiter,print.level=prlevel,
       ucuts=ucutsp,pr=prp,A=A,M=M)
    klobj=rvineKLss(out$estimate,ucutsp,A,prp,iprint=iprint)
    if(iprint) cat("KLss=",klobj$KLss,"\n")
    parmat[,ii]=out$estimate
    mprobmat[,ii]=prp
    vineprmat[,ii]=klobj$vinepr
    vKLdiv[ii]=klobj$KLdiv
    vKLss[ii]=klobj$KLss
  }
  list(mordprobitprmat=mprobmat, parmat=parmat,
    vKLdiv=vKLdiv, vKLss=vKLss, rvineprmat=vineprmat)
}

#============================================================
