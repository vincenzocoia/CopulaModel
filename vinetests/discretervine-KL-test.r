# KL divergence of discrete R-vine with multivariate ordinal probit

library(CopulaModel)
library(combinat)
library(mvtnorm)

# wrapper function 1
# perm = permutation of 1:d
# ucuts = (ncateg+1) x d matrix of ordinal cutpoints, first row approx 0
#                   last row approx 1
# rmat = dxd correlation matrix
# A = dxd vine array with 1:d on diagonal
# pr = probability vector from multivariate ordinal probit
# iprint = print flag for intermediate results
# prlevel = print.level for nlm()
# Output: KL sample size is printed within the function
wrap1r=function(perm,ucuts,rmat,A,pr,iprint=F,prlevel=1)
{ d=ncol(ucuts)
  tem=varray2M(A); M=tem$mxarray
  ncateg=nrow(ucuts)-1
  # permute all
  ucutsp=ucuts[,perm]
  rmatp=rmat[perm,perm]
  nn=length(pr)
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
  if(iprint) print(stpar)
  #out0=nlm(rvineKLfn,stpar,hessian=F,iterlim=50,print.level=prlevel,ucuts=ucutsp,pr=prp,A=A)
  out=nlm(f90rvineKL,stpar,hessian=F,iterlim=50,print.level=prlevel,ucuts=ucutsp,pr=prp,A=A,M=M)
  klobj=rvineKLss(out$estimate,ucutsp,A,prp,iprint=iprint)
  cat("KLss=",klobj$KLss,"\n")
  invisible(0)
}

# wrapper function 2
# d = dimension
# ncateg = number of ordinal cutpoints, between 3 to 9
# rmat = dxd correlation matrix
# A = dxd vine array with 1:d on diagonal
# iprint = print flag for intermediate results
# prlevel = print.level for nlm()
# Output: KL sample size is printed within the function
#   enumeration through perms of vine array A
#   calling routine can enumerate through different A when d>=5
wrap2r=function(d,ncateg,rmat,A,iprint=F,prlevel=1)
{ # create random cutpoints, each in (1:9)*0.1
  tem=varray2M(A); M=tem$mxarray
  ucuts=matrix(0,ncateg-1,d)
  for(j in 1:d)
  { ucuts[1,j]=floor(runif(1,0,10-ncateg))+1
    for(k in 2:(ncateg-1))
    { ucuts[k,j]=floor(runif(1,ucuts[k-1,j],9-ncateg+k))+1 }
  }
  ucuts=ucuts/10
  cat("\ncut points in U(0,1) scale\n")
  print(ucuts)
  cat("rmat\n")
  print(rmat)
  zcuts=qnorm(ucuts)
  zcuts=rbind(rep(-6,d),zcuts,rep(6,d))
  ucuts=rbind(rep(0.00001,d),ucuts,rep(.99999,d)) # for vine
  pr=pmfmordprobit(zcuts,rmat5,iprint=F)
  permd=permn(d)
  permd=t(array(unlist(permd), dim = c(d,prod(1:d))))
  ii=(permd[,d-1]<permd[,d])
  permd=permd[ii,]
  nn=length(pr)
  for(ii in 1:nrow(permd))
  { perm=permd[ii,]
    cat("\nii= ", ii, " perm=", perm,"\n")
    rmatp=rmat[perm,perm]
    ucutsp=ucuts[,perm]
    pcmat=cor2pcor.rvine(rmatp,A)$pctree 
    stpar=NULL
    for(ell in 1:(d-1)) stpar=c(stpar,pcmat[ell,(ell+1):d])
    # need to permute pr also
    if(iprint) print(stpar)
    prp=pr
    for(i in 1:nn)
    { jj=d2v(d,ncateg,i-1)
      jj2=jj[perm]
      i2=v2d(jj2-1,ncateg)
      prp[i2+1]=pr[i]
    }
    #out=nlm(rvineKLfn,stpar,hessian=F,iterlim=50,print.level=prlevel,
    #   ucuts=ucutsp,pr=prp,A=A)
    out=nlm(f90rvineKL,stpar,hessian=F,iterlim=50,print.level=prlevel,
       ucuts=ucutsp,pr=prp,A=A,M=M)
    if(ii>4) prlevel=0
    klobj=rvineKLss(out$estimate,ucutsp,A,prp,iprint=iprint)
    cat("KLss=",klobj$KLss,"\n")
  }
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  invisible(0)
}


set.seed(12345)
d=5
ucuts5=matrix(c(.4,.5,.4,.3,.2, .7,.8,.6,.6,.6),2,5,byrow=T)
rh=0.5
rmat5=toeplitz(rh^(0:4))
zcuts=qnorm(ucuts5)
zcuts=rbind(rep(-6,d),zcuts,rep(6,d))
ucuts5=rbind(rep(0.00001,d),ucuts5,rep(.99999,d)) # for vine
pr=pmfmordprobit(zcuts,rmat5,iprint=F)
#print(pr)
cat("sum of multivariate probit probabilities=", sum(pr),"\n")
rh1=.6; rh2=.5
acf10=ar2acf(rh1,rh2,10)
ar2mat5=toeplitz(acf10[1:5])
D5=Dvinearray(5,iNO=T)
perm5=permn(5)
perm5=t(array(unlist(perm5), dim = c(5,120)))
ii=(perm5[,4]<perm5[,5])
perm5b=perm5[ii,]

cat("\nd=5, enumerate through perms of Dvine\n")
#for(ii in 1:nrow(perm5b))
for(ii in 1:5)
{ perm=perm5b[ii,]
  cat("\nperm=",perm,"\n")
  wrap1r(perm,ucuts5,rmat5,D5,pr,iprint=F)
}
cat("\n============================================================\n")

# allow for some randomness in ucuts with wrap2r
set.seed(123)
for(isim in 1:2)
{ cat("\nrandom ucuts with AR1, isim=", isim,"\n")
  for(bnum in 0:7) 
  { cat("bnum=", bnum,"\n")
    A5=vnum2array(5,bnum) 
    wrap2r(5,3,rmat5,A5,iprint=F) 
    cat("\n============================================================\n")
  }
}
for(isim in 1:2)
{ cat("\nrandom ucuts with AR2, isim=", isim,"\n")
  for(bnum in 0:7) 
  { cat("bnum=", bnum,"\n")
    A5=vnum2array(5,bnum) 
    wrap2r(5,3,ar2mat5,A5,iprint=F) 
    cat("\n============================================================\n")
  }
}

