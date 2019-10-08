# Functions for multivariate ordinal response including
# polychoric correlation for various input formats

# decimal to vector
# d = dimension, 
# ncat = #categories, 
# ii = non-negative integer
# izero = T if categories start at 0
# izero = F if categories start at ,
# Output: vector of size d, each element in 1..ncat or 0..(ncat-1)
d2v=function(d,ncat,ii,izero=F)
{ tt=ii
  jj=rep(0,d)
  for(i in seq(d,1,-1))
  { jj[i]=tt%%ncat; tt=floor(tt/ncat); }
  if(!izero) jj=jj+1
  jj
}

# vector to decimal, 
# jj = vector with elements in 0:(ncat-1), d=length of vector
# ncat = number of ordinal categories
# Output:  decimal representation of vector when all possible
#    d-vector of 0:(ncat-1) are put in lexicographical order
v2d=function(jj,ncat)
{ d=length(jj)
  pw=1; s=0
  for(i in seq(d,1,-1))
  { s=s+pw*jj[i]; pw=pw*ncat; }
  s
}

# decimal to binary vector
# d = dimension
# ii = integer ii  in 0 to 2^d-1
# Output: binary d-vector jj corresponding to integer ii  
d2b=function(d,ii)
{ s=ii; jj=rep(0,d)
  for(i in d:1)
  { jj[i]=s%%2; s=floor(s/2) }
  jj
}


# This function should only be used if ncat^d is not too large.
# odat = nxd matrix of ordinal data with categories 0:(ncat-1) or 1:ncat;
#   in the latter case, data are converted to 0:(ncat-1).
# ncat = number of ordinal categories
# Output: data set with frequencies of distinct observed vectors.
ordinal2fr=function(odat,ncat)
{ n=nrow(odat)
  d=ncol(odat)
  if(min(odat)==1) odat=odat-1 # relabel 0:(ncat-1)
  indx=rep(0,n)
  for(i in 1:n) { indx[i]=v2d(odat[i,],ncat) }
  tbl=table(indx)
  aa=names(tbl)
  aa=as.integer(aa)
  fr=as.integer(tbl)
  nn=length(fr)
  odatfr=matrix(0,nn,d+1)
  odatfr[,d+1]=fr
  for(i in 1:nn) { odatfr[i,1:d]=d2v(d,ncat,aa[i],izero=T) }
  odatfr
} 

# Calculate the matrix with the cutpoints in the uniform(0,1) scale
# odat = data matrix of ordinal response where a row represents an
#  individual's response and a column represents an item or variable.
# Number of ordinal categories for each variable is ncat
#  and the ordinal categories are coded as 0,...,(ncat-1). or 1:ncat
#  in a contiguous manner; in the latter case, data are converted to 0:(ncat-1)
# Currently supported: items have the same number of categories.
# Output:
# A matrix with the cutpoints (without the boundary cutpoints) in the
# uniform scale where the number of rows corresponds to an ordinal category
# and each column represents an item.
unifcuts=function(odat)
{ n=nrow(odat)
  d=ncol(odat)
  ncat=length(table(odat))
  if(min(odat)==1) odat=odat-1 # relabel 0:(ncat-1)
  a=matrix(NA,ncat-1,d)
  for(j in 1:d) 
  { fr=table(c(odat[,j],0:(ncat-1)))  # augment data in case categ is missing
    pmf=(fr-1)/n  # 0 is possible for some categories
    cdf=cumsum(pmf)
    a[,j]=cdf[-ncat]
  }
  a
}

# Version for multivariate ordinal with no predictors/covariates
# odatfr = ordinal data matrix, nrec x (d+1), d=#items, nrec=num of records
#    column d+1 has frequencies of the d-vectors
# ncat = nc=number of ordinal categories
#    categories are in 0,1,...,ncat-1
# zcuts = (ncat+1)xd matrix in N(0,1) scale, obtained via unifcuts and qnorm 
# iprint = flag to print observed vs expected table assuming bivariate Gaussian
# prlevel = print.level in nlm()
# Outputs : polychoric correlation matrix based on two-stage estimate and 
#  an indicator if the 2-stage correlation matrix estimate is positive definite
polychoric=function(odatfr,zcuts,iprint=F,prlevel=0)
{ d=ncol(odatfr)-1
  fr=odatfr[,d+1]
  dat=odatfr[,1:d]
  nc=nrow(zcuts)-1
  sampsize=sum(fr)
  kk=0
  dd=d*(d-1)/2
  nllkvec=rep(0,dd)
  rhovec=rep(0,dd)
  summ=array(0,c(nc,2*nc,dd))
  for(j2 in 2:d)
  { for(j1 in 1:(j2-1))
    { if(iprint) cat("\nitems:", j1, j2, "\n")
      kk=kk+1
      y1=dat[,j1]; y2=dat[,j2]
      # addition in case of zero counts for some categories
      y1=c(y1,0:(nc-1))
      y2=c(y2,0:(nc-1))
      fradd=c(fr,rep(1,nc))
      bcount=tapply(fradd,list(y1,y2),sum)
      bcount[is.na(bcount)]=0
      #print(bcount)
      diag(bcount)=diag(bcount)-1
      if(iprint) print(bcount)
      out=nlm(bprobitnllk,p=.5,hessian=T,iterlim=50,print.level=prlevel,
        zcuts=zcuts,bfr=bcount,jj1=j1,jj2=j2)
      if(iprint) 
      { rhose=1/sqrt(out$hess)
        cat("nllk min, estimate and SE:",out$min, out$estimate, rhose,"\n") 
      }
      nllkvec[kk]=out$min
      rhovec[kk]=out$estimate
      # model-based expected table 
      nc1=nc+1
      z1=matrix(zcuts[,j1],nc1,nc1)
      z2=matrix(zcuts[,j2],nc1,nc1,byrow=T)
      cdf2= pbnorm(z1,z2,out$estimate)
      pmf=apply(cdf2,2,diff)
      pmf2=apply(t(pmf),2,diff)
      pmf2=t(pmf2)
      pmf2[pmf2<=0]=1.e-15
      ecount=pmf2*sampsize
      if(iprint) { cat("observed and expected\n"); print(cbind(bcount,ecount))}
      summ[,,kk]=cbind(bcount,ecount)
    }
  }
  #cat("\npolychloric correlation matrix\n")
  rmat=corvec2mat(rhovec)
  iposdef=isposdef(rmat)
  list(polych=rmat,zcuts=zcuts,iposdef=iposdef)
}

# This function is used by polychoric().
# pbnorm is used to get the bivariate cdf, 
#    then row and column differencing are used to get the bivariate pmf
# rho = correlation parameter 
# zcuts = (ncat+1)xd matrix in N(0,1) scale, 
#    ncat=nc=#categories for each variable 
# bfr = vector of bivariate frequencies
# jj1, jj2 = indices of 2 variables
# Output: negative log-likelihood
bprobitnllk=function(rho,zcuts,bfr,jj1,jj2)
{ nc=nrow(zcuts)-1
  nlk=0.; 
  if(rho<=-1. | rho>=1.) { return(1.e10) }
  nc1=nc+1
  z1=matrix(zcuts[,jj1],nc1,nc1)
  z2=matrix(zcuts[,jj2],nc1,nc1,byrow=T)
  cdf2= pbnorm(z1,z2,rho)
  pmf=apply(cdf2,2,diff)
  pmf2=apply(t(pmf),2,diff)
  pmf2=t(pmf2)
  pmf2[pmf2<=0]=1.e-15
  nlk=-sum(bfr*log(pmf2))
  nlk
}

#============================================================

# This function is used by polychoric.wPred()
# assume ordprobit.univar and mord2uu have been used to get zzdat.
# zzdat = nxd matrix with corners of rectangle for each vector observation 
# rho = latent correlation parameter
# jj1, jj2 = indices of 2 variables
# Output: negative log-likelihood
bprobitwPrednllk=function(rho,zzdat,jj1,jj2)
{ if(rho<=-1. | rho>=1.) { return(1.e10) }
  n=nrow(zzdat); d=ncol(zzdat)/2
  nlk=0; 
  for(i in 1:n)
  { zcdf1=zzdat[i,d+jj1]; zcdf2=zzdat[i,d+jj2];
    zcdf1m=zzdat[i,jj1]; zcdf2m=zzdat[i,jj2];
    bpmf= pbnorm(zcdf1,zcdf2,rho) - pbnorm(zcdf1m,zcdf2,rho) -
           pbnorm(zcdf1,zcdf2m,rho) + pbnorm(zcdf1m,zcdf2m,rho)
    if(bpmf<=0.) bpmf=1.e-10
    nlk=nlk-log(bpmf)
  }
  nlk
}

# polychoric correlation matrix estimate for ordinal response with predictors, 
# assume ordprobit.univar and mord2uu have been used to get zzdat
# zzdat = nxd matrix with corners of rectangle for each vector observation 
# iprint = flag for printing intermediate results
# prlevel = print.level for nlm
# Outputs : polychoric correlation matrix  based on two-stage estimate and 
#  an indicator if the 2-stage correlation matrix estimate is positive definite
polychoric.wPred=function(zzdat,iprint=F,prlevel=0)
{ n=nrow(zzdat); d=ncol(zzdat)/2
  kk=0
  dd=d*(d-1)/2
  nllkvec=rep(0,dd)
  rhovec=rep(0,dd)
  for(j1 in 2:d)
  { for(j2 in 1:(j1-1))
    { if(iprint) cat("\nvariables:", j1, j2, "\n")
      kk=kk+1
      out=nlm(bprobitwPrednllk,p=.5,hessian=F,iterlim=50,print.level=prlevel,
        zzdat=zzdat,jj1=j1,jj2=j2)
      if(iprint) { cat("nllk min and estimate :",out$min, out$estimate,"\n") }
      nllkvec[kk]=out$min
      rhovec[kk]=out$estimate
    }
  }
  rmat=corvec2mat(rhovec)
  iposdef=isposdef(rmat)
  list(polych=rmat,iposdef=iposdef)
}

# polychoric correlation from a bivariate table/array
# bivtab = bivariate ordinal, which will be labeled as 0:ncat[j] 
#                         or 1:ncat[j] for j=1,2
# iprint = flag to print observed vs expected table assuming bivariate Gaussian
# prlevel = print.level for nlm()
# Output: 
#   polychoric correlation and its SE (with empirical univariate cutpoints)
polychoric.bivtab=function(bivtab,iprint=F,prlevel=0)
{ nc1=nrow(bivtab)
  nc2=ncol(bivtab)
  n=sum(bivtab)
  fr=apply(bivtab,1,sum)
  pmf=fr/n  
  cdf=cumsum(pmf)
  cut1=c(0,cdf)
  fr=apply(bivtab,2,sum)
  pmf=fr/n  
  cdf=cumsum(pmf)
  cut2=c(0,cdf)
  zcut1=qnorm(cut1); zcut1[1]=-6; zcut1[nc1+1]=6
  zcut2=qnorm(cut2); zcut2[1]=-6; zcut2[nc2+1]=6
  bprobnllk= function (rho) 
  { nlk=0
    if (rho<= -1 | rho>=1) { return(1e+10) }
    z1=matrix(zcut1,nc1+1,nc2+1)
    z2=matrix(zcut2,nc1+1,nc2+1, byrow=T)
    cdf2=pbnorm(z1,z2,rho)
    pmf=apply(cdf2,2,diff)
    pmf2=apply(t(pmf),2,diff)
    pmf2=t(pmf2)
    pmf2[pmf2<=0]=1e-15
    nlk=-sum(bivtab*log(pmf2))
    nlk
  }
  out=nlm(bprobnllk,p=.5,hessian=T,iterlim=50,print.level=prlevel)
  rho=out$estimate
  SE=sqrt(1/out$hessian)
  # model-based expected table 
  z1=matrix(zcut1,nc1+1,nc2+1)
  z2=matrix(zcut2,nc1+1,nc2+1, byrow = T)
  cdf2= pbnorm(z1,z2,out$estimate)
  pmf=apply(cdf2,2,diff)
  pmf2=apply(t(pmf),2,diff)
  pmf2=t(pmf2)
  pmf2[pmf2<=0]=1.e-15
  ecount=pmf2*n
  if(iprint) { cat("observed and expected\n"); print(cbind(bivtab,ecount)) }
  c(rho,SE)
}

# polychoric correlation from data set of ordinal variables
# This function could fail(?) if some ordinal category has zero counts
#   for one or more of the ordinal variables
#   (in this case, use polychoric() with some preprocessing)
# odat = ordinal data set with d columns
# iprint = flag to print observed vs expected table assuming bivariate Gaussian
#  and SEs for the polychoric correlations for each bivariate margin
# prlevel = print.level for nlm()
# Outputs : 
#  polychoric correlation matrix and
#  an indicator if the 2-stage correlation matrix estimate is positive definite
polychoric0=function(odat,iprint=F,prlevel=0)
{ d=ncol(odat)
  rmat=matrix(1,d,d)
  for(j2 in 2:d)
  { for(j1 in 1:(j2-1))
    { atab=table(odat[,j1],odat[,j2])
      if(iprint) cat("\n",j1,j2,"\n")
      out=polychoric.bivtab(atab,iprint=iprint,prlevel=prlevel)
      if(iprint) cat("polychoric and SE:",  out,"\n")
      rmat[j1,j2]=out[1]; rmat[j2,j1]=out[1]
    }
  }
  iposdef=isposdef(rmat)
  list(polych=rmat, iposdef=iposdef)
}

