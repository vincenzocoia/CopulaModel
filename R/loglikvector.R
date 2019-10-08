# Functions for Vuong's procedure, and vectors of log-likelihood 
# for each observation (evaluated at MLE or IFM estimate) 

# Vuong's procedure for 2 competing models
# llkv1,llkv2 = two loglikelihood vectors of the same length as data set
#   computed from a common data set with 2 different models
# Output: interval from Vuong's procedure
#  negative interval means that model 1 is better
#  positive interval means that model 2 is better
#  interval that includes 0 implies models not significantly different
vuongllkr=function(llkv1,llkv2)
{ n=length(llkv1)
  ddif=llkv2-llkv1
  Del12=mean(ddif)
  sig12=sd(ddif)
  moe=1.96*sig12/sqrt(n)
  print(c(Del12,sig12,moe))
  c(Del12-moe,Del12+moe)
}

# Vuong's procedure for 2 competing models
# version 2 with Schwarz/BIC correction of dim(cpar)*log(n)/(2*n)
# llkv1,llkv2 = two loglikelihood vectors of the same length as data set
#   computed from a common data set with 2 different models
# dim1, dim2 = parameter vector dimensions for models 1 and 2
# Output: interval from Vuong's procedure with Schwarz correction
#  negative interval means that model 1 is better
#  positive interval means that model 2 is better
#  interval that includes 0 implies models not significantly different
vuong2llkr=function(llkv1,llkv2,dim1,dim2)
{ n=length(llkv1)
  ddif=llkv2-llkv1
  Del12=mean(ddif)
  sig12=sd(ddif)
  moe=1.96*sig12/sqrt(n)
  cat("Vuong Delta, sigma, moe: ", Del12,sig12,moe,"\n")
  Del12=Del12-(dim2-dim1)*log(n)/2/n # Schwarz correction
  cat("Vuong/Schwarz Delta, sigma, moe: ", Del12,sig12,moe,"\n")
  c(Del12-moe,Del12+moe)
}


# Clarke's procedure for 2 competing models (not as good)
# llkv1,llkv2 = two loglikelihood vectors of the same length as data set
#   computed from a common data set with 2 different models
# Output: interval from Clarke's procedure
#  interval below 0.5 means that model 1 is better
#  interval above 0.5 means that model 2 is better
#clarkellkr= function(llkv1,llkv2,iprint=FALSE)
#{ n=length(llkv1)
#  ddif=llkv2-llkv1
#  phat=mean(ddif>=0)
#  sig12=sqrt(phat*(1-phat))
#  moe=1.96*sig12/sqrt(n)
#  if(iprint) cat("phat,sig12,moe", phat,sig12,moe,"\n")
#  c(phat-moe,phat+moe)
#}

#============================================================
# For repeated measures discrete data : assume conversion to u1vec, u2vec
#  the corners of rectangles

# log-likelihood vectors

# param = parameter for mrectpr()
# uudat = nx(2d) matrix with corners of rectangle on U(0,1) scale
#  on each observation, based on sample size n,
#    d = clsize = cluster size, 
# mrectpr = function for multivariate rectangle probability
# Output: llkv vector for discrete with pr=rectangle probability
#      length nrow(uudat)
mdiscretellkv=function(param,uudat,mrectpr)
{ d=ncol(uudat)/2
  n=nrow(uudat)
  llkv=rep(0,n)
  jj1=1:d; jj2=(d+1):(2*d)
  for(i in 1:n)
  { u1=uudat[i,jj1]; u2=uudat[i,jj2]
    pr=mrectpr(u1,u2,param)
    if(is.na(pr) | pr<=0.)  pr=1.e-15
    llkv[i]=log(pr)
  }
  llkv  
}

# convert to N(0,1) scale for exchangeable discretized mvn 
# param = exchangeable correlation parameter
# zzdat = nx(2d) matrix with corners of rectangle on N(0,1) scale
#   d = #response variable or cluster size, 
# Output: llkv vector for discretrized positive exchangeable MVN,
#      length nrow(zzdat)
emvndiscretellkv=function(param,zzdat)
{ d=ncol(zzdat)/2
  rh=param
  n=nrow(zzdat)
  llkv=rep(0,n)
  jj1=1:d; jj2=(d+1):(2*d)
  for(i in 1:n)
  { z1=zzdat[i,jj1]; z2=zzdat[i,jj2]
    pr=exchmvn(z1,z2,rh)
    if(is.na(pr) | pr<=0.)  pr=1.e-15
    llkv[i]=log(pr)
  }
  llkv  
}

# param = parameter of bivariate copula cdf pcop,
#   current version: scalar parameter for each pair-copula
# uudat = nx(2d) matrix with corners of rectangle on U(0,1) scale
# A = vine array with 1:d on diagonal
# pcopnames = string vector with names of copula cdfs (length ntrunc)
# iprint = print flag for printing in rvinepmf.discrete
# Output: llkv vector with discrete rectangle probabilities, length nrow(uudat)
rvinediscretellkv=function(param,uudat,A,pcopnames,iprint=FALSE)
{ ntrunc=length(pcopnames) # assume >=1
  d=ncol(uudat)/2
  n=nrow(uudat)
  parmat=matrix(0,d,d)
  ii=0
  for(ell in 1:ntrunc)
  { parmat[ell,(ell+1):d]=param[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  out=varray2M(A) # maximal array
  M=out$mxarray
  llkv=rep(0,n)
  jj1=1:d; jj2=(d+1):(2*d)
  for(i in 1:n)
  { pr=rvinepmf.discrete(parmat,uudat[i,jj1],uudat[i,jj2],A,M,pcopnames,iprint)
    if(pr<=0) pr=1.e-15
    llkv[i]=log(pr)
  }
  llkv
}

# ir1factpmf() is in ir1fact.R
# ir2factpmf() is in ir2fact.R

#============================================================

# Continuous response, for R-vine, see rvinenllk.R for 
#  rvinellkv1.trunc() and rvinellkv.trunc

#============================================================

# Interface to f90 code for 1-factor, 2-factor, bi-factor, nested-factor
# param = parameter vector with MLE
# udat = nxd matrix of uniform scores
# strmodel = one of "1factor", "2factor", "bifactor", "nestedfactor" 
# copname = something like "frank", "bb1", "bb1frank", "bb1frk", "t", "tapprox"
# nq = number of quadrature points per dimension
# nu = degree of freedom parameter for copname= "t" or "tapprox"; 
#    it is scalar or has dimension 2 depending on 1 or more factors.
# Output: llkv = vector for log-likelihood with dimension=nrow(udat)
strfactllkv=function(param,udat,strmodel,copname,nq,grsize=0,nu=0,ipdf=1)
{ n=nrow(udat)
  d=ncol(udat)
  gl=gausslegendre(nq)
  llkv=rep(0,n)
  if(length(nu)>1) { nu1=nu[1]; nu2=nu[2] }
  else { nu1=0; nu2=0 }
  if(strmodel=="1factor") { f90cop=f90cop1nllk }
  else if(strmodel=="2factor") { f90cop=f90cop2nllk }
  else if(strmodel=="bifactor") { f90cop=f90str2nllk }
  else if(strmodel=="nestedfactor") { f90cop=f90str1nllk }
  for(i in 1:n)
  { dstr=list(data=matrix(udat[i,],1,d), copname=copname, quad=gl,repar=0,
      nu=nu,nu1=nu1,nu2=nu2, grsize=grsize, pdf=ipdf)
    llkv[i]=-f90cop(param,dstr,iprfn=F)$fnval
  }
  llkv
}

# llk vector for mvt with 2-factor, bi-factor, nested-factor structures
# modification of mvtbifct()
# param = MLE, vector of length np = 2*(d_1+...+d_G); 
# tdata = nxd matrix of t(df) scores; 
# grsize = vector with group sizes: d_1,d_2,...,d_G (G=mgrp)
# df = degree of freedom parameter for multivariate Student t
# if full = F, reduced model is used: Vj = rh2j*V0 + sqrt(1-rh2j^2)*epsj; 
#    Zij = rh1ij*Vj + sqrt(1-rh1ij^2)*epsij
# if full = F, in param: values rh2j (Vj and V0) go first and 
#    then rh1ij (Z_{ij} and Vj: 1st group, then second group etc) 
# if full = F, all parameters are unconditional
# Output: llkv = vector for log-likelihood with dimension=nrow(udat)
mvtbifactllkv=function(param,tdata,grsize,df,full=TRUE)
{ d=ncol(tdata);
  n=nrow(tdata);
  if(full==TRUE)  # bi-factor model
  { mgrp=length(param)/2;
    rh1=param[1:mgrp]; rh2=param[(mgrp+1):(2*mgrp)];
  }
  if(full==FALSE) # nested factor model
  { mgrp=length(grsize);
    ind=0;
    a2=param[(mgrp+1):(mgrp+d)]; a1=rep(0,d);
    for (jg in 1:mgrp)
    { ind1=ind+1;  ind2=ind+grsize[jg];
      a1[ind1:ind2]=param[jg]; 
      ind=ind+grsize[jg];
    }
    rh11=a1*a2;
    rh21=a2*sqrt(1-a1^2)/sqrt(1-rh11^2);
    rh1=rh11; rh2=rh21;
  }
  #nllk = bifcttlik(tdata,grsize,rh1,rh2,df);
  sigmat=bifct(grsize,rh1,rh2); 
  dett=sigmat$fctdet;
  invt=sigmat$fctinv; 
  nllkv=rep(0,n);
  for (i in 1:n)
  { tem=tdata[i,]; 
    if(df<=200)
    { const=lgamma(.5*(d+df))+(d-1)*lgamma(.5*df)-d*lgamma(.5*(1+df));
      nllkv[i]=-const+0.5*(log(dett)+(d+df)*log(1+t(tem)%*%invt%*%tem/df)-
         (1+df)*sum(log(1+tem^2/df)));
    }
    else if(df>200)
    { nllkv[i]=0.5*(log(dett)+t(tem)%*%invt%*%tem-sum(tem^2)); }
  }
  -nllkv
} 

# tri-factor extension : extra argument sbgrsize
# param = MLE, vector of length np = 3*(d_1+...+d_G); 
# tdata = nxd matrix of t(df) scores;
# grsize = vector with group sizes: d_1,d_2,...,d_G 
# sbgrsize = vector with subgroup sizes: 
#  the subgroup sizes are k_{11},...,k_{1d_1}, ... ,k_{G1},...,k_{Gd_G}
#  k_{j1}+...+k_{jd_j}=d_j
# df = degree of freedom parameter for multivariate Student t
# if full = F, reduced model is used
# if full = F, in param: values for (Vk and V0) go first and 
#    then (V_{jk} and Vk: 1st group, then second group etc)
#    then (Z_{ijk} and V_{jk}: 1st subgroup, then second subgroup etc) 
# if full = F, all parameters are unconditional
# Output: llkv = vector for log-likelihood with dimension=nrow(udat)
mvttrifactllkv=function(param,tdata,grsize,sbgrsize,df,full=TRUE)
{ d=ncol(tdata);
  n=nrow(tdata);
  if(full==TRUE)  # tri-factor model
  { mgrp=length(param)/3;
    rh1=param[1:mgrp]; rh2=param[(mgrp+1):(2*mgrp)];
    rh3=param[(2*mgrp+1):(3*mgrp)];
  }
  if(full==FALSE) # nested factor model
  { mgrp=length(grsize);
    msbgrp=length(sbgrsize)
    ind=0; indsb=0;
    a3=param[(mgrp+msbgrp+1):(mgrp+msbgrp+d)]; a2=rep(0,d); a1=rep(0,d);
    for (jsbg in 1:msbgrp)
    { indsb1=indsb+1;  indsb2=indsb+sbgrsize[jsbg];
      a2[indsb1:indsb2]=param[jsbg+mgrp]; 
      indsb=indsb+sbgrsize[jsbg]; 
    }
    for (jg in 1:mgrp)
    { ind1=ind + 1; ind2=ind+grsize[jg];
      a1[ind1:ind2]=param[jg];
      ind=ind+grsize[jg]; 
    }
    rh11=a1*a2*a3;
    rh21=a2*a3*sqrt(1-a1^2)/sqrt(1-rh11^2);
    rh31=a3*sqrt(1-a2^2)/sqrt(1-rh11^2)/sqrt(1-rh21^2);
    rh1=rh11; rh2=rh21; rh3=rh31;
  }
  #nllk = trifcttlik(tdata,grsize,sbgrsize,rh1,rh2,rh3,df)
  sigmat=trifct(grsize,sbgrsize,rh1,rh2,rh3); 
  dett=sigmat$fctdet;
  invt=sigmat$fctinv; 
  nllkv=rep(0,n);
  for (i in 1:n)
  { tem=tdata[i,]; 
    if(df<=200)
    { const=lgamma(.5*(d+df))+(d-1)*lgamma(.5*df)-d*lgamma(.5*(1+df));
      nllkv[i]=nllkv[i]-const+0.5*(log(dett)+(d+df)*
         log(1+t(tem)%*%invt%*%tem/df)-(1+df)*sum(log(1+tem^2/df)));
    }
    else if(df>200)
    { nllkv[i]=nllkv[i]+0.5*(log(dett)+t(tem)%*%invt%*%tem-sum(tem^2)); }
  }
  -nllkv
} 

