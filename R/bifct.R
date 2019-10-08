# Functions for bi-factor and tri-factor correlation matrices
# and maximum likelihood for multivariate t with these structures
# Code written by Pavel Krupskii

# grsize = vector with group sizes: d_1,d_2,...,d_G
# rh1 = vector of correlation with global latent variable
# rh2 = vector of partial correlation with group latent variable given global;
#  for variable linked to group g
# Outputs: fctmat = correlation matrix; fctdet = det(fctmat); 
#   fctinv = solve(fctmat)
bifct=function(grsize,rh1,rh2)
{ mgrp=length(grsize); # mgrp=G=#groups
  d=sum(grsize);
  a2=rh2*sqrt(1-rh1^2);  # loadings for factor 2
  A=rep(0,mgrp); Aj=rep(0,d);
  eta=rep(0,mgrp); etaa=rep(0,d);
  dzeta1=rh1^2/(1-rh2^2)/(1-rh1^2); 
  dzeta2=a2/(1-rh2^2)/(1-rh1^2);
  st=0;
  fctdiag=matrix(0,nrow=d,ncol=d);

  for(j in 1:mgrp)
  { ind=(st+1):(st+grsize[j]);
    fctdiag[ind,ind]=1;
    A[j]     = 1-grsize[j]+sum(1/(1-rh2[ind]^2));
    Aj[ind]  = A[j];
    eta[j]   = sum((dzeta2[ind])*(rh1[ind]));  
    etaa[ind]= eta[j]/A[j];     
    st=st+grsize[j];
  }
  fctmat=outer(rh1,rh1)+outer(a2,a2)*fctdiag;
  diag(fctmat)=1;
  xi=dzeta1/rh1-etaa*dzeta2;
  A0=1-sum(eta^2/A)+sum(dzeta1);
  fctdet=A0*prod(A)*prod((1-rh1^2)*(1-rh2^2));
  fctinv= -(1/A0)*outer(xi,xi) - (1/Aj)*outer(dzeta2,dzeta2)*fctdiag;
  diag(fctinv)= 1/(1-rh1^2)/(1-rh2^2) - dzeta2^2/Aj - xi^2/A0;
  list(fctmat=fctmat,fctdet=fctdet,fctinv=fctinv)
}  
    
#version 2, using the inverse and det of a smaller matrix
# see documentation for bifct()
bifct2=function(grsize,rh1,rh2)
{ matp=rh1;
  a2=rh2*sqrt(1-rh1^2);  # loadings for factor 2
  mgrp=length(grsize);
  dvar=sum(grsize);
  st=0;
  for(jg in 1:mgrp)
  { tem=rep(0,dvar);
    tem.ind=(st+1):(st+grsize[jg]);
    tem[tem.ind]=a2[tem.ind];
    matp=cbind(matp,tem);
    st=st+grsize[jg];
  }
  vrho=1-rh1^2-a2^2;
  invrho=1/vrho;
  invrho=diag(invrho);
  prmat=invrho%*%matp;
  #print(t(matp)%*%prmat)
  nmat=t(matp)%*%prmat+diag(rep(1,mgrp+1));  
  fctmat=matp%*%t(matp) + diag(vrho);
  fctdet=det(nmat)*prod(vrho);
  fctinv=invrho - prmat%*%solve(nmat,t(prmat));
  list(fctmat=fctmat,fctdet=fctdet,fctinv=fctinv)
}


# model with 3 nested factors, using the inverse and det of a smaller matrix
# grsize = vector with group sizes: d_1,d_2,...,d_G (G=mgrp)
#  the groups are k_{11},...,k_{1d_1}, ... ,k_{G1},...,k_{Gd_G}
# sbgrsize = vector with subgroup sizes: 
#  the subgroup sizes are k_{11},...,k_{1d_1}, ... ,k_{G1},...,k_{Gd_G}
#  k_{j1}+...+k_{jd_j}=d_j
# rh1 = vector of correlation with global latent variable
# rh2 = vector of partial correlation with group latent variable given global;
#  for variable linked to group g
# rh3 = vector of partial correlation with subgroup latent variable given
#    global and group latent; ith variable linked to subgroup sg
# Outputs: fctmat = correlation matrix; fctdet = det(fctmat); 
#   fctinv = solve(fctmat)
trifct=function(grsize,sbgrsize,rh1,rh2,rh3)
{ matp=rh1;
  a2=rh2*sqrt(1-rh1^2);                # loadings for factor 2
  a3=rh3*sqrt(1-rh1^2)*sqrt(1-rh2^2);  # loadings for factor 3
  mgrp=length(grsize);
  msbgrp=length(sbgrsize);
  dvar=sum(grsize);
  st=0; st2=0; jsbg=1;
  for(jg in 1:mgrp)
  { tem=rep(0,dvar); 
    tem.ind=(st+1):(st+grsize[jg]);
    tem[tem.ind]=a2[tem.ind];
    matp=cbind(matp,tem);
    while (st2< st+grsize[jg])
    { tem2=rep(0,dvar); 
      tem2.ind=(st2+1):(st2+sbgrsize[jsbg]); 
      tem2[tem2.ind]=a3[tem2.ind];
      matp=cbind(matp,tem2);
      st2=st2+sbgrsize[jsbg]; 
      jsbg=jsbg+1;
    }
    st=st+grsize[jg];     
  }
  vrho=1-rh1^2-a2^2-a3^2;
  invrho=1/vrho;
  invrho=diag(invrho);
  prmat=invrho%*%matp;
  #print(t(matp)%*%prmat)
  nmat=t(matp)%*%prmat+diag(rep(1,msbgrp+mgrp+1));  
  fctmat=matp%*%t(matp)+diag(vrho);
  fctdet=det(nmat)*prod(vrho);
  fctinv=invrho- prmat%*%solve(nmat,t(prmat));
  list(fctmat=fctmat,fctdet=fctdet,fctinv=fctinv)
}

#============================================================

# Functions below have been replaced by faster versions in f90
#  with gradients of the negative log-likelihood.
# This file is kept for possible back-checks.
  
# if grsize=c(1,...,1), get the same results as in 1-factor model for rh1;
# (and rh2 cannot be identified)
# if grsize=c(ncol(data)), get the same results as in 2-factor model

# data =  nxd matrix of t(df) scores
# grsize = vector with group sizes: d_1,d_2,...,d_G
# rh1 = vector of correlation with global latent variable
# rh2 = vector of partial correlation with group latent variable given global;
#  for variable linked to group g
# df = degree of freedom parameter for multivariate Student t
# if df>200, the multivariate Gaussian distribution is used .
# Output: negative log-likelihood for bifactor t model;
bifcttnllk= function(data,grsize,rh1,rh2,df)
{ n=nrow(data);
  d=ncol(data);
  sigmat=bifct(grsize,rh1,rh2); 
  dett=sigmat$fctdet;
  invt=sigmat$fctinv; 
  nllk=0;
  for(i in 1:n)
  { tem=data[i,]; 
    if(df<=200)
    { const=lgamma(.5*(d+df))+(d-1)*lgamma(.5*df)-d*lgamma(.5*(1+df));
      nllk=nllk-const+0.5*(log(dett)+(d+df)*log(1+t(tem)%*%invt%*%tem/df)-
         (1+df)*sum(log(1+tem^2/df)));
    }
    if(df>200)
    { nllk=nllk+0.5*(log(dett)+t(tem)%*%invt%*%tem-sum(tem^2)); }
  }
  nllk
}  


# negative log-likelihood for trifactor t model; 
# data = nxd matrix of t(df) scores
# grsize = vector with group sizes: d_1,d_2,...,d_G 
# sbgrsize = vector with subgroup sizes: 
#  the subgroup sizes are k_{11},...,k_{1d_1}, ... ,k_{G1},...,k_{Gd_G}
#  k_{j1}+...+k_{jd_j}=d_j
# rh1 = vector of correlation with global latent variable
# rh2 = vector of partial correlation with group latent variable given global;
#  for variable linked to group g
# rh3 = vector of partial correlation with subgroup latent variable given
#    global and group latent; ith variable linked to subgroup sg
# df = degree of freedom parameter for multivariate Student t
# if df>200, the multivariate Gaussian distribution is used .
# Output: negative log-likelihood for trifactor t model;
trifcttnllk= function(data,grsize,sbgrsize,rh1,rh2,rh3,df)
{ n=nrow(data);
  d=ncol(data);
  sigmat=trifct(grsize,sbgrsize,rh1,rh2,rh3); 
  dett=sigmat$fctdet;
  invt=sigmat$fctinv; 
  nllk=0;
  for(i in 1:n)
  { tem=data[i,]; 
    if(df<=200)
    { const=lgamma(.5*(d+df))+(d-1)*lgamma(.5*df)-d*lgamma(.5*(1+df));
      nllk=nllk-const+0.5*(log(dett)+(d+df)*log(1+t(tem)%*%invt%*%tem/df)-
        (1+df)*sum(log(1+tem^2/df)));
    }
    if(df>200)
    { nllk=nllk+0.5*(log(dett)+t(tem)%*%invt%*%tem-sum(tem^2)); }
  }
  nllk
}          

# tdata = nxd matrix of t(df) scores; 
# start = vector of starting values for *partial* correlations
# start is a vector of length 2*(d_1+...+d_G); 
#  start = c(start.rh1, start.rh2)
#  starting values for rh1 (Z_{ij} and V_0) go first,
#    then starting values for rh2 (Z_{ij} and Vj given V0)
# grsize = vector with group sizes: d_1,d_2,...,d_G (G=mgrp)
# df = degree of freedom parameter for multivariate Student t
# prlevel = print.level for nlm()
# mxiter = max number of iterations for nlm()
# if full = F, reduced model is used: Vj = rh2j*V0 + sqrt(1-rh2j^2)*epsj; 
#    Zij = rh1ij*Vj + sqrt(1-rh1ij^2)*epsij
# if full = F, in param: values rh2j (Vj and V0) go first and 
#    then rh1ij (Z_{ij} and Vj: 1st group, then second group etc) 
# if full = F, all parameters are unconditional
# Output: mle object (output of nlm)
mvtbifct= function(tdata,start,grsize,df,prlevel=0,full=T,mxiter=100)
{ d=ncol(tdata);
  n=nrow(tdata);
  err=0;
  nllkfn= function(param)
  { if(max(abs(param))> 0.999) { return(1.e10) }
    if(full==T)  # bifactor model
    { mgrp=length(param)/2;
      rh1=param[1:mgrp]; rh2=param[(mgrp+1):(2*mgrp)];
    }
    if(full==F) # nested factor model
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
    nllk=bifcttnllk(tdata,grsize,rh1,rh2,df);
    return(nllk)
  } 

  mle= try(nlm(nllkfn,p=start,hessian=F,iterlim=mxiter,print.level=prlevel),T)
  if(is.list(mle)==F) { print("mvtbifct() has crashed"); return(NULL); }
  pmle=mle$estimate
  cat("MLE: \n")
  if(full==T) pmle=matrix(pmle,ncol=2);
  print(pmle)
  #cat("SEs: \n")
  #if (max(abs(pmle)) < 0.995 & rcond(mle$hessian) > 1e-12) 
  #{ acov=solve(mle$hessian); } 
  #else {print("NA")}
  cat("nllk:  \n"); print(mle$minimum)
  mle
}  


# tdata = nxd matrix of t(df) scores; 
# start = vector of starting values for *partial* correlations
# start is a vector of length np = 3*(d_1+...+d_G) if full = T; 
# start is a vector of length np = d_1+...+d_G + mgrp + msbgrp if full = F
# where G=mgrp = number of groups, msbgrp = number of subgroups 
#  start = c(start.rh1, start.rh2, start.rh3)
# starting values for rh1 (Z_{ikg} and V_0) go first,
#    then starting values for rh2 (Z_{ikg} and Vg given V0)
#    then starting values for rh3 (Z_{ikg} and Vkg given V0, Vk)
# grsize = vector with group sizes: d_1,d_2,...,d_G 
# sbgrsize = vector with subgroup sizes: 
#  the subgroup sizes are k_{11},...,k_{1d_1}, ... ,k_{G1},...,k_{Gd_G}
#  k_{g1}+...+k_{gd_g}=d_g
# df = degree of freedom parameter for multivariate Student t
# prlevel = print.level for nlm()
# if full = F, reduced model is used
# if full = F, in param: values for (Vg and V0) go first and 
#    then (V_{kg} and Vg: 1st group, then second group etc)
#    then (Z_{ikg} and V_{kg}: 1st subgroup, then second subgroup etc) 
# if full = F, all parameters are unconditional
# repar = T: transform correlation parameters: r0 = 2*r/(1+r^2), -1 <= r0 <= 1
# mxiter = max number of iterations for nlm()
# Output: mle object (output of nlm)
mvttrifct= function(tdata,start,grsize,sbgrsize,df,prlevel=0,full=T,
  repar=F,mxiter=100)
{ d=ncol(tdata);
  n=nrow(tdata);
  nllkfn= function(param)
  { if(repar==T) { param=2*param/(1+param^2); }
    if(max(abs(param))> 0.999 & repar==F) { return(1.e10) }
    if(full==T)  # trifactor model
    { mgrp=length(param)/3;
      rh1=param[1:mgrp]; rh2=param[(mgrp+1):(2*mgrp)]; 
      rh3=param[(2*mgrp+1):(3*mgrp)];
    }
    if(full==F) # nested factor model
    { mgrp=length(grsize);
      msbgrp=length(sbgrsize);
      ind=0; indsb=0;
      a3=param[(mgrp+msbgrp+1):(mgrp+msbgrp+d)]; a2=rep(0,d); a1=rep(0,d);
      for (jsbg in 1:msbgrp)
      { indsb1=indsb+1;  indsb2=indsb+sbgrsize[jsbg];
        a2[indsb1:indsb2]=param[jsbg+mgrp]; 
        indsb=indsb+sbgrsize[jsbg]; 
      }
      for (jg in 1:mgrp)
      { ind1=ind+1; ind2=ind+grsize[jg];
        a1[ind1:ind2]=param[jg];
        ind=ind+grsize[jg]; 
      }
      rh11= a1*a2*a3;
      rh21= a2*a3*sqrt(1-a1^2)/sqrt(1-rh11^2);
      rh31= a3*sqrt(1-a2^2)/sqrt(1-rh11^2)/sqrt(1-rh21^2);
      rh1=rh11; rh2=rh21; rh3=rh31;
    }
    nllk=trifcttnllk(tdata,grsize,sbgrsize,rh1,rh2,rh3,df);
    return(nllk)
  } 

  mle=nlm(nllkfn,p=start,hessian=T,iterlim=mxiter,print.level=prlevel)
  pmle=mle$estimate
  if(repar==T) { pmle=2*pmle/(1+pmle^2); se.adj=2*(1-pmle^2)/(1+pmle^2)^2; }
  cat("MLE: \n")
  if(full==T) pmle=matrix(pmle,ncol=3);
  print(pmle)
  cat("SEs: \n")
  if (max(abs(pmle))<0.995 & rcond(mle$hessian)>1e-12) 
  { acov=solve(mle$hessian); 
    if(repar==F) print(sqrt(diag(acov))) 
    if(repar==T) print(se.adj*sqrt(diag(acov)))
  } 
  else {print("NA")}
  cat("nllk:  \n"); print(mle$minimum)
  mle
}  

