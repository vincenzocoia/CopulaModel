# Multivariate Gaussian and t with 1-factor and 2-factor structures
# Code written by Pavel Krupskii
# This has been replaced by newer, faster versions with gradients of the
# negative log-likelihood. This file is kept for possible back-checks.

# MVN 1-factor model
# z = d-vector in N(0,1) scale
# rhvec = vector of correlations of length d
# Output: log copula density
mvn1fact= function(z,rhvec)
{ rhtr=rhvec/(1-rhvec^2);
  m0=prod(1-rhvec^2);
  m1=1+sum(rhvec*rhtr);
  det=m0*m1;
  sigmainv= -outer(rhtr,rhtr)/m1;
  diag(sigmainv)= 1/(1-rhvec^2) - rhtr^2/m1;
  func= -0.5*(t(z)%*%sigmainv%*%z - sum(z^2) + log(det));
  func
}

# MVt model with 1-factor correlation structure
# z = d-vector in t(df) scale
# rhvec = vector of correlations of length d
# df = degree of freedom parameter for multivariate Student t
# Output: log copula density
mvt1fact= function(z,rhvec,df)
{ d=sum(z+1-z);
  const= lgamma(.5*(d + df))+(d-1)*lgamma(.5*df)-d*lgamma(.5*(1+df));
  rhtr=rhvec/(1-rhvec^2);
  m0=prod(1-rhvec^2);
  m1=1+sum(rhvec*rhtr);
  det=m0*m1;
  sigmainv= -outer(rhtr,rhtr)/m1;
  diag(sigmainv)= 1/(1-rhvec^2) - rhtr^2/m1;
  func= -0.5*log(det) - 0.5*(d+df)*log(1 + t(z)%*%sigmainv%*%z/df);
  func= const + func + 0.5*(1+df)*sum(log(1+z^2/df));
  func
}

# MVN copula density with 2 factors 
# z = d-vector in N(0,1) scale
# rhvec = vector of correlations/partial correlations of length 2*d;
#    the second d components are partial correlations
# Output: log copula density
mvn2fact= function(z,rhvec)
{ d=length(z);    
  rh1=rhvec[1:d];
  rh2=rhvec[(d+1):(2*d)]*sqrt(1-rh1^2); 
  dl=1-rh1^2-rh2^2;
  a11=1+sum(rh1^2/dl);
  a22=1+sum(rh2^2/dl);
  a12=sum(rh1*rh2/dl);
  m0=a11*a22-a12^2;
  det=prod(dl)*m0;
  sigma.inv= -(a11*outer(rh2/dl,rh2/dl) -a12*(outer(rh1/dl,rh2/dl) + 
    outer(rh2/dl,rh1/dl)) + a22*outer(rh1/dl,rh1/dl))/m0;
  diag(sigma.inv)= diag(sigma.inv) + 1/dl;
  func= -0.5*(t(z)%*%sigma.inv%*%z - sum(z^2) + log(det));
  func
}

# MVt copula density with 2 factors
# z = d-vector in t(df) scale
# rhvec = vector of correlations/partial correlations of length 2*d;
#    the second d components are partial correlations
# df = degree of freedom parameter for multivariate Student t
# Output: log copula density
mvt2fact= function(z,rhvec,df)
{ d=length(z);
  rh1=rhvec[1:d];
  rh2=rhvec[(d+1):(2*d)]*sqrt(1-rh1^2); 
  const= lgamma(.5*(d+df))+(d-1)*lgamma(.5*df) - d*lgamma(.5*(1+df));
  dl=1-rh1^2-rh2^2;
  a11=1+sum(rh1^2/dl);
  a22=1+sum(rh2^2/dl);
  a12=sum(rh1*rh2/dl);
  m0=a11*a22-a12^2;
  det=prod(dl)*m0;
  sigma.inv= -(a11*outer(rh2/dl,rh2/dl) -a12*(outer(rh1/dl,rh2/dl) + 
    outer(rh2/dl,rh1/dl)) + a22*outer(rh1/dl,rh1/dl))/m0;
  diag(sigma.inv)= diag(sigma.inv) + 1/dl;
  func= -0.5*log(det) - 0.5*(d+df)*log(1 + t(z)%*%sigma.inv%*%z/df);
  func= const + func + 0.5*(1+df)*sum(log(1+z^2/df));
  func
}


# 1-factor MVt/MVN model
# start = starting point of correlation parameters 
# df = degree of freedom parameter for multivariate Student t
# tdat = nxd matrix of t(df) scores
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# if df>200 MVN model is used
# Output: nlm object with $estimate, $hessian, $minimum
ml1mvtfact= function(start,df,tdat,prlevel=0,mxiter=100)
{ nllkfn1= function(param)
  { d=ncol(tdat)  
    n=nrow(tdat)
    rhvec=param
    nlb=sum(param<=-1 | param>=1)
    if(nlb>0) { return(1.e10) }
    nllk=0
    for(i in 1:n)
    { tvec=tdat[i,];
      if(df<200)  { integl=mvt1fact(tvec,rhvec,df); } 
      if(df>=200) { integl=mvn1fact(tvec,rhvec); }
      nllk=nllk-integl;
    }
    return(nllk);
  }
  
  mle=nlm(nllkfn1,p=start,hessian=T,iterlim=mxiter,print.level=prlevel)
  if(prlevel>0)
  { cat("MLE: \n")
    print(mle$estimate)
    cat("SEs: \n")
    acov=solve(mle$hessian)
    print(sqrt(diag(acov)))
    cat("nllk:  \n")
    print(mle$minimum)
  }
  mle
}  

# 2-factor MVt/MVN model
# start = param= starting point of correlations/partial correlations
# df = degree of freedom parameter for multivariate Student t
# ifixed = vector of length(param) of True/False, such that
#        ifixed[i]=T iff param[i] is fixed at the given value
# tdat = nxd matrix of t(df) scores
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# if df>200, MVN model is used
# Output: nlm object with $estimate, $hessian, $minimum
ml2mvtfact= function(start,df,ifixed,tdat,prlevel=0,mxiter=100)
{ nllkfn2= function(param)
  { d=ncol(tdat)  
    n=nrow(tdat)
    rhvec=param
    rhvec[ifixed]=start[ifixed]
    nlb=sum(param<=-1 | param>=1)
    if(nlb>0) { return(1.e10) }
    nllk=0
    for(i in 1:n)
    { tvec=tdat[i,];
      if(df<200)  {integl=mvt2fact(tvec,rhvec,df);} 
      if(df>=200) {integl=mvn2fact(tvec,rhvec);}
      nllk=nllk-integl;
    }
    return(nllk);
  }
  
  mle=nlm(nllkfn2,p=start,hessian=T,iterlim=mxiter,print.level=prlevel)
  if(prlevel>0)
  { m0 = mle$estimate;
    cat("MLE: \n")
    print(mle$estimate)
    cat("SEs: \n")
    acov=solve(mle$hessian)
    print(sqrt(diag(acov)))
    cat("nllk:  \n")
    print(mle$minimum)
  }
  mle
}
