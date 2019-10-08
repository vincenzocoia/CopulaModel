# Copula models for count response data, including discretized MVN
# diagnostics from bivariate margins

# rho = latent correlation parameter
# param = univariate parameter, common to the two margins; assumed
#   estimated in a previous step
# ucdf = function for the univariate cdf
# xdat1 = covariate matrix for margin 1
# xdat2 = covariate matrix for margin 2
#   covariates can be varying over repeated measurements 
#     that is, xdat1, xdat2 can be different 
# y1,y2 = integer-valued responses for variables 1 and 2
# Output : negative log-likelihood
latentBVNnllk1=function(rho,param,ucdf,xdat1,xdat2,y1,y2)
{ if(rho<=-1. | rho>=1.) { return(1.e10) }
  n=length(y1)
  if(is.vector(xdat1)) xdat1=matrix(xdat1,n,1)
  if(is.vector(xdat2)) xdat2=matrix(xdat2,n,1)
  nlk=0; 
  for(i in 1:n)
  { ucdf1=ucdf(y1[i],param,xdat1[i,])
    ucdf2=ucdf(y2[i],param,xdat2[i,])
    ucdf1m=ucdf(y1[i]-1,param,xdat1[i,])
    ucdf2m=ucdf(y2[i]-1,param,xdat2[i,])
    zcdf1=qnorm(ucdf1); zcdf2=qnorm(ucdf2)
    zcdf1m=qnorm(ucdf1m); zcdf2m=qnorm(ucdf2m)
    if(zcdf1m<=-6) zcdf1m=-6
    if(zcdf2m<=-6) zcdf2m=-6
    bpmf= pbnorm(zcdf1,zcdf2,rho) - pbnorm(zcdf1m,zcdf2,rho) -
           pbnorm(zcdf1,zcdf2m,rho) + pbnorm(zcdf1m,zcdf2m,rho)
    if(bpmf<=0.) bpmf=1.e-10
    nlk=nlk-log(bpmf)
  }
  nlk
}

# rho = latent correlation parameter
# par1 = univariate parameter of margin 1; assumed estimated in a previous step
# par2 = univariate parameter of margin 1; assumed estimated in a previous step
# cdf1, cdf2 = functions for univariate cdfs of margins 1, 2
# xdat1, xdat2 = covariate vectors of margins 1,2 
# y1,y2 = integer-valued responses for variables 1 and 2
# Output : negative log-likelihood
latentBVNnllk2=function(rho,par1,par2,cdf1,cdf2,xdat1,xdat2,y1,y2)
{ if(rho<=-1. | rho>=1.) { return(1.e10) }
  n=length(y1)
  if(is.vector(xdat1)) xdat1=matrix(xdat1,n,1)
  if(is.vector(xdat2)) xdat2=matrix(xdat2,n,1)
  nlk=0
  for(i in 1:n)
  { ucdf1=cdf1(y1[i],par1,xdat1[i,])
    ucdf2=cdf2(y2[i],par2,xdat2[i,])
    ucdf1m=cdf1(y1[i]-1,par1,xdat1[i,])
    ucdf2m=cdf2(y2[i]-1,par2,xdat2[i,])
    zcdf1=qnorm(ucdf1); zcdf2=qnorm(ucdf2)
    zcdf1m=qnorm(ucdf1m); zcdf2m=qnorm(ucdf2m)
    if(zcdf1m<=-6) zcdf1m=-6
    if(zcdf2m<=-6) zcdf2m=-6
    bpmf= pbnorm(zcdf1,zcdf2,rho) - pbnorm(zcdf1m,zcdf2,rho) -
           pbnorm(zcdf1,zcdf2m,rho) + pbnorm(zcdf1m,zcdf2m,rho)
    if(bpmf<=0.) bpmf=1.e-10
    nlk=nlk-log(bpmf)
  }
  nlk
}

# independent estimating equations for univariate parameter vector
#  when param is common to the different margins (repeated measurements)
# upmf is the model
# ydat is (n*d)x1 vector of counts
# xdat is (n*d)xq matrix, n=#subjects, d=#repeated measurements, q=#covariates
# LB is lower bound vector for param
# UB is upper bound vector for param
# Output : negative of sum of univariate log-likelihoods
ieenllk=function(param,upmf,ydat,xdat,LB,UB)
{ if(any(param<=LB) | any(param>=UB)) return(1.e10)
  nrec=length(ydat)
  nlk=0
  for(i in 1:nrec)
  { y=ydat[i]; 
    if(is.matrix(xdat)) { x=xdat[i,] } else { x=xdat[i] }
    pr=upmf(y,param,x)
    if(is.nan(pr)) pr=1.e-10
    if(pr<=0) pr=1.e-10
    nlk=nlk-log(pr)   
  }
  nlk
}

# Common univariate parameter for different margins.
# Fit univariate model and assess Exp vs Obs for each margin,
# estimate latent correlation parameter separately for each pair,
# assess bivariate Exp vs Obs.
# Inputs: 
# ydat = (n*d)x1 vector of counts
# xdat = (n*d)xq matrix, n=#subjects, d=#repeated measurements, q=#covariates
#   xdat has d=nrep consecutive rows for each subject/unit
# nrep = d is number of repeated measurements per subject
# upmf,ucdf  = functions for model pmf,cdf
# upmfcdf = function that computes the pmf,cdf up to the value mx
# mx = bound used for Expected vs Observed tables in univariate/bivariate
# ustart = starting parameter point for univariate model
# prlevel = print.level for nlm optimization
# LB = lower bound vector for parameter
# UB = upper bound vector for parameter
# Outputs:
# uparam = univariate parameter estimate 
# length(uparam)>=ncol(xdat)+1; other parameters are not betas 
# rhvec = vector of latent correlation estimates; order is 12,13,23,14,24,...
# E1arr, O1arr = expected and observed counts for univariate
# E2arr, O2arr = expected and observed counts for bivariate
MVNlatent1=function(ydat,xdat,nrep,upmf,ucdf,upmfcdf,mx,ustart,LB,UB,prlevel=0)
{ fit=nlm(ieenllk,p=ustart,hessian=T,print.level=prlevel,upmf=upmf,
     ydat=ydat,xdat=xdat,LB=LB,UB=UB)
  uparam=fit$estimate
  mx1=mx+1
  nrec=length(ydat); n=nrec/nrep # sample size
  if(is.vector(xdat)) xdat=matrix(xdat,nrec,1)
  nc=ncol(xdat); nc1=nc+1
  d=nrep; dd=d*(d-1)/2
  rhvec=rep(0,dd)
  O1arr=matrix(0,mx1,d); E1arr=matrix(0,mx1,d)
  for(j in 1:d)
  { cat("\nUnivariate E vs O, variable ", j, "\n")
    jj=seq(j,nrec,nrep); yy=ydat[jj]; x=xdat[jj,];
    ys=yy;
    ys[ys>=mx]=mx
    ys=c(ys,0:mx)
    Obs1=table(ys)-1
    Exp1=rep(0,mx1)
    for(i in 1:n)
    { if(is.matrix(x)) { out1=upmfcdf(mx,uparam,x[i,]) }
      else { out1=upmfcdf(mx,uparam,x[i]) }
      cdf1=out1[,3]
      cdf1[mx1]=1
      cdf1=c(0,cdf1)
      pmf1=diff(cdf1)
      Exp1=Exp1+pmf1
    }
    E1arr[,j]=Exp1
    O1arr[,j]=Obs1
    cat(Exp1,"\n")
    print(c(Obs1))
    cat("--------------------------------------------------\n")
  }
  O2arr=array(0,c(mx1,mx1,dd)); E2arr=array(0,c(mx1,mx1,dd))
  ii=0
  for(k in 2:d)
  { kk=seq(k,nrec,nrep); y2=ydat[kk]
    for(j in 1:(k-1))
    { jj=seq(j,nrec,nrep); y1=ydat[jj]
      cat("\nBivariate E vs O for ", j,k, "\n")
      rhstart=cor(y1,y2)
      biv=nlm(latentBVNnllk1,p=rhstart,hessian=T,print.level=prlevel,
        param=uparam,ucdf=ucdf,xdat1=xdat[jj,],xdat2=xdat[kk,],y1=y1,y2=y2)
      #print(biv$estimate)
      ii=ii+1; rhvec[ii]=biv$estimate
      y1s=y1; y2s=y2;
      y1s[y1s>=mx]=mx
      y2s[y2s>=mx]=mx
      y1s=c(y1s,0:mx)
      y2s=c(y2s,0:mx)
      Obs=table(y1s,y2s)
      diag(Obs)=diag(Obs)-1
      Exp=matrix(0,mx1,mx1)
      for(i in 1:n)
      { x1=xdat[jj[i],]; x2=xdat[kk[i],]
        out1=upmfcdf(mx,uparam,x1)
        out2=upmfcdf(mx,uparam,x2)
        cdf1=out1[,3]; cdf2=out2[,3]
        zz1=qnorm(cdf1[1:mx]); zz1=c(-6,zz1,6)
        zz2=qnorm(cdf2[1:mx]); zz2=c(-6,zz2,6)
        z1=matrix(zz1,mx1+1,mx1+1)
        z2=matrix(zz2,mx1+1,mx1+1,byrow=T) 
        bcdf=pbnorm(z1,z2,rhvec[ii])
        pmf=apply(bcdf,2,diff)
        pmf2=apply(t(pmf),2,diff)
        pmf2=t(pmf2)
        Exp=Exp+pmf2
      }
      print(Exp)
      print(Obs)
      E2arr[,,ii]=Exp
      O2arr[,,ii]=Obs
      cat("\n============================================================\n")
    }
  }
  list(uparam=uparam,rhvec=rhvec,E1arr=E1arr,O1arr=O1arr,E2arr=E2arr,O2arr=O2arr)
}

# Different univariate parameter for different margins.
# Fit univariate model and assess Exp vs Obs for each margin,
# estimate latent correlation parameter separately for each pair,
# assess bivariate Exp vs Obs.
# Inputs: 
# ydat = (n*d)x1 vector of counts
# xdat = (n*d)xq matrix, n=#subjects, d=#repeated measurements, q=#covariates
#   xdat has d=nrep consecutive rows for each subject/unit
# nrep = d is number of repeated measurements per subject
# unllks = vector of strings for the univariate negative log-likelihoods
# ucdfs = vector of strings for cdf functions for margins 1,...,d 
# upmfcdfs = vector of strings for pmf/cdf functions which compute up to mx
# mx = bound used for Expected vs Observed tables in univariate/bivariate
# ustart = starting parameter point for univariate model
# prlevel = print.level for nlm optimization
# LB = lower bound vector for parameter
# UB = upper bound vector for parameter
# Outputs:
# uparam = matrix with univariate parameter estimates for all univariate margin 
# rhvec = vector of latent correlation estimates; order is 12,13,23,14,24,...
# E1arr, O1arr = expected and observed counts for univariate
# E2arr, O2arr = expected and observed counts for bivariate
MVNlatent2=function(ydat,xdat,nrep,unllks,upmfs,ucdfs,upmfcdfs,mx,ustart,LB,UB,
  prlevel=0)
{ mx1=mx+1
  np=length(ustart)
  nrec=length(ydat); n=nrec/nrep # sample size
  if(is.vector(xdat)) xdat=matrix(xdat,nrec,1)
  nc=ncol(xdat); nc1=nc+1
  d=nrep; dd=d*(d-1)/2
  parmat=matrix(0,np,d)
  rhvec=rep(0,dd)
  O1arr=matrix(0,mx1,d); E1arr=matrix(0,mx1,d)
  for(j in 1:d)
  { jj=seq(j,nrec,nrep); yy=ydat[jj]; x=xdat[jj,];
    unllk=match.fun(unllks[j])
    fit=nlm(unllk,p=ustart,hessian=T,print.level=prlevel,xdat=x,y=yy) 
    parmat[,j]=fit$estimate
    cat("\nUnivariate E vs O, variable ", j, "\n")
    ys=yy;
    ys[ys>=mx]=mx
    ys=c(ys,0:mx)
    Obs1=table(ys)-1
    Exp1=rep(0,mx1)
    upmfcdf=match.fun(upmfcdfs[j])
    for(i in 1:n)
    { if(is.matrix(x)) { out1=upmfcdf(mx,parmat[,j],x[i,]) }
      else { out1=upmfcdf(mx,parmat[,j],x[i]) }
      # out1=upmfcdf(mx,parmat[,j],x[i,])
      cdf1=out1[,3]
      cdf1[mx1]=1
      cdf1=c(0,cdf1)
      pmf1=diff(cdf1)
      Exp1=Exp1+pmf1
    }
    E1arr[,j]=Exp1
    O1arr[,j]=Obs1
    cat(Exp1,"\n")
    print(c(Obs1))
    cat("--------------------------------------------------\n")
  }
  O2arr=array(0,c(mx1,mx1,dd)); E2arr=array(0,c(mx1,mx1,dd))
  ii=0
  for(k in 2:d)
  { kk=seq(k,nrec,nrep); y2=ydat[kk]
    ucdf2=match.fun(ucdfs[k])
    upmfcdf2=match.fun(upmfcdfs[k])
    for(j in 1:(k-1))
    { jj=seq(j,nrec,nrep); y1=ydat[jj]
      ucdf1=match.fun(ucdfs[j])
      upmfcdf1=match.fun(upmfcdfs[j])
      cat("\nBivariate E vs O for ", j,k, "\n")
      rhstart=cor(y1,y2)
      biv=nlm(latentBVNnllk2,p=rhstart,hessian=T,print.level=prlevel,
        par1=parmat[,j],par2=parmat[,k],cdf1=ucdf1,cdf2=ucdf2,
        xdat1=xdat[jj,],xdat2=xdat[kk,],y1=y1,y2=y2)
      #print(biv$estimate)
      ii=ii+1; rhvec[ii]=biv$estimate
      y1s=y1; y2s=y2;
      y1s[y1s>=mx]=mx
      y2s[y2s>=mx]=mx
      y1s=c(y1s,0:mx)
      y2s=c(y2s,0:mx)
      Obs=table(y1s,y2s)
      diag(Obs)=diag(Obs)-1
      Exp=matrix(0,mx1,mx1)
      for(i in 1:n)
      { x1=xdat[jj[i],]; x2=xdat[kk[i],]
        out1=upmfcdf1(mx,parmat[,j],x1)
        out2=upmfcdf2(mx,parmat[,k],x2)
        cdf1=out1[,3]; cdf2=out2[,3]
        zz1=qnorm(cdf1[1:mx]); zz1=c(-6,zz1,6)
        zz2=qnorm(cdf2[1:mx]); zz2=c(-6,zz2,6)
        z1=matrix(zz1,mx1+1,mx1+1)
        z2=matrix(zz2,mx1+1,mx1+1,byrow=T) 
        bcdf=pbnorm(z1,z2,rhvec[ii])
        pmf=apply(bcdf,2,diff)
        pmf2=apply(t(pmf),2,diff)
        pmf2=t(pmf2)
        Exp=Exp+pmf2
      }
      print(Exp)
      print(Obs)
      E2arr[,,ii]=Exp
      O2arr[,,ii]=Obs
      cat("\n============================================================\n")
    }
  }
  list(uparam=parmat,rhvec=rhvec,E1arr=E1arr,O1arr=O1arr,E2arr=E2arr,O2arr=O2arr)
}

# Assume same model for each univariate margin.
# ydat = (n*d)x1 vector of counts
# xdat = (n*d)xq matrix, n=#subjects, d=#repeated measurements, q=#covariates
# nrep = d = consecutive rows for each subject/unit in xdat
# mx = bound used for Expected vs Observed tables in univariate/bivariate
# uparam = univariate parameter estimate 
# upmfcdfs = vector of strings for function pmf/cdf which computes up to mx
# cpar = vector of copula parameters for tree 1, 
#     currently assumes scalar for each pair copula
# A = dxd vine array, with 1:d on diagonal
# pcop = copula cdf with scalar parameter for tree 1 of vine
#   length(uparam)>=ncol(xdat)+1; other parameters are not betas 
# Outputs:
# E2arr, O2arr = expected and observed counts for bivariate margins, tree 1
vinebivEvsO1=function(ydat,xdat,nrep,mx,uparam,upmfcdf,cpar,A,pcop)
{ nrec=nrow(xdat)
  d=nrep
  n=nrec/d
  mx1=mx+1
  # bivariate margins in tree 1 are (A[1,j],j) for j=2,..,d
  O2arr=array(0,c(mx1,mx1,d-1)); E2arr=array(0,c(mx1,mx1,d-1))
  for(j in 2:d)
  { cat("\nUnivariate E vs O, variables ", A[1,j], j, "\n")
    jj=seq(j,nrec,nrep); y1=ydat[jj]
    kk=seq(A[1,j],nrec,nrep); y2=ydat[kk]
    y1s=y1; y2s=y2;
    y1s[y1s>=mx]=mx
    y2s[y2s>=mx]=mx
    y1s=c(y1s,0:mx)
    y2s=c(y2s,0:mx)
    Obs=table(y1s,y2s)
    diag(Obs)=diag(Obs)-1
    Exp=matrix(0,mx1,mx1)
    for(i in 1:n)
    { x1=xdat[jj[i],]; x2=xdat[kk[i],]
      out1=upmfcdf(mx,uparam,x1)
      out2=upmfcdf(mx,uparam,x2)
      cdf1=out1[,3]; cdf2=out2[,3]
      uu1=cdf1[1:mx]; uu1=c(0.000001,uu1,0.999999)  
      # problem with boundary for Gumbel and other copulas
      uu2=cdf2[1:mx]; uu2=c(0.000001,uu2,0.999999)
      u1=matrix(uu1,mx1+1,mx1+1)
      u2=matrix(uu2,mx1+1,mx1+1,byrow=T) 
      # assume pcop can take matrix input
      bcdf=pcop(u1,u2,cpar[j-1])
      pmf=apply(bcdf,2,diff)
      pmf2=apply(t(pmf),2,diff)
      pmf2=t(pmf2)
      Exp=Exp+pmf2
    }
    print(Exp)
    print(Obs)
    E2arr[,,j-1]=Exp
    O2arr[,,j-1]=Obs
    cat("\n============================================================\n")
  }
  list(E2arr=E2arr,O2arr=O2arr)
}

