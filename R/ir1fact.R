# Functions for 1-factor copula for item response (ordinal) response
# Updated so that ordinal categories can be 1:ncat or 0:(ncat-1)
#  for log-likelihoods, where ncat = number of ordinal categories.

# simulation of 1-factor model for item response
# ucuts = (ncat-1)xd matrix, increasing in each column and bounded in (0,1)
# n = sample size
# parobj1 = vector of length d or a dxm matrix if the copula has m parameters 
# qcond1 = function for the conditional quantile function of the copula
# copname1 = name of linking copula for factor 1
# Output: d-dimensional random sample of ordinal categories,
#   each value in 0,1,...,(ncat-1)
sim1irfact=function(ucuts,n,parobj1,qcond1,copname1,ivect=F)
{ # should do some checks for matching dim d of ucuts and param
  copname=tolower(copname1)
  yy=sim1fact(n,parobj1,qcond1,copname1,ivect=ivect)
  d=ncol(yy)
  for(j in 1:d)
  { y=cut(yy[,j],c(0,ucuts[,j],1))
    yy[,j]=as.integer(y)-1
  }
  yy
}

# testing
# ucuts=matrix(c(.3,.6,.4,.7,.5,.8),2,3)
# param=c(1,1.5,2)
# set.seed(123)
# ydat=sim1irfact(ucuts,n=5,param,qcond1=qcondfrk,copname1="frank")
# print(ydat)
# for(j in 1:length(param)) print(table(ydat[,j]))

#============================================================

# input:
# dstruct = list with $dat, $cutp, $quad (quadrature points and nodes)
# nq = number of quadrature points for Gauss-Legendre
# dat = data matrix where the number of rows corresponds to an
#   individual's response and each column represents an item
#   Number of ordinal categories for each item, coded as 0,...,(ncat-1)
#                                                     or 1,...,ncat.
#   Currently supported are items that have the same number of categories.
# param = vector with the copula parameters
# cutp = matrix with the cutpoints (without the boundary cutpoints) in the
#   uniform scale, where the number of rows is ncat-1,
#   and each column represents an item.
# pcondcop = function for conditional copula cdf C_{2|1}
#   Assume the same bivariate copula family for all links
# Output: vector with the 1-factor probabilities for all the data points
# This function requires R library abind.
ir1factpmf=function(param,dstruct,pcondcop)
{ dat=dstruct$dat
  if(min(dat)==1) dat=dat-1 # relabel 0:(ncat-1)
  gl=dstruct$quad
  cutp=dstruct$cutp
  nq=length(gl$nodes)
  ncat=nrow(cutp)+1
  d=ncol(cutp)
  n=nrow(dat)
  fden=array(NA,dim=c(nq,ncat,d))
  condcdf=array(NA,dim=c(nq,ncat-1,d))
  for(j in 1:d)
  { for(k in 1:(ncat-1)) condcdf[,k,j]=pcondcop(cutp[k,j],gl$nodes,param[j]) }
  arr0=array(0,dim=c(nq,1,d))
  arr1=array(1,dim=c(nq,1,d))
  condcdf=abind(arr0,condcdf,arr1,along=2)
  for(j in 1:d)
  { for(k in 1:ncat) fden[,k,j]=condcdf[,k+1,j]-condcdf[,k,j] }
  fproduct=matrix(NA,n,nq)
  for(i in 1:n)
  { fproduct.i=1
    for(j in 1:d)
    { temp=fden[,dat[i,j]+1,j]
      fproduct.i=fproduct.i*temp
    }
    fproduct[i,]=fproduct.i
  }
  probs=fproduct%*%gl$w
}

# wrapper function for 1-factor copula model for discrete data 
# nq = number of quadrature points
# start = starting point (d-vector)
# ydata = nxd matrix of discrete data in {0,...,ncat-1} or {1,...,ncat} 
# pcond = function for copula conditional cdf, linking to latent variable 
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# Output: MLE as nlm object
ml1irfact=function(nq,start,ydata,pcond,LB=0,UB=50,ihess=F,prlevel=0,mxiter=50)
{ gl=gausslegendre(nq)
  cutp=unifcuts(ydata)
  if(min(ydata)==1) ydata=ydata-1 # relabel 0:(ncat-1)
  dstruct=list(dat=ydata,quad=gl,cutp=cutp)
  FMlik1= function(theta,dstruct,pcondcop)
  { if(sum(theta<=LB)>0) return(1.e10)
    if(sum(theta>=UB)>0) return(1.e10)
    tem=ir1factpmf(theta,dstruct,pcondcop)
    if(any(tem<=0) || any(is.nan(tem))) return(1.e10)
    loglik=sum(log(tem))
    -loglik
  }
  mle=nlm(FMlik1,p=start,dstruct=dstruct,pcondcop=pcond,hessian=T,
     iterlim=mxiter,print.level=prlevel)
  if(prlevel>0)
  { cat("nllk: \n")
    print(mle$minimum)
    cat("MLE: \n")
    print(mle$estimate)
    if(ihess) 
    { iposdef=isposdef(mle$hessian)
      if(iposdef)
      { cat("SEs: \n")
        acov=solve(mle$hessian)
        SEs=sqrt(diag(acov))
        print(SEs)
        mle$SE=SEs
      }
      else cat("Hessian not positive definite\n")
    }
  }
  mle
}


# wrapper for 1-factor copula model, calls to f90 code
# nq = number of quadrature points
# start = starting point (d-vector)
# ydata = nxd matrix of discrete data in {0,...,ncat-1} or {1,...,ncat} 
# copname = name of model, e.g., "gumbel" or "t"
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# return MLE as nlm object
f90ml1irfact=function(nq,start,ydata,copname,LB=0,UB=40,ihess=F,prlevel=0,
  mxiter=50,nu1=3)
{ gl=gausslegendre(nq)
  ucutp=unifcuts(ydata)
  ncat=nrow(ucutp)+1
  d=ncol(ucutp);
  n=nrow(ydata);
  ydata=as.integer(ydata)  # this becomes a vector
  if(min(ydata)==1) ydata=ydata-1 # relabel 0:(ncat-1)
  copname=tolower(copname)
  wl=gl$weights
  xl=gl$nodes
  FMlik1f90= function(param)
  { if(sum(param<=LB)>0) return(1.e10)
    if(sum(param>=UB)>0) return(1.e10)
    npar=length(param)
    if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("irgum1fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp),
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="normal" | copname=="gaussian")
    { #zcutp=qnorm(ucutp)
      #zl=qnorm(xl)
      out= .Fortran("irgau1fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="t")
    { # convert to tcutp and tl for cutpoints and quad pointd
      #tl=qt(xl,nu)
      out= .Fortran("irt1fact",
        as.integer(npar), as.double(param), as.double(nu1),
        as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else { cat("copname not available\n"); return (NA); }
    nllk=out$nllk
    if(is.nan(nllk) || is.infinite(nllk) ) { return(1.e10) }
    if(any(is.nan(out$grad)) || any(is.infinite(out$grad)) ) { return(1.e10) }
    attr(nllk,"gradient") =out$grad  # for nlm with analytic gradient
    nllk
  }
  mle=nlm(FMlik1f90,p=start,hessian=ihess,iterlim=mxiter,print.level=prlevel,
     check.analyticals = F)
  if(prlevel>0)
  { cat("nllk: \n")
    print(mle$minimum)
    cat("MLE: \n")
    print(mle$estimate)
    if(ihess) 
    { iposdef=isposdef(mle$hessian)
      if(iposdef)
      { cat("SEs: \n")
        acov=solve(mle$hessian)
        SEs=sqrt(diag(acov))
        print(SEs)
        mle$SE=SEs
      }
      else cat("Hessian not positive definite\n")
    }
  }
  mle
}

#============================================================

# negative log-likelihoods and derivatives computed in f90,
# also function for computing Fisher information via f90

# 1-factor copula model for item response;
#  implemented are Gaussian, t (fixed nu) and Gumbel for linking copulas
#    to the latent variable
# version for input to pdhessmin and pdhessminb
# param = parameter vector
# dstruct = list with data set $dat, $cutp for cutpoints in (0,1),
#           copula name $copname,
#           $quad is list with quadrature weights and nodes, 
#           $nu if copname="t"
# iprfn = indicator for printing of function and gradient (within NR iterations)
# Output: nllk, grad, hess
ir1nllk=function(param,dstruct,iprfn=F)
{ ydata=dstruct$dat
  copname=dstruct$copname
  copname=tolower(copname)
  if(copname=="t") nu=dstruct$nu 
  gl=dstruct$quad
  ucutp=dstruct$cutp
  ncat=nrow(ucutp)+1
  d=ncol(ucutp);
  n=nrow(ydata);
  ydata=as.integer(ydata)
  wl=gl$weights
  xl=gl$nodes
  nq=length(xl)
  npar=length(param)
  # check for copula type and call different f90 routines
  if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("irgum1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="normal" | copname=="gaussian")
  { #zcutp=qnorm(ucutp)
    #zl=qnorm(xl)
    out= .Fortran("irgau1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t") # need to add f90 code
  { # convert to tcutp and tl for cutpoints and quad pointd
    #tl=qt(xl,nu)
    #tcutp=qnorm(ucutp)
    out= .Fortran("irt1fact",
      as.integer(npar), as.double(param), as.double(nu),
      as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else { cat("copname not available\n"); return (NA); }
  if(iprfn) print(cbind(param,out$grad))
  list(fnval=out$nllk, grad=out$grad, hess=matrix(out$hess,npar,npar)) 
}

# 1-factor copula model for item response with calculations in f90
# param = vector of parameters for 1-factor copula model
# ucutp = (ncat-1)xd matrix of uniform cutpoints
#    ncat=#categories, d=#items
# nq = number of quadrature points for Gauss-Legendre
# copname = copula name: gumbel, normal, t, ...
# dfdefault = df default value for t model
# nn = nominal sample size for inverse of Fisher information
# Output: Fisher information matrix
f90irfisherinfo1=function(param,ucutp,nq=15,copname,nu1=3,nn=1000)
{ gl=gausslegendre(nq)
  ncat=nrow(ucutp)+1
  d=ncol(ucutp)
  wl=gl$weights
  xl=gl$nodes
  nq=length(xl)
  npar=length(param)
  copname=tolower(copname)
  # check for copula type and call different f90 routines
  if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("irgum1finfo",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="normal" | copname=="gaussian")
  { #zcutp=qnorm(ucutp)
    #zl=qnorm(xl)
    out= .Fortran("irgau1finfo",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t")
  { # convert to tcutp and tl for cutpoints and quad pointd
    #nu=dfdefault
    #tl=qt(xl,nu)
    out= .Fortran("irt1finfo",
      as.integer(npar), as.double(param), as.double(nu1),
      as.integer(d), as.integer(ncat), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else { cat("copname not available\n"); return (NA); }
  finfo=matrix(out$finfo,npar,npar) 
  acov=solve(finfo)
  SEvec=sqrt(diag(acov/nn))
  list(finfo=finfo, SE=SEvec)
}

