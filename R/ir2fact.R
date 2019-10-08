# Functions for 2-factor copula for item response (ordinal) response
# Updated so that ordinal categories can be 1:ncat or 0:(ncat-1)
#  for log-likelihoods, where ncat = number of ordinal categories.

# simulation of 2-factor model for item response
# ucuts = (ncat-1)xd matrix, increasing in each column and bounded in (0,1)
# n = sample size
# parobj1 = d-dimensional vector of dependence parameters for factor 1 
#         or dxm matrix where m is #parameters for the biv copula
# parobj2 =  d-dimensional vector or dxm matrix for factor 2 
# qcond1 = function for the conditional quantile function of factor 1
# qcond2 = function for the conditional quantile function of factor 2
# copname1 = abbreviated name of copula for factor 1
# copname2 = abbreviated name of copula for factor 2
# Output: d-dimensional random sample of ordinal categories,
#   each value in 0,1,...,(ncat-1)
sim2irfact=function(ucuts,n,parobj1,parobj2,qcond1,qcond2,copname1="",
  copname2="",ivect=F)
{ # should do some checks for matching dim d of ucuts and param
  copname1=tolower(copname1)
  copname2=tolower(copname2)
  yy=sim2fact(n,parobj1,parobj2,qcond1,qcond2,copname1,copname2,ivect=ivect)
  d=ncol(yy)
  for(j in 1:d)
  { y=cut(yy[,j],c(0,ucuts[,j],1))
    yy[,j]=as.integer(y)-1
  }
  yy
}

#============================================================

# purpose: calculates the pmf for the 2-factor model 
# we allow two different copula families, one for factor 1 and one for factor 2.
# input:
# dstruct = list with $dat, $cutp, $quad (quadrature points and nodes)
# nq = number of quadrature points for Gauss-Legendre
# dat = data matrix where the number of rows corresponds to an
#   individual's response and each column represents an item
#   Number of ordinal categories for each item, coded as 0,...,(ncat-1)
#                                                     or 1,...,,ncat.
#   Currently supported are items that have the same number of categories.
# param = vector that combines theta,delta where
#   theta: a vector with the copula parameters at the 1st factor
#   delta: a vector with the copula parameters at the 2nd factor
# cutp = matrix with the cutpoints (without the boundary cutpoints) in the
#   uniform scale where the number of rows corresponds to an ordinal category
#   and each column represents an item.
# pcondcop1 = function with the conditional copula cdf C_{2|1} for factor 1
# pcondcop2 = function with the conditional copula cdf C_{2|1} for factor 2
# Output: vector with the 2-factor probabilities for all the data points
# This function requires R library abind.
ir2factpmf=function(param,dstruct,pcondcop1,pcondcop2)
{ dat=dstruct$dat
  if(min(dat)==1) dat=dat-1 # relabel 0:(ncat-1)
  gl=dstruct$quad
  cutp=dstruct$cutp
  nq=length(gl$nodes)
  ncat=nrow(cutp)+1
  d=ncol(cutp)
  n=nrow(dat)
  npar=length(param)
  theta=param[1:d]; delta=param[(d+1):npar]
  condcdf=array(NA,dim=c(nq,ncat-1,d))
  for(j in 1:d)
  { for(k in 1:(ncat-1)) condcdf[,k,j]=pcondcop1(cutp[k,j],gl$nodes,theta[j]) }
  condcdf2=array(NA,dim=c(nq,nq,ncat-1,d))
  for(j in 1:d)
  { for(k in 1:(ncat-1))
    { for(u in 1:nq) condcdf2[,u,k,j]=pcondcop2(condcdf[u,k,j],gl$nodes,delta[j]) }
  }
  arr0=array(0,dim=c(nq,nq,1,d))
  arr1=array(1,dim=c(nq,nq,1,d))
  condcdf2=abind(arr0,condcdf2,arr1,along=3)
  fden2=array(NA,dim=c(nq,nq,ncat,d))
  for(j in 1:d)
  { for(k in 1:ncat) fden2[,,k,j]=condcdf2[,,k+1,j]-condcdf2[,,k,j] }
  fproduct=matrix(NA,n,nq)
  for(i in 1:n)
  { fproduct.i=1
    for(j in 1:d)
    { temp=fden2[,,dat[i,j]+1,j]
      fproduct.i=fproduct.i*temp
    }
    fproduct[i,]=as.vector(fproduct.i%*%gl$w)
  }
  fproduct%*%gl$w
}

# wrapper function for 2-factor copula model for discrete data 
# nq = number of quadrature points
# start = starting point (d-vector)
# ydata = a nxd matrix of discrete data in {0,...,ncat-1} or {1,...,ncat} 
# pcond = a bivariate copula conditional cdf
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# Output: MLE as nlm object
ml2irfact=function(nq,start,ydata,pcond1,pcond2,LB=0,UB=40,ihess=F,prlevel=0,mxiter=50)
{ gl=gausslegendre(nq)
  cutp=unifcuts(ydata)
  if(min(ydata)==1) ydata=ydata-1 # relabel 0:(ncat-1)
  dstruct=list(dat=ydata,quad=gl,cutp=cutp)
  FMlik2= function(param,dstruct,pcondcop1,pcondcop2)
  { if(sum(param<=LB)>0) return(1.e10)
    if(sum(param>=UB)>0) return(1.e10)
    tem=ir2factpmf(param,dstruct,pcondcop1,pcondcop2)
    if(any(tem<=0) || any(is.nan(tem))) return(1.e10)
    loglik=sum(log(tem))
    -loglik
  }
  mle=nlm(FMlik2,p=start,dstruct=dstruct,pcondcop1=pcond1,pcondcop2=pcond2,
      hessian=ihess,iterlim=mxiter,print.level=prlevel)
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


# wrapper for 2-factor copula model, calls to f90 code
# nq = number of quadrature points
# start = starting point (d-vector)
# ydata = a nxd matrix of discrete data in {0,...,ncat-1} or {1,...,ncat} 
# copname = name of model, e.g., "gumbel" or "t"
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# Output: MLE as nlm object
f90ml2irfact=function(nq,start,ydata,copname,LB=0,UB=40,ihess=F,prlevel=0,
  mxiter=50,nu1=5,nu2=5)
{ gl=gausslegendre(nq)
  ucutp=unifcuts(ydata)
  ncat=nrow(ucutp)+1
  d=ncol(ucutp);
  n=nrow(ydata);
  ydata=as.integer(ydata) # this becomes a vector
  if(min(ydata)==1) ydata=ydata-1 # relabel 0:(ncat-1)
  copname=tolower(copname)
  wl=gl$weights
  xl=gl$nodes
  FMlik2f90= function(param)
  { if(sum(param<=LB)>0) return(1.e10)
    if(sum(param>=UB)>0) return(1.e10)
    npar=length(param)
    if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("irgum2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp),
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="normal" | copname=="gaussian")  
    { out= .Fortran("irgau2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="t") 
    { out= .Fortran("irt2fact",
        as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
        as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="tgumbel") 
    { out= .Fortran("irtgum2fact",
        as.integer(npar), as.double(param), as.double(nu1),
        as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumbelt") 
    { out= .Fortran("irgumt2fact",
        as.integer(npar), as.double(param), as.double(nu2),
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
    attr(nllk,"gradient")=out$grad  # for nlm with analytic gradient
    nllk
  }
  mle=nlm(FMlik2f90,p=start,hessian=ihess,iterlim=mxiter,print.level=prlevel,
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

# wrapper for 2-factor copula model, calls to f90 code, some fixed parameters 
# nq = number of quadrature points
# start = starting point (d-vector)
# ifixed = vector of length(param) of True/False, such that
#        ifixed[i]=T iff param[i] is fixed at the given value start[j]
# ydata = a nxd matrix of discrete data in {0,...,K-1} 
# copname = name of the model, e.g., "gumbel" or  "t"
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# Output: MLE as nlm object
f90ml2irfactb=function(nq,start,ifixed,ydata,copname,LB=0,UB=40,ihess=F,
   prlevel=0,mxiter=50,nu1=5,nu2=5)
{ gl=gausslegendre(nq)
  ucutp=unifcuts(ydata)
  ncat=nrow(ucutp)+1
  d=ncol(ucutp);
  n=nrow(ydata);
  ydata=as.integer(ydata)
  copname=tolower(copname)
  wl=gl$weights
  xl=gl$nodes
  FMlik2f90= function(param)
  { param0=start0; param0[!ifixed]=param
    if(sum(param0<=LB)>0) return(1.e10)
    if(sum(param0>=UB)>0) return(1.e10)
    np0=length(param0)
    npar=length(param)
    if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("irgum2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp),
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="normal" | copname=="gaussian")  
    { out= .Fortran("irgau2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="t") 
    { out= .Fortran("irt2fact",
        as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
        as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="tgumbel") 
    { out= .Fortran("irtgum2fact",
        as.integer(npar), as.double(param), as.double(nu1),
        as.integer(d), as.integer(ncat),
        as.integer(n), as.integer(ydata), as.double(ucutp), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumbelt") 
    { out= .Fortran("irgumt2fact",
        as.integer(npar), as.double(param), as.double(nu2),
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
    attr(nllk,"gradient")=out$grad[!ifixed]  # for nlm with analytic gradient
    nllk
  }
  start0=start  # in order to keep the fixed components
  mle=nlm(FMlik2f90,p=start[!ifixed],hessian=ihess,iterlim=mxiter,
    print.level=prlevel,check.analyticals = F)
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

# 2-factor copula model for item response
#  implemented are Gaussian, t (fixed nu) and Gumbel for linking copulas
#    to both latent variables
#   also Gumbel/t and t/Gumbel for links to latent variables 1 and 2 resp,
# version for input to pdhessmin and pdhessminb
# param = parameter vector
# dstruct = list with  data set $dat, $cutp for cutpoints in (0,1),
#           copula name $copname,
#           $quad is list with quadrature weights and nodes, 
#           $nu if copname="t"
# iprfn = print flag for function and gradient (within NR iterations)
# Output: nllk, grad, hess
ir2nllk=function(param,dstruct,iprfn=F)
{ ydata=dstruct$dat
  copname=dstruct$copname
  copname=tolower(copname)
  # if(copname=="t") get dfdefault from dstruct
  if(copname=="t") { nu1=dstruct$nu1; nu2=dstruct$nu2 }
  else if(copname=="tgumbel") { nu1=dstruct$nu1 }
  else if(copname=="gumbelt") { nu2=dstruct$nu2 }
  gl=dstruct$quad
  ucutp=dstruct$cutp
  ncat=nrow(ucutp)+1
  d=ncol(ucutp);
  n=nrow(ydata);
  ydata=as.integer(ydata) # this becomes a vector
  wl=gl$weights
  xl=gl$nodes
  nq=length(xl)
  npar=length(param)
  # check for copula type and call different f90 routines
  if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("irgum2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="normal" | copname=="gaussian") # 
  { #zcutp=qnorm(ucutp)
    #zl=qnorm(xl)
    out= .Fortran("irgau2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t") # get nu1, nu2 from dstruct
  { #cat(npar,nu1,nu2,d,ncat,n,nq,"\n")
    #print(ucutp)
    #print(dim(ydata))
    #print(ydata[1:2])
    out= .Fortran("irt2fact",
      as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
      as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  # add also mixed copulas for factors 1,2
  else if(copname=="tgumbel") # get nu1 from dstruct
  { out= .Fortran("irtgum2fact",
      as.integer(npar), as.double(param), as.double(nu1),
      as.integer(d), as.integer(ncat),
      as.integer(n), as.integer(ydata), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbelt") # get nu2 from dstruct
  { out= .Fortran("irgumt2fact",
      as.integer(npar), as.double(param), as.double(nu2),
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

# 2-factor copula model for item response with calculations in f90
# param = vector of parameters for 2-factor copula model
# ucutp = (ncat-1)xd matrix of uniform cutpoints
#    ncat=#categories, d=#items
# copname = copula name: gumbel, normal, t, ...
# nq = number of quadrature points
# ifixed = vector of length(parobj1)+length(parobj2) indicating which
#    parameters for fixed (e.g., conditional independence for factor 2)
# dfdefault1, dfdefault2 = df default values for t model
# nn = nominal sample size for inverse of Fisher information
# Output: Fisher information matrix
f90irfisherinfo2=function(param,ucutp,copname,nq=15,ifixed,nu1=3,nu2=3,nn=1000)
{ gl=gausslegendre(nq)
  ncat=nrow(ucutp)+1
  d=ncol(ucutp)
  wl=gl$weights
  xl=gl$nodes
  nq=length(xl)
  npar=length(param)
  copname=tolower(copname)
  if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("irgum2finfo",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="normal" | copname=="gaussian")
  { #zcutp=qnorm(ucutp)
    #zl=qnorm(xl)
    out= .Fortran("irgau2finfo",
      as.integer(npar), as.double(param), as.integer(d), as.integer(ncat),
      as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t")
  { # convert to tcutp and tl for cutpoints and quad pointd
    #nu1=dfdefault1
    #nu2=dfdefault2
    #tl1=qt(xl,nu1)
    #tl2=qt(xl,nu2)
    out= .Fortran("irt2finfo",
      as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
      as.integer(d), as.integer(ncat), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbelt")
  { out= .Fortran("irgumt2finfo",
      as.integer(npar), as.double(param), as.double(nu2),
      as.integer(d), as.integer(ncat), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="tgumbel")
  { out= .Fortran("irtgum2finfo",
      as.integer(npar), as.double(param), as.double(nu1),
      as.integer(d), as.integer(ncat), as.double(ucutp), 
      as.integer(nq), as.double(wl), as.double(xl), 
      finfo=as.double(rep(0,npar*npar))  )
  }
  else { cat("copname not available\n"); return (NA); }
  finfo= matrix(out$finfo,npar,npar) 
  finfo=finfo[!ifixed,!ifixed]
  acov=solve(finfo)
  SEvec=sqrt(diag(acov/nn))
  list(finfo=finfo, SE=SEvec)
}

