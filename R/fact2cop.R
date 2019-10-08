# density and log-likelihood functions for 2-factor copula
# gradient computed from f90 is passed to nlm()

#coptab=c("frank","mtcj","mtcjr","fgm")
#cop1name="frank"
#copname="gumbel"
#cop2name="fgm"
#copname %in% coptab
#cop1name %in% coptab && cop2name %in% coptab


# wrapper function for 2-factor copula model, calls to f90 code
# nq = number of quadrature points
# start = starting point (d-vector or dxm matrix)
# udata = nxd matrix of uniform scores
# copname = string for the model : e.g., "frank", "gumbel" or "gumfrank"
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# repar = 1 for Gumbel copula with cpar -> 1+param^2
# repar = 2 for second parameter of BB1 copula with delta -> 1+param^2
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# nu1,nu2 = degree of freedom parameters for factors 1,2 if copname ="t"
# Output: MLE as nlm object (estimate, Hessian, SEs, nllk)
f90ml2fact=function(nq,start,udata,copname,LB=0,UB=40,repar=0,ihess=F,
   prlevel=0,mxiter=100,nu1=3,nu2=3)
{ copname=tolower(copname);
  gl=gausslegendre(nq)
  wl=gl$weights
  xl=gl$nodes
  nllkfn2= function(param)
  { d=ncol(udata)  
    n=nrow(udata)
    npar=length(param)
    if(repar==1) { pr0=param; param=1+param^2; }
    else if(repar==2)
    { pr0=param; ind0=1:(2*d); param[ind0]=param[ind0]^2+rep(c(0,1),d); }
    nlb=sum(param<=LB | param>=UB)
    if(nlb>0) { return(1.e10) }
    if(copname=="frank" | copname=="frk")
    { out= .Fortran("frk2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("gum2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    if(copname=="gumfrank" | copname=="gumfrk")
    { out= .Fortran("gumfrk2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    if(copname=="bb1frank" | copname=="bb1frk")
    { out= .Fortran("bb1frk2fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    if(copname=="t")
    { tl1=qt(xl,nu1)   # nu1 is a scalar
      tl2=qt(xl,nu2)   # nu2 is a scalar
      tdata=qt(udata,nu1) # transform data and nodes to t(nu1) scale
      out= .Fortran("t2fact",
        as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
        as.integer(d), as.integer(n), as.double(tdata), 
        as.integer(nq), as.double(wl), as.double(tl1), as.double(tl2),
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }

    nllk=out$nllk
    if(is.nan(nllk) || is.infinite(nllk) ) { return(1.e10) }
    if(any(is.nan(out$grad)) || any(is.infinite(out$grad)) ) { return(1.e10) }
    gr0=out$grad;
    if(repar==1||repar==2) { gr0[ind0]=2*pr0[ind0]*gr0[ind0]; }
    attr(nllk,"gradient")=gr0  # for nlm with analytic gradient
    nllk
  }
  mle=nlm(nllkfn2,p=start,hessian=ihess,iterlim=mxiter,print.level=prlevel,
     check.analyticals=F)
  if(prlevel>0)
  { cat("nllk: \n")
    print(mle$minimum)
    cat("MLE: \n")
    if(repar<1) { print(mle$estimate); }
    else if(repar==1) { print(mle$estimate^2+1); }
    else if(repar==2) 
    { tem = mle$estimate; d = length(tem)/3; 
      ind0 = 1:(2*d); tem[ind0]=tem[ind0]^2+rep(c(0,1),d); 
      print(tem); 
    }
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
# start = starting point (d-vector or dxm matrix)
# ifixed = vector of length(param) of True/False, such that
#        ifixed[i]=T iff param[i] is fixed at the given value start[j]
# udata = nxd matrix of uniform scores
# copname = string for the model : e.g., "frank", "gumbel" or "gumfrank"
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# nu1,nu2 = degree of freedom parameters for factors 1,2 if copname ="t"
# Output: MLE as nlm object (estimate, Hessian, SEs, nllk)
f90ml2factb=function(nq,start,ifixed,udata,copname,LB=0,UB=40,ihess=F,
   prlevel=0,mxiter=100,nu1=3,nu2=3)
{ copname=tolower(copname);
  gl=gausslegendre(nq)
  wl=gl$weights
  xl=gl$nodes
  nllkfn2= function(param)
  { d=ncol(udata)  
    n=nrow(udata)
    param0=start0; param0[!ifixed]=param  # full length parameter for f90
    np0=length(param0)
    nlb=sum(param0<=LB | param0>=UB)
    if(nlb>0) { return(1.e10) }
    if(copname=="frank" | copname=="frk")
    { out= .Fortran("frk2fact",
        as.integer(np0), as.double(param0), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,np0)),
        hess=as.double(rep(0,np0*np0))  )
    }
    if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("gum2fact",
        as.integer(np0), as.double(param0), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,np0)),
        hess=as.double(rep(0,np0*np0))  )
    }
    if(copname=="gumfrank" | copname=="gumfrk")
    { out= .Fortran("gumfrk2fact",
        as.integer(np0), as.double(param0), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,np0)),
        hess=as.double(rep(0,np0*np0))  )
    }
    if(copname=="bb1frank" | copname=="bb1frk")
    { out= .Fortran("bb1frk2fact",
        as.integer(np0), as.double(param0), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,np0)),
        hess=as.double(rep(0,np0*np0))  )
    }
    if(copname=="t")
    { tl1=qt(xl,nu1)   # nu1 is a scalar
      tl2=qt(xl,nu2)   # nu2 is a scalar
      tdata=qt(udata,nu1) # transform data and nodes to t(nu1) scale
      out= .Fortran("t2fact",
        as.integer(np0), as.double(param0), as.double(nu1), as.double(nu2),
        as.integer(d), as.integer(n), as.double(tdata), 
        as.integer(nq), as.double(wl), as.double(tl1), as.double(tl2),
        nllk=as.double(0.),grad=as.double(rep(0,np0)),
        hess=as.double(rep(0,np0*np0))  )
    }

    nllk=out$nllk
    attr(nllk,"gradient")=out$grad[!ifixed]  # for nlm with analytic gradient
    nllk
  }
  start0=start  # in order to keep the fixed components
  mle=nlm(nllkfn2,p=start[!ifixed],hessian=ihess,iterlim=mxiter,
     print.level=prlevel,check.analyticals=F)
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

# alternative names
ml2fact=f90ml2fact
ml2factb=f90ml2factb

#============================================================


# 2-factor copula, negative log-likelihood and derivatives computed in f90
# currently for Frank, Gumbel and Gumbel/Frank (latter means Gumbel for 
#   factor 1 and Frank for factor 2) and BB1/Frank
#  param = parameter vector
#  dstruct = list with  data set $data, copula name $copname,
#           $quad is list with quadrature weights and nodes, 
#           $repar is code for reparametrization (for Gumbel, BB1)
#  iprfn = print flag for function and gradient (within NR iterations)
# Output: nllk, grad, hess
f90cop2nllk=function(param,dstruct,iprfn=F)
{ udata=dstruct$data
  copname=dstruct$copname
  copname=tolower(copname)
  gl=dstruct$quad
  repar=dstruct$repar
  if(copname=="t" | copname=="tapprox") 
  { nu1=dstruct$nu1; nu2=dstruct$nu2; }
  d=ncol(udata);
  n=nrow(udata);
  wl=gl$weights
  xl=gl$nodes
  nq=length(xl)
  if(repar==2) { ind1=1:(2*d); pr0=param[ind1]; param[ind1]=(param[ind1])^2+rep(c(0,1),d); } 
  else if(repar==1) { pr0=param; param=param^2+1; }
  npar=length(param)
  if(copname=="frank" | copname=="frk")
  { out= .Fortran("frk2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="lfrank" | copname=="lfrk")
  { out= .Fortran("lfrk2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("gum2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumfrank" | copname=="gumfrk")
  { out= .Fortran("gumfrk2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="bb1frank" | copname=="bb1frk")
  { out= .Fortran("bb1frk2fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t")
  { tl1=qt(xl,nu1)   # nu1 is a scalar
    tl2=qt(xl,nu2)   # nu2 is a scalar
    tdata=qt(udata,nu1) # transform data and nodes to t(nu1) scale
    out= .Fortran("t2fact",
      as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
      as.integer(d), as.integer(n), as.double(tdata), 
      as.integer(nq), as.double(wl), as.double(tl1), as.double(tl2), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="tapprox")
  { tl1=qt(xl,nu1)   # nu1 is a scalar
    tl2=qt(xl,nu2)   # nu2 is a scalar
    tdata=qt(udata,nu1) # transform data and nodes to t(nu1) scale
    nu1plusone=nu1+1 # for conditional
    # can replace pp by something better later
    pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
    nipol=length(pp)
    qq=qt(pp,nu1plusone)
    pder=pcderiv(qq,pp)  # get derivs for the interpolation
    # pass qq,pp,pder to f90 and make call to pchev in fortran
    out= .Fortran("t2ipol",
      as.integer(npar), as.double(param), as.double(nu1), as.double(nu2),
      as.integer(d), as.integer(n), as.double(tdata), 
      as.integer(nq), as.double(wl), as.double(tl1), as.double(tl2),
      as.integer(nipol), as.double(qq), as.double(pp), as.double(pder), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else { cat("copname not available\n"); return (NA); }
  if(iprfn) print(cbind(param,out$grad))
  nllk=out$nllk; hess=matrix(out$hess,npar,npar); grad=out$grad;

  if(repar==1)  
  { tem=diag(2*pr0);  # jacobian
    hess=tem%*%hess%*%tem + 2*diag(grad)  
    grad=2*pr0*grad; # pr0 is transformed parameter
  }
  else if(repar==2)  
  { tem=diag(c(2*pr0,rep(1,d)));  # jacobian
    hess=tem%*%hess%*%tem + 2*diag(c(grad[ind1],rep(0,d)))  
    grad[ind1]=2*pr0*grad[ind1]; # pr0 is transformed parameter
  }

  list(fnval=nllk, grad=grad, hess=hess) 
}

