# density and log-likelihood functions for 1-factor copula

# This function works for the 1-parameter bivariate copula families;
# it works for m-parameter bivariate if
#   dcop accepts input of form :  scalar, d-vector , dxm matrix
# u0 = latent variable, uvec=vector of dimension d, dcop = copula density
# uvec = vector of length d, components in (0,1)
# dcop = name of function of bivariate copula density
# param = d-vector or mxd matrix, parameter of dcop 
# Output: integrand for 1-factor copula density
d1factcop= function(u0,uvec,dcop,param)
{ # dcop(u0,uvec,param) is a vector
  #if(is.matrix(param)) param=t(param)  # for BB1?
  tem=dcop(u0,uvec,param)
  prod(tem)  # product
}

# wrapper function for 1-factor copula model 
# nq = number of quadrature points
# start = starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)
# udata = nxd matrix of uniform scores
# dcop = name of function for a bivariate copula density
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# prlevel = printlevel for nlm()
# mxiter = maximum number of iteration for nlm()
# Output: nlm object with minimum, estimate, hessian
ml1fact=function(nq,start,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
{ start=c(start)
  np=length(start)
  gl=gausslegendre(nq)
  wl=gl$weights
  xl=gl$nodes
  nllkfn1= function(param)
  { d=ncol(udata)  
    n=nrow(udata)
    nlb=sum(param<=LB | param>=UB)
    if(nlb>0) { return(1.e10) }
    nllk=0
    if(np==d) par0=param
    else { par0=matrix(param,nrow=d,byrow=T) }  # to handle BB1 
    for(i in 1:n)
    { uvec=udata[i,]
      integl=0; 
      for(iq in 1:nq)
      { fvalgl=d1factcop(xl[iq],uvec,dcop,par0)
        integl=integl+wl[iq]*fvalgl
      }
      nllk=nllk-log(integl);
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
    cat("nllk: \n")
    print(mle$minimum)
  }
  mle
}  


# wrapper function for 1-factor copula model (some parameters can be fixed)
# nq = number of quadrature points
# start = starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)
# ifixed = vector of length(param) of True/False, such that
#        ifixed[i]=T iff param[i] is fixed at the given value start[j]
# udata = nxd matrix of uniform scores
# dcop = name of function for a bivariate copula density
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# prlevel = printlevel for nlm()
# mxiter = maximum number of iteration for nlm()
# Output: nlm object with minimum, estimate, hessian
ml1factb=function(nq,start,ifixed,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
{ start=c(start)
  np=length(start)
  gl=gausslegendre(nq)
  wl=gl$weights
  xl=gl$nodes
  nllkfn1= function(param)
  { d=ncol(udata)  
    n=nrow(udata)
    param0=start0; param0[!ifixed]=param  # full length parameter for d1factcop
    np0=length(param0)
    nlb=sum(param0<=LB | param0>=UB)
    if(nlb>0) { return(1.e10) }
    nllk=0
    if(np==d) par0=param0
    else { par0=matrix(param0,nrow=d,byrow=T) }  # to handle BB1 
    for(i in 1:n)
    { uvec=udata[i,]
      integl=0; 
      for(iq in 1:nq)
      { fvalgl=d1factcop(xl[iq],uvec,dcop,par0)
        integl=integl+wl[iq]*fvalgl
      }
      nllk=nllk-log(integl);
    }
    return(nllk);
  }
  start0=start # in order to keep the fixed components
  mle=nlm(nllkfn1,p=start[!ifixed],hessian=T,iterlim=mxiter,print.level=prlevel)
  if(prlevel>0)
  { cat("MLE: \n")
    print(mle$estimate)
    cat("SEs: \n")
    acov=solve(mle$hessian)
    print(sqrt(diag(acov)))
    cat("nllk: \n")
    print(mle$minimum)
  }
  mle
}  

#============================================================

# version for input to pdhessmin and pdhessminb
# for BB1, param is 2d-vector with th1,de1,th2,de2,...
#  param = parameter vector
#  dstruct = list with  data set $data, copula name $copname,
#           $quad is list with quadrature weights and nodes, 
#           $repar is code for reparametrization (for Gumbel, BB1)
#  iprfn = print flag for function and gradient (within NR iterations)
# Output: nllk, grad, hess
f90cop1nllk=function(param,dstruct,iprfn=F)
{ udata=dstruct$data
  copname=dstruct$copname
  copname=tolower(copname)
  gl=dstruct$quad
  repar=dstruct$repar
  if(copname=="t") nu=dstruct$nu # assumed scalar
  d=ncol(udata);
  n=nrow(udata);
  wl=gl$weights
  xl=gl$nodes
  nq=length(xl)
  npar=length(param)
  if(repar==2) { pr0=param; param=param^2+rep(c(0,1),d); } 
  if(repar==1) { pr0=param; param=param^2+1; }
  # check for copula type and call different f90 routines
  if(copname=="frank" | copname=="frk")
  { out= .Fortran("frk1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="lfrank" | copname=="lfrk")
  { out= .Fortran("lfrk1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("gum1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t")
  { tl=qt(xl,nu)   # nu is scalar
    tdata=qt(udata,nu) # transform data and nodes to t(nu) scale
    out= .Fortran("t1fact",
      as.integer(npar), as.double(param), as.double(nu),
      as.integer(d), as.integer(n), as.double(tdata), 
      as.integer(nq), as.double(wl), as.double(tl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  # assume param(2xd) has been converted to a column vector
  else if(copname=="bb1")
  { out= .Fortran("bb11fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else { cat("copname not available\n"); return (NA); }
  if(iprfn) print(cbind(param,out$grad))
  nllk=out$nllk; hess=matrix(out$hess,npar,npar); grad=out$grad;
  if(repar>0) 
  { tem=diag(2*pr0);  # jacobian
    hess=tem%*%hess%*%tem + 2*diag(grad)  
    grad=2*pr0*grad; # pr0 is transformed parameter
  }
  list(fnval=nllk, grad=grad, hess=hess)
}

# version using nlm()
# wrapper function for 1-factor copula model 
# nq = number of quadrature points
# start = starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)
# udata = nxd matrix of uniform scores
# copname = name of copula family such as "gumbel", "frank", "bb1", "t"
# LB = lower bound on parameters (scalar or same dimension as start) 
# UB = upper bound on parameters (scalar or same dimension as start)
# ihess = flag for hessian option in nlm()
# prlevel = printlevel for nlm()
# mxiter = max number of iterations for nlm()
# nu = degree of freedom parameter if copname ="t"
# Output: MLE as nlm object (estimate, Hessian, SEs, nllk)
f90ml1fact=function(nq,start,udata,copname,LB=0,UB=40,ihess=F,prlevel=0,
  mxiter=100,nu=3)
{ # repar hasn't been added to this version
  copname=tolower(copname)
  gl=gausslegendre(nq)
  wl=gl$weights
  xl=gl$nodes
  nllkfn1= function(param)
  { d=ncol(udata)  
    n=nrow(udata)
    npar=length(param)
    nlb=sum(param<=LB | param>=UB)
    if(nlb>0) { return(1.e10) }
    if(copname=="frank" | copname=="frk")
    { out= .Fortran("frk1fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("gum1fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="t")
    { tl=qt(xl,nu)   # nu is scalar
      tdata = qt(udata,nu) # transform data and nodes to t(nu) scale
      out= .Fortran("t1fact",
        as.integer(npar), as.double(param), as.double(nu),
        as.integer(d), as.integer(n), as.double(tdata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    # assume param(2xd) has been converted to a column vector
    else if(copname=="bb1")
    { out= .Fortran("bb11fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
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
  mle=nlm(nllkfn1,p=start,hessian=ihess,iterlim=mxiter,print.level=prlevel,
      check.analyticals=F)
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



