# Code written by Pavel Krupskii

# structured factor copula models : functions to input to pdhessmin
# negative log-likelihoods and derivatives computed in f90

# nested structured factor copula model 
# negative log-likelihoods and derivatives computed in f90
# version for input to pdhessmin
# inputs
#  param = parameter vector 
#  dstruct = list with  data set $data, copula name $copname,
#           $quad is list with quadrature weights and nodes, 
#           $repar is code for reparametrization (for Gumbel, BB1)
#           $grsize is a vector with group sizes
#  iprfn =indicator for printing of function and gradient (within NR iterations)
# V_j is the latent variable for group j
# V_0 is the global/common latent variable
# parameters for copulas linking U_{ij} and V_j go *at the end* (i's with j=1 then j=2 etc)
# parameters for copulas linking V_j and V_0 go *first* (j=1,2 etc) 
# Outputs: nllk, grad, hess

f90str1nllk=function(param,dstruct,iprfn=F)
{ udata=dstruct$data
  copname=dstruct$copname
  copname=tolower(copname)
  gl=dstruct$quad
  repar=dstruct$repar
  grsize=dstruct$grsize # vector of group sizes
  mgrp=length(grsize); # total number of groups
  dvar=ncol(udata); # total number of variables
  n=nrow(udata);
  wl=gl$weights
  xl=gl$nodes
  if(copname=="t") 
  { nu=dstruct$nu; # assumed a vector of length 2
    tl1=qt(xl,nu[1]);
    tl2=qt(xl,nu[2]);
    tl=cbind(tl1,tl2);
  }
  else if(copname=="tbb1" | copname == "tgum" | copname == "tgumbel") 
  { nu=dstruct$nu; # assumed a scalar
    nu=nu[1];
    tl=qt(xl,nu);
  }
  nq=length(xl)
  npar=length(param)
  if(repar==2) 
  { pr0=param; param=param^2+1; 
    param[(1+mgrp):(dvar+mgrp)] = param[(1+mgrp):(dvar+mgrp)] - 1;
  } 
  else if(repar==1) { pr0=param; param=param^2+1; }
  # check for copula type and call different f90 routines
  if(copname=="frank" | copname=="frk") 
  { out= .Fortran("strfrk1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbel" | copname=="gum")
  { out= .Fortran("strgum1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="frkgum" | copname=="frankgumbel")
  { out= .Fortran("strfrkgum1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="frkbb1" | copname=="frankbb1")
  { out= .Fortran("strfrkbb1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbb1" | copname=="gumbelbb1")
  { out= .Fortran("strgumbb1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="tbb1")
  { out= .Fortran("strtbb1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), as.double(tl),
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="tgum" | copname=="tgumbel")
  { out= .Fortran("strtgum1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), as.double(tl),
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t")
  { out= .Fortran("strt1",
      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
      as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
      as.integer(nq), as.double(wl), as.double(tl), 
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

# old version
#  param = parameter vector (same order as simnestfact()
# V_j is the latent variable for group j
# V_0 is the global latent variable
#f90strnllk= function(param, dstruct, iprfn=F)
#{ udata=dstruct$data
#  copname=dstruct$copname
#  copname=tolower(copname)
#  gl=dstruct$quad
#  repar=dstruct$repar
#  grsize=dstruct$grsize # vector of group sizes
#  mgrp=length(grsize); # total number of groups
#  dvar=ncol(udata); # total number of variables
#  n=nrow(udata);
#  wl=gl$weights
#  xl=gl$nodes
#  if(copname=="t") 
#  { nu=dstruct$nu; # assumed a vector of length 2
#    tl1=qt(xl,nu[1]);
#    tl2=qt(xl,nu[2]);
#    tl=cbind(tl1,tl2);
#  }
#  nq=length(xl)
#  npar=length(param)
#  #if(repar==2) { pr0=param; param=param^2+rep(c(0,1),mgrp); } 
#  if(repar==1) { pr0=param; param=param^2+1; }
#  # check for copula type and call different f90 routines
#  if(copname=="frank") 
#  { out= .Fortran("strfrk",
#      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
#      as.integer(dvar), as.integer(grsize),as.double(udata), 
#      as.integer(nq), as.double(wl), as.double(xl), 
#      nllk=as.double(0.),grad=as.double(rep(0,npar)),
#      hess=as.double(rep(0,npar*npar))  )
#  }
#  else if(copname=="gumbel")
#  { out= .Fortran("strgum",
#      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
#      as.integer(dvar), as.integer(grsize),as.double(udata), 
#      as.integer(nq), as.double(wl), as.double(xl), 
#      nllk=as.double(0.),grad=as.double(rep(0,npar)),
#      hess=as.double(rep(0,npar*npar))  )
#  }
#  else if(copname=="t")
#  { out= .Fortran("strt",
#      as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
#      as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
#      as.integer(nq), as.double(wl), as.double(tl), 
#      nllk=as.double(0.),grad=as.double(rep(0,npar)),
#      hess=as.double(rep(0,npar*npar))  )
#  }
#  else { cat("copname not available\n"); return (NA); }
#  if(iprfn) print(cbind(param,out$grad))
#  nllk=out$nllk; hess=matrix(out$hess,npar,npar); grad=out$grad;
#  if(repar>0) 
#  { tem=diag(2*pr0);  # jacobian
#    hess=tem%*%hess%*%tem + 2*diag(grad)  
#    grad=2*pr0*grad; # pr0 is transformed parameter
#  }
#  list(fnval=nllk, grad=grad, hess=hess)
#}


# bi-factor structured factor copula model 
# negative log-likelihoods and derivatives computed in f90
# version for input to pdhessmin
# inputs
#  param = parameter vector
#  dstruct = list with  data set $data, copula name $copname,
#           $quad is list with quadrature weights and nodes, 
#           $repar is code for reparametrization (for Gumbel, BB1)
#           $grsize is a vector with group sizes
#  iprfn =indicator for printing of function and gradient (within NR iterations)
# if dstruct$pdf == 1 the function evaluates nllk only 
#       (and returns zero gradient and hessian)
# V_0 is the global latent variable that loads on all variables
# V_j is a latent variable that loads only for variables in group j
# parameters for copulas linking U_{ij} and V_0 go first
# parameters for copulas linking U_{ij} and V_g (j in group g) given V_0 go at the end 
#  (by group g=1,2,..,mgrp etc)
# Outputs: nllk, grad, hess

f90str2nllk=function(param,dstruct,iprfn=F)
{ udata=dstruct$data
  copname=dstruct$copname
  copname=tolower(copname)
  gl=dstruct$quad
  repar=dstruct$repar
  grsize=dstruct$grsize # vector of group sizes
  mgrp=length(grsize); # total number of groups
  dvar=ncol(udata); # total number of variables
  n=nrow(udata);
  wl=gl$weights
  xl=gl$nodes
  ipdf=dstruct$pdf # ipdf=1 for pdf only, 0 otherwise
  if(copname=="t" | copname=="tapprox") 
  { nu=dstruct$nu; # assumed a vector of length 2
    tl1=qt(xl,nu[1]);
    tl2=qt(xl,nu[2]);
    tl=cbind(tl1,tl2);
  }
  nq=length(xl)
  npar=length(param)
  #if(repar==2) { pr0=param; param=param^2+rep(c(0,1),mgrp); } 
  if(repar==1) { pr0=param; param=param^2+1; }
  
  if(ipdf==0)  # derivatives
  { # check for copula type and call different f90 routines
    if(copname=="frank" | copname=="frk") 
    { out= .Fortran("strfrk2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("strgum2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumfrk" | copname=="gumbelfrank")
    { out= .Fortran("strgumfrk2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="bb1frk" | copname=="bb1frank")
    { out= .Fortran("strbb1frk2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="bb1gum" | copname=="bb1gumbel")
    { out= .Fortran("strbb1gum2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }    
    else if(copname=="t")
    { out= .Fortran("strt2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    # add option of copname=="tapprox") using strt1ipol
    else if(copname=="tapprox")
    { nu1plusone=nu[1]+1 # for conditional
      # can replace pp by something better later
      pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
      nipol=length(pp)
      qq=qt(pp,nu1plusone)
      pder=pcderiv(qq,pp)  # get derivs for the interpolation
      # pass qq,pp,pder to f90 and make call to pchev in fortran
      out= .Fortran("strt2ipol",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        as.integer(nipol), as.double(qq), as.double(pp), as.double(pder), 
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
  }

  if(ipdf==1) # no derivatives
  { # check for copula type and call different f90 routines
    if(copname=="frank" | copname=="frk") 
    { out= .Fortran("strfrk2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="gumbel" | copname=="gum")
    { out= .Fortran("strgum2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="gumfrk" | copname=="gumbelfrank")
    { out= .Fortran("strgumfrk2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="bb1frk" | copname=="bb1frank")
    { out= .Fortran("strbb1frk2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="bb1gum" | copname=="bb1gumbel")
    { out= .Fortran("strbb1gum2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="t")
    { out= .Fortran("strt2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        nllk=as.double(0.) )
    }
    else if(copname=="tapprox")
    { nu1plusone=nu[1]+1 # for conditional
      # can replace pp by something better later
      pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
      nipol=length(pp)
      qq=qt(pp,nu1plusone)
      pder=pcderiv(qq,pp)  # get derivs for the interpolation
      # pass qq,pp,pder to f90 and make call to pchev in fortran
      out= .Fortran("strt2ipolnllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        as.integer(nipol), as.double(qq), as.double(pp), as.double(pder), 
        nllk=as.double(0.) )
    }
    else { cat("copname not available\n"); return (NA); }
    nllk=out$nllk; grad=0; hess=0;
  }
  list(fnval=nllk, grad=grad, hess=hess)
}

