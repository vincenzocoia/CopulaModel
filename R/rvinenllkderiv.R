# D-vine and R-vine code for nllk with gradient wrt parameters 
# More general are rvinenllkder1.trunc() and rvinenllkder2.trunc()

# parvec = vector of length dd=(d*(d-1)-(d-ntrunc-1)*(d-ntrunc))/2; 
#     1-parameter family for each pair-copula
# udat = nxd matrix with uniform scores 
# logdcopdernames = vector of logdcopder (by tree), length ntrunc
# pconddernames = vector of pcondder (by tree), length ntrunc
# logdcopder = vector of names of log copula density with derivatives 
#    with respect to u,v,cpar of log(dcop) 
# pcondder = vector of names of conditional cdfs 
#   (only one needed for permutation symmetric pair-copulas)
#   with derivatives with respect to u,v,cpar of  pcond(v|u;cpar)
# logdcopder,pcondder: output order of a function is 
#        value, partial u, partial v, partial cpar
# LB,UB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: nllk and grad of nllk as attribute (for nlm optimization)
dvinenllkder1.trunc=function(parvec,udat,logdcopdernames,pconddernames,
  LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) 
  { nllk=1.e10; nllkder=rep(1,length(parvec)) 
    attr(nllk,"gradient")=nllkder
    return(nllk)
  }
  ntrunc=length(logdcopdernames)
  d=ncol(udat)  
  ii=0; th=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { th[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  d1=d-1; n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); 
  vder=array(0,c(n,d,4)); vpder=array(0,c(n,d,4)); dcopder=array(0,c(n,d,4));
  dd=(d*(d-1))/2; grad=rep(0,dd)
  # vectors to store derivs for recursions, forward and backward terms
  # mapping from parameter in A array to grad
  qq=0; Q=matrix(0,d,d)
  for(i in 1:(d-1))
  { for(j in (i+1):d) 
    { qq=qq+1; Q[i,j]=qq; Q[j,i]=qq }
  }
  # matrix for derivs of different pcond wrt param
  fder=array(0,c(n,dd,d)); bder=array(0,c(n,dd,d))
  nllk=0
  # tree 1
  logdcopder=match.fun(logdcopdernames[1])
  for(j in 2:d) 
  { dcopder[,j,]=logdcopder(udat[,j-1],udat[,j],th[1,j])
    nllk=nllk-sum(dcopder[,j,1])
    qq=Q[1,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
  }
  # tree 2
  if(ntrunc>=2)
  { pcondder=match.fun(pconddernames[1])
    for(j in 2:d1) 
    { vpder[,j,]=pcondder(udat[,j-1],udat[,j],th[1,j]) 
      vp[,j]=vpder[,j,1]
      qq=Q[1,j]; bder[,qq,j]=vpder[,j,4]
    }
    for(j in 3:d) 
    { vder[,j,]=pcondder(udat[,j],udat[,j-1],th[1,j]) 
      v[,j]=vder[,j,1]
      qq=Q[1,j]; fder[,qq,j]=vder[,j,4]
    }
    #cat("fder\n"); print(fder)
    #cat("bder\n"); print(bder)
    logdcopder=match.fun(logdcopdernames[2])
    for(j in 3:d) 
    { dcopder[,j,]=logdcopder(vp[,j-1],v[,j],th[2,j])
      nllk=nllk-sum(dcopder[,j,1])
      temmat=dcopder[,j,2]*bder[,,j-1]+dcopder[,j,3]*fder[,,j]
      grad=grad-apply(temmat,2,sum)
      qq=Q[2,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
    }
    w=v; wp=vp; wfder=fder; wbder=bder
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcondder=match.fun(pconddernames[ell-1]) 
      logdcopder=match.fun(logdcopdernames[ell])
      for(j in ell:d1) 
      { vpder[,j,]=pcondder(wp[,j-1],w[,j],th[ell-1,j]) 
        vp[,j]=vpder[,j,1]
        bder[,,j]=vpder[,j,2]*wfder[,,j]+ vpder[,j,3]*wbder[,,j-1]
        qq=Q[ell-1,j]; bder[,qq,j]=vpder[,j,4]
      }
      for(j in (ell+1):d) 
      { vder[,j,]=pcondder(w[,j],wp[,j-1],th[ell-1,j])  
        v[,j]=vder[,j,1]
        fder[,,j]=vder[,j,2]*wbder[,,j-1]+ vder[,j,3]*wfder[,,j]
        qq=Q[ell-1,j]; fder[,qq,j]=vder[,j,4]
      }
      #cat("fder\n"); print(fder)
      #cat("bder\n"); print(bder)
      for(j in (ell+1):d) 
      { dcopder[,j,]=logdcopder(vp[,j-1],v[,j],th[ell,j])
        nllk=nllk-sum(dcopder[,j,1])
        temmat=dcopder[,j,2]*bder[,,j-1]+dcopder[,j,3]*fder[,,j]
        grad=grad-apply(temmat,2,sum)
        qq=Q[ell,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
      }
      w=v; wp=vp; wfder=fder; wbder=bder
    }
  }
  dd=(d*(d-1)-(d-ntrunc-1)*(d-ntrunc))/2; # should match length(parvec)
  attr(nllk,"gradient")=grad[1:dd]
  nllk
}

#============================================================


# versions of nllk with gradient for R-vine
# rvinenllkder1.trunc : 1-dimensional parameter for each edge
# rvinenllkder2.trunc : 2-dimensional parameter for each edge in tree 1
#                      and 1-dimensional parameter for remaining edges
#   rvinenllkder1.trunc can be used with bivariate t copulas with
#   fixed shape parameter using logdbvtcop.deriv, pcondbvtcop.deriv with 
#    dfdefault=nu1.
#    add logdbvtcop2.deriv, pcondbvt2.deriv with dfdefault2 etc if
#    shape parameters are different for different trees

# parvec = vector of length dd=(d*(d-1)-(d-ntrunc-1)*(d-ntrunc))/2; 
#     1-parameter family for each pair-copula
# udat = nxd matrix with uniform scores 
# A = dxd vine array with 1:d on diagonal
# logdcopdernames = vector of logdcopder (by tree), length ntrunc
# pconddernames = vector of pcondder (by tree), length ntrunc
# logdcopder = vector of names of log copula density with derivatives 
#    with respect to u,v,cpar of log(dcop) 
# pcondder = vector of names of conditional cdfs 
#   (only one needed for permutation symmetric pair-copulas)
#   with derivatives with respect to u,v,cpar of  pcond(v|u;cpar)
# logdcopder,pcondder: output order of a function is 
#        value, partial u, partial v, partial cpar
# LB,UB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: nllk and grad of nllk as attribute (for nlm optimization)
rvinenllkder1.trunc=function(parvec,udat,A,logdcopdernames,pconddernames,
  LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) 
  { nllk=1.e10; nllkder=rep(1,length(parvec)) 
    attr(nllk,"gradient")=nllkder
    return(nllk)
  }
  ntrunc=length(logdcopdernames)
  d=ncol(A)  
  ii=0; th=matrix(0,d,d)
  for(ell in 1:ntrunc)
  { th[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
    ii=ii+(d-ell)
  }
  out=varray2M(A)
  M=out$mxarray; icomp=out$icomp
  d1=d-1; n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d)
  vder=array(0,c(n,d,4)); vpder=array(0,c(n,d,4)); dcopder=array(0,c(n,d,4));
  dd=(d*(d-1))/2; grad=rep(0,dd)
  # vectors to store derivs for recursions, forward and backward terms
  # mapping from parameter in A array to grad
  qq=0; Q=matrix(0,d,d)
  for(i in 1:(d-1))
  { for(j in (i+1):d) 
    { qq=qq+1; Q[i,j]=qq; Q[j,i]=qq }
  }
  # matrix for derivs of different pcond wrt param
  fder=array(0,c(n,dd,d)); bder=array(0,c(n,dd,d)); sder=array(0,c(n,dd,d))
  nllk=0
  # tree 1
  logdcopder=match.fun(logdcopdernames[1])
  for(j in 2:d) 
  { dcopder[,j,]=logdcopder(udat[,A[1,j]],udat[,j],th[1,j])
    nllk=nllk-sum(dcopder[,j,1])
    qq=Q[1,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
  }
  # tree 2
  if(ntrunc>=2)
  { pcondder=match.fun(pconddernames[1])
    for(j in 2:d) 
    { if(icomp[1,j]==1)  # need not be commented out
      { vpder[,j,]=pcondder(udat[,A[1,j]],udat[,j],th[1,j]) 
        vp[,j]=vpder[,j,1]
        qq=Q[1,j]; bder[,qq,j]=vpder[,j,4]
      }
    }
    for(j in 2:d) 
    { vder[,j,]=pcondder(udat[,j],udat[,A[1,j]],th[1,j])
      v[,j]=vder[,j,1]
      qq=Q[1,j]; fder[,qq,j]=vder[,j,4]
    }
    # convert to loop
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) sder[,,j]=bder[,,M[2,j]] else sder[,,j]=fder[,,A[2,j]] }
    logdcopder=match.fun(logdcopdernames[2])
    for(j in 3:d) 
    { dcopder[,j,]=logdcopder(s[,j],v[,j],th[2,j])
      nllk=nllk-sum(dcopder[,j,1])
      # how to convert this?
      temmat=dcopder[,j,2]*sder[,,j]+dcopder[,j,3]*fder[,,j]
      grad=grad-apply(temmat,2,sum)
      qq=Q[2,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
    }
    w=v; wp=vp; wfder=fder; wbder=bder
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcondder=match.fun(pconddernames[ell-1]) 
      logdcopder=match.fun(logdcopdernames[ell])
      for(j in ell:d) 
      { if(icomp[ell-1,j]==1) # need not be commented out
        { vpder[,j,]=pcondder(s[,j],w[,j],th[ell-1,j])
          vp[,j]=vpder[,j,1]
          bder[,,j]=vpder[,j,2]*wfder[,,j]+ vpder[,j,3]*sder[,,j]
          qq=Q[ell-1,j]; bder[,qq,j]=vpder[,j,4]
        }
      }
      for(j in ell:d) 
      { vder[,j,]=pcondder(w[,j],s[,j],th[ell-1,j])
        v[,j]=vder[,j,1]
        fder[,,j]=vder[,j,2]*sder[,,j]+ vder[,j,3]*wfder[,,j]
        qq=Q[ell-1,j]; fder[,qq,j]=vder[,j,4]
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) { sder[,,j]=bder[,,M[ell,j]] }
        else { sder[,,j]=fder[,,A[ell,j]] } 
      }
      for(j in (ell+1):d) 
      { dcopder[,j,]=logdcopder(s[,j],v[,j],th[ell,j])
        nllk=nllk-sum(dcopder[,j,1])
        # how to convert this?
        temmat=dcopder[,j,2]*sder[,,j]+dcopder[,j,3]*fder[,,j]
        grad=grad-apply(temmat,2,sum)
        qq=Q[ell,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
      }
      w=v; wp=vp; wfder=fder; wbder=bder
    }
  }
  dd=(d*(d-1)-(d-ntrunc-1)*(d-ntrunc))/2; # should match length(parvec)
  attr(nllk,"gradient")=grad[1:dd]
  nllk
}

# parvec is vector of length dd=(d*(d-1)-(d-ntrunc-1)*(d-ntrunc))/2 + (d-1); 
#     2-parameter for each pair-copula for tree 1
#     1-parameter for each pair-copula for remaining trees
# udat = nxd matrix with uniform scores 
# A = dxd vine array with 1:d on diagonal
# logdcopdernames = vector of logdcopder (by tree), length ntrunc
# pconddernames = vector of pcondder (by tree), length ntrunc
# pcondder = vector of names of conditional cdfs 
#   (only one needed for permutation symmetric pair-copulas)
#   with derivatives with respect to u,v,cpar of  pcond(v|u;cpar)
# logdcopder,pcondder: output order of a function is 
#        value, partial u, partial v, partial cpar
# LB,UB = lower/upper bound for parvec: constants or vectors of length(parvec)
# Output: nllk and grad of nllk as attribute (for nlm optimization)
rvinenllkder2.trunc=function(parvec,udat,A,logdcopdernames,pconddernames,
  LB=0,UB=10)
{ if(any(parvec<=LB) | any(parvec>=UB)) 
  { nllk=1.e10; nllkder=rep(1,length(parvec)) 
    attr(nllk,"gradient")=nllkder
    return(nllk)
  }
  ntrunc=length(logdcopdernames)
  d=ncol(A)  
  ii=0; th=matrix(0,d,d)
  th[1,2:d]=parvec[seq(1,2*d-3,2)] 
  th[2:d,1]=parvec[seq(2,2*d-2,2)]  # store second parameter in lower triangle
  ii=2*(d-1)
  if(ntrunc>=2)
  { for(ell in 2:ntrunc)
    { th[ell,(ell+1):d]=parvec[ii+(1:(d-ell))] 
      ii=ii+(d-ell)
    }
  }
  out=varray2M(A)
  M=out$mxarray; icomp=out$icomp
  d1=d-1; n=nrow(udat)
  v=matrix(0,n,d); vp=matrix(0,n,d); s=matrix(0,n,d)
  vder=array(0,c(n,d,5)); vpder=array(0,c(n,d,5)); dcopder=array(0,c(n,d,5));
  dd=(d*(d-1))/2 + d1; grad=rep(0,dd)
  # vectors to store derivs for recursions, forward and backward terms
  # mapping from parameter in A array to grad
  qq=1; Q=matrix(0,d,d)
  for(i in 1:(d-1))
  { if(i==1) { incr=2 } else { incr=1 }
    for(j in (i+1):d) 
    { Q[i,j]=qq; Q[j,i]=qq; qq=qq+incr;  }
  }
  # matrix for derivs of different pcond wrt param
  fder=array(0,c(n,dd,d)); bder=array(0,c(n,dd,d)); sder=array(0,c(n,dd,d))
  nllk=0
  # tree 1
  logdcopder=match.fun(logdcopdernames[1])
  for(j in 2:d) 
  { dcopder[,j,]=logdcopder(udat[,A[1,j]],udat[,j],c(th[1,j],th[j,1]))
    nllk=nllk-sum(dcopder[,j,1])
    qq=Q[1,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
    grad[qq+1]=grad[qq+1]-sum(dcopder[,j,5])  # second parameter tree1
  }
  # tree 2
  if(ntrunc>=2)
  { pcondder=match.fun(pconddernames[1])
    for(j in 2:d) 
    { if(icomp[1,j]==1)  # need not be commented out
      { vpder[,j,]=pcondder(udat[,A[1,j]],udat[,j],c(th[1,j],th[j,1])) 
        vp[,j]=vpder[,j,1]
        qq=Q[1,j]; bder[,qq,j]=vpder[,j,4]
        bder[,qq+1,j]=vpder[,j,5]  # second parameter tree1
      }
    }
    for(j in 2:d) 
    { vder[,j,]=pcondder(udat[,j],udat[,A[1,j]],c(th[1,j],th[j,1]))
      v[,j]=vder[,j,1]
      qq=Q[1,j]; fder[,qq,j]=vder[,j,4]
      fder[,qq+1,j]=vder[,j,5]  # second parameter tree1
    }
    # convert to loop
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) s[,j]=vp[,M[2,j]] else s[,j]=v[,A[2,j]] }
    for(j in 3:d) 
    { if(A[2,j]<M[2,j]) sder[,,j]=bder[,,M[2,j]] else sder[,,j]=fder[,,A[2,j]] }
    logdcopder=match.fun(logdcopdernames[2])
    for(j in 3:d) 
    { dcopder[,j,1:4]=logdcopder(s[,j],v[,j],th[2,j])
      nllk=nllk-sum(dcopder[,j,1])
      temmat=dcopder[,j,2]*sder[,,j]+dcopder[,j,3]*fder[,,j]
      grad=grad-apply(temmat,2,sum)
      qq=Q[2,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
    }
    w=v; wp=vp; wfder=fder; wbder=bder
  }
  # remaining trees
  if(ntrunc>=3)
  { for(ell in 3:ntrunc)
    { pcondder=match.fun(pconddernames[ell-1]) 
      logdcopder=match.fun(logdcopdernames[ell])
      for(j in ell:d) 
      { if(icomp[ell-1,j]==1) # need not be commented out
        { vpder[,j,1:4]=pcondder(s[,j],w[,j],th[ell-1,j])
          vp[,j]=vpder[,j,1]
          bder[,,j]=vpder[,j,2]*wfder[,,j]+ vpder[,j,3]*sder[,,j]
          qq=Q[ell-1,j]; bder[,qq,j]=vpder[,j,4]
        }
      }
      for(j in ell:d) 
      { vder[,j,1:4]=pcondder(w[,j],s[,j],th[ell-1,j])
        v[,j]=vder[,j,1]
        fder[,,j]=vder[,j,2]*sder[,,j]+ vder[,j,3]*wfder[,,j]
        qq=Q[ell-1,j]; fder[,qq,j]=vder[,j,4]
      }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) s[,j]=vp[,M[ell,j]] else s[,j]=v[,A[ell,j]] }
      for(j in (ell+1):d) 
      { if(A[ell,j]<M[ell,j]) { sder[,,j]=bder[,,M[ell,j]] }
        else { sder[,,j]=fder[,,A[ell,j]] } 
      }
      for(j in (ell+1):d) 
      { dcopder[,j,1:4]=logdcopder(s[,j],v[,j],th[ell,j])
        nllk=nllk-sum(dcopder[,j,1])
        temmat=dcopder[,j,2]*sder[,,j]+dcopder[,j,3]*fder[,,j]
        grad=grad-apply(temmat,2,sum)
        qq=Q[ell,j]; grad[qq]=grad[qq]-sum(dcopder[,j,4])
      }
      w=v; wp=vp; wfder=fder; wbder=bder
    }
  }
  dd=(d*(d-1)-(d-ntrunc-1)*(d-ntrunc))/2 + d1; # should match length(parvec)
  attr(nllk,"gradient")=grad[1:dd]
  nllk
}



