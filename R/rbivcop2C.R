# functions for random sample from bivariate copulas with links to C
# 2-parameter copula familes  bb1-bb10

# n = sample size
# cpar = copula parameter 
# type = "qcond" (based on C_{2|1}^{-1}) or "mix" (mixture representation)
# icheck = print flag for summaries
# Ouput: nx2 matrix of bivariate data with uniform(0,1) margins 

# Also can use the multivariate versions
#   rmbb1, rmbb2, rmbb3, rmbb6, rmbb7, rmbb10 with d=2

# cpar with th>0, de>1, 
# type = "qcond" or "mix"
rbb1=function(n,cpar,type="qcond",icheck=F)
{ if(type=="qcond")
  { uvec=rep(0,n)
    vvec=rep(0,n)
    out= .C("rbb1", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
    uu=cbind(out$uvec,out$vvec)
  }
  else # use stochastic representation gamma mix of power of Gumbel
  { th1=1/cpar[1]; de=cpar[2]
    r=rgamma(n,th1)
    xx=rmgum(n,2,de)  
    uu=(-log(xx))/r
    uu=(1+uu)^(-th1)
  }
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>0, de>0, 
# type = "qcond" or "mix"
rbb2=function(n,cpar,type="qcond",icheck=F)
{ if(type=="qcond")
  { uvec=rep(0,n)
    vvec=rep(0,n)
    out= .C("rbb2", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
    uu=cbind(out$uvec,out$vvec)
  }
  else # use stochastic representation gamma mix of power of MTCJ
  { th1=1/cpar[1]; de=cpar[2]
    r=rgamma(n,th1)
    xx=rmtcj(n,de/r)
    uu=(-log(xx))/r
    uu=(1+uu)^(-th1)
    uu
  }
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>1, de>0, 
# type = "qcond" or "mix"
rbb3=function(n,cpar,type="qcond",icheck=F)
{ if(type=="qcond")
  { uvec=rep(0,n)
    vvec=rep(0,n)
    out= .C("rbb3", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
    uu=cbind(out$uvec,out$vvec)
  }
  else # use stochastic representation postable mix of power of MTCJ
  { th1=1/cpar[1]; de=cpar[2]
    r=rpostable(n,th1)
    xx=rmtcj(n,de/r)
    uu=(-log(xx))/r
    uu=exp(-uu^th1)
    uu
  }
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>1, de>1, 
# type = "qcond" or "mix"
rbb6=function(n,cpar,type="qcond",icheck=F)
{ if(type=="qcond")
  { uvec=rep(0,n)
    vvec=rep(0,n)
    out= .C("rbb6", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
    uu=cbind(out$uvec,out$vvec)
  }
  else # use stochastic representation Sibuya mix of power of Gumbel
  { th1=1/cpar[1]; de=cpar[2]
    r=rsibuya(n,th1)
    xx=rgum(n,de)
    uu=1-(1-xx^(1/r))^th1
    uu
  }
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>1, de>0, 
# type = "qcond" or "mix"
rbb7=function(n,cpar,type="qcond",icheck=F)
{ if(type=="qcond")
  { uvec=rep(0,n)
    vvec=rep(0,n)
    out= .C("rbb7", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
    uu=cbind(out$uvec,out$vvec)
  }
  else # use stochastic representation Sibuya mix of power of MTCJ
  { th1=1/cpar[1]; de=cpar[2]
    r=rsibuya(n,th1)
    xx=rmtcj(n,de/r)
    uu=1-(1-xx^(1/r))^th1
    uu
  }
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>1, ga>0
rbb9=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rbb9", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  uu=cbind(out$uvec,out$vvec)
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>0, de>0
rbb4=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rbb4", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  uu=cbind(out$uvec,out$vvec)
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with th>1, de>0
rbb5=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rbb5", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  uu=cbind(out$uvec,out$vvec)
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

# cpar with vth>1, 0<de<=1
rbb8=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rbb8", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  uu=cbind(out$uvec,out$vvec)
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}


# cpar with th>0, 0<ppi<=1, 
rbb10=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rbb10", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  uu=cbind(out$uvec,out$vvec)
  # add "mix" option?
  if(icheck) 
  { cat("\ncpar=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    cat("correlation = ", cor(uu[,1],uu[,2]),"\n")
  }
  uu
}

