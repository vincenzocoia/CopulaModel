# functions for random sample from bivariate copulas with links to C
# 1-parameter copula familes with conditional cdf not in closed form  
# pla = Plackett
# joe = Joe/B5
# gum = Gumbel
# gal = Galambos
# hr = Huelser-Reiss

# n = sample size
# cpar = copula parameter 
# icheck = print flag for summaries
# Output: nx2 matrix of bivariate data with uniform(0,1) margins

# bivariate Plackett, cpar>0
rpla=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rpla", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(out$uvec); s2=mean(out$vvec)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(out$uvec*out$vvec)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(out$uvec,out$vvec)
}

# bivariate Joe/B5, cpar>1
rjoe=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rjoeb5", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(out$uvec); s2=mean(out$vvec)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(out$uvec*out$vvec)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(out$uvec,out$vvec)
}

# bivariate Gumbel, cpar>1
rgum=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rgum", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(out$uvec); s2=mean(out$vvec)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(out$uvec*out$vvec)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(out$uvec,out$vvec)
}

# reflected Gumbel, cpar>1
rgumr=function(n,cpar)
{ out=rgum(n,cpar)
  1-out
}

# bivariate Galambos, cpar>0
rgal=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rgal", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(out$uvec); s2=mean(out$vvec)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(out$uvec*out$vvec)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(out$uvec,out$vvec)
}

# bivariate Huesler-Reiss, cpar>0
rhr=function(n,cpar,icheck=F)
{ uvec=rep(0,n)
  vvec=rep(0,n)
  out= .C("rhr", as.integer(n), as.double(cpar), 
           uvec=as.double(uvec), vvec=as.double(vvec))
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(out$uvec); s2=mean(out$vvec)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(out$uvec*out$vvec)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(out$uvec,out$vvec)
}
