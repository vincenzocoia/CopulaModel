# Functions for random sample bivariate copulas 
# 1-parameter copula familes with conditional cdf in closed form
# frk = Frank
# mtcj = Mardia-Takahasi-Clayton-Cook-Johnson

# n = sample size
# cpar = copula parameter 
# icheck = print flag for summaries
# Output: nx2 matrix of bivariate data with uniform(0,1) margins 

# bivariate Frank, cpar >0 or <0
rfrk=function(n,cpar,icheck=F)
{ u1=runif(n)
  p=runif(n)
  u2=qcondfrk(p,u1,cpar)
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(u1); s2=mean(u2)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(u1*u2)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(u1,u2)
}

# bivariate MTCJ, cpar>0 
rmtcj=function(n,cpar,icheck=F)
{ u1=runif(n)
  p=runif(n)
  u2=qcondmtcj(p,u1,cpar)
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(u1); s2=mean(u2)
    cat("average u =", s1,"\n")
    cat("average v =", s2,"\n")
    s12=mean(u1*u2)
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  cbind(u1,u2)
}

