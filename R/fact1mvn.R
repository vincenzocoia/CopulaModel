# rectangle probability for 1-factor multivariate normal

# lb = vector of lower limits of integral/probability  
# ub = vector of upper limits of integral/probability 
# load1 = loading vector for factor 1;
#   means are 0 and variances are 1;
#   the correlation matrix is: rmat=outer(load1,load1), diag(rmat)=1.
# eps = tolerance for numerical integration 
# B = upper limit of integration and negative of lower limit
# Output: rectangle probability
fact1mvn=function(lb,ub,load1,eps=1.e-6,B=6)
{ sqload1=sqrt(1-load1^2)
  d=length(load1) # same as length of lb, ub
  integfn= function(z)
  { nn=length(z)
    tem=rep(1,nn)
    for(j in 1:d)
    { utmp=(ub[j]-load1[j]*z)/sqload1[j]
      ltmp=(lb[j]-load1[j]*z)/sqload1[j]
      tem=tem*(pnorm(utmp)-pnorm(ltmp))
    }
    #print(tem)
    tem*dnorm(z)
  }
  out=integrate(integfn,-B,B,rel.tol=eps)
  out$value
}

