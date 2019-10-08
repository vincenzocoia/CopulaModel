# Functions for Kullback-Leibler and Jeffreys divergences and
#  Kullback-Leibler sample sizes for two 1-parameter bivariate copula families

# rho = BVN parameter (BVN copula is copula1)
# dcop2 = name of function for bivariate copula density 2
# param2 = copula parameter
# copname2 = copula name
# UB = upper bound for integration (after converting from U(0,1) to N(0,1)
# iprint = print flag for intermediate results
# Outputs: vector with KLcop1true,KLcop2true,Jeffreys,sampsize1,sampsize2
KLcopvsbvn=function(rho,dcop2,param2,copname2="bivcop",UB=7,iprint=F)
{ cv=qnorm(0.95) 
  if(iprint) cat("rho=", rho, "  param=", param2, " copula=", copname2, "\n")

  fcop= function(x,param)
  { u=pnorm(x[1]); v=pnorm(x[2]);
    pdf=dcop2(u,v,param)*dnorm(x[1])*dnorm(x[2])
    pdf
  }
  
  # Jeffreys divergence, integrand cannot have other arguments
  kl1d= function(x)
  { c1=dbvn(x,rho)
    c2=fcop(x,param2)
    c1*log(c1/c2)
  }
  kl2d= function(x)
  { c1=dbvn(x,rho)
    c2=fcop(x,param2)
    c2*log(c2/c1)
  }
  kl1d2= function(x)
  { c1=dbvn(x,rho)
    c2=fcop(x,param2)
    c1*(log(c1/c2))^2
  }
  kl2d2= function(x)
  { c1=dbvn(x,rho)
    c2=fcop(x,param2)
    c2*(log(c2/c1))^2
  }
  # next 3 lines for testing only
  dens= function(x) { fcop(x,param2) }
  out=adaptIntegrate(dens,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  cat("check that this is 1: ",out$integral,"\n")

  Del1=adaptIntegrate(kl1d,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(Del1$integral)
  Del2=adaptIntegrate(kl2d,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(Del2$integral)
  sig1=adaptIntegrate(kl1d2,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(sig1$integral)
  sig2=adaptIntegrate(kl2d2,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(sig2$integral)
  Del1=Del1$integral
  Del2=Del2$integral
  if(iprint) cat("Del1=", Del1, " Del2=", Del2, " Jeffreys=", Del1+Del2,"\n")
  sig1=sig1$integral
  sig2=sig2$integral
  sig1=sqrt(sig1-Del1^2)
  sig2=sqrt(sig2-Del2^2)
  if(iprint) cat("sig1=", sig1, " sig2=", sig2,"\n")
  n1=cv*sig1/Del1; n2=cv*sig2/Del2;
  n1=n1*n1; n2=n2*n2;
  if(iprint) cat(" sample sizes to discriminate with probability 0.95: ", n1, " ", n2,"\n\n")
  out=c(Del1,Del2,Del1+Del2,n1,n2)
  #names(out)=c("KLbvntrue","KLcoptrue","Jeffreys","sampsize1","sampsize2")
  names(out)=c("KLcop1true","KLcop2true","Jeffreys","sampsize1","sampsize2")
  out
}

# ============================================================

# KL for two different copula densities
# KL divergence family c2 relative to c1 : 
# when parameter of c1 is de1  \int c1 * log c1/c2 = Delta
# copname1 = copula name for first copula
# copname2 = copula name for second copula
# param1 = copula1 parameter
# param2 = copula2 parameter
# dcop1 = name of function for bivariate copula density 1
# dcop2 = name of function for bivariate copula density 2
# UB = upper bound for integration (after converting from U(0,1) to N(0,1)
# iprint = print flag for intermediate results
# Outputs: vector with KLcop1true,KLcop2true,Jeffreys,sampsize1,sampsize2
KLcopvscop=function(copname1="cop1",param1,dcop1,copname2="cop2",param2,dcop2,UB=7,iprint=F)
{ cv=qnorm(0.95) 
  if(iprint) cat("\nparam1=", param1, " copula=", copname1, "\n")
  if(iprint) cat("param2=", param2, " copula=", copname2, "\n")
  fcop1= function(x,param1)
  { u=pnorm(x[1]); v=pnorm(x[2]);
    tem=dcop1(u,v,param1)
    if(is.nan(tem)) return(0)
    if(is.infinite(tem)) return(0)  # might not always be the correct thing
    pdf=ifelse(tem<=0,0,tem*dnorm(x[1])*dnorm(x[2]))
    #pdf=dcop1(u,v,param1)*dnorm(x[1])*dnorm(x[2])
    pdf
  }
  fcop2= function(x,param2)
  { u=pnorm(x[1]); v=pnorm(x[2]);
    tem=dcop2(u,v,param2)
    #print(c(x,param2,tem))
    if(is.nan(tem)) return(0)
    if(is.infinite(tem)) return(0)  # might not always be the correct thing 
    pdf=ifelse(tem<=0,0,tem*dnorm(x[1])*dnorm(x[2]))
    #pdf=dcop2(u,v,param2)*dnorm(x[1])*dnorm(x[2])
    pdf
  }
  
  # Jeffreys divergence, integrand cannot have other arguments
  kl1d= function(x)
  { c1=fcop1(x,param1)
    c2=fcop2(x,param2)
    #if(c1<0| c2<0) print(c(x,c1,c2,c1/c2))
    #c1*log(c1/c2)
    if(c1<=0) { out=0 } else if(c2<=0) { out=0. } else { out=c1*log(c1/c2) }
    out
  }
  kl2d= function(x)
  { c1=fcop1(x,param1)
    c2=fcop2(x,param2)
    #c2*log(c2/c1)
    #print(c(x,c1,c2,c1/c2))
    # based on output, set integrand to 0 when negative pdf
    if(c2<=0) { out=0 } else if(c1<=0) { out=0. } else { out=c2*log(c2/c1) }
    out

  }
  kl1d2= function(x)
  { c1=fcop1(x,param1)
    c2=fcop2(x,param2)
    #c1*(log(c1/c2))^2
    if(c1<=0) { out=0 } else if(c2<=0) { out=0. } else { out=c1*(log(c1/c2))^2 }
    out
  }
  kl2d2= function(x)
  { c1=fcop1(x,param1)
    c2=fcop2(x,param2)
    #c2*(log(c2/c1))^2
    if(c2<=0) { out=0 } else if(c1<=0) { out=0. } else { out=c2*(log(c2/c1))^2 }
    out
  }
  # next 6 lines for testing only
  dens1= function(x) { fcop1(x,param1) }
  dens2= function(x) { fcop2(x,param2) }
  out=adaptIntegrate(dens1,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  if(iprint) cat("check that this is 1: ",out$integral,"\n")
  out=adaptIntegrate(dens2,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  if(iprint) cat("check that this is 1: ",out$integral,"\n")
  Del1=adaptIntegrate(kl1d,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(Del1$integral)
  Del2=adaptIntegrate(kl2d,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(Del2$integral)
  sig1=adaptIntegrate(kl1d2,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(sig1$integral)
  sig2=adaptIntegrate(kl2d2,c(-UB,-UB),c(UB,UB),tol=1.e-4)
  #print(sig2$integral)
  Del1=Del1$integral
  Del2=Del2$integral
  if(iprint) cat("Del1=", Del1, " Del2=", Del2, " Jeffreys=", Del1+Del2,"\n")
  sig1=sig1$integral
  sig2=sig2$integral
  sig1=sqrt(sig1-Del1^2)
  sig2=sqrt(sig2-Del2^2)
  if(iprint) cat("sig1=", sig1, " sig2=", sig2,"\n")
  n1=cv*sig1/Del1; n2=cv*sig2/Del2;
  n1=n1*n1; n2=n2*n2;
  if(iprint) cat(" sample sizes to discriminate with probability 0.95: ", n1, " ", n2,"\n")
  out=c(Del1,Del2,Del1+Del2,n1,n2)
  names(out)=c("KLcop1true","KLcop2true","Jeffreys","sampsize1","sampsize2")
  out
}

