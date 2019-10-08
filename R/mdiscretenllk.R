# Functions for negative log-likelihoods of multivariate discrete models based 
# on copula and exchangeanle multivariate normal dependence

# cpar = parameter of multivariate copula cdf pmcop
# uudat has dimension nx(2d) with corners of rectangle 
#     d = #response variable or cluster size, 
# pmcop = function for multivariate copula
# LB = vector of lower bounds
# UB = vector of upper bounds
# Output: nllk for discrete with pr=rectangle probability
mdiscretenllk=function(cpar,uudat,pmcop,LB,UB)
{ d=ncol(uudat)/2
  if(any(cpar<=LB) | any(cpar>=UB)) return(1.e10)
  n=nrow(uudat)
  nllk=0
  jj1=1:d; jj2=(d+1):(2*d)
  for(i in 1:n)
  { u1=uudat[i,jj1]; u2=uudat[i,jj2]
    pr=rectmult(u1,u2,cpar,pmcop)
    if(is.na(pr) | pr<=0.)  pr=1.e-15
    nllk=nllk-log(pr)
  }
  nllk  
}

# Exchangeable MVN rectangle probability 
# cpar = positive correlation parameter
# zzdat has dimension nx(2d) with corners of rectangle in N(0,1) scale
#   d = #response variables or cluster size, 
# Output: nllk for discretrized positive exchangeable MVN 
emvndiscretenllk=function(cpar,zzdat)
{ d=ncol(zzdat)/2
  if(cpar<=0. | cpar>=1.) return(1.e10)
  rho=cpar
  n=nrow(zzdat)
  nllk=0
  jj1=1:d; jj2=(d+1):(2*d)
  for(i in 1:n)
  { z1=zzdat[i,jj1]; z2=zzdat[i,jj2]
    pr=exchmvn(z1,z2,rho)
    if(is.na(pr) | pr<=0.)  pr=1.e-15
    nllk=nllk-log(pr)
  }
  nllk  
}


