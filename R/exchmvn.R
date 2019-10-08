# Functions for positive exchangeable MVN and derivative with respect
#    to a_k, b_k or rho,
# Extracted from mprobit package

# lb = vector of lower limits of integral/probability  
# ub = vector of upper limits of integral/probability 
# rho = correlation (positive constant over pairs)
# mu = mean vector 
# scale = standard deviation 
# eps = tolerance for numerical integration 
# Output: rectangle probability
exchmvn=function(lb,ub,rho, mu=0,scale=1,eps=1.e-06)
{ #eps=1.e-6 is recommended value because of bound on Romberg iterations
  if(rho<0 || rho>=1) stop("0<=rho<1")
  m=length(ub)
  if(m!=length(lb)) stop("lengths of lb and ub must be the same")
  tem=scale
  a=(lb-mu)/tem
  b=(ub-mu)/tem
  out=.C("r_exchmvn",
    as.integer(m), as.double(a), as.double(b), as.double(rho),
    as.double(eps), pr = as.double(eps))
  out$pr
}

# partial derivative of exchmvn with respect to lower or upper limit
# lb = vector of lower limits of integral/probability  
# ub = vector of upper limits of integral/probability 
# rho = correlation (positive constant over pairs)
# eps = tolerance for numerical integration 
# k = argument of lb/ub for derivative, 
# ksign = 1 for upper limit ub[k], -1 for lower limit lb[k]
# Output: derivative of exchmvn with respect to lb[k] or ub[k]
exchmvn.deriv.margin=function(lb,ub,rho,k,ksign, eps=1.e-06)
{ #eps=1.e-6 is recommended value because of bound on Romberg iterations
  if(rho<0 || rho>=1) stop("0<=rho<1")
  m=length(lb)
  if(m!=length(ub)) stop("lengths of lb and ub must be the same")
  if(k<1) k=1
  if(k>m) k=m
  if(abs(ksign)!=1) ksign=1
  out=.C("r_emvnd",
    as.integer(m), as.double(lb), as.double(ub), as.double(rho),
    as.integer(k), as.integer(ksign),
    as.double(eps), deriv = as.double(eps))
  out$deriv
}

# lb = vector of lower limits of integral/probability  
# ub = vector of upper limits of integral/probability 
# rho = correlation (positive constant over pairs)
# Output: derivative of exchmvn with respect to rho
exchmvn.deriv.rho=function(lb,ub,rho, eps=1.e-06)
{ #eps=1.e-6 is recommended value because of bound on Romberg iterations
  if(rho<0 || rho>=1) stop("0<=rho<1")
  m=length(lb)
  if(m!=length(ub)) stop("lengths of lb and ub must be the same")
  out=.C("r_emvndrh",
    as.integer(m), as.double(lb), as.double(ub), as.double(rho),
    as.double(eps), deriv = as.double(eps))
  out$deriv
}

# test case 
# m=5
# a=rep(-1,m)
# b=rep(2,m)
# rho=.6
# print(exchmvn(a,b,rho))

# answer should be 0.5437533059

# names of functions made longer and self-explanatory

#b[m]=1.5
#print(emvnd(a,b,rho,1,-1))
#print(emvnd(a,b,rho,1,1))
#print(emvnd(a,b,rho,m,-1))
#print(emvnd(a,b,rho,m,1))
#print(emvndrh(a,b,rho))

# answers should be 
#integ:     dera1=-0.084819, derb1=0.025738
#integ:     dera1=-0.085593, derb1=0.093674
#integ.   : derrh=0.417630
