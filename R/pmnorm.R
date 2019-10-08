# bivariate and trivariate Gaussian/normal rectangle probabilities 
#  using code of Schervish (1984).
# lb = vector of lower limits of integral/probability  
# ub = vector of upper limits of integral/probability 
# mu = mean vector 
# sigma = covariance matrix 
# eps = tolerance for accuracy
# Outputs: 
#   $prob = rectangle probability
#   $bound = error bound
#   $err = error code
pmnorm=function(lb,ub,mu,sigma, eps=1.e-05)
{ n=length(lb)
  if(n!=length(ub))
  stop("lengths of lb and ub must be the same")
  tem=sqrt(diag(sigma))
  a=(lb-mu)/tem
  b=(ub-mu)/tem
  tem=diag(1/tem)
  corr=tem%*%sigma%*%tem
  sig=NULL
  for(i in (2:n)) { sig=c(sig, corr[i, 1:(i-1)]) }
  inf = 1*(b> -9) + 2*(a<9) - 1
  
  bound=.5*eps
  out= .C("mvnscher",
    as.double(b), as.double(a), as.double(sig),
    as.double(eps), as.integer(n), as.integer(inf),
    prob=as.double(eps), perr=as.double(bound), ifault=as.integer(n))
  out=list(pr=out$prob, bound=out$perr, err=out$ifault)
  out
}
