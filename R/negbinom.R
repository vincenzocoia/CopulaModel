# functions for negative binomial regression models 

# dnbinom(x, size, prob, mu, log = FALSE)
# size: target for number of successful trials
# prob: probability of success in each trial
# mu: alternative parametrization via mean

# NB(theta,xi) has pmf 
# f(y;theta,xi)= { Gamma(theta+y) xi^y \over Gamma(theta) y! (1+xi)^(theta+y) }

# param = b[0] b[1] ... b[nc] theta, nc=ncol(xdat)
#    length(param)=ncol(xdat)+2
#    mu(x)= exp(b[0]+b^T x)), xi(x)=mu(x)/theta; theta fixed
#    sigma2/mu=1/p=1+xi, p=probability parameter
# y = vector of non-negative integers, length(y)=nrow(xdat)
# xdat = matrix of covariates (without intercept)
# Output: negative log-likelihood
nb2nllk=function(param,y,xdat)
{ np=length(param)
  gam=param[np]
  if(gam<=0) return(1.e10)
  bvec=param[2:(np-1)]
  b0=param[1]
  if(is.vector(xdat)) { mu=exp(b0+xdat*bvec) }
  else { mu=exp(b0+xdat%*%bvec) }
  #theta=1/gam (old parametrization)
  theta=gam
  pr=dnbinom(y,size=theta,mu=mu)
  nllk=-sum(log(pr))
  nllk
}

# param=b[0] b[1] ... b[nc] gam, nc=ncol(xdat)
#    length(param)=ncol(xdat)+2
#    mu(x)= exp(b[0]+b^T x)), theta(x)=mu(x)/gam; xi=gam and p=1/(1+xi) fixed 
#    sigma2/mu=1/p=1+xi, p=probability parameter
# y = vector of non-negative integers, length(y)=nrow(xdat)
# xdat = matrix of covariates (without intercept)
# Output: negative log-likelihood
nb1nllk=function(param,y,xdat)
{ np=length(param)
  gam=param[np]
  if(gam<=0) return(1.e10)
  bvec=param[2:(np-1)]
  b0=param[1]
  if(is.vector(xdat)) { mu=exp(b0+xdat*bvec) }
  else { mu=exp(b0+xdat%*%bvec) }
  theta=mu/gam
  pr=dnbinom(y,size=theta,prob=1/(1+gam))
  nllk=-sum(log(pr))
  nllk
}

# This code is OK if counts are not large.
# NB probabilities from 0 to ub, using cumprod
# ub = upper bound value for pmf and cdf values
# theta = convolution parameter (size in dnbinom)
# p = probability parameter in (0,1), p=1/(1+xi)
# Output: 
#  3-column matrix with 0:ub in column 1, pmf in column 2, cdf in column 3
nbpmfcdf=function(ub,theta,p)
{ qq=1-p
  p0=p^theta
  if(ub<=0) return(matrix(c(0,p0,p0),1,3))
  pr=cumprod(c(p0,qq*(theta+1:ub-1)/(1:ub)))
  cdf=cumsum(pr)
  cbind(0:ub,pr,cdf)
} 

# mu(x)= exp(b[0]+b^T x), theta(x)=mu(x)/gam, p=1/(1+xi)
# ub = upper bound value for pmf and cdf values
# param =  b0,bvec,gam=xi, 
#    length(param)=length(x)+2
# x = vector of covariates
# Output: 
#  3-column matrix with 0:ub in column 1, pmf in column 2, cdf in column 3
nb1pmfcdf=function(ub,param,x)
{ np=length(param)
  gam=param[np]
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  theta=mu/gam; p=1/(1+gam)
  qq=1-p
  p0=p^theta
  if(ub<=0) return(matrix(c(0,p0,p0),1,3))
  pr=cumprod(c(p0,qq*(theta+1:ub-1)/(1:ub)))
  cdf=cumsum(pr)
  cbind(0:ub,pr,cdf)
}

# mu(x)= exp(b[0]+b^T x), xi(x)=mu(x)/gam, p=1/(1+xi)
# ub = upper bound value for pmf and cdf values
# param = b0,bvec,gam=theta 
# x is a vector, length(param)=length(x)+2
#    length(param)=length(x)+2
# Output: 
#  3-column matrix with 0:ub in column 1, pmf in column 2, cdf in column 3
nb2pmfcdf=function(ub,param,x)
{ np=length(param)
  gam=param[np]
  bvec=param[2:(np-1)]
  b0=param[1]
  #theta=1/gam (old parametrization)
  theta=gam
  mu=exp(b0+sum(bvec*x)); xi=mu/gam; p=1/(1+xi)
  qq=1-p
  p0=p^theta
  if(ub<=0) return(matrix(c(0,p0,p0),1,3))
  pr=cumprod(c(p0,qq*(theta+1:ub-1)/(1:ub)))
  cdf=cumsum(pr)
  cbind(0:ub,pr,cdf)
}

# mu(x)= exp(b[0]+b^T x), theta(x)=mu(x)/gam, p=1/(1+xi)
# y = non-negative integer
# param = vector (b0, bvec, gam=xi) with xi>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: cdf values
nb1cdf=function(y,param,x)
{ np=length(param)
  gam=param[np]
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  theta=mu/gam
  cdf=pnbinom(y,size=theta,prob=1/(1+gam))
  cdf
}

# mu(x)= exp(b[0]+b^T x), xi(x)=mu(x)/gam, p=1/(1+xi)
# y = non-negative integer
# param = vector (b0, bvec, gam=theta) with theta>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: cdf values
nb2cdf=function(y,param,x)
{ np=length(param)
  gam=param[np]
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  #theta=1/gam (old parametrization)
  theta=gam
  cdf=pnbinom(y,size=theta,mu=mu)
  cdf
}

# mu(x)= exp(b[0]+b^T x), theta(x)=mu(x)/gam, p=1/(1+xi)
# y = non-negative integer
# param = vector (b0, bvec, gam=xi) with xi>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: pmf values
nb1pmf=function(y,param,x)
{ np=length(param)
  gam=param[np]
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  theta=mu/gam
  pmf=dnbinom(y,size=theta,prob=1/(1+gam))
  pmf
}

# mu(x)= exp(b[0]+b^T x), xi(x)=mu(x)/gam, p=1/(1+xi)
# y = non-negative integer
# param = vector (b0, bvec, gam=theta) with theta>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: pmf values
nb2pmf=function(y,param,x)
{ np=length(param)
  gam=param[np]
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  #theta=1/gam (old parametrization)
  theta=gam
  pmf=dnbinom(y,size=theta,mu=mu)
  pmf
}


