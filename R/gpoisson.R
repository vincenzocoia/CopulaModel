# functions for generalized Poisson regression models

# f(y;theta,vrh)= theta*(theta+vrh*y)^(y-1) * exp(-theta-vrh*y)/y!
#            = theta^y * (1+vrh*y/theta)^(y-1) * exp(-theta-vrh*y)/y!

# theta=th>0, 0<=vrh<1; vrh=0 for Poisson, y in {0,1,2,...}
# xi=(overdispersion index minus one)=(1-vrh)^{-2}-1 or
# 1-vrh=sqrt(1/(1+xi))


# Generalized Poisson pmf
# y = non-negative integer vector or scalar
# Output: pmf values
dgpois=function(y,param)
{ if(is.matrix(param)) { th=param[,1]; vrh=param[,2] }
  else { th=param[1]; vrh=param[2] }
  tem=(y-1)*log(th+vrh*y)-th-vrh*y+log(th)-lgamma(y+1)
  exp(tem)
}

# Generalized Poisson cdf
# y = non-negative integer
# param=(theta,vrh) where theta>0, 0<=vrh<1; 
#   vrh=0 for Poisson, 
# Output: cdf values
pgpois=function(y,param)
{ if(y<0) return(0)
  th=param[1]; vrh=param[2]
  p0=exp(-th)
  if(y==0) return(p0) 
  ii=1:y
  tem=exp(-vrh)*th/ii
  pp=cumprod(c(p0,tem))
  mul=c(1,(1+vrh*ii/th)^(ii-1))
  pr=pp*mul
  sum(pr)
}

# This will work if mean of distribution is not too large;
#  that is, this function is meant for small counts.
# Generalized Poisson pmf and cdf from 0 to ub
# ub = upper bound value for pmf and cdf values
# theta = convolution parameter
# vrh = other parameter 
# Output: 
#  3-column matrix with 0:ub in column 1, pmf in column 2, cdf in column 3
gpoispmfcdf=function(ub,theta,vrh)
{ p0=exp(-theta)
  if(ub<=0) return(matrix(c(0,p0,p0),1,3))
  ii=1:ub
  tem=exp(-vrh)*theta/ii
  pp=cumprod(c(p0,tem))
  mul=c(1,(1+vrh*ii/theta)^(ii-1))
  pr=pp*mul
  cdf=cumsum(pr)
  cbind(0:ub,pr,cdf)
}

#============================================================

# reparametrization for GP regression
# theta=convolution parameter, vrh is second parameter linked to overdispersion
# mu=theta/(1-vrh), sigma2=theta/(1-vrh)^3,
# sigma2/mu=1/(1-vrh)^2 = 1+xi
# xi=1/(1-vrh)^2-1 >0,  vrh=1-sqrt(1/(1+xi))
# GP1: vrh,xi is fixed and theta is function of covariates through mu
# GP2: theta is fixed and vrh is function of covariates through mu

# mu(x)= exp(b[0]+b^T x), theta(x)=mu(x)*(1-vrh), 1-vrh=sqrt(1/(1+xi))
# y = non-negative integer, 
# param = vector  (b0, bvec, xi) with xi>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: cdf values
gp1cdf=function(y,param,x)
{ np=length(param)
  xi=param[np]; vrh=1-sqrt(1/(1+xi))
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  theta=mu*(1-vrh)
  cdf=pgpois(y,c(theta,vrh))
  cdf
}

# mu(x)= exp(b[0]+b^T x), 1-vrh(x)=theta/mu(x), 1-vrh=sqrt(1/(1+xi))
# y = non-negative integer, 
# param = vector  (b0, bvec, theta) with theta>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: cdf values
gp2cdf=function(y,param,x)
{ np=length(param)
  theta=param[np]; 
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  vrh=1-theta/mu
  cdf=pgpois(y,c(theta,vrh))
  cdf
}

# mu(x)= exp(b[0]+b^T x), theta(x)=mu(x)*(1-vrh), 1-vrh=sqrt(1/(1+xi))
# y = non-negative integer, 
# param = vector  (b0, bvec, xi) with xi>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: pmf values
gp1pmf=function(y,param,x)
{ np=length(param)
  xi=param[np]; vrh=1-sqrt(1/(1+xi))
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  theta=mu*(1-vrh)
  pmf=dgpois(y,c(theta,vrh))
  pmf
}

# mu(x)= exp(b[0]+b^T x), 1-vrh(x)=theta/mu(x), 1-vrh=sqrt(1/(1+xi))
# y = non-negative integer, 
# param = vector  (b0, bvec, theta) with theta>0 
# x = covariate vector, 
#   length(param)=length(x)+2
# Output: pmf values
gp2pmf=function(y,param,x)
{ np=length(param)
  theta=param[np]; 
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  vrh=1-theta/mu
  pmf=dgpois(y,c(theta,vrh))
  pmf
}

# mu(x)= exp(b[0]+b^T x)), 1-vrh(x)=theta/mu(x), 1-vrh=sqrt(1/(1+xi))
# param = b[0] b[1] ... b[nc] theta>0, nc=ncol(xdat)
#    length(param)=ncol(xdat)+2
# y = vector of non-negative integers, length(y)=nrow(xdat)
# xdat = matrix of covariates (without intercept)
# Output: negative log-likelihood
gp2nllk=function(param,y,xdat)
{ np=length(param)
  theta=param[np]
  if(theta<=0) return(1.e10)
  bvec=param[2:(np-1)]
  b0=param[1]
  if(is.vector(xdat)) { mu=exp(b0+xdat*bvec) }
  else { mu=exp(b0+xdat%*%bvec) }
  vrh=1-theta/mu
  logpr= (y-1)*log(theta+vrh*y)-theta-vrh*y+log(theta)-lgamma(y+1)
  nllk=-sum(logpr)
  nllk
}

# mu(x)= exp(b[0]+b^T x)), theta(x)=mu(x)*(1-vrh), 1-vrh=sqrt(1/(1+xi))
# param = b[0] b[1] ... b[nc] xi>=0, nc=ncol(xdat)
#    length(param)=ncol(xdat)+2
# y = vector of non-negative integers, length(y)=nrow(xdat)
# xdat = matrix of covariates (without intercept)
# Output: negative log-likelihood
gp1nllk=function(param,y,xdat)
{ np=length(param)
  xi=param[np]; 
  if(xi<0) return(1.e10)
  vrh=1-sqrt(1/(1+xi))
  bvec=param[2:(np-1)]
  b0=param[1]
  if(is.vector(xdat)) { mu=exp(b0+xdat*bvec) }
  else { mu=exp(b0+xdat%*%bvec) }
  theta=mu*(1-vrh)
  logpr= (y-1)*log(theta+vrh*y)-theta-vrh*y+log(theta)-lgamma(y+1)
  nllk=-sum(logpr)
  nllk
}

# mu(x)= exp(b[0]+b^T x)), theta(x)=mu(x)*(1-vrh), 1-vrh=sqrt(1/(1+xi))
# ub = upper bound value for pmf and cdf values
# param = b[0] b[1] ... b[nc] xi>=0, 
#    length(param)=length(x)+2
# x = vector of covariates
# Output: 
#  3-column matrix with 0:ub in column 1, pmf in column 2, cdf in column 3
gp1pmfcdf=function(ub,param,x)
{ np=length(param)
  xi=param[np]; vrh=1-sqrt(1/(1+xi))
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  theta=mu*(1-vrh)
  gpoispmfcdf(ub,theta,vrh)
}

# mu(x)= exp(b[0]+b^T x)), 1-vrh(x)=theta/mu(x), 1-vrh=sqrt(1/(1+xi))
# ub = upper bound value for pmf and cdf values
# param = b[0] b[1] ... b[nc] theta>=0, 
#    length(param)=length(x)+2
# x = vector of covariates
# Output: 
#  3-column matrix with 0:ub in column 1, pmf in column 2, cdf in column 3
gp2pmfcdf=function(ub,param,x)
{ np=length(param)
  theta=param[np]; 
  bvec=param[2:(np-1)]
  b0=param[1]
  mu=exp(b0+sum(bvec*x))
  vrh=1-theta/mu
  gpoispmfcdf(ub,theta,vrh)
}

#============================================================

# checks
#theta=1.5
#vrh=.4
#out=gpoispmfcdf(15,theta,vrh)
#print(out)
#pmf=dgpois(0:15,c(theta,vrh))
#print(pmf)
#cdf=cumsum(pmf)
#print(cdf)
#print(cdf-out[,3])
#cdfx=pgpois(15,c(theta,vrh))
#print(cdfx)
#pmf2=dgpois(4,cbind(c(2,2,2),c(0,.5,1)))
#print(pmf2)

#xdat=matrix(c(1:3,0,2,2),2,3,byrow=T)
#param=c(.1,.2,.3,.4,2)
#y=c(3,0)
#nllk1=gp1nllk(param,y,xdat)
#print(nllk1)
#print(gp1nllk(param,y[1],matrix(xdat[1,],1,3)))
#print(gp1nllk(param,y[2],matrix(xdat[2,],1,3)))
#tem1=gp1pmfcdf(5,param,xdat[1,])
#print(tem1)
#tem2=gp1pmfcdf(5,param,xdat[2,])
#print(tem2)
#print(log(tem1[4,2])+log(tem2[1,2])) # negative of nllk1

#nllk2=gp2nllk(param,y,xdat)
#print(nllk2)
#print(gp2nllk(param,y[1],matrix(xdat[1,],1,3)))
#print(gp2nllk(param,y[2],matrix(xdat[2,],1,3)))
#tem1=gp2pmfcdf(5,param,xdat[1,])
#print(tem1)
#tem2=gp2pmfcdf(5,param,xdat[2,])
#print(tem2)
#print(log(tem1[4,2])+log(tem2[1,2])) # negative of nllk2
