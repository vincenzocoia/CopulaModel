# functions for IG(m=mu,varsigma=lambda)
# eta=sqrt(lambda)=sqrt(varsigma) = convolution parameter
# varsigma abbreviated to vsi below

# Seshadri (1993) The Inverse Gaussian Distribution. Clarendon Press.
#  IG0(mu,vsi) with m=mu>0, vsi=varsigma>0 means
#  f(x;m,vsi)= {sqrt{vsi}\over sqrt{2\pi x^3}}
#   exp{ -{vsi\over 2m^2} {(x-m)^2\over x}},  x>0.
# In this parametrization m=ze*eta is the mean and vsi=eta^2
# (since then m^2/vsi=ze^2).
# On page 81, Theorem 2.8:
# F(x;m,vsi) is the cdf; then
# F(x;m,vsi)=F(x/m;1,vsi/m) and
#  F(x;m,\vsi) = \cases{
#   0.5 G(a) + 0.5 exp{2\vsi m^{-1}}  G(a+4\vsi m^{-1}) ,  0<= x<= m 
#   1-0.5 G(a) + 0.5 exp{2\vsi m^{-1}}  G(a+4\vsi m^{-1}), #   m< x< \oo  }
# where a=\vsi(x-m)^2/[m^2x] and 
#  G(w)=\int_w^\oo (2\pi t)^{-1/2} exp{-t/2} dt = 2*Phi(-w^{1/2})
#   =2*Phibar(w^{1/2});
# the latter via the transform to z^2 since G(w)=P(W>w)
# where W~ chi^2_1 or Gamma(1/2,rate=2).


# x = positive value or vector
# mu = mean parameter m 
# vsi = second parameter 
# (mu,vsi) <-> (eta,ze);  eta=convolution parameter 
#   ze is like a scale parameter that can be set to 1 for the copula
#   eta=sqrt(vsi), ze=mu/eta; mu=ze*eta, vsi=eta^2
# Invariance property is pIG(x,mu,vsi)=pIG(x/mu,1,vsi/m)
# Output: cdf of inverse Gaussian
pIG=function(x,mu,vsi)
{ a=vsi*(x-mu)^2/mu^2/x
  Ga=2*pnorm(-sqrt(a))
  Ga4=2*pnorm(-sqrt(a+4*vsi/mu))
  tem=.5*Ga
  ii=(x>mu)
  tem[ii]=1-tem[ii]
  cdf=tem+.5*Ga4*exp(2*vsi/mu)
  cdf
}

# quantile function based on qinvgauss in statmod package 
# p = value in (0,1): maybe can be a vector
# mu = mean parameter m 
# vsi = second parameter 
#  mu. vsi are scalars (or vectors with same length as p)
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# mxstep = bound on step size for Newton-Raphson iterations
# iprint = print flag for iterations
# Output: quantile function of inverse Gaussian
qIG=function(p,mu,vsi,mxiter=10,eps=1.e-6,mxstep=5,iprint=F) 
{ thi=vsi/mu; sqthi=sqrt(thi)
  u=qnorm(p)
  r1=1 + u/sqthi + u^2/(2*thi) + u^3/(8*thi*sqthi) # statmod, can be negative
  r1[r1<0]=0.1
  x=r1 # starting point
  diff=1; iter=0
  while(iter<mxiter & max(abs(diff))>eps)
  { h=pIG(x,1,thi)-p
    diff=h/dIG(x,1,thi)
    x=x-diff
    iter=iter+1
    while(max(abs(diff))>mxstep | min(x)<=0) { diff=diff/2; x=x+diff }
    if(iprint) cat("iter=", iter, " x=", x, ", diff=",h, "\n")
  }
  x*mu
}


# x = positive value or vector
# mu = mean parameter = m (Seshadri's book)
# vsi = second parameter = lambda (Seshadri's book)
# Output: quantile function of inverse Gaussian
dIG=function(x,mu,vsi)
{ tem=vsi*(x-mu)^2/mu^2/x
  pdf=exp(-.5*tem)*sqrt(vsi/(2*pi*x))/x
  pdf
}

# transformation function for IG
# x = positive value
gforIG= function(x,mu,vsi)
{ sqrt(vsi)*(mu-x)/(mu*sqrt(x)) }

# inverse of gforIG, 
# y = positive value
ginvforIG= function(y,mu,vsi)
{ tem=2*mu*vsi+(mu*y)^2
  discr=mu^2*y^2 +4*mu*vsi
  #print(discr)
  #root1=(tem-sqrt(discr))/(2*vsi)
  #root2=(tem+sqrt(discr))/(2*vsi)
  (tem-mu*y*sqrt(discr))/(2*vsi)
}

# n = simulation sample size
# mu = mean parameter = m (Seshadri's book)
# vsi = second parameter = lambda (Seshadri's book)
#  eta=sqrt(vsi) is convolution parameter 
# Output: random sample from inverse Gaussian
rIG=function(n,mu,vsi)
{ z=rnorm(n)
  t2=rt(n,2)
  eta=sqrt(mu/(2*vsi))
  y=ifelse(eta*z>t2,z,-z)
  ginvforIG(y,mu,vsi)
}

# checks for gforIG
#m=3
#lm=2
#x=1:10
#y=gforIG(x,m,lm)
#tem=ginvforIG(y,m,lm)
#print(y)
#print(x)
#print(tem)

# check random sample
#set.seed(123)
#cat("m=", m, " lm=", lm, "\n")
#n=10000
#x=rIG(n,m,lm)
#print(summary(x))
#print(var(x))
#cat("theoretical values mean=", m, " var=", m^3/lm,"\n")

#============================================================

# inverse Gaussian convolution 1-factor model 
# see Section 4.28 of dmwc

# random IG convolution 1-factor model
# let marginal parameters be etavec=thetavec+theta0
# th=sqrt(lambda)=sqrt(varsigma) = convolution parameter
# zeta=ze is fixed non-convolution parameter = m/th, m=mu=mean=ze*th
# LT is exp{(lm/mu)*[1-sqrt(1+2*mu^2*s/lm)]}
#   = exp((th/ze)*[1-sqrt(1+2*ze^2*s)]}
# IG(th0,ze) convolute IG(thj,ze) leads to IG(th0+thj,ze)
# can take ze=1 for use with copula

# n = sample size
# th0 = scalar for convolution parameter of the shared/common component
# thvec = vector of convolution parameters of individual components, length d
# zeta = scale parameter that doesn't appear in the copula
rIGconv=function(n,th0,thvec,ze=1)
{ #set.seed(seed)
  d=length(thvec)
  y=matrix(0,n,d)
  vsivec=thvec^2; vsi0=th0^2  
  mvec=ze*thvec; m0=ze*th0    # mean
  for(i in 1:n)
  { yindiv=rIG(d,mvec,vsivec)
    y[i,]=yindiv+rIG(1,m0,vsi0)
    # marginal distribution has etavec=th0+thvec, ze in transformed parameters
    #  and vsi=etavec^2, muvec=ze*etavec  in original parameters 
  }
  y
}

# cdf of X1=Z0+Z1, ... Xd=Z0+Zd using adaptive integration
# Z0~IG(theta0), Z1~IG(theta1), .. Zd~IG(thetad) independent
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
#   mu=theta, vsi=theta^2 so that mean=var=theta
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: cdf of inverse Gaussian 1-factor model
pmIGfact=function(xvec,th0,thvec,zero=0)
{ d=length(xvec)
  Finteg= function(z) 
  { intg=dIG(z,th0,th0^2)
    for(j in 1:d) intg=intg*pIG(xvec[j]-z,thvec[j],thvec[j]^2) 
    intg
  }
  out=integrate(Finteg,zero,min(xvec)-zero)
  out$value
}


# cdf of X1=Z0+Z1, ... Xd=Z0+Zd using Gauss-Legendre quadrature
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# Output: cdf of inverse Gaussian 1-factor model
pmIGfact.gl=function(xvec,th0,thvec,gl)
{ d=length(xvec)
  xmin=min(xvec)
  z=gl$nodes*xmin
  w=gl$weights
  intg=dIG(z,th0,th0^2)
  for(j in 1:d) intg=intg*pIG(xvec[j]-z,thvec[j],thvec[j]^2) 
  cdf=xmin*sum(w*intg)
  cdf
}

# pdf of X1=Z0+Z1, ... Xd=Z0+Zd using adaptive integration
# Z0~IG(theta0), Z1~IG(theta1), .. Zd~IG(thetad) independent
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
#   mu=theta, vsi=theta^2 so that mean=var=theta
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: pdf of inverse Gaussian 1-factor model
dmIGfact=function(xvec,th0,thvec,zero=0)
{ d=length(xvec)
  Finteg= function(z) 
  { intg=dIG(z,th0,th0^2)
    for(j in 1:d) intg=intg*dIG(xvec[j]-z,thvec[j],thvec[j]^2) 
    intg
  }
  out=integrate(Finteg,zero,min(xvec)-zero)
  out$value
}

# pdf of X1=Z0+Z1, ... Xd=Z0+Zd using Gauss-Legendre quadrature
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# Output: pdf of inverse Gaussian 1-factor model
dmIGfact.gl=function(xvec,th0,thvec,gl)
{ d=length(xvec)
  xmin=min(xvec)
  z=gl$nodes*xmin
  w=gl$weights
  intg=dIG(z,th0,th0^2)
  for(j in 1:d) intg=intg*dIG(xvec[j]-z,thvec[j],thvec[j]^2) 
  cdf=xmin*sum(w*intg)
  cdf
}

# bivariate marginal density of inverse Gaussian factor model
# x1, x2 = positive values
# th0 = theta0 = positive parameter of the common component
# th1, th2 = parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# Output: pdf of bivariate margin of inverse Gaussian 1-factor model
dbIGfact=function(x1,x2,th0,th1,th2,zero=0)
{ Finteg= function(z) 
  { intg=dIG(z,th0,th0^2)
    intg=intg*dIG(x1-z,th1,th1^2) *dIG(x2-z,th2,th2^2)
    intg
  }
  out=integrate(Finteg,zero,min(x1,x2)-zero)
  out$value
}

#============================================================

# for copula need qIG 

# uvec = d-vector with values in (0,1)
# param = vector of dimension d+1 with th0, thvec
#   see parametrizations in rIGconv() above
#   marginal IG distribution has etavec=th0+thvec, ze in transformed parameters
#   and vsi=etavec^2, muvec=ze*etavec  in original parameters 
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# Output: copula cdf for pmIGfact, 
pmIGfcop.gl=function(uvec,param,gl)
{ th0=param[1]; thvec=param[-1]; 
  # assume ze=1
  vsivec=thvec^2; vsi0=th0^2  
  mvec=thvec; m0=th0    # mean
  etavec=th0+thvec
  xvec=qIG(uvec,etavec,etavec^2)  
  pmIGfact.gl(xvec,th0,thvec,gl)
}

# copula pdf for pmIGfact, 
# uvec = d-vector with values in (0,1)
# param = vector of dimension d+1 with th0, thvec
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# Output: copula pdf for pmIGfact, 
dmIGfcop.gl=function(uvec,param,gl)
{ th0=param[1]; thvec=param[-1]; 
  # assume ze=1
  vsivec=thvec^2; vsi0=th0^2  
  mvec=thvec; m0=th0    # mean
  etavec=th0+thvec
  xvec=qIG(uvec,etavec,etavec^2)
  numer=dmIGfact.gl(xvec,th0,thvec,gl)
  denom=prod(dIG(xvec,etavec,etavec^2))
  numer/denom
}

