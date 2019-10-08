# functions for gamma convolution factor model

# cdf of X1=Z0+Z1, ... Xd=Z0+Zd using adaptive integration
# Z0~Gamma(theta0), Z1~Gamma(theta1), .. Zd~Gamma(thetad) independent
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: cdf of gamma 1-factor model
pmgamfact=function(xvec,th0,thvec,zero=0)
{ d=length(xvec)
  Finteg= function(z) 
  { intg=dgamma(z,th0)
    for(j in 1:d) intg=intg*pgamma(xvec[j]-z,thvec[j]) 
    intg
  }
  out=integrate(Finteg,zero,min(xvec)-zero)
  out$value
}

# pdf of X1=Z0+Z1, ... Xd=Z0+Zd using adaptive integration
# Z0~Gamma(theta0), Z1~Gamma(theta1) Zd~Gamma(thetad) independent
# If some parameters <1, the density could be infinite
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: pdf of gamma 1-factor model (where it is finite)
dmgamfact=function(xvec,th0,thvec,zero=0)
{ d=length(xvec)
  finteg= function(z) 
  { intg=dgamma(z,th0)
    for(j in 1:d) intg=intg*dgamma(xvec[j]-z,thvec[j]) 
    intg
  }
  out=integrate(finteg,zero,min(xvec)-zero)
  out$value
}

# pdf of X1=Z0+Z1, ... Xd=Z0+Zd using Gauss-Legendre quadrature
# xvec = d-vector of positive values
# th0 = theta0 = positive parameter of the common component
# thvec = d-vector of parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# Output: pdf of gamma 1-factor model (where it is finite)
dmgamfact.gl=function(xvec,th0,thvec,gl)
{ d=length(xvec)
  xmin=min(xvec)
  z=gl$nodes*xmin
  w=gl$weights
  intg=dgamma(z,th0)
  for(j in 1:d) intg=intg*dgamma(xvec[j]-z,thvec[j]) 
  pdf=xmin*sum(w*intg)
  pdf
}


# bivariate cdf of X1=Z0+Z1, X2=Z0+Z2
# Z0~Gamma(theta0), Z1~Gamma(theta1) Z2~Gamma(theta2) independent
# x1, x2 = positive values
# th0 = theta0 = positive parameter of the common component
# th1, th2 = parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: cdf of bivariate margin of gamma 1-factor model
pbgamfact=function(x1,x2,th0,th1,th2,zero=0)
{ Finteg= function(z) pgamma(x1-z,th1) * pgamma(x2-z,th2) * dgamma(z,th0)
  out=integrate(Finteg,zero,min(x1,x2)-zero)
  out$value
}

# bivariate pdf of X1=Z0+Z1, X2=Z0+Z2 using adaptive integration
# Z0~Gamma(theta0), Z1~Gamma(theta1) Z2~Gamma(theta2) independent
# x1, x2 = positive values
# th0 = theta0 = positive parameter of the common component
# th1, th2 = parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: pdf of bivariate margin of gamma 1-factor model
dbgamfact=function(x1,x2,th0,th1,th2,zero=0)
{ finteg= function(z) dgamma(x1-z,th1) * dgamma(x2-z,th2) * dgamma(z,th0)
  out=integrate(finteg,zero,min(x1,x2)-zero)
  out$value
}

# bivariate pdf of X1=Z0+Z1, X2=Z0+Z2 using Gauss-Legendre
# Z0~Gamma(theta0), Z1~Gamma(theta1) Z2~Gamma(theta2) independent
# x1, x2 = positive values
# th0 = theta0 = positive parameter of the common component
# th1, th2 = parameters for the individual components
# gl = Gauss-Legendre quadrature object with $nodes and $weights 
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: pdf of bivariate margin of gamma 1-factor model
dbgamfact.gl=function(x1,x2,th0,th1,th2,gl)
{ xmin=min(x1,x2)
  z=gl$nodes*xmin
  w=gl$weights
  intg=dgamma(z,th0)
  intg=intg*dgamma(x1-z,th1)*dgamma(x2-z,th2)
  pdf=xmin*sum(w*intg)
  pdf
}

# 0<u<1, 0<v<1
# param = (th0,th1,th2)
# Output: copula cdf for pbgamfact, 
pbgamfcop=function(u,v,param)
{ th0=param[1]; th1=param[2]; th2=param[3];
  x1=qgamma(u,th0+th1)
  x2=qgamma(v,th0+th2)
  pbgamfact(x1,x2,th0,th1,th2)
}

# 0<u<1, 0<v<1
# param = (th0,th1,th2)
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: copula pdf for pbgamfact
dbgamfcop=function(u,v,param,zero=0)
{ th0=param[1]; th1=param[2]; th2=param[3];
  x1=qgamma(u,th0+th1)
  x2=qgamma(v,th0+th2)
  numer=dbgamfact(x1,x2,th0,th1,th2,zero)
  denom=dgamma(x1,th0+th1)*dgamma(x2,th0+th2)
  numer/denom
}

# specify gldefault as global variable if needed 
#   (for some functions with dcop as an argument),
#  for example, gldefault=gausslegendre(25)
# 0<u<1, 0<v<1
# param = (th0,th1,th2)
# gl = Gauss-Legendre object with components $nodes and $weights
# Output: copula pdf for pbgamfact with Gauss-Legendre
dbgamfcop.gl=function(u,v,param,gl=gldefault)
{ th0=param[1]; th1=param[2]; th2=param[3];
  x1=qgamma(u,th0+th1)
  x2=qgamma(v,th0+th2)
  numer=dbgamfact.gl(x1,x2,th0,th1,th2,gl)
  denom=dgamma(x1,th0+th1)*dgamma(x2,th0+th2)
  numer/denom
}

# 0<u<1, 0<v<1
# param = (th0,th1,th2)
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: copula conditional cdf of pbgamfact given variable 1: C_{2|1}(v|u)
pcondbgamfcop21=function(v,u,param,zero=0)
{ th0=param[1]; th1=param[2]; th2=param[3];
  x1=qgamma(u,th0+th1)
  x2=qgamma(v,th0+th2)
  gamfact1= function(x1,x2,th0,th1,th2)
  { finteg= function(z) dgamma(x1-z,th1) * pgamma(x2-z,th2) * dgamma(z,th0)
    out=integrate(finteg,zero,min(x1,x2)-zero)
    out$value
  }
  numer=gamfact1(x1,x2,th0,th1,th2)
  denom=dgamma(x1,th0+th1) 
  numer/denom
}

# 0<u<1, 0<v<1
# param = (th0,th1,th2)
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: copula conditional cdf of pbgamfact given variable 2: C_{1|2}(u|v)
pcondbgamfcop12=function(u,v,param,zero=0)
{ th0=param[1]; th1=param[2]; th2=param[3];
  x1=qgamma(u,th0+th1)
  x2=qgamma(v,th0+th2)
  gamfact2= function(x1,x2,th0,th1,th2)
  { finteg= function(z) pgamma(x1-z,th1) * dgamma(x2-z,th2) * dgamma(z,th0)
    out=integrate(finteg,zero,min(x1,x2)-zero)
    out$value
  }
  numer=gamfact2(x1,x2,th0,th1,th2)
  denom=dgamma(x2,th0+th2)
  numer/denom
}

# uvec = d-vector with values in (0,1)
# param = vector of dimension d+1 with th0, thvec
# Output: copula cdf for pmgamfact
pmgamfcop=function(uvec,param)
{ th0=param[1]; thvec=param[-1]; 
  xvec=qgamma(uvec,th0+thvec)
  pmgamfact(xvec,th0,thvec)
}

# uvec = d-vector with values in (0,1)
# param = vector of dimension d+1 with th0, thvec
# zero = 0 or something like 0.0001 as tolerance for integration
# Output: copula pdf for pmgamfact
dmgamfcop=function(uvec,param,zero=0)
{ th0=param[1]; thvec=param[-1]; 
  xvec=qgamma(uvec,th0+thvec)
  numer=dmgamfact(xvec,th0,thvec,zero)
  denom=prod(dgamma(xvec,th0+thvec))
  numer/denom
}

# copula pdf for pmgamfact with Gauss-Legendre quadrature
# uvec = d-vector with values in (0,1)
# param = vector of dimension d+1 with th0, thvec
# gl = Gauss-Legendre object with components $nodes and $weights
# Output: copula pdf for pmgamfact
dmgamfcop.gl=function(uvec,param,gl)
{ th0=param[1]; thvec=param[-1]; 
  xvec=qgamma(uvec,th0+thvec)
  numer=dmgamfact.gl(xvec,th0,thvec,gl)
  denom=prod(dgamma(xvec,th0+thvec))
  numer/denom
}



# random gamma convolution 1-factor model
# n = sample size
# th0 = scalar for shape parameter of the shared/common component
# thvec = vector of shape parameters of individual components, length d
rgammaconv=function(n,th0,thvec)
{ d=length(thvec)
  y=matrix(0,n,d)
  for(i in 1:n)
  { yindiv=rgamma(d,thvec,1)
    y[i,]=yindiv+rgamma(1,th0,1)
  }
  y
}

