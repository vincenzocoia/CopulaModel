# bivariate copula of trivariate 1-factor copula model given variable 3,
#  and conditional Spearman rho
# Examples uses Gumbel and Frank

library(CopulaModel)
# gldefault can be set globally

# trivariate margin of 1-factor copula, 
# uu = 2-column matrix of values in (0,1)
# u3 = scalar in (0,1)
# param = parameter of copula
# pcondcop = function for conditional cdf of bivariate linking copula
# gl = Gauss-Legendre object (nodes and weights)
# Output: C_{123}(uu[,1],uu[,2],u3)
ptfact1cop= function(uu,u3,param,pcondcop,gl=gldefault)
{ #gl <- gausslegendre(nq);
  wl <- gl$weights; xl <- gl$nodes; #nq=length(wl)
  if(is.matrix(param)) { par.1=param[1,]; par.2=param[2,]; par.3=param[3,] } 
  else { par.1=param[1]; par.2=param[2]; par.3=param[3] }          
  if(is.vector(uu)) uu=matrix(uu,nrow=1)
  m=nrow(uu) 
  tcdf = rep(0,m); 
  for(i in 1:m) 
  { tem=pcondcop(uu[i,1],xl,par.1)*pcondcop(uu[i,2],xl,par.2)*pcondcop(u3,xl,par.3)
    tcdf[i]=sum(wl*tem)
  }
  tcdf
}

# conditional 12|3 for ptfact1cop, 
# uu = 2-column matrix of values in (0,1)
# u3 = scalar in (0,1)
# param = parameter of copula
# pcondcop = function for conditional cdf of bivariate linking copula
# dcop = function for density of bivariate linking copula
# gl = Gauss-Legendre object (nodes and weights)
# Output: bivariate conditional cdf C_{12|3}(uu|u3)
#   this is not a copula when u3 is fixed
pcondtfact12g3= function(uu,u3,param,pcondcop,dcop,gl=gldefault)
{ #gl <- gausslegendre(nq);
  wl <- gl$weights; xl <- gl$nodes; #nq=length(wl)
  if(is.matrix(param)) { par.1=param[1,]; par.2=param[2,]; par.3=param[3,] } 
  else { par.1=param[1]; par.2=param[2]; par.3=param[3] }          
  if(is.vector(uu)) uu=matrix(uu,nrow=1)
  m=nrow(uu) 
  bcdf = rep(0,m); 
  for(i in 1:m) 
  { tem=dcop(u3,xl,par.3)*pcondcop(uu[i,1],xl,par.1)*pcondcop(uu[i,2],xl,par.2)
    # denominator is 1 
    bcdf[i]=sum(wl*tem)
  }
  bcdf
}

# conditional cdf 1|3 for ptfact1cop, 
# uu = vector of values in (0,1)
# u3 = scalar in (0,1)
# param = parameter of copula
# pcondcop = function for conditional cdf of bivariate linking copula
# dcop = function for density of bivariate linking copula
# gl = Gauss-Legendre object (nodes and weights)
# Output: univariate conditional cdf C_{1|3}(uu|u3)
pcondtfact1g3 = function(uu,u3,param,pcondcop,dcop,gl=gldefault)
{ #gl <- gausslegendre(nq);
  wl <- gl$weights; xl <- gl$nodes; #nq=length(wl)
  if(is.matrix(param)) { par.1=param[1,]; par.3=param[2,]; } 
  else { par.1=param[1]; par.3=param[2]; }          
  m=length(uu) 
  ucdf = rep(0,m); 
  for(i in 1:m) 
  { tem=dcop(u3,xl,par.3)*pcondcop(uu[i],xl,par.1)
    ucdf[i]=sum(wl*tem)
  }
  ucdf
}

# (1,3) bivariate marginal density for ptfact1cop, uu is vector, u3 is scalar 
# uu = vector of values in (0,1)
# u3 = scalar in (0,1)
# param = parameter of copula
# dcop = function for density of bivariate linking copula
# gl = Gauss-Legendre object (nodes and weights)
# Output: univariate conditional pdf c_{1|3}(uu|u3)
dtfact1g3 = function(uu,u3,param,dcop,gl=gldefault)
{ #gl <- gausslegendre(nq);
  wl <- gl$weights; xl <- gl$nodes; 
  if(is.matrix(param)) { par.1=param[1,]; par.3=param[2,]; } 
  else { par.1=param[1]; par.3=param[2]; }          
  m=length(uu) 
  bpdf = rep(0,m); 
  for(i in 1:m) 
  { tem=dcop(u3,xl,par.3)*dcop(uu[i],xl,par.1)
    bpdf[i]=sum(wl*tem)
  }
  bpdf
}

# conditional quantile cdf 1|3 for ptfact1cop, pp is vector, u3 is scalar 
# pp = vector of values in (0,1)
# u3 = scalar in (0,1)
# param = parameter of copula
# xcdf = vector of cdf values
# yuu = vector of matching quantiles
#    xcdf,yy are used for interpolation
# qderiv = vector of derivatives from pcderiv
# Output: univariate conditional quantile C_{1|3}^{-1}(uu|u3)
qcondtfact1g3 = function(pp,u3,param,xcdf,yuu,qderiv)
{ if(is.matrix(param)) { par.1=param[1,]; par.3=param[2,]; } 
  else { par.1=param[1]; par.3=param[2]; }          
  #deriv=pcderiv(xcdf,yuu)
  #print(cbind(deriv,qderiv))
  qq=pcinterpolate(xcdf,yuu,qderiv,pp)
  if(is.matrix(qq)) { return(qq[,1]) } else { return(qq[1]) }
}

# conditional 12|3 for ptfact1cop; copula of C_{12|3} with u3 fixed
# uu = 2-column matrix of values in (0,1)
# u3 = scalar in (0,1)
# param = parameter of copula
# pcondcop = function for conditional cdf of bivariate linking copula
# dcop = function for density of bivariate linking copula
# gl = Gauss-Legendre object (nodes and weights)
# Output: copula of C_{12|3} with u3 fixed
ptfact12g3cop = function(u1,u2,param,u3,pcondcop,dcop,gl=gldefault,
   yuu,mcdf1,qderiv1,mcdf2,qderiv2)
{ #gl <- gausslegendre(nq);
  wl <- gl$weights; xl <- gl$nodes; #nq=length(wl)
  if(is.matrix(param)) { par.1=param[1,]; par.2=param[2,]; par.3=param[3,] } 
  else { par.1=param[1]; par.2=param[2]; par.3=param[3] }          
  x1=pcinterpolate(mcdf1,yuu,qderiv1,u1)  # F_{1|3}^{-1}(u1) in element 1
  x2=pcinterpolate(mcdf2,yuu,qderiv2,u2)  # F_{2|3}^{-1}(u2)
  tem=dcop(u3,xl,par.3)*pcondcop(x1[1],xl,par.1)*pcondcop(x2[1],xl,par.2)
  bcdf=sum(wl*tem)
  bcdf
}

#============================================================

# given u3 and param, get univariate margin C_{1|3}, C_{2|3} 
# and their derivatives, then use pcinterpolate to get quantiles for copula
gldefault=gausslegendre(25)

u3=.4
uvec=seq(.05,.95,.05)
param=c(2,1.5,1.5)
mcdf1=pcondtfact1g3(uvec,u3,param[c(1,3)],pcondgum,dgum,gl=gldefault)
mcdf1=c(0,mcdf1,1)
mcdf2=pcondtfact1g3(uvec,u3,param[c(2,3)],pcondgum,dgum,gl=gldefault)
mcdf2=c(0,mcdf2,1)
uvec=c(0.000001,uvec,0.999999)
#cat("density of (1,3) and (2,3) margins\n")
mpdf1=dtfact1g3(uvec,u3,param[c(1,3)],dgum,gl=gldefault)
mpdf2=dtfact1g3(uvec,u3,param[c(2,3)],dgum,gl=gldefault)

# univariate quantile function
cat("\nconditional cdf of (1|3) and (2|3) margins\n")
pp=seq(.1,.9,.1)
# what is derivative at endpoints?
qderiv1=1/mpdf1
qq1=qcondtfact1g3(pp,u3,param[c(1,3)],mcdf1,uvec,qderiv1)
print(qq1)
qq2=qcondtfact1g3(pp,u3,param[c(2,3)],mcdf2,uvec,1/mpdf2)
print(qq2)

cat("\ncdf of (1,2|3) margin\n")
for(u1 in seq(.1,.9,.1))
{ copcdf=ptfact12g3cop(u1,u1+.05,param,u3,pcondgum,dgum,gl=gldefault,
   uvec,mcdf1,1/mpdf1,mcdf2,1/mpdf2)
  print(copcdf)
}

# copula of (1,2|3) with U3=u3 given
ptemcop=function(u1,u2,param)
{ ptfact12g3cop(u1,u2,param,u3,pcondgum,dgum,gl=gldefault,
   uvec,mcdf1,1/mpdf1,mcdf2,1/mpdf2)
}

spear=rhoS(param,cop=ptemcop,zero=.0001)
#print(spear)
# 0.2412805 for u3=0.4

# set up Spearman rho as a function of u3
# write this into a wrapper??
param=c(2,1.5,1.5)
cat("\nGumbel : ", param,"\n")
cat("Spearman's rho conditional on different u3\n")
for(u3 in seq(.1,.9,.1))
{ uvec=seq(.05,.95,.05)
  uvec=c(0.0001,uvec,0.9999) # problems with .00001 and .000001 for some u3
  mcdf1=pcondtfact1g3(uvec,u3,param[c(1,3)],pcondgum,dgum,gl=gldefault)
  mcdf2=pcondtfact1g3(uvec,u3,param[c(2,3)],pcondgum,dgum,gl=gldefault)
  mpdf1=dtfact1g3(uvec,u3,param[c(1,3)],dgum,gl=gldefault)
  mpdf2=dtfact1g3(uvec,u3,param[c(2,3)],dgum,gl=gldefault)
  ptemcop=function(u1,u2,param)
  { ptfact12g3cop(u1,u2,param,u3,pcondgum,dgum,gl=gldefault,
     uvec,mcdf1,1/mpdf1,mcdf2,1/mpdf2)
  }
  spear=rhoS(param,cop=ptemcop,zero=.0001)
  print(c(u3,spear))
}

#frk.b2cpar(c(.5,.6,.7))
param=c(4.875023,6.562560,9.101941)
cat("\nFrank : ", param,"\n")
cat("Spearman's rho conditional on different u3\n")
for(u3 in seq(.1,.9,.1))
{ uvec=seq(.05,.95,.05)
  uvec=c(0.000001,uvec,0.999999) # 
  mcdf1=pcondtfact1g3(uvec,u3,param[c(1,3)],pcondfrk,dfrk,gl=gldefault)
  mcdf2=pcondtfact1g3(uvec,u3,param[c(2,3)],pcondfrk,dfrk,gl=gldefault)
  mpdf1=dtfact1g3(uvec,u3,param[c(1,3)],dfrk,gl=gldefault)
  mpdf2=dtfact1g3(uvec,u3,param[c(2,3)],dfrk,gl=gldefault)
  ptemcop=function(u1,u2,param)
  { ptfact12g3cop(u1,u2,param,u3,pcondfrk,dfrk,gl=gldefault,
     uvec,mcdf1,1/mpdf1,mcdf2,1/mpdf2)
  }
  spear=rhoS(param,cop=ptemcop,zero=.000001)
  print(c(u3,spear))
}

# convert into wrappers and try more parameter vectors to see patterns
# also increase nq
