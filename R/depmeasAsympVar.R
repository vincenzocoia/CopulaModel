# asymptotic variance of bivariate dependence measures tau, rhoS, beta
# for bivariate copula with parameter cpar

# cpar = copula parameter 
# pcop = bivariate copula cdf
# pcond12 = conditional cdf C_{1|2}, 
# pcond21 = conditional cdf C_{2|1}, 
#    above two are the same if C is permutation symmetric
# zero = 0 or something like 1.e-6
# tol = accuracy for 2-dimensional integration on [zero,1-zero]^2 
# Output: asymptotic variance of Kendall's tau
ktau.avar=function(cpar,pcop,pcond12,pcond21,zero=0,tol=1.e-5)
{ tauavarfn= function(u) 
  { u1=u[1]; u2=u[2];
    tem=2*pcop(u1,u2,cpar)+1-u1-u2
    tem=tem*(2*pcond12(u1,u2,cpar)-1)*pcond21(u2,u1,cpar)
  }
  lb=zero; ub=1-zero
  e123=adaptIntegrate(tauavarfn,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
  e123=1/3-2*e123$integ
  tau=ktau(cpar,icond=T,pcond12,pcond21,zero=0,tol=tol)
  avar=16*(e123-(1+tau)^2/4)
  avar
}

# cpar = copula parameter 
# pcop = bivariate copula cdf
# Output: asymptotic variance of Blomqvist's beta 
blomqvist.avar=function(cpar,pcop)
{ blom=4*pcop(0.5,0.5,cpar)-1
  avar=(1-blom^2)
  avar
}

# cpar = copula parameter 
# pcop = bivariate copula cdf
# pcond12 = conditional cdf C_{1|2}, 
# pcond21 = conditional cdf C_{2|1}, 
#    above two are the same if C is permutation symmetric
# nq = number of quadrature points for Gauss-Legendre in inner integral
# zero = 0 or something like 1.e-6
# tol = accuracy for 2-dimensional integration on [zero,1-zero]^2 
# Output: asymptotic variance of Spearman's rho 
rhoS.avar=function(cpar,pcop,pcond12,pcond21,nq=25,zero=0,tol=1.e-5)
{ gl=gausslegendre(nq)
  wl=gl$weights
  xl=gl$nodes
  rhavarfn= function(u) 
  { u1=u[1]; u2=u[2];
    g1=sum(wl*pcop(u1,xl,cpar))+.5-u1
    g2=sum(wl*pcop(xl,u2,cpar))+.5-u2
    g2p=sum(wl*pcond12(xl,u2,cpar))-1
    tem1=(u1+g1)^2
    tem2=u1*(u1*u2+g1+g2+u2*g2p) + (g1+g2)*g2p
    tem2=tem2*2*pcond21(u2,u1,cpar)
    tem1-tem2
  }
  lb=zero; ub=1-zero
  e123=adaptIntegrate(rhavarfn,lowerLimit=c(lb,lb),upperLimit=c(ub,ub),absError=tol)
  e123=e123$integ
  rh=rhoS(cpar,cop=pcop)
  avar=144*e123-9*(3+rh)^2
  avar
}

