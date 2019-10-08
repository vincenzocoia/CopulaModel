# bivariate asymmetric Gumbel (or bilogistic)
# 0<u<1, 0<v<1 for all functions
# cpar = copula parameter = (ze,eta) with 
# ze, eta each in (0,1)

# B(w)=A(w,1-w) ; A is the exponent of the bivariate exponential distribution
# B(w;ze,eta)=w z^{1-ze}+(1-w)(1-z)^{1-eta},  0< w< 1,
# z=z(w;ze,eta) is the root of the equation 
#  (1-ze) w (1-z)^eta - (1-eta) (1-w) z^ze = 0
# ww = vector of values in (0,1)
# cpar = (ze,eta) in (0,1)^2
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for intermediate results
# Output: B and its first/second order derivatives at ww values
Basymgum=function(ww,cpar,mxiter=30,eps=1.e-7,iprint=F)
{ ze=cpar[1]; eta=cpar[2]
  ze1=1-ze; eta1=1-eta
  # solve equation 
  w=ww; w1=1-w
  iter=0; dif=1.;
  z=1./(1.+(w1/w)^(2./(ze+eta)))  # starting guess based on ze=eta
  if(iprint) cat("\nw=", w, "\n")
  while(iter<mxiter & max(abs(dif))>eps)
  { za=z^(ze); z1=1.-z; zb=z1^(eta);
    g= w*ze1*zb- w1*eta1*za; 
    gp= -w*eta*ze1*zb/z1 -w1*ze*eta1*za/z;
    dif=g/gp; iter=iter+1; z=z-dif;
    while(min(z)<=0. | max(z)>=1.) {dif=dif/2.; z=z+dif;}
    if(iprint) cat(iter,w,z,dif,"\n");
  }
  if(iter>=mxiter) cat("did not converge\n")
  if(iprint) cat(w,z,"\n");
  bfn=w*z/za+w1*(z1/zb)
  bw=z/za-z1/zb
  bz=w*ze1/za-w1*eta1/zb
  zp=-(ze1*zb+eta1*za)/gp
  bder= bw+bz*zp
  zpp=-(-ze1*eta*zb/z1+eta1*ze*za/z)*zp/gp
  bwz=ze1/za+eta1/zb
  bder2=bwz*zp+bz*zpp
  list(Bfn=bfn,Bder=bder,Bder2=bder2)
}

# copula cdf for asymmetric Gumbel
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (0,1)^2
# Output: cdf
pasymgum=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Aval=ss*Basymgum(w,cpar)$Bfn
  cdf=exp(-Aval)
  cdf
}

# C_{2|1} for asymmetric Gumbel
# v = scalar or vector of values in (0,1)
# u = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (0,1)^2
# Output: conditional cdf
pcondasymgum21=function(v,u,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Bout=Basymgum(w,cpar)
  Aval=ss*Bout$Bfn
  cdf=exp(-Aval)
  Aderx=Bout$Bfn+(1-w)*Bout$Bder
  ccdf=cdf*Aderx/u
  ccdf
}

# C_{1|2} for asymmetric Gumbel
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (0,1)^2
# Output: conditional cdf
pcondasymgum12=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Bout=Basymgum(w,cpar)
  Aval=ss*Bout$Bfn
  cdf=exp(-Aval)
  Adery=Bout$Bfn-w*Bout$Bder
  ccdf=cdf*Adery/v
  ccdf
}

# copula density for asymmetric Gumbel
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (0,1)^2
# Output: copula density
dasymgum=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Bout=Basymgum(w,cpar)
  Aval=ss*Bout$Bfn
  cdf=exp(-Aval)
  Aderx=Bout$Bfn+(1-w)*Bout$Bder
  Adery=Bout$Bfn-w*Bout$Bder
  # need second deriv of B at w
  Aderxy=-w*(1-w)*Bout$Bder2/ss
  pdf=cdf*(Aderx*Adery-Aderxy)/u/v
  pdf
}

# Kendall's tau for asymmetric Gumbel 
# cpar = copula parameter (ze,eta) in (0,1)^2
asymgum.cpar2tau=function(cpar)
{ tauasym= function(w)
  { Bout=Basymgum(w,cpar)
    B=Bout$Bfn
    Bp=Bout$Bder
    w1=1-w
    ((2*w-1)*Bp*B + w*w1*Bp*Bp)/(B*B)
  }
  tem=integrate(tauasym,0.0000001,.9999999, rel.tol=1e-06)
  tau=tem$value
  tau
}

# Spearman's rho for asymmetric Gumbel 
# cpar = copula parameter (ze,eta) in (0,1)^2
asymgum.cpar2rhoS=function(cpar)
{ spasym= function(w)
  { Bout=Basymgum(w,cpar)
    B=Bout$Bfn
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(spasym,0.0000001,.9999999, rel.tol=1e-06)
  rho=12*tem$value-3
  rho
}

