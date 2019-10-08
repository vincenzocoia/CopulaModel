# bivariate asymmetric Galambos (or negative bilogistic)
# 0<u<1, 0<v<1 for all functions
# cpar = copula parameter = (ze,eta) with 
# ze, eta each in (-oo,0), that is, negative

# B(w)=A(w,1-w) ; A is the exponent of the bivariate exponential distribution
# B(w;ze,eta)=1-w u^{1-ze}-(1-w)(1-u)^{1-eta}
# u=u(w;ze,eta) is the root of the equation 
#  (1-ze) w u^{-ze}-(1-eta)(1-w)(1-u)^{-eta}=0
# ww = vector of values in (0,1)
# cpar = (ze,eta) in (-oo,0)^2
# mxiter = maximum number of iterations
# eps = tolerance for convergence
# iprint = print flag for intermediate results
# Output: B and its first/second order derivatives at ww values
Basymgal=function(ww,cpar,mxiter=30,eps=1.e-7,iprint=F)
{ ze=cpar[1]; eta=cpar[2]
  ze1=1-ze; eta1=1-eta
  # solve equation 
  w=ww; w1=1-w
  iter=0; dif=1.;
  z=(w1/w)^(-2./(ze+eta))  # t in notes
  z=z/(1.+z)  # starting guess based on ze=eta
  if(iprint) cat("\nw=", w, "\n")
  while(iter<mxiter & max(abs(dif))>eps)
  {  za=z^(-ze); z1=1.-z; zb=z1^(-eta);
     g= w*ze1*za- w1*eta1*zb; 
     gp= -w*ze*ze1*za/z -w1*eta*eta1*zb/z1;
     dif=g/gp; iter=iter+1; z=z-dif;
     while(min(z)<=0. | max(z)>=1.) {dif=dif/2.; z=z+dif;}
     if(iprint) cat(iter,w,z,dif,"\n");
  }
  if(iter>=mxiter) cat("did not converge\n")
  if(iprint) cat(w,z,"\n");
  bfn=1-w*(z^ze1)-(1-w)*(z1^eta1)
  bw=-z*za+z1*zb
  bz=-w*ze1*za+w1*eta1*zb
  zp=-(ze1*za+eta1*zb)/gp
  bder= bw+bz*zp
  zpp=-(-ze1*ze*za/z-eta1*eta*zb/z1)*zp/gp
  bwz=-ze1*za-eta1*zb
  bder2=bwz*zp+bz*zpp
  list(Bfn=bfn,Bder=bder,Bder2=bder2)
}

# copula cdf for asymmetric Galambos
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (-oo,0)^2
# Output: cdf
pasymgal=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Aval=ss*Basymgal(w,cpar)$Bfn
  cdf=exp(-Aval)
  cdf
}

# C_{2|1} for asymmetric Galambos
# v = scalar or vector of values in (0,1)
# u = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (-oo,0)^2
# Output: conditional cdf
pcondasymgal21=function(v,u,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Bout=Basymgal(w,cpar)
  Aval=ss*Bout$Bfn
  cdf=exp(-Aval)
  Aderx=Bout$Bfn+(1-w)*Bout$Bder
  ccdf=cdf*Aderx/u
  ccdf
}

# C_{1|2} for asymmetric Galambos
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (-oo,0)^2
# Output: conditional cdf
pcondasymgal12=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Bout=Basymgal(w,cpar)
  Aval=ss*Bout$Bfn
  cdf=exp(-Aval)
  Adery=Bout$Bfn-w*Bout$Bder
  ccdf=cdf*Adery/v
  ccdf
}

# copula density for asymmetric Galambos
# u = scalar or vector of values in (0,1)
# v = scalar or vector of values in (0,1)
# cpar = copula parameter (ze,eta) in (-oo,0)^2
# Output: copula density
dasymgal=function(u,v,cpar)
{ x=-log(u)
  y=-log(v)
  ss=x+y; w=x/ss
  Bout=Basymgal(w,cpar)
  Aval=ss*Bout$Bfn
  cdf=exp(-Aval)
  Aderx=Bout$Bfn+(1-w)*Bout$Bder
  Adery=Bout$Bfn-w*Bout$Bder
  # need second deriv of B at w
  Aderxy=-w*(1-w)*Bout$Bder2/ss
  pdf=cdf*(Aderx*Adery-Aderxy)/u/v
  pdf
}

# Kendall's tau for asymmetric Galambos
# cpar = copula parameter (ze,eta) in (-oo,0)^2
asymgal.cpar2tau=function(cpar)
{ tauasym= function(w)
  { Bout=Basymgal(w,cpar)
    B=Bout$Bfn
    Bp=Bout$Bder
    w1=1-w
    ((2*w-1)*Bp*B + w*w1*Bp*Bp)/(B*B)
  }
  tem=integrate(tauasym,0.0000001,.9999999, rel.tol=1e-06)
  tau=tem$value
  tau
}

# Spearman's rho for asymmetric Galambos
# cpar = copula parameter (ze,eta) in (-oo,0)^2
asymgal.cpar2rhoS=function(cpar)
{ spasym= function(w)
  { Bout=Basymgal(w,cpar)
    B=Bout$Bfn
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(spasym,0.0000001,.9999999, rel.tol=1e-06)
  rho=12*tem$value-3
  rho
}

