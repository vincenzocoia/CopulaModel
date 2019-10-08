# Archimedean copula with (normalized) integrated Mittag-Leffler LT:
# incomplete beta, integrated LT of BB1, 2-parameter family
# 2 parameters vth=1/ga>0 and de>1 

# numerically decreases in concordance as ga increases
# reparametrize vth=1/ga to get increasing in concordance 
# numerically increases in concordance as de increases

# For several functions below:
# s = positive real
# parameters ga>0, de>1 with ze>de, ga=ze-de

# Output: integrated beta LT, 
ibetalt= function(s,ga,de)
{ sde=s^(1/de)
  w=sde/(1+sde)
  1-pbeta(w,de,ga)
}

# con=beta(de,ze-de) in calling routine
# Output: first derivative integrated beta LT, 
ibetalt1= function(s,ga,de,con)
{ s1=s^(1/de)
  ze=de+ga
  -(1+s1)^(-ze) /(de*con)
}


# con=beta(de,ze-de) in calling routine
# Output: second derivative integrated beta LT, 
ibetalt2= function(s,ga,de,con)
{ s1=s^(1/de)
  ze=de+ga
  (s1/s)*(1+s1)^(-ze-1) * ze/(de*de*con)
}

# tt with value in (0,1)
# Output: inverse of integrated beta LT
ibetaltinv= function(tt,ga,de)
{ tem=1-qbeta(1-tt,de,ga)
  (1./tem-1)^de
}

# Archimedean copula with 2 parameters
# u,v = vectors or scalars in (0,1)
# cpar = (vth,de) where vth=1/ga>0, de>1 with ze>de>1, ga=ze-de>0
# Output: copula cdf
pimitlefA=function(u,v,cpar)
{ if(is.matrix(cpar)) { ga=1/cpar[,1]; de=cpar[,2] }
  else { ga=1/cpar[1]; de=cpar[2] }
  tem=ibetaltinv(u,ga,de)+ibetaltinv(v,ga,de)
  ibetalt(tem,ga,de)
}

# conditional cdf of Archimedean copula based on ibetalt
# v,u = vectors or scalars in (0,1)
# cpar = (vth,de) where vth=1/ga>0, de>1 with ze>de>1, ga=ze-de>0
# Output: conditional cdf
pcondimitlefA=function(v,u,cpar)
{ if(is.matrix(cpar)) { ga=1/cpar[,1]; de=cpar[,2] }
  else { ga=1/cpar[1]; de=cpar[2] }
  tem1=ibetaltinv(u,ga,de)
  tem2=ibetaltinv(v,ga,de)
  con=beta(de,ga)
  ibetalt1(tem1+tem2,ga,de,con)/ibetalt1(tem1,ga,de,con)
}

# density of Archimedean copula based on ibetalt
# u,v = vectors or scalars in (0,1)
# cpar = (vth,de) where vth=1/ga>0, de>1 with ze>de>1, ga=ze-de>0
# Output: copula density
dimitlefA=function(u,v,cpar)
{ if(is.matrix(cpar)) { ga=1/cpar[,1]; de=cpar[,2] }
  else { ga=1/cpar[1]; de=cpar[2] }
  tem1=ibetaltinv(u,ga,de)
  tem2=ibetaltinv(v,ga,de)
  con=beta(de,ga)
  ibetalt2(tem1+tem2,ga,de,con)/(ibetalt1(tem1,ga,de,con)*ibetalt1(tem2,ga,de,con))
}

# inverse of conditional cdf C_{2|1}^{-1}(p|u;cpar)
# vectorized version in C would be much faster
# p,u = vectors or scalars in (0,1)
# cpar = (vth,de) where vth=1/ga>0, de>1 with ze>de>1, ga=ze-de>0
# Output: quantile function
qcondimitlefA=function(p,u,cpar, eps=1.e-8,mxiter=30,iprint=F)
{ iter=0; diff=1.;
  v=u;
  while(iter<mxiter & max(abs(diff))>eps)
  { num=pcondimitlefA(v,u,cpar)-p; 
    den=dimitlefA(u,v,cpar);
    diff=num/den;
    v=v-diff;
    while(min(v)<0. | max(v)>1.) { diff=diff/2.; v=v+diff;}
    iter=iter+1;
    if(iprint) cat(iter, diff, v,"\n")
  }
  if(iter>mxiter) cat("*** did not converge\n")
  # might need to add bisection method when cpar implies tau>0.9
  v
}


# cpar = (vth,de) where vth=1/ga>0, de>1 with ze>de>1, ga=ze-de>0
# Output: Kendall's tau
imitlefA.cpar2tau=function(cpar)
{ if(is.matrix(cpar)) { ga=1/cpar[,1]; de=cpar[,2] }
  else { ga=1/cpar[1]; de=cpar[2] }
  con=beta(de,ga)
  con2=beta(2*de,2*ga)
  1-4*con2/de/con^2
}


#============================================================

# reflected/survival copula has skewness to upper tail
# cpar = (vth,de) where vth=1/ga>0, de>1 with ze>de>1, ga=ze-de>0
pimitlefAr=function(u,v,cpar)
{ u+v-1+pimitlefA(1-u,1-v,cpar) }

pcondimitlefAr=function(v,u,cpar)
{ 1-pcondimitlefA(1-v,1-u,cpar) }

dimitlefAr=function(u,v,cpar)
{ dimitlefA(1-u,1-v,cpar) }

