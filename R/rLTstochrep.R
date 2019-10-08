
# LT of Archimedean and mix-maxid copulas
# Simulation of the random variables with the given LT
# Simulation of multivariate Archimedean based on stochastic representation

# n = simulation sample size in all functions

#============================================================

# 1-parameter LT and corresponding Archimedean copulas

# generate positive stable random variable (Chambers, Mallows and Stuck 1986)
# n = simulation sample size 
# alp = 1/cpar in (0,1)
# Output: random sample of size n with positive stable distribution
rpostable=function(n,alp)
{ v=runif(n,0,pi)
  a1=1/(1-alp); rat=alp/(1-alp)
  h= sin((1-alp)*v)* (sin(alp*v))^rat / (sin(v))^a1
  w=rexp(n)
  x=(h/w)^(1/rat)
  x
}

# generate Sibuya random variable (Devroye 1993)
# n = simulation sample size 
# alp = 1/cpar in (0,1)
# Output: random sample of size n with Sibuya distribution
rsibuya=function(n,alp)
{ mu=rexp(n)*rgamma(n,1-alp)/rgamma(n,alp)
  x=1+rpois(n,mu)
  x
}

# generate log series random variable (Kemp 1981)
# Kemp's algorithm1 for log series rv
# this works for cpar>37 with cpar as parameter and not alp
# n = simulation sample size 
# cpar = copula parameter of Frank with positive dependence (cpar>0)
# Output: random sample of size n with logarithmic series distribution
rlogseries=function(n,cpar)
{ #cpar=-log(1-alp)
  alp=1-exp(-cpar)  # rounds to 1 for cpar>37.4
  x=rep(1,n)
  v=runif(n)
  ii=(v<alp)
  u=runif(n)
  tem=exp(-cpar*u[ii])
  x[ii]=floor(1+log(v[ii])/log(1-tem))
  # when cpar exceeds 37, this can become -Inf
  # because 1-exp(-cpar*u[ii]) rounds off as 1
  isinf=is.infinite(x) 
  # occurs where u large enough, e.g., >0.8 for cpar=45
  x[isinf]=floor(1-log(v[isinf])/tem[isinf])   # first order approx
  #print(cbind(u[isinf],v[isinf],x[isinf]))
  x
}

#============================================================

# Version in standalone R
# multivariate frank = log series Archimedean
# n = simulation sample size 
# d = dimension
# cpar = copula parameter >0
# Output: nxd matrix with U(0,1) margins
rmfrk0=function(n,d,cpar)
{ cpar0=exp(-cpar)
  alp=1-cpar0
  r=rlogseries(n,cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r
    # r infinite => tem=0 => next line is infinite ?
    tem2= -log(1-alp*exp(-tem))/cpar
    isinf=is.infinite(tem2)
    #print(cbind(tem2[isinf],tem[isinf]))
    # Inf occurs when tem[isinf] close enough to 0
    tem2[isinf]= -log(cpar0+alp*tem[isinf])/cpar
    uu[,j]=tem2
  }
  uu
}

# Version with link to C
# n = simulation sample size 
# d = dimension
# cpar = copula parameter >0
# icheck =T to print out means and correlation
# Output: nxd matrix with U(0,1) margins
rmfrk=function(n,d,cpar,icheck=F)
{ uu=rep(0,n*d)
  out= .C("rmfrk", as.integer(n), as.integer(d), as.double(cpar), 
           uu=as.double(uu))
  uu=matrix(out$uu,n,d)
  if(icheck) 
  { cat("\ncopula parameter=",cpar,"\n")
    s1=mean(uu[,1]); s2=mean(uu[,2])
    cat("average u1 =", s1,"\n")
    cat("average u2 =", s2,"\n")
    s12=mean(uu[,1]*uu[,2])
    cat("correlation = ", 12*(s12-s1*s2),"\n")
  }
  uu
}


# multivariate MTCJ = Archimedean with gamma LT
# n = simulation sample size 
# d = dimension
# cpar = copula parameter >0
# Output: nxd matrix with U(0,1) margins
rmmtcj=function(n,d,cpar) 
{ r=rgamma(n,1/cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r
    uu[,j]=(1+tem)^(-1/cpar)
  }
  uu
}


# multivariate Gumbel = Archimedean with positive stable LT
# n = simulation sample size 
# d = dimension
# cpar = copula parameter >1
# Output: nxd matrix with U(0,1) margins
rmgum=function(n,d,cpar)
{ alp=1/cpar
  r=rpostable(n,alp)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r
    uu[,j]=exp(-tem^alp)
  }
  uu
}

# multivariate Joe = Archimedean with Sibuya LT
# n = simulation sample size 
# d = dimension
# cpar = copula parameter >1
# Output: nxd matrix with U(0,1) margins
rmjoe=function(n,d,cpar)
{ alp=1/cpar
  r=rsibuya(n,alp)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=runif(n)^(1/r)
    uu[,j]=1-(1-tem)^alp
  }
  uu
}

#============================================================
# 2-parameter LTs and BB copulas

# gamma stopped positive stable or Mittag-Leffler LT 
# n = simulation sample size 
# param = (th,de) with th>0, de>=1
# Output: random sample of size n with variable with Mittag-Leffler LT
rmitlef=function(n,param)
{ th=param[1]; de=param[2]; alp=1/de
  tt=rgamma(n,1/th)
  # next few lines same as rpostable(n,alp)
  v=runif(n,0,pi)
  a1=1/(1-alp); rat=alp/(1-alp)
  h= sin((1-alp)*v)* (sin(alp*v))^rat / (sin(v))^a1
  w=rexp(n)
  x=(h/w)^(1/rat)
  # next line same as rpostable(n,alp) * (tt^de)
  x*(tt^de)
}
  
# LTF gamma stopped gamma
# n = simulation sample size 
# param = (th,de) with th>0, de0
# Output: random sample of size n with variable with gamma stopped gamma LT
rgammaSgamma=function(n,param)
{ th=param[1]; de=param[2]; 
  tt=rgamma(n,1/th)
  z=rgamma(n,tt/de)
  z
}

# LTG positive stable stopped gamma
# n = simulation sample size 
# param = (th,de) with th>1, de>0
# Output: random sample of size n with variable with postable stopped gamma LT
rpostableSgamma=function(n,param)
{ th=param[1]; de=param[2]; 
  tt=rpostable(n,1/th)
  z=rgamma(n,tt/de)
  z
}

# LTH Sibuya stopped positive stable
# n = simulation sample size 
# param = (th,de) with th>1, de>=1
# Output: random sample of size n with variable with Sibuya stopped postable LT
rsibuyaSpostable=function(n,param)
{ th=param[1]; de=param[2];  alp=1/de
  tt=rsibuya(n,1/th)
  z=rpostable(n,alp) * (tt^de)
  z
}

# LTI Sibuya stopped gamma
# n = simulation sample size 
# param = (th,de) with th>1, de>0
# Output: random sample of size n with variable with Sibuya stopped gamma LT
rsibuyaSgamma=function(n,param)
{ th=param[1]; de=param[2];  alp=1/de
  tt=rsibuya(n,1/th)
  z=rgamma(n,tt/de)
  z
}

#============================================================

# n = simulation sample size 
# d = dimension
# cpar = (th,de), th>0, de>=1
# Output: multivariate BB1 random sample 
rmbb1=function(n,d,cpar)
{ th1=1/cpar[1]; # th1=1/th
  de1=1/cpar[2]
  r=rmitlef(n,cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r # next LTpsi(tem)
    uu[,j]=(1+tem^de1)^(-th1)
  }
  uu
}

# n = simulation sample size 
# d = dimension
# cpar = (th,de), th>0, de>0
# Output: multivariate BB2 random sample 
rmbb2=function(n,d,cpar)
{ th1=1/cpar[1]; # th1=1/th
  de1=1/cpar[2]
  r=rgammaSgamma(n,cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r # next LTpsi(tem)
    uu[,j]=(1+de1*log(1+tem))^(-th1)
  }
  uu
}

# n = simulation sample size 
# d = dimension
# cpar = (th,de), th>1, de>0
# Output: multivariate BB3 random sample 
rmbb3=function(n,d,cpar)
{ th1=1/cpar[1]; # th1=1/th
  de1=1/cpar[2]
  r=rpostableSgamma(n,cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r # next LTpsi(tem)
    uu[,j]=exp(-(de1*log(1+tem))^th1)
  }
  uu
}


# n = simulation sample size 
# d = dimension
# cpar = (th,de), th>1, de>=1
# Output: multivariate BB6 random sample 
rmbb6=function(n,d,cpar)
{ th1=1/cpar[1]; # th1=1/th
  de1=1/cpar[2]
  r=rsibuyaSpostable(n,cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r # next LTpsi(tem)
    uu[,j]=1-(1-exp(-tem^de1))^th1
  }
  uu
}

# n = simulation sample size 
# d = dimension
# cpar = (th,de), th>1, de>0
# Output: multivariate BB7 random sample 
rmbb7=function(n,d,cpar)
{ th1=1/cpar[1]; # th1=1/th
  de1=1/cpar[2]
  r=rsibuyaSgamma(n,cpar)
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=rexp(n)/r # next LTpsi(tem)
    uu[,j]=1-(1-(1+tem)^(-de1))^th1
  }
  uu
}

#============================================================

# n = simulation sample size 
# d = dimension
# cpar = (th,ppi), th>0, 0<ppi<1
# Output: multivariate shifted NB Archimedean random sample 
rmbb10=function(n,d,cpar)
{ th1=1/cpar[1]; ppi=cpar[2]
  r=rnbinom(n,size=th1,prob=1-ppi)+th1
  uu=matrix(0,n,d)
  for(j in 1:d) 
  { tem=runif(n)^(1/r)
    tem2=(1-ppi)*tem/(1-ppi*tem)
    uu[,j]=tem2^th1
  }
  uu
}

