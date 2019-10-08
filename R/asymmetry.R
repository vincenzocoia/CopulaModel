# Functions for (a) measures for bivariate reflection and permutation asymmetry 
# and (b) generation from bivariate copulas with extreme asymmetry

# A reference for rbreflasym() and the measures of reflection asymmetry is:
# Rosco, J-F and Joe H (2014). Statistical Papers, v 54, 709-726,

# A bivariate copula with extreme tail asymmetry: upper tail dependence =1
#  lower tail order is infinite
# n = sample size, 
# param = parameter in (0,1)
# Output: nx2 matrix with U(0,1) margins
rbreflasym=function(n,param=.25)
{ u=runif(n)
  v=u
  a=1-param
  i1=(u<a)
  v[i1]=a-u[i1]
  cbind(u,v)
}

# Rotation of copula in rbreflasym so that U-V is positively skewed
# n = sample size, 
# param = parameter in (0,1)
# Output: nx2 matrix with U(0,1) margins
rbpermasym=function(n,param=.25)
{ u=runif(n)
  v=1-u
  a=1-param
  i1=(u<a)
  v[i1]=u[i1]+param
  cbind(u,v)
}

# n = sample size
# cpar = (de,p)
#   de>1, 0<=p<1, p=0  for no skewing
#   p=0 for independence
# Output: nx2 matrix with U(0,1) margins
rbasymgum1=function(n,cpar)
{ de=cpar[1]; p=cpar[2]
  uu=rgum(n,de)
  x1=-log(uu[,1])
  x2=rexp(n)
  x=pmin(x1/p,x2/(1-p))
  uu[,1]=exp(-x)
  uu
}

# bivariate Gumbel with one skew parameter 0<p<1
# stochastic representation is in function rbasymgum1
# 0<u<1, 0<v<1: could be vectors
# cpar = (de,p)
#   de>1, 0<=p<1, p=0  for no skewing
#   p=0 for independence
# Output: bivariate cdf
pbasymgum1=function(u,v,cpar)
{ de=cpar[1]; p=cpar[2]
  if(de>50)  # assume M-O
  { tem=pmin(u,v^p); cdf=tem*(v^(1-p)) }
  else
  { x=-log(u); y=-log(v)
    tem=(y^de+(x*p)^de)^(1/de) + (1-p)*x
    cdf=exp(-tem)
  }
  cdf
}

# bivariate Marshall-Olkin skew 1-parameter subfamily, 
#  limit as Gumbel parameter de->oo
# n = sample size, 
# p = parameter in (0,1)
#   p=0 for independence, p=1 for comonotonicity
# Output: nx2 matrix with U(0,1) margins
rbMO1=function(n,p)
{ y=rexp(n)
  x2=rexp(n)
  x=pmin(y/p,x2/(1-p))
  xy=cbind(x,y)
  exp(-xy)
}

# skewness for reflection asymmetry based on third moment
# max value is 27/256
# uu = nx2 matrix, each column with U(0,1) margins
# Output: skewness coefficient and estimated SE
skewrefl=function(uu)
{ udiff=uu[,1]+uu[,2]-1
  sk=udiff^3
  mn=mean(sk); vv=mean(sk^2)-mn^2
  se=sqrt(vv/length(udiff))
  c(mn,se)
}

# skewness for permutation asymmetry based on third moment
# max value is 27/256
# uu = nx2 matrix, each column with U(0,1) margins
# Output: skewness coefficient and estimated SE
skewperm=function(uu)
{ udiff=uu[,1]-uu[,2]
  sk=udiff^3
  mn=mean(sk); vv=mean(sk^2)-mn^2
  se=sqrt(vv/length(udiff))
  c(mn,se)
}

# skewness for permutation asymmetry based on extreme quantiles 
# uu = nx2 matrix, each column with U(0,1) margins
# p = quantile value between 0.5 and 1
# nrep = number of replications for bootstrap SE
# Output: skewness coefficient and estimated SE
qskewperm=function(uu,p=.05,nrep=100)
{ udiff=uu[,1]-uu[,2]
  tem=quantile(udiff,c(p,.5,1-p),type=8)
  zeta=(tem[3]-2*tem[2]+tem[1])/(tem[3]-tem[1])
  # bootstrap is better than jackknife
  n=length(udiff); 
  zrep=rep(0,nrep)
  for(ij in 1:nrep)
  { ii=sample(n,replace=T); utem=udiff[ii]
    tem1=quantile(utem,c(p,.5,1-p),type=8)
    zeta1=(tem1[3]-2*tem1[2]+tem1[1])/(tem1[3]-tem1[1])
    zrep[ij]=zeta1
  }
  out=c(zeta,mean(zrep),sqrt(sum((zrep-zeta)^2)/(nrep-1)))
  names(out)=c("qskewperm","qavg","SE")
  out
}

# skewness for reflection asymmetry based on extreme quantiles
# uu = nx2 matrix, each column with U(0,1) margins
# p = quantile value between 0.5 and 1
# nrep = number of replications for bootstrap SE
# Output: skewness coefficient and estimated SE
qskewrefl=function(uu,p=.05,nrep=100)
{ udiff=uu[,1]+uu[,2]-1
  tem=quantile(udiff,c(p,.5,1-p),type=8)
  zeta=(tem[3]-2*tem[2]+tem[1])/(tem[3]-tem[1])
  n=length(udiff); 
  zrep=rep(0,nrep)
  for(ij in 1:nrep)
  { ii=sample(n,replace=T); utem=udiff[ii]
    tem1=quantile(utem,c(p,.5,1-p),type=8)
    zeta1=(tem1[3]-2*tem1[2]+tem1[1])/(tem1[3]-tem1[1])
    zrep[ij]=zeta1
  }
  out=c(zeta,mean(zrep),sqrt(sum((zrep-zeta)^2)/(nrep-1)))
  names(out)=c("qskewrefl","qavg","SE")
  out
}

