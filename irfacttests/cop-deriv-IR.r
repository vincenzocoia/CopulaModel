################################################################################
#  various copula cdfs derivatives                         
#  written by Aristidis Nikoloulopoulos for Psychometrika publication
################################################################################

# v = vector or value in (0,1)
# u = vector or value in (0,1)
# rho = correlation parameter
# Output: C_{2|1}(v|u} for bivariate normal copula
pcondnor=function(v,u,rho)
{ val=pnorm((qnorm(v)-rho*qnorm(u))/sqrt(1-rho^2))
  val[ v<=0 | u<=0 | u>=1 ] = 0
  val[ v==1 ] = 1
  val
}


# v = vector or value in (0,1)
# u = vector or value in (0,1)
# rho = correlation parameter
# df = degree of freedom parameter (set dfdefault value before calling)
# Output: c(u,v) for Student t copula density
dtcop=function(u,v,rho,df=dfdefault)
{ tem1=cbind(u,v)
  tem2=qt(tem1,df)
  tem3=dt(tem2,df)
  x=tem2[,1]
  y=tem2[,2]
  z=tem3[,1]
  w=tem3[,2]
  x2=x*x
  y2=y*y
  r2=rho*rho
  val=df/2*(1+(x2+y2-2*rho*x*y)/(df*(1-r2)))^(-(df+2)/2)/(df*pi*z*w*sqrt(1-r2))
  val[ u<=0 | v<=0 ] = 0
  val[ u>=1 | v>=1 ] = 0
  val
}


# v = vector or value in (0,1)
# u = vector or value in (0,1)
# rho = correlation parameter
# Output: derivative of C_{2|1}(v|u} with respect to the correlation parameter
# for bivariate normal copula
pconddotbvncop=function(u,v,rho)
{ qnormv=qnorm(v)
  tem1=qnorm(u)-rho*qnormv
  tem=1-rho^2
  tem2=sqrt(tem)
  tem3=(-qnormv*tem2 + tem1*rho/tem2)/tem
  val=dnorm(tem1/tem2)*tem3
  val[ u<=0 | v<=0 ] = 0
  val[ u>=1 | v>=1 ] = 0
  val
}

# v = vector or value in (0,1)
# u = vector or value in (0,1)
# rho = correlation parameter
# df = degree of freedom parameter (set dfdefault value before calling)
# Output: derivative of C_{2|1}(v|u} with respect to the correlation parameter
# for Student t copula
pconddott=function(u,v,rho,df=dfdefault)
{ qtv=qt(v,df)
  tem1=qt(u,df)-rho*qtv
  tem=(1-rho^2)*(df+qtv^2)/(df+1)
  tem2=sqrt(tem)
  tem3=-rho*(df+qtv^2)/(df+1)/tem2
  tem4=(-qtv*tem2 - tem1*tem3)/tem
  val=dt(tem1/tem2,df+1)*tem4
  val[ u<=0 | v<=0 ] = 0
  val[ u>=1 | v>=1 ] = 0
  val
}

# v = vector or value in (0,1)
# u = vector or value in (0,1)
# cpar = copula parameter
# Output: derivative of C_{2|1}(v|u} with respect to the copula parameter
# for Gumbel copula
pconddotgum=function(u,v,cpar)
{ logu=-log(u)
  loglogu=log(logu)
  logv=-log(v)
  loglogv=log(logv)
  logva=logv^(cpar - 1)
  logva1=logv^cpar
  logua1=logu^cpar
  temp=logua1 + logva1
  ia=1/cpar
  ia1=ia-1
  ia2=ia^2
  iv=(1/v)
  exptemp=exp(-temp^ia)
  temp2=(logua1 * loglogu + logva1 * loglogv)
  temp3=(logva * (cpar * iv))
  temp4=temp^ia1
  val=exptemp * ((temp^(ia -2) * (ia1*temp2) - 
     temp4*(log((logua1 + logva1)) * ia2)) * (ia *  temp3) +  
     temp4 *(ia *(logva*loglogv * (cpar*iv) + logva*iv) - ia2*temp3)) - 
     exptemp * (temp4 * (ia*temp2) - temp^ia * (log(temp) *ia2)) * 
     (temp4 * (ia * temp3))
  val[ u<=0 | v<=0 ] = 0
  val[ u>=1 | v>=1 ] = 0
  val
}

# v = vector or value in (0,1)
# u = vector or value in (0,1)
# cpar = copula parameter
# Output: derivative of C_{2|1}(v|u} with respect to the copula parameter
# for survival Gumbel copula
pconddotgumr=function(u,v,cpar)
{ u=1-u
  v=1-v
  val=-pconddotgum(u,v,cpar)
  val[ u<=0 | v<=0 ] = 0
  val[ u>=1 | v>=1 ] = 0
  val
}

