
# functions for trivariate mixture of max-id copulas, used for
#  Markov order2 time series
# param = (LTpar,cpar12) where
#  LTpar is parameter of a Laplace transform
#  cpar12 is parameter of a bivariate max-id copula

# positive stable LT and pfrk for H12,H32
# uu = 3-dimensional vector with values in (0,1)
# param = (LTpar,cpar12) 
# Output: trivariate cdf
pmxid3ps=function(uu,param)
{ LTpar=param[1]; cpar12=param[2]
  ss=.5*(-log(uu))^(LTpar)
  vv=exp(-ss)
  H12=pfrk(vv[1],vv[2],cpar12)
  H32=pfrk(vv[3],vv[2],cpar12)
  tem=-log(H12)-log(H32)+ss[1]+ss[3]
  exp(-tem^(1/LTpar))
}

# reflected pmxid3ps
pmxidr3ps=function(uu,param)
{ vv=1-uu; LTpar=param[1]; #cpar12=param[2]
  1-sum(vv)+pmxid2ps(vv[1],vv[2],param)
    pmxid2ps(vv[2],vv[3],param)+ pgum(vv[1],vv[3],LTpar) - pmxid3ps(vv,param)
}

# (1,2) bivariate margin of pmxid2ps
# argument u,v to follow pattern of other bivariate copulas
# 0<u,v<1 could be vectors
# param = (LTpar,cpar12) 
# Output: bivariate cdf
pmxid2ps=function(u,v,param)
{ LTpar=param[1]; cpar12=param[2]
  su=.5*(-log(u))^(LTpar)
  sv=.5*(-log(v))^(LTpar)
  eu=exp(-su); ev=exp(-sv)
  H12=pfrk(eu,ev,cpar12)
  tem=-log(H12)+su+sv
  exp(-tem^(1/LTpar))
}

# reflected pmxid2ps
pmxidr2ps=function(u,v,param)
{ u+v-1+pmxid2ps(1-u,1-v,param) }

# log series LT and pfrk for H12,H32
# uu = 3-dimensional vector with values in (0,1)
# param = (LTpar,cpar12) 
# Output: trivariate cdf
pmxid3ls=function(uu,param)
{ LTpar=param[1]; cpar12=param[2]
  ss=sqrt(1-exp(-uu*LTpar))
  sde=sqrt(1-exp(-LTpar))
  vv=ss/sde
  H12=pfrk(vv[1],vv[2],cpar12)
  H32=pfrk(vv[3],vv[2],cpar12)
  tem=1-H12*H32*ss[1]*ss[3]
  -(1/LTpar)*log(tem)
}

# (1,2) bivariate margin of pmxid2ls
# argument u,v to follow pattern of other bivariate copulas
# 0<u,v<1 could be vectors
# param = (LTpar,cpar12) 
# Output: bivariate cdf
pmxid2ls=function(u,v,param)
{ LTpar=param[1]; cpar12=param[2]
  su=sqrt(1-exp(-u*LTpar))
  sv=sqrt(1-exp(-v*LTpar))
  #ss=sqrt(1-exp(-uu*LTpar))
  sde=sqrt(1-exp(-LTpar))
  tu=su/sde; tv=sv/sde
  H12=pfrk(tu,tv,cpar12)
  tem=1-H12*su*sv
  -(1/LTpar)*log(tem)
}

#============================================================

# positive stable LT and pgum for H12,H32 (example of MM1 in mmdc)
# uu = 3-dimensional vector with values in (0,1)
# param = (LTpar,cpar12) 
# Output: trivariate cdf
ptmm1=function(uu,param)
{ LTpar=param[1]; cpar12=param[2]
  z=(-log(uu))
  ss=.5*(z^LTpar)
  tem12=(ss[1]^cpar12+ss[2]^cpar12)^(1/cpar12)
  tem32=(ss[3]^cpar12+ss[2]^cpar12)^(1/cpar12)
  tem=tem12+tem32+ss[1]+ss[3]
  exp(-tem^(1/LTpar))
}

# (1,2) bivariate margin of ptmm1
# 0<u,v<1 could be vectors
# param = (LTpar,cpar12) 
# Output: bivariate cdf
pbmm1=function(u,v,param)
{ LTpar=param[1]; cpar12=param[2]
  zu=(-log(u)); zv=(-log(v))
  su=.5*(zu^LTpar); sv=.5*(zv^LTpar)
  tem12=(su^cpar12+sv^cpar12)^(1/cpar12)
  tem=tem12+su+sv
  exp(-tem^(1/LTpar))
}

# conditionals C_{2|1} and C_{3|12} in order to check Kendall tau 
# and for simulation

# conditional of pbmm1 with respect to u
# 0<u,v<1 could be vectors
# param = (LTpar,cpar12) 
# Output: conditional (univariate) cdf
pcondbmm1=function(v,u,param)
{ LTpar=param[1]; cpar12=param[2]
  zu=(-log(u)); zv=(-log(v))
  su=.5*(zu^LTpar); sv=.5*(zv^LTpar)
  supow=su^cpar12
  sm=supow+sv^cpar12
  tem12=sm^(1/cpar12)
  tem=tem12+su+sv
  tempow=tem^(1/LTpar)
  cdf=exp(-tempow)
  deru= (tem12/sm)*supow/zu + su/zu
  ccdf=cdf*(tempow/tem)*deru/u
  ccdf
}

# density of pbmm1 
# 0<u,v<1 could be vectors
# param = (LTpar,cpar12) 
# Output: bivariate pdf
dbmm1=function(u,v,param)
{ LTpar=param[1]; cpar12=param[2]
  zu=(-log(u)); zv=(-log(v))
  su=.5*(zu^LTpar); sv=.5*(zv^LTpar)
  supow=su^cpar12; svpow=sv^cpar12
  sm=supow+svpow
  tem12=sm^(1/cpar12)
  tem=tem12+su+sv
  tempow=tem^(1/LTpar)
  cdf=exp(-tempow)
  deru= (tem12/sm)*supow/zu + su/zu
  derv= (tem12/sm)*svpow/zv + sv/zv
  pdf=cdf*( deru*derv * (tempow + LTpar-1) +
        tem * (tem12/sm^2) * LTpar*(cpar12-1)*supow*svpow/zu/zv )
  pdf=pdf *(tempow/tem^2) /u/v
}

#cpar=c(1.5,1.2)
#u=.3
#u=.8
#v=seq(.4,.9,.1)
#chkcopderiv(u,v,cpar,bcdf=pbmm1,pcond=pcondbmm1,bpdf=dbmm1,str="bivMM1",eps=1.e-5)

#============================================================

