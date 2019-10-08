
# GARCH filter function; this requires library(fGarch)

# logret is Nxd matrix of log returns, d=number of assets
# ar=T for GARCH(1,1)-AR(1)
# ar=F for GARCH(1,1)
# m1:m2 is range of cases for subset  n=nr=m2-m1+1 = number of cases in subset
# iprint=T to print GARCH estimates in intermediate calculations

# Output: 
#  filter= GARCH filtered data (nxd) 
#  uscore= empirical uniform scores (nxd) 
#  zscore= empirical normal scores (nxd) 
#  uscmodel= model-based uniform scores (nxd) 
#  zscmodel= model-based normal scores (nxd) 
#  sigmat= matrix of estimated volatilities (nxd) 
#  coef = 6xd with GARCH parameters mu, ar1, omega, alpha1, beta1, df=nu
#          if ar=T
#  coef = 5xd with GARCH parameters mu, omega, alpha1, beta1, df=nu
#          if ar=F

# empirical uniform and normal scores are based on (rank-.5)/n 
# model-based  uniform and normal scores are based on distribution of
#  the innovations

gfiltersubset=function(lgret,ar,m1,m2,iprint=F)
{ nr=m2-m1+1
  subdat=lgret[m1:m2,]
  d=ncol(subdat)
  filter=matrix(0,nr,d) 
  uscore=matrix(0,nr,d)
  zscore=matrix(0,nr,d)
  uscmodel=matrix(0,nr,d)
  zscmodel=matrix(0,nr,d)
  sigmat=matrix(0,nr,d)
  coef=matrix(0,6,d)
  if(!ar) coef = matrix(0,5,d);
  qn=qnorm(((1:nr)-0.5)/nr)
  for (j in 1:d)
  { if(ar) 
    { tem=garchFit(~arma(1,0)+garch(1,1), data=subdat[,j], cond.dist="std",
           trace=F,desc="") 
    }
    else
    { tem=garchFit(~garch(1,1),data=subdat[,j],cond.dist="std",trace=F,desc="") 
    }
    # use std = Student t innovations rather than normal
    if(iprint) { cat("variable ",j,"\n"); print(coef(tem)) }
    filter[,j]=tem@residuals/ tem@sigma.t;
    sigmat[,j]=tem@sigma.t;
    #print(tem@fit$params$params)
    if(ar) { coef[,j]=tem@fit$params$params[c(1:4,6,9)] }
    else { coef[,j]=tem@fit$params$params[c(1:3,5,8)] }
    jj=rank(filter[,j]);
    uscore[,j]=(jj-.5)/nr;
    zscore[,j]=qn[jj];
    uscmodel[,j]=pstd(filter[,j]); # pstd is in fGarch package
    zscmodel[,j]=qnorm(uscmodel[,j])
  }
  list(filter=filter,uscore=uscore,zscore=zscore,
    uscmodel=uscmodel,zscmodel=zscmodel,sigmat=sigmat,coef=coef)
}

