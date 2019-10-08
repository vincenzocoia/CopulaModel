# cdf for bivariate normal, vectorized inputs

# R interface,
# link to C code which is translation of Donnelly (1973), Comm ACM.
# If z1,z2,rho are compatible, the function will probably work,
# z1 = (vector of) variable 1
# z2 = (vector of) variable 2
# rho = correlation in (-1,1)
# Output
#   bivariate cdfs at (z1,z2,rho)
#   -1 for any illegal or NA/NaN inputs
#   default output is a vector
pbnorm=function(z1,z2,rho,icheck=F)
{ 
  #if(!is.loaded("pbnorm"))  dyn.load("./pbnorm.so")
  n1=length(z1)
  n2=length(z2)
  nr=length(rho)
  # check for matrices etc
  imat=is.matrix(z1)
  if(imat) { ir=nrow(z1); ic=ncol(z1) }
  nn=max(n1,n2,nr)
  #if(imat) print(c(ir,ic))
  if(n1<nn & n1>1) return(NA)
  if(n2<nn & n2>1) return(NA)
  if(nr<nn & nr>1) return(NA)
  if(nn>1)
  { if(n1==1) z1=rep(z1,nn)
    if(n2==1) z2=rep(z2,nn)
    if(nr==1) rho=rep(rho,nn)
  }
  if(icheck)
  { # handling infinites and NaN
    z1[is.infinite(z1) & z1>0]=20
    z2[is.infinite(z2) & z2>0]=20
    z1[is.infinite(z1) & z1<0]=-20
    z2[is.infinite(z2) & z2<0]=-20
    rho[is.nan(z1)]=-2; z1[is.nan(z1)]=0
    rho[is.nan(z2)]=-2; z2[is.nan(z2)]=0
    rho[is.nan(rho)]=-2
    rho[is.na(z1)]=-2; z1[is.na(z1)]=0
    rho[is.na(z2)]=-2; z2[is.na(z2)]=0
    rho[is.na(rho)]=-2
  }
  out= .C("pbnorm", as.integer(nn), as.double(z1), as.double(z2),
         as.double(rho), bcdf=as.double(rep(0,nn)) )
  #print(out$bcdf)
  #if(imat) { bcdf=matrix(out$bcdf,ir,ic) } else { bcdf=out$bcdf }
  bcdf=out$bcdf
  if(imat) bcdf=matrix(bcdf,ir,ic)
  bcdf
}

