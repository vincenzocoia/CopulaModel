# cdf for bivariate t, vectorized inputs

# R interface for bivariate t cdf, nu=df is integer-valued only
# link to fortran source code in the mvtnorm R package source
# If z1,z2,rho, nu are compatible, the function will probably work,
# z1 = (vector of) variable 1
# z2 = (vector of) variable 2
# rho = correlation in (-1,1)
# nu = degree of freedom
# Output
#   bivariate cdfs at (z1,z2,rho,nu)
#   -1 for any illegal or NA/NaN inputs
#   default output is a vector
pbvt=function(z1,z2,param,icheck=F)
{
  #if(!is.loaded("pbvt"))  dyn.load("./pbvt.so")
  if(is.matrix(param)) { rho=param[,1]; nu=param[,2] }
  else { rho=param[1]; nu=param[2] }
  n1=length(z1)
  n2=length(z2)
  nr=length(rho)
  ndf=length(nu)
  # check for matrices etc
  imat=is.matrix(z1)
  if(imat) { ir=nrow(z1); ic=ncol(z1) }
  nn=max(n1,n2,nr,ndf)
  #if(imat) print(c(ir,ic))
  if(n1<nn & n1>1) return(NA)
  if(n2<nn & n2>1) return(NA)
  if(nr<nn & nr>1) return(NA)
  if(ndf<nn & ndf>1) return(NA)
  if(nn>1)
  { if(n1==1) z1=rep(z1,nn)
    if(n2==1) z2=rep(z2,nn)
    if(nr==1) rho=rep(rho,nn)
    if(ndf==1) nu=rep(floor(nu),nn)
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
    nu[is.nan(nu)]=-2
    nu[is.nan(z1)]=-2; nu[is.nan(z2)]=-2; 
    rho[is.na(z1)]=-2; z1[is.na(z1)]=0
    rho[is.na(z2)]=-2; z2[is.na(z2)]=0
    rho[is.na(rho)]=-2
    nu[is.na(nu)]=-2
    nu[is.na(z1)]=-2; nu[is.na(z2)]=-2; 
    #print(z1)
    #print(z2)
    #print(rho)
  }
  out= .Fortran("pbvt", as.integer(nn), as.double(z1), as.double(z2),
         as.double(rho), as.integer(nu), bcdf=as.double(rep(0,nn)) )
  #print(out$bcdf)
  if(imat) { bcdf=matrix(out$bcdf,ir,ic) } else bcdf=out$bcdf
  bcdf
}

