# random generation from copula GARCH with factor model: links to C

# copula GARCH with 1-factor model
# n = simulation sample size
# garchpar = 6xd matrix rows are mu,ar1,omega,alpha1,beta1,nu 
# cpar = copula parameter vector
# sigma0 = d-vector with starting values for conditional SDs (GARCH output) 
# copcode = copula code for 1-factor copula
#  Currently copcode options are
#  bivariate Gaussian: 1
#  bivariate t: 2
#  Gumbel: 3
#  reflected Gumbel: -3
#  Frank: 5
#  BB1: 9
# Output: dxn logret matrix and nx1 portfolio return vector
#  To get nxd logret matrix: 
#   transpose in this code, or make small changes in the C code
rgarch1fact=function(n,garchpar,cpar,sigma0,copcode)
{ d=ncol(garchpar)
  lgret=matrix(0,d,n)
  portfret=rep(0,n)
  out= .C("rgarch1fact", as.integer(d), as.integer(n), as.double(garchpar[1,]),
    as.double(garchpar[2,]), as.double(garchpar[3,]), as.double(garchpar[4,]),
    as.double(garchpar[5,]), as.double(garchpar[6,]), as.double(sigma0),
    as.double(cpar), as.integer(copcode), 
    lgret=as.double(lgret), portfret=as.double(portfret) )
  list(logret=matrix(out$lgret,d,n),portfret=out$portfret)
}

# copula GARCH with 2-factor model
# n = simulation sample size
# garchpar = 6xd matrix rows are mu,ar1,omega,alpha1,beta1,nu 
# cpar = copula parameter vector
# sigma0 = d-vector with starting values for conditional SDs (GARCH output) 
# copcode = copula code for 1-factor copula
#  Currently copcode options are
#  bivariate Gaussian: 1
#  bivariate t: 2
#  Gumbel: 3
#  reflected Gumbel: -3
#  Frank: 5
#  BB1/Frank (BB1 for factor 1 and Frank for factor 2): 9
# Output: dxn logret matrix and nx1 portfolio return vector
#  To get nxd logret matrix: 
#   transpose in this code, or make small changes in the C code
rgarch2fact=function(n,garchpar,cpar,sigma0,copcode)
{ d=ncol(garchpar)
  lgret=matrix(0,d,n)
  portfret=rep(0,n)
  out= .C("rgarch2fact", as.integer(d), as.integer(n), as.double(garchpar[1,]),
    as.double(garchpar[2,]), as.double(garchpar[3,]), as.double(garchpar[4,]),
    as.double(garchpar[5,]), as.double(garchpar[6,]), as.double(sigma0),
    as.double(cpar), as.integer(copcode), 
    lgret=as.double(lgret), portfret=as.double(portfret) )
  list(logret=matrix(out$lgret,d,n),portfret=out$portfret)
}

# copula GARCH with bi-factor model
# n = simulation sample size
# grsize = vector of group sizes for m groups with sum(grsize)=d 
# garchpar = 6xd matrix rows are mu,ar1,omega,alpha1,beta1,nu 
# cpar = copula parameter vector
# sigma0 = d-vector with starting values for conditional SDs (GARCH output) 
# copcode = copula code for 1-factor copula
#  Currently copcode options are
#  bivariate Gaussian: 1
#  bivariate t: 2
#  Gumbel: 3
#  reflected Gumbel: -3
#  Frank: 5
#  BB1/Frank (BB1 for global factor and Frank for group factor): 9
# Output: dxn logret matrix and nx1 portfolio return vector
#  To get nxd logret matrix: 
#   transpose in this code, or make small changes in the C code
rgarchbifact=function(n,grsize,garchpar,cpar,sigma0,copcode)
{ d=ncol(garchpar)
  m=length(grsize)  # sum(m) should equal d
  lgret=matrix(0,d,n)
  portfret=rep(0,n)
  out= .C("rgarchbifact", as.integer(m), as.integer(grsize), as.integer(n), 
    as.double(garchpar[1,]), as.double(garchpar[2,]), as.double(garchpar[3,]), 
    as.double(garchpar[4,]), as.double(garchpar[5,]), as.double(garchpar[6,]), 
    as.double(sigma0), as.double(cpar), as.integer(copcode), 
    lgret=as.double(lgret), portfret=as.double(portfret) )
  list(logret=matrix(out$lgret,d,n),portfret=out$portfret)
}

# copula GARCH with nested-factor model
# n = simulation sample size
# grsize = vector of group sizes for m groups with sum(grsize)=d 
# garchpar = 6xd matrix rows are mu,ar1,omega,alpha1,beta1,nu 
# cpar = copula parameter vector
# sigma0 = d-vector with starting values for conditional SDs (GARCH output) 
# copcode = copula code for 1-factor copula
#  Currently copcode options are
#  bivariate Gaussian: 1
#  bivariate t: 2
#  Gumbel: 3
#  reflected Gumbel: -3
#  Frank: 5
#  Gumbel/BB1 (Gumbel for global factor and BB1 for group factors): 11
# Output: dxn logret matrix and nx1 portfolio return vector
#  To get nxd logret matrix: 
#   transpose in this code, or make small changes in the C code
rgarchnestfact=function(n,grsize,garchpar,cpar,sigma0,copcode)
{ d=ncol(garchpar)
  m=length(grsize)  # sum(m) should equal d
  lgret=matrix(0,d,n)
  portfret=rep(0,n)
  out= .C("rgarchnestfact", as.integer(m), as.integer(grsize), as.integer(n), 
    as.double(garchpar[1,]), as.double(garchpar[2,]), as.double(garchpar[3,]), 
    as.double(garchpar[4,]), as.double(garchpar[5,]), as.double(garchpar[6,]), 
    as.double(sigma0), as.double(cpar), as.integer(copcode), 
    lgret=as.double(lgret), portfret=as.double(portfret) )
  list(logret=matrix(out$lgret,d,n),portfret=out$portfret)
}

#============================================================

# copula GARCH with bi-factor mvt model and interpolation for pt()
# better to create pp,qq,pder outside this function?
# copcode=2 is the only choice; the arguments are same as rgarchbifact()
rgarchbifactmvt=function(n,grsize,garchpar,cpar,sigma0,copcode=2)
{ d=ncol(garchpar)
  m=length(grsize)  # sum(m) should equal d
  lgret=matrix(0,d,n)
  portfret=rep(0,n)
  pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
  nipol=length(pp) # number of points for interpolation
  df=cpar[2*d+1]  # mvt df parameter 
  #cat("df=",df,"\n")
  qq=qt(pp,df)
  pder=pcderiv(qq,pp)
  #print(qq)
  out= .C("rgarchbifactmvt", as.integer(m), as.integer(grsize), as.integer(n), 
    as.double(garchpar[1,]), as.double(garchpar[2,]), as.double(garchpar[3,]), 
    as.double(garchpar[4,]), as.double(garchpar[5,]), as.double(garchpar[6,]), 
    as.double(sigma0), as.double(cpar), 
    as.integer(nipol), as.double(qq), as.double(pp), as.double(pder),
    lgret=as.double(lgret), portfret=as.double(portfret) )
  list(logret=matrix(out$lgret,d,n),portfret=out$portfret)
}

#============================================================

# copula GARCH with nested-factor mvt model and interpolation for pt()
# copcode=2 is the only choice; the arguments are same as rgarchnestfact()
rgarchnestfactmvt=function(n,grsize,garchpar,cpar,sigma0,copcode=2)
{ d=ncol(garchpar)
  m=length(grsize)  # sum(m) should equal d
  lgret=matrix(0,d,n)
  portfret=rep(0,n)
  pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
  nipol=length(pp) # number of points for interpolation
  df=cpar[m+d+1]  # mvt df parameter 
  #cat("df=",df,"\n")
  qq=qt(pp,df)
  pder=pcderiv(qq,pp)
  #print(qq)
  out=.C("rgarchnestfactmvt", as.integer(m), as.integer(grsize), as.integer(n), 
    as.double(garchpar[1,]), as.double(garchpar[2,]), as.double(garchpar[3,]), 
    as.double(garchpar[4,]), as.double(garchpar[5,]), as.double(garchpar[6,]), 
    as.double(sigma0), as.double(cpar), 
    as.integer(nipol), as.double(qq), as.double(pp), as.double(pder),
    lgret=as.double(lgret), portfret=as.double(portfret) )
  list(logret=matrix(out$lgret,d,n),portfret=out$portfret)
}
