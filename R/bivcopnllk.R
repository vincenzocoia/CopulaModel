# MLE for bivariate copula model, IFM and full log-likelihood

# second stage of IFM for dependence parameter
# param = copula parameter,
# udat = nx2 matrix of uniform scores 
# logdcop = function with log of copula density 
# ivect = T if logdcop can be used with vectorized inputs, ivect=F otherwise
# LB = lower bound (vector) for copula parameters 
# UB = lower bound (vector) for copula parameters 
# Output: negative log-likelihood
bivcopnllk=function(param,udat,logdcop,ivect=T,LB=0,UB=1000)
{ if(any(param<=LB) | any(param>=UB)) return(1.e10)
  if(ivect) { nllk=-sum(logdcop(udat[,1],udat[,2],param)) }
  else 
  { nllk=0; n=nrow(udat)
    for(i in 1:n) nllk=nllk-logdcop(udat[i,1],udat[i,2],param) 
  }
  if(is.nan(nllk) | is.infinite(nllk)) nllk=1.e10
  nllk
}

# ML of univariate + dependence parameters
# param = (margin1,margin2,copula) parameter vector
# xdat = data (untransformed)
# logdcop = function for log of copula density 
# logpdf1 = function for log of first univariate marginal density 
# cdf1 = function for first univariate marginal density
# np1 = dimension of parameter vector for cdf1
# logpdf2 = function for log of second univariate marginal density 
# cdf2 = function for second univariate marginal density
# np2 = dimension of parameter vector for cdf2
#   length(param)-np1-np2 is number of copula parameters
# ivect = T if logdcop can be used with vectorized inputs, ivect=F otherwise
# LB = lower bound vector for parameters 
# UB = lower bound vector for parameters 
# Output: negative log-likelihood
bivmodnllk=function(param,xdat,logdcop,logpdf1,cdf1,np1,logpdf2,cdf2,np2,
  ivect=T,LB,UB)
{ np=length(param)
  par1=param[1:np1]; par2=param[(np1+1):(np1+np2)]
  if(any(param<=LB) | any(param>=UB)) return(1.e10)  
  udat1=cdf1(xdat[,1],par1)
  lpdf1=logpdf1(xdat[,1],par1)
  udat2=cdf2(xdat[,2],par2)
  lpdf2=logpdf2(xdat[,2],par2)
  par3=param[(np1+np2+1):np]
  #lbpdf=logdcop(udat1,udat2,par3)
  if(ivect) { lbpdf=logdcop(udat1,udat2,par3) }
  else 
  { n=nrow(xdat); lbpdf=rep(0,n)
    for(i in 1:n) lbpdf[i]=logdcop(udat1[i],udat2[i],par3) 
  }
  nllk=-sum(lbpdf)-sum(lpdf1)-sum(lpdf2)
  if(is.nan(nllk) | is.infinite(nllk)) nllk=1.e10
  nllk
}

#============================================================

# copula nllk when univariate quantile function is computed using 
#  monotone interpolation for table lookup
# bivariate copula is  C(u,v,cpar)=F_{12}(F^{-1}(u),F^{-1}(v))
# and density of F is updf

# second stage of IFM for dependence parameter
# param = copula parameter,
# udat = nx2 matrix of uniform scores 
# logbpdf = function with log of density of F_{12} (vectorizable)
# logupdf = function with log of univariate density updf (vectorizable)
# uquant = function for F^{-1}
# iunivar = indices such that par1=param[iunivar] is parameter of F
# ppvec = quantiles to use for interpolation, 
#         there is a default if it is input as NULL
# LB = lower bound (vector) for copula parameters 
# UB = lower bound (vector) for copula parameters 
# Output: negative log-likelihood
bivcopnllk.ipol=function(param,udat,logbpdf,logupdf,uquant,iunivar,
  ppvec=NULL,LB,UB)
{ if(any(param<=LB) | any(param>=UB)) return(1.e10)
  par1=param[iunivar]
  if(length(ppvec)==0) ppvec=c(0.0001,0.001,seq(.01,.99,.02),.999,.9999)
  # univariate quantiles for copula
  qqvec=uquant(ppvec,par1)
  ppder=pcderiv(ppvec,qqvec)
  x1=pcinterpolate(ppvec,qqvec,ppder,udat[,1])[,1]
  x2=pcinterpolate(ppvec,qqvec,ppder,udat[,2])[,1]
  lgnumer=sum(logbpdf(x1,x2,param)) 
  lgdenom=sum(logupdf(x1,par1))+sum(logupdf(x2,par1))
  -lgnumer+lgdenom
}

