
# Multivariate ordinal response with no predictors/covariates
# Analogue of polychoric() for bivariate margins when 
#  bivariate Gaussian copula is replaced by another bivariate copula 
# odatfr = ordinal data matrix, nrec x (d+1), d=#items, nrec=num of records
#    column d+1 has frequencies of the d-vectors
#   ncateg = number of ordinal categories
#   categories are in 0,1,...,ncateg-1
# ucuts = (ncateg+1)xd matrix of cutpoints in U(0,1) scale, obtained via unifcuts 
# pcop = bivariate copula family
# cparstart = vector of starting point; length d*(d-1)/2 times dimension of
#       parameter for pcop
# LB = lower bound on copula parameter of pcop
# UB = upper bound on copula parameter of pcop
# iprint = flag for printing of intermediate results
# prlevel = print.level for nlm()
# Outputs:
#   $nllkvec = vector of length d*(d-1)/2 of bivariate nllk
#   $cparvec = vector of length d*(d-1)/2 times the dim of parameter in pcop, 
#       with bivariate MLEs for each pair
#   $summary = summary array of observed vs model-expected counts
#       dimension is ncateg x (2*ncateg) x d*(d-1)/2
ordinal.bivcop=function(odatfr,ucuts,pcop,cparstart,LB=0,UB=10,
   iprint=F,prlevel=0)
{ d=ncol(odatfr)-1
  fr=odatfr[,d+1]
  dat=odatfr[,1:d]
  nc=nrow(ucuts)-1  # same as ncateg in documentation
  sampsize=sum(fr)
  kk=0
  dd=d*(d-1)/2 
  npar=length(cparstart)/dd
  #cat("npar=",npar,"\n")
  if(npar==0 | npar!=floor(npar)) 
  { cat("length of cparstart should be d*(d-1)/2 times dim of parameter for pcop\n");
    return(0)
  } 
  nllkvec=rep(0,dd)
  cparvec=rep(0,npar)
  summ=array(0,c(nc,2*nc,dd))
  for(j2 in 2:d)
  { for(j1 in 1:(j2-1))
    { if(iprint) cat("\nitems:", j1, j2, "\n")
      kk=kk+1
      y1=dat[,j1]; y2=dat[,j2]
      # addition in case of zero counts for some categories
      y1=c(y1,0:(nc-1))
      y2=c(y2,0:(nc-1))
      fradd=c(fr,rep(1,nc))
      bcount=tapply(fradd,list(y1,y2),sum)
      bcount[is.na(bcount)]=0
      diag(bcount)=diag(bcount)-1
      #if(iprint) print(bcount)
      kindex=((kk-1)*npar+1):(kk*npar)
      out=nlm(bivcopOrdinalnllk,p=cparstart[kindex],hessian=T,iterlim=50,
        print.level=prlevel, ucuts=ucuts,bfr=bcount,jj1=j1,jj2=j2,
        pcop=pcop,LB=LB,UB=UB)
      if(iprint) 
      { cparse=sqrt(diag(solve(out$hess)))
        cat("nllk min ",out$min,"\n") 
        cat("estimate ",out$estimate,"\n") 
        cat("SE ", cparse,"\n") 
      }
      nllkvec[kk]=out$min
      cparvec[kindex]=out$estimate
      # model-based expected bivariate table 
      nc1=nc+1
      u1=matrix(ucuts[,j1],nc1,nc1)
      u2=matrix(ucuts[,j2],nc1,nc1,byrow=T)
      cdf2= pcop(u1,u2,out$estimate)
      pmf=apply(cdf2,2,diff)
      pmf2=apply(t(pmf),2,diff)
      pmf2=t(pmf2)
      pmf2[pmf2<=0]=1.e-15
      ecount=pmf2*sampsize
      if(iprint) { cat("observed and expected\n"); print(cbind(bcount,ecount))}
      summ[,,kk]=cbind(bcount,ecount)
    }
  }
  list(nllkvec=nllkvec,cparvec=cparvec,summary=summ)
}

# This function is used by ordinal.bivcop;
# pcop is used to compute the bivariate cdf, 
#  then row and column differencing is used to get the bivariate pmf
# cpar = copula parameter 
# ucuts = (ncateg+1)xd  matrix of cutpoints on U(0,1) scale
# bfr = vector of bivariate frequencies
# jj1, jj2 = indices of 2 variables (between 1 and d)
# pcop = bivariate copula family
# LB = lower bound on copula parameter
# UB = upper bound on copula parameter
# Output: negative log-likelihood
bivcopOrdinalnllk=function(cpar,ucuts,bfr,jj1,jj2,pcop,LB=0,UB=10)
{ nc=nrow(ucuts)-1
  nlk=0.; 
  if(any(cpar<=LB) | any(cpar>=UB)) { return(1.e10) }
  nc1=nc+1
  u1=matrix(ucuts[,jj1],nc1,nc1)
  u2=matrix(ucuts[,jj2],nc1,nc1,byrow=T)
  cdf2= pcop(u1,u2,cpar)
  pmf=apply(cdf2,2,diff)
  pmf2=apply(t(pmf),2,diff)
  pmf2=t(pmf2)
  pmf2[pmf2<=0]=1.e-15
  nlk=-sum(bfr*log(pmf2))
  nlk
}

