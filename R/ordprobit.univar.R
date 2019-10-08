# Maximum likelihood for ordinal probit with covariates
# modified from mprobit R package

# ML for ordinal probit model, using modified Newton-Raphson
# x = (n*npred) matrix of predictors/covariates
# y = nx1 vector of ordinal categories in 1,...,ncateg or 0...ncateg-1
# iprint = print flag for NR iterations
# mxiter = maximum number of NR iterations
# toler = tolerance for convergence
ordprobit.univar=function(x,y,iprint=F,mxiter=20,toler=1.e-6)
{ if(!is.vector(x))
  { if(nrow(x)!=length(y)) stop("x, y not same length") }     
  else if(length(x)!=length(y)) { stop("x, y not same length") }
  if(is.vector(x)) x=as.matrix(x)
  n=length(y)
  # assume y in 1,...,ncateg or 0...ncateg-1
  ncateg=length(unique(y))
  if(min(y)==0) y=y+1 # convert to 1...ncateg
  npred=ncol(x)
  np=ncateg-1+npred # number of parameters
  # centering of x so that can use start of 0 for beta
  xmn=apply(x,2,"mean")
  xc=scale(x,center=xmn,scale=F)
  # starting point for NR
  cum=(1:(ncateg-1))
  cutp=rep(0,ncateg-1)
  for(k in cum)
  { pr=sum(y<=k)
    if (pr==0) { pr=1 }
    cutp[k]=qnorm(pr/n)
  }
  b=rep(0,npred)
  # loop
  mxdif=1
  iter=0
  while(iter<mxiter & mxdif>toler)
  { tem=xc%*%b
    cutb=c(-10,cutp,10)
    ub=tem+cutb[y+1]; lb=tem+cutb[y]
    ucdf=pnorm(ub); lcdf=pnorm(lb)
    updf=dnorm(ub); lpdf=dnorm(lb)
    # score vector
    dbeta=rep(0,npred); dcut=rep(0,ncateg+1)
    # Hessian matrix
    d2beta=matrix(0,npred,npred)
    d2bcut=matrix(0,npred,ncateg+1)
    d2cut=matrix(0,ncateg+1,ncateg+1)
    for(i in 1:n)
    { uderi=-updf[i]*ub[i]
      lderi=-lpdf[i]*lb[i]
      xx=xc[i,]
      pri=ucdf[i]-lcdf[i]; prderi=updf[i]-lpdf[i]
      dbeta=dbeta+xx*prderi/pri 
      k=y[i]
      dcut[k+1]=dcut[k+1]+updf[i]/pri
      dcut[k]=dcut[k]-lpdf[i]/pri
      pr2=pri^2
      d2beta=d2beta+ outer(xx,xx)*((uderi-lderi)*pri-prderi^2)/pr2 
      d2bcut[,k+1]=d2bcut[,k+1]+xx*(uderi*pri-updf[i]*prderi)/pr2
      d2bcut[,k]=d2bcut[,k]+xx*(-lderi*pri+lpdf[i]*prderi)/pr2
      d2cut[k+1,k+1]=d2cut[k+1,k+1]+ (uderi*pri-updf[i]^2)/pr2
      d2cut[k,k]=d2cut[k,k]+ (-lderi*pri-lpdf[i]^2)/pr2
      tem2=updf[i]*lpdf[i]/pr2
      d2cut[k,k+1]=d2cut[k,k+1]+tem2
      d2cut[k+1,k]=d2cut[k+1,k]+tem2
    }  
    sc=c(dcut[2:ncateg],dbeta)
    if(npred==1) d2bcut=matrix(c(d2bcut[,2:ncateg]),npred,ncateg-1)
    else d2bcut=d2bcut[,2:ncateg]
    d2cut=d2cut[2:ncateg,2:ncateg]
    h=cbind(d2cut,t(d2bcut))
    h=rbind(h,cbind(d2bcut,d2beta))
    dif=solve(h,sc)
    mxdif=max(abs(dif))
    cutp=cutp-dif[1:(ncateg-1)]
    b=b-dif[ncateg:np]
    # modification for cutp out of order
    chk=cutp[-1]-cutp[1:(ncateg-2)]
    ibad=sum(chk<=0)
    while(ibad>0)
    { dif=dif/2
      mxdif=mxdif/2
      cutp=cutp+dif[1:(ncateg-1)]
      b=b+dif[ncateg:np]
      chk=cutp[-1]-cutp[1:(ncateg-2)]
      ibad=sum(chk<=0)
    }
    iter=iter+1
    if(iprint)
    { cat("iter=",iter,", (with centered x's) cutp=", cutp, ", b=",b,"\n")
      cat("         scorevec=", sc,"\n\n")
    }
  }

  if(iter>=mxiter) cat("*** did not converge, check with iprint=T\n")
  # cutpoints with original x
  for(j in 1:npred)
  { cutp=cutp-b[j]*xmn[j] }
  if(iprint) cat("(with original x's) cutp=", cutp,"\n")

  # Hessian with original x's, repeat of previous code with x instead of xc
  tem=x%*%b
  cutb=c(-10,cutp,10)
  ub=tem+cutb[y+1]; lb=tem+cutb[y]
  ucdf=pnorm(ub); lcdf=pnorm(lb)
  updf=dnorm(ub); lpdf=dnorm(lb)
  nllk=0
  dbeta=rep(0,npred); dcut=rep(0,ncateg+1)
  d2beta=matrix(0,npred,npred); d2bcut=matrix(0,npred,ncateg+1)
  d2cut=matrix(0,ncateg+1,ncateg+1)
  for(i in 1:n)
  { uderi=-updf[i]*ub[i]; lderi=-lpdf[i]*lb[i]
    xx=x[i,]
    pri=ucdf[i]-lcdf[i]
    nllk=nllk-log(pri)
    prderi=updf[i]-lpdf[i]
    dbeta=dbeta+xx*prderi/pri 
    k=y[i]
    dcut[k+1]=dcut[k+1]+updf[i]/pri
    dcut[k]=dcut[k]-lpdf[i]/pri
    pr2=pri^2
    d2beta=d2beta+ outer(xx,xx)*((uderi-lderi)*pri-prderi^2)/pr2 
    d2bcut[,k+1]=d2bcut[,k+1]+xx*(uderi*pri-updf[i]*prderi)/pr2
    d2bcut[,k]=d2bcut[,k]+xx*(-lderi*pri+lpdf[i]*prderi)/pr2
    d2cut[k+1,k+1]=d2cut[k+1,k+1]+ (uderi*pri-updf[i]^2)/pr2
    d2cut[k,k]=d2cut[k,k]+ (-lderi*pri-lpdf[i]^2)/pr2
    tem2=updf[i]*lpdf[i]/pr2
    d2cut[k,k+1]=d2cut[k,k+1]+tem2
    d2cut[k+1,k]=d2cut[k+1,k]+tem2
  }  
  sc=c(dcut[2:ncateg],dbeta)
  if(npred==1) d2bcut=matrix(c(d2bcut[,2:ncateg]),npred,ncateg-1)
  if(npred>1) d2bcut=d2bcut[,2:ncateg]
  d2cut=d2cut[2:ncateg,2:ncateg]
  h=cbind(d2cut,t(d2bcut))
  h=rbind(h,cbind(d2bcut,d2beta))
  h=-h
  covm=solve(h)
  if(iprint)
  { cat("nllk= ", nllk,"\n")
    cat("cutpts= ", cutp,"\n")
    cat("beta= ", b,"\n")
    cat("SEs : ",sqrt(diag(covm)),"\n\n")
  }
  list(negloglik=nllk, cutpts=cutp, beta=b, cov=covm)
}

# multivariate ordinal (repeated measures) to uudat
# yvec = integer-valued nn-vector with values in 1...ncateg or 0..(ncateg-1)
# xmat = nn*npred matrix,
#   nn = nrep*ncl, ncl = #clusters, 
# nrep = #repeated measures per cluster
# b0cut = vector of cutpoints
# bvec = vector of regression coefficients
# Output: uudat matrix with dimension nx(2d) with corners of rectangle 
#   uu1[1] ... uu1[d] uu2[1] ... uu2[d]
mord2uu=function(xmat,yvec,nrep,b0cut,bvec)
{ d=nrep
  nn=length(yvec)
  if(min(yvec)==0) yvec=yvec+1  # convert to 1...ncateg
  ncl=nn/d  # number of clusters
  b0=c(-6,b0cut,6)
  if(is.matrix(xmat)) { tem=xmat%*%bvec }  # bvec is regression coef vector
  else { tem=xmat*bvec }
  low=b0[yvec]+tem
  upp=b0[yvec+1]+tem
  low=matrix(low,ncl,d,byrow=T)
  upp=matrix(upp,ncl,d,byrow=T)
  zzdat=cbind(low,upp)
  uudat=pnorm(zzdat)
  list(uudat=uudat,zzdat=zzdat)
}

#load("ord1.RData")
#xvec=c(t(xx))
#yvec=c(t(yy))

#out=ordprobit.univar(xvec,yvec,iprint=T)
#iter= 1 , (with centered x's) cutp= -0.5253297 0.5063273 , b= 0.3619286 
#         scorevec= -1.776357e-13 2.966516e-13 79.70742 
#
#iter= 2 , (with centered x's) cutp= -0.5258449 0.5068241 , b= 0.3642923 
#         scorevec= -0.3445046 0.3425231 0.4906575 
#
#iter= 3 , (with centered x's) cutp= -0.5258453 0.5068245 , b= 0.3642926 
#         scorevec= -0.0003125367 0.0003139111 3.675881e-05 
#
#(with original x's) cutp= -0.5285182 0.5041516 
#nllk=  859.4584 
#cutpts=  -0.5285182 0.5041516 
#beta=  0.3642926 
#SEs :  0.04696766 0.04678269 0.06786274

#names(out)
#[1] "negloglik" "cutpts"    "beta"      "cov" 
#out$cutpts
#[1] -0.5285182  0.5041516
# out$beta
#[1] 0.3642926

