# Fisher information for 1-factor and 2-factor item response models
# code written Aristidis Nikoloulopoulos

# d = dimension, 
# K = #categories, 
# ii = non-negative integer
# Output: vector of size d, each element in 1..K or 0..(K-1)
#d2v=function(d,K,ii,izero=F)
#{ 
#  tt=ii
#  jj=rep(0,d)
#  for(i in seq(d,1,-1))
#  { jj[i]=tt%%K; tt=floor(tt/K); }
#  if(!izero) jj=jj+1
#  jj
#}


# Fisher information matrix for 1-factor item response (IR) model
#   this works if ncat^d doesn't overflow
# nq = # quadature points
# ucuts = (ncat-1)xd matrix cut points on uniform scale
#   ncat=#categories, d=#items
# theta = vector of parameters for factor 1
# pcondcop = function for conditional cdf C_{2|1} for factor1
# pconddotcop = function for derivative of conditional cdf for factor1 wrt theta
# dcop = function for copula pdf for factor1
# Output: Fisher information matrix
irfisherinfo1=function(nq,ucuts,theta,pcondcop,pconddotcop,dcop)
{ gl=gausslegendre(nq)
  d=ncol(ucuts)
  ncat=nrow(ucuts)+1
  nvect=ncat^d
  condcdf=array(1,dim=c(nq,ncat-1,d))
  conddotcdf=array(0,dim=c(nq,ncat-1,d))
  for(j in 1:d)
  { for(k in 1:(ncat-1))
    { condcdf[,k,j]=pcondcop(ucuts[k,j],gl$nodes,theta[j])
      conddotcdf[,k,j]=pconddotcop(ucuts[k,j],gl$nodes,theta[j])
    }
  }
  arr0=array(0,dim=c(nq,1,d))
  arr1=array(1,dim=c(nq,1,d))
  condcdf=abind(arr0,condcdf,arr1,along=2)
  conddotcdf=abind(arr0,conddotcdf,arr0,along=2)
  fden=array(NA,dim=c(nq,ncat,d))
  fdotden=array(NA,dim=c(nq,ncat,d))
  for(j in 1:d)
  { for(k in 1:ncat)
    { fden[,k,j]=condcdf[,k+1,j]-condcdf[,k,j]
      fdotden[,k,j]=conddotcdf[,k+1,j]-conddotcdf[,k,j]
    }
  }
  der=matrix(NA,nvect,d)
  pr=rep(0,nvect)
  s=0
  for(i in 1:nvect)
  { y=d2v(d,ncat,i-1)
    for(jj in 1:d)
    { fproduct.i=1
      for(j in (1:d)[-jj])
      { temp=fden[,y[j],j]
        fproduct.i=fproduct.i*temp
      }
      der[i,jj]=(fproduct.i*fdotden[,y[jj],jj])%*%gl$w
    }
    deri=der[i,]
    pr[i]=c((fproduct.i*fden[,y[d],d])%*%gl$w)
    s=s+deri%*%t(deri)/pr[i]
  }
  #print(sum(pr)) # 1
  #print(apply(der,2,sum)) # rep(0,d)
  s
}

# Fisher information matrix for 2-factor IR model
#   this works if ncat^d doesn't overflow
# nq = # quadature points
# ucuts = (ncat-1)xd matrix cut points on uniform scale
#   ncat=#categories, d=#items
# theta = vector of parameters for factor 1
# delta = vector of parameters for factor 2
# pcondcop1 = function for conditional cdf C_{2|1} for factor1
# pcondcop2 = function for conditional cdf C_{2|1} for factor2
# pconddotcop1 = function for derivative of conditional cdf for factor1 wrt theta
# pconddotcop2 = function for derivative of conditional cdf for factor2 wrt delta
# dcop1 = function for copula pdf for factor1
# dcop2 = function for copula pdf for factor2
# iindep = TRUE if 2*d-1 linking copulas (last variable has conditional indep)
# iindep = FALSE if 2*d linking copulas
# Output: Fisher information matrix
irfisherinfo2=function(nq,ucuts,theta,delta,pcondcop1,pcondcop2,pconddotcop1,
  pconddotcop2,dcop1,dcop2,iindep=F)
{ gl=gausslegendre(nq)
  d=ncol(ucuts)
  ncat=nrow(ucuts)+1
  nvect=ncat^d
  condcdf=array(NA,dim=c(nq,ncat-1,d))
  conddotcdf=array(NA,dim=c(nq,ncat-1,d))
  cden=array(NA,dim=c(nq,ncat-1,d))
  for(j in 1:d)
  { for(k in 1:(ncat-1))
    { condcdf[,k,j]=pcondcop1(ucuts[k,j],gl$nodes,theta[j])
      conddotcdf[,k,j]=pconddotcop1(ucuts[k,j],gl$nodes,theta[j])
      cden[,k,j]=dcop1(ucuts[k,j],gl$n,theta[j])
    }
  }
  condcdf2=array(NA,dim=c(nq,nq,ncat-1,d))
  conddotcdf2=array(NA,dim=c(nq,nq,ncat-1,d))
  cden3=array(NA,dim=c(nq,nq,ncat-1,d))
  for(j in 1:d)
  { for(k in 1:(ncat-1))
  { for(u in 1:nq)
  { condcdf2[,u,k,j]=pcondcop2(condcdf[u,k,j],gl$nodes,delta[j])
    temp=dcop2(condcdf[u,k,j],gl$nodes,delta[j])
    cden3[,u,k,j]=temp*conddotcdf[u,k,j]
    conddotcdf2[,u,k,j]=pconddotcop2(condcdf[u,k,j],gl$nodes,delta[j])
  }}}
  arr0=array(0,dim=c(nq,nq,1,d))
  arr1=array(1,dim=c(nq,nq,1,d))
  condcdf2=abind(arr0,condcdf2,arr1,along=3)
  conddotcdf2=abind(arr0,conddotcdf2,arr0,along=3)
  cden3=abind(arr0,cden3,arr0,along=3)
  fden2=array(NA,dim=c(nq,nq,ncat,d))
  fbarden2=array(NA,dim=c(nq,nq,ncat,d))
  fdotden2=array(NA,dim=c(nq,nq,ncat,d))
  for(j in 1:d)
  { for(k in 1:ncat)
    { fden2[,,k,j]=condcdf2[,,k+1,j]-condcdf2[,,k,j]
      fdotden2[,,k,j]=conddotcdf2[,,k+1,j]-conddotcdf2[,,k,j]
      fbarden2[,,k,j]=cden3[,,k+1,j]-cden3[,,k,j]
    }
  }
  der.theta=matrix(NA,nvect,d)
  der.delta=matrix(NA,nvect,d)
  s=0
  for(i in 1:nvect)
  { y=d2v(d,ncat,i-1)
    for(jj in 1:d)
    { fproduct.i=1
      for(j in (1:d)[-jj])
      { temp=fden2[,,y[j],j]
        fproduct.i=fproduct.i*temp
      }
      der.theta[i,jj]=as.vector((fproduct.i*fbarden2[,,y[jj],jj])%*%gl$w)%*%gl$w
      der.delta[i,jj]=as.vector((fproduct.i*fdotden2[,,y[jj],jj])%*%gl$w)%*%gl$w
    }
    if(iindep) { der=c(der.theta[i,],der.delta[i,-d]) }
    else { der=c(der.theta[i,],der.delta[i,]) }
    prob=as.vector((fproduct.i*fden2[,,y[d],d])%*%gl$w)%*%gl$w
    s=s+der%*%t(der)/prob[1,1]
  }
  s
}
 
