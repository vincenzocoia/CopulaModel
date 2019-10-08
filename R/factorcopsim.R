# Simulation for 1-factor and 2-factor copula models

# n = sample size, parobj has d-dimensional vector of dependence parameters  
#         or dx(npar) matrix where npar is #parameters for the pair-copula
# parobj1 = d-dimensional vector of dependence parameters for factor 1 
#         or dxm matrix where m is #parameters for the pair-copula
# qcond1 = the function for the conditional quantile function of the copula
# copname1 = string for name of copula
# ivect = T if qcond function has vectorized inputs
# Output: d-dimensional random sample with U(0,1) margins
sim1fact=function(n,parobj1,qcond1,copname1,ivect=F)
{ vv=runif(n) # latent
  if(is.vector(parobj1)) d=length(parobj1)
  else d=nrow(parobj1)   # assuming matrix here
  copname=tolower(copname1)
  yy=matrix(0,n,d)
  if(copname=="frank" | copname=="frk" | copname=="mtcj" |  copname=="mtcjr" | 
     copname=="fgm" | ivect)
  { # can vectorize because qcond1 has closed form or has been vectorized
    for(j in 1:d)
    { qq=runif(n)
      if(is.vector(parobj1)) th=parobj1[j]
      else th=parobj1[j,]
      yy[,j]=qcond1(qq,vv,th)   
    }
  }
  else # cannot vectorize
  { for(j in 1:d)
    { if(is.vector(parobj1)) th=parobj1[j]
      else th=parobj1[j,]
      for(i in 1:n)
      { v1=vv[i]
        qq=runif(1)
        yy[i,j]=qcond1(qq,v1,th)
      }
    }
  }
  yy
}


# n = sample size, 
# parobj1 = d-dimensional vector of dependence parameters for factor 1 
#         or dxm matrix where m is #parameters for the pair-copula
# parobj2 =  d-dimensional vector or dxm matrix for factor 2 
# qcond1 = function for the conditional quantile function of factor 1
# qcond2 = function for the conditional quantile function of factor 2
# copname1 = string for name of copula for factor 1
# copname2 = string for name of copula for factor 2
# ivect = T if qcond functions have vectorized inputs
# Output : a d-dimensional random sample with U(0,1) margins
sim2fact=function(n,parobj1,parobj2,qcond1,qcond2,copname1="a1",copname2="a2",
  ivect=F)
{ copname1=tolower(copname1); copname2=tolower(copname2)
  coptab=c("frank","frk","mtcj","mtcjr","fgm")
  vv1=runif(n)
  vv2=runif(n)
  if(is.vector(parobj1)) d=length(parobj1)
  else d=nrow(parobj1)   # assuming matrix here
  yy=matrix(0,n,d)
  if(copname1 %in% coptab && copname2 %in% coptab | ivect)
  { # can vectorize because qcond has closed form or has been vectorized
    for(j in 1:d)
    { qq=runif(n)
      tem=qcond2(qq,vv2,parobj2[j])   # parobj1 is vector in these cases
      yy[,j]=qcond1(tem,vv1,parobj1[j])   # parobj2 is vector in these cases
    }
  }
  else # cannot vectorize 
  { for(j in 1:d)
    { if(is.vector(parobj1)) th1=parobj1[j] else th1=parobj1[j,]
      if(is.vector(parobj2)) th2=parobj2[j] else th2=parobj2[j,]
      for(i in 1:n)
      { v1=vv1[i]; v2=vv2[i]
        qq=runif(1)
        tem=qcond2(qq,v2,th2)
        yy[i,j]=qcond1(tem,v1,th1)
      }
    }
  }
  yy
}

