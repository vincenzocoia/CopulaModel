# multivariate t with p-factor, bi-factor, tri-factor correlation structure
# R interface to tpfactnllk/tbifactnllk/ttrifactnllk nllk+grad 

# rhvec = vector of length d*p with partial corr representation of loadings
# tdata =  nxd data set of t(df)-scores
# df = degree of freedom parameter >0
#   add ifixed option??
# Output: negative log-likelihood of copula for mvt p-factor model
mvtpfactnllk=function(rhvec,tdata,df)
{ d=ncol(tdata)
  n=nrow(tdata)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  p=length(rhvec)/d
  if(df< 300)
  { out= .Fortran("tpfactnllk",
      as.integer(d), as.integer(p),  as.integer(n), as.double(df),
      as.double(rhvec), as.double(tdata),
      nllk=as.double(0), grad=as.double(rep(0,d*p))  )
  }
  else # df>=300
  { Robs=cor(tdata)  
    out= .Fortran("pfactnllk",
      as.integer(d), as.integer(p), as.double(rhvec), as.double(Robs),
      as.integer(n),
      nllk=as.double(0), grad=as.double(rep(0,d*p))  )
    lgdenom = -.5*sum(tdata^2)-n*d*.5*log(2*pi)
    out$nllk = out$nllk + lgdenom
  }
  #nllk=out$nllk; 
  #attr(nllk,"gradient") = out$grad;
  #nllk
  list(nllk=out$nllk,lgrad=out$grad)
}

# rhvec = vector of length d*2 for partial corr representation of loadings,
#   first d correlations with common factor, then
#   partial correlations with group factor given common factor
# grsize = vector of group sizes for bi-factor model
# tdata =  nxd data set of t(df)-scores
# df = degree of freedom parameter >0
# Output: negative log-likelihood of copula for mvt bi-factor model
mvtbifactnllk=function(rhvec,grsize,tdata,df)
{ d=ncol(tdata)
  n=nrow(tdata)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  mgrp=length(grsize)
  if(df< 300)
  { out= .Fortran("tbifactnllk",
      as.integer(d), as.integer(mgrp), as.integer(grsize),
      as.integer(n), as.double(df),  as.double(rhvec), as.double(tdata), 
      nllk=as.double(0), grad=as.double(rep(0,d*2))  )
  }
  else # df>=300
  { Robs=cor(tdata)
    out= .Fortran("bifactnllk",
      as.integer(d), as.integer(mgrp), as.double(rhvec), as.integer(grsize),
      as.double(Robs), as.integer(n),
      nllk=as.double(0), grad=as.double(rep(0,d*2))  )
    lgdenom = -.5*sum(tdata^2)-n*d*.5*log(2*pi)
    out$nllk = out$nllk + lgdenom
  }
  #nllk=out$nllk; 
  #attr(nllk,"gradient") = out$grad;
  #nllk
  list(nllk=out$nllk,lgrad=out$grad)
}

#============================================================
# front ends to the nllk functions, original version written by Pavel Krupskii

# MLE in a MVt model with a p-factor correlation structure
# tdata = nxd matrix of t-scores
# start = vector of length p*d with starting values
# pfact = number of factors (such as 1,2,3,... )
# df = degrees of freedom parameter >0
# prlevel = print.level for nlm()
# mxiter = maximum number of iterations for nlm()
# Outputs: nlm object with ($code,$estimate,$gradient,$iterations,$minimum) 
# Note the minimum nllk can be the same for different parameter vectors
#   because of invariance of the loading matrix to rotations
mvtpfact=function(tdata,start,pfact,df,prlevel=0,mxiter=100)
{ d=ncol(tdata);
  n=nrow(tdata);
  if(length(start) != pfact*d) { cat("incorrect number of parameters in start\n"); return(NULL); }
  nllkfn= function(param)
  { if(max(abs(param))>0.999) { return(1.e10) }
    out=mvtpfactnllk(param,tdata,df);
    nllk=out$nllk; 
    attr(nllk,"gradient")=out$lgrad;
    nllk
  } 

  # *** should hessian=T be set so that hessian is returned?
  mle= try(nlm(nllkfn,p=start,hessian=F,iterlim=mxiter,print.level=prlevel,
       check.analyticals = F),T) 
  if(is.list(mle)==F) { cat("mvtpfact() has failed\n"); return(NULL); }
  pmle=mle$estimate
  cat("MLE: \n")
  pmle=matrix(pmle,ncol=pfact);
  print(pmle)
  #cat("SEs: \n")
  #if (max(abs(pmle)) < 0.995 & rcond(mle$hessian) > 1e-12) 
  #{ acov=solve(mle$hessian); } 
  #else {print("NA")}
  cat("nllk:  \n"); cat(mle$minimum,"\n")
  mle
}  

# MLE in a MVt model with a bi-factor correlation structure
# tdata = nxd matrix of t-scores
# start = vector of length 2*d with starting values of partial correlations
#    values for correlations of observed Z_{gj} and common latent V_0 go first,
#    then partial correlations of Z_{gj} and V_g given V_0 (j in group g)
# grsize = vector of group sizes for bi-factor model
# df = degrees of freedom parameter >0
# prlevel = print.level for nlm()
# full = T for bi-factor
# full = F for nested-factor (reduced model with fewer parameters)
#     V_g = rh2_g*V_0 + sqrt(1-rh2_g^2)*eps_g; 
#     Z_{gj} = rh1_{gj}*V_g + sqrt(1-rh1_{gj}^2)*eps_{gj} (j in group g)
# mxiter = maximum number of iterations for nlm()
# Outputs: nlm object with ($code,$estimate,$gradient,$iterations,$minimum) 
# Note the minimum nllk can be the same for different parameter vectors
#   if some group size values are 1 or 2.
mvtbifact=function(tdata,start,grsize,df,prlevel=0,full=T,mxiter=100)
{ d=ncol(tdata);
  n=nrow(tdata);
  if(full==T & length(start)!=2*d) 
  { cat("start should have length 2*d for bi-factor\n"); return(0) }
  if(full==F & length(start)!=d+length(grsize)) 
  { cat("start should have length d+length(grsize) for nested-factor\n"); 
    return(0) 
  }
  #err=0;
  nllkfn= function(param)
  { if(max(abs(param))>0.999) { return(1.e10) }
    param0 = param;
    if(full==F) # nested factor model
    { mgrp=length(grsize);
      ind=0;
      a2=param[(mgrp+1):(mgrp+d)]; a1=rep(0,d); 
      for (jg in 1:mgrp)
      { ind1=ind+1;  ind2=ind+grsize[jg]; 
        a1[ind1:ind2]=param[jg];  
        ind=ind+grsize[jg];  
      }
      rh1=a1*a2; 
      rh2=a2*sqrt(1-a1^2)/sqrt(1-rh1^2); 
      param0=c(rh1,rh2); 
    }
    #out = bifcttlik.vec(tdata,grsize,param0,df);
    out=mvtbifactnllk(param0,grsize,tdata,df);
    nllk=out$nllk; 
    lgrd=out$lgrad;
    if(full==T) attr(nllk,"gradient")=lgrd;
    if(full==F)
    { lgrd0=rep(0,mgrp+d); 
      lgrad1=lgrd[1:d]; lgrad2=lgrd[(d+1):(2*d)];
      b1=(1-(a1*a2)^2)^1.5; b2=sqrt(1-a1^2);
      tem1= lgrad1*a2-lgrad2*(1-a2^2)*a1*a2/b1/b2;
      ind0=0; 
      for (jg in 1:mgrp) 
      { ind1=1+ind0; ind2=grsize[jg]+ind0;  
        lgrd0[jg]=sum(tem1[ind1:ind2]); 
        ind0=ind0+grsize[jg]; 
      }   
      lgrd0[(mgrp+1):(mgrp+d)]=lgrad1*a1+lgrad2*b2/b1;
      attr(nllk,"gradient")=lgrd0;
    }
    nllk
  } 

  # *** should hessian=T be set so that hessian is returned?
  mle= try(nlm(nllkfn,p=start,hessian=F,iterlim=mxiter,print.level=prlevel,
     check.analyticals = F),T)
  if(is.list(mle)==F) { cat("mvtbifact() has failed\n"); return(NULL); }
  bifmle=mle$estimate
  cat("MLE: \n")
  if(full==T) bifmle=matrix(bifmle,ncol=2);
  print(bifmle)
  #cat("SEs: \n")
  #if (max(abs(pmle)) < 0.995 & rcond(mle$hessian) > 1e-12) 
  #{ acov=solve(mle$hessian); } 
  #else {print("NA")}
  cat("nllk:  \n"); print(mle$minimum)
  mle
}  

#============================================================

# tri-factor

# rhvec = vector of length d*3 for partial corr representation of loadings,
#   first d correlations with common factor, then
#   partial correlations with group factor given common factor, then
#   partial correlations with subgroup factor given common and group factors
# grsize = vector of group sizes 
# sbgrsize = vector of subgroup sizes 
# tdata =  nxd data set of t(df)-scores
# df = degree of freedom parameter >0
# Output: negative log-likelihood of copula for mvt bi-factor model
mvttrifactnllk=function(rhvec,grsize,sbgrsize,tdata,df)
{ d=ncol(tdata)
  n=nrow(tdata)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  mgrp=length(grsize)
  msbgrp=length(sbgrsize)
  if(df<300)
  { out= .Fortran("ttrifactnllk",
      as.integer(d), as.integer(mgrp), as.integer(msbgrp), 
      as.integer(grsize), as.integer(sbgrsize),
      as.integer(n), as.double(df),  as.double(rhvec), as.double(tdata), 
      nllk=as.double(0), grad=as.double(rep(0,d*3))  )
  }
  else # df>=300
  { Robs=cor(tdata)
    out= .Fortran("trifactnllk",
      as.integer(d), as.integer(mgrp), as.integer(msbgrp),
      as.double(rhvec), as.integer(grsize), as.integer(sbgrsize),
      as.double(Robs), as.integer(n),
      nllk=as.double(0), grad=as.double(rep(0,d*3))  )
    lgdenom = -.5*sum(tdata^2)-n*d*.5*log(2*pi)
    out$nllk = out$nllk + lgdenom
  }

  #nllk=out$nllk; 
  #attr(nllk,"gradient") = out$grad;
  #nllk
  list(nllk=out$nllk,lgrad=out$grad)
}


# MLE in a MVt model with a tri-factor or 3-nested correlation structure
# tdata = nxd matrix of t-scores
# start = vector of length 3*d with starting values of partial correlations
#    values for correlations of observed Z_{gj} and common latent V_0 go first,
#    then partial correlations of Z_{gj} and V_g given V_0 (j in group g)
#    then partial correlations of Z_{ghj} and V_{gh} given V_0,V_g 
#                    (j in group g, subgroup h)
# grsize = vector of group sizes 
# sbgrsize = vector of subgroup sizes 
# df = degrees of freedom parameter >0
# prlevel = print.level for nlm()
# full = T for tri-factor
# full = F for 3-nested-factor (reduced model with fewer parameters)
#     V_g = rh2_g*V_0 + sqrt(1-rh2_g^2)*eps_g; 
#     Z_{gj} = rh1_{gj}*V_g + sqrt(1-rh1_{gj}^2)*eps_{gj} (j in group g)
#     Z_{ghj}
# mxiter = maximum number of iterations for nlm()
# Outputs: nlm object with ($code,$estimate,$gradient,$iterations,$minimum) 
# Note the minimum nllk can be the same for different parameter vectors
#   if some group/subgroup size values are 1 or 2.
mvttrifact=function(tdata,start,grsize,sbgrsize,df,prlevel=0,full=T,mxiter=150)
{ 
  d=ncol(tdata);
  n=nrow(tdata);
  if(full==T & length(start)!=3*d) 
  { cat("start should have length 3*d for tri-factor\n"); return(0) }
  if(!subgr.consistent(grsize,sbgrsize))
  { cat("grsize and sbgrsize not consistent\n"); return(-1.e10) }
  nllkfn= function(param)
  { if(max(abs(param))>0.999) { return(1.e10) }
    param0 = param;
    if(full==F) # nested factor model
    { par0=nest2condcor(param0,grsize,sbgrsize); 
      th11=par0$th11; th22=par0$th22; th33=par0$th33;
      th1c=par0$th1c; th2c=par0$th2c; th3c=par0$th3c;
      param0=c(th1c,th2c,th3c); 
    }
    out=mvttrifactnllk(param0,grsize,sbgrsize,tdata,df);
    nllk=out$nllk; 
    lgrd=out$lgrad;
    if(full==T) attr(nllk,"gradient")=lgrd;
    if(full==F) 
    { nstlgrd=grad2nestgrad(th11,th22,th33,th1c,th2c,th3c,lgrd,grsize,sbgrsize);
      attr(nllk,"gradient")=nstlgrd; 
    }
    nllk
  } 
  # should hessian=T be set so that hessian is returned
  mle= try(nlm(nllkfn,p=start,hessian=F,iterlim=mxiter,print.level=prlevel,
     check.analyticals = F),T)
  if(is.list(mle)==F) { cat("mvttrifact() has failed\n"); return(NULL); }
  trifmle=mle$estimate
  cat("MLE: \n")
  if(full==T) trifmle=matrix(trifmle,ncol=3);
  print(trifmle)
  #cat("SEs: \n")
  #if (max(abs(pmle)) < 0.995 & rcond(mle$hessian) > 1e-12) 
  #{ acov=solve(mle$hessian); } 
  #else {print("NA")}
  cat("nllk:  \n"); print(mle$minimum)
  mle
}  

