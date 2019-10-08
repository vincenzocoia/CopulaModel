# Functions for simulation from structured copula models
# Code written by Pavel Krupskii

# Simulate data for 2-nested copula model (groups) [old version]
# if cop = 1, data have standard normal marginals
# if cop > 1, data have uniform marginals
# nn = sample size; 
# grsize = vector of group sizes
# cop = 1: normal, 3: gumbel, 5: frank
# param = vector of parameters, length is d+mgrp : starting with mgrp?
# the order in param is the same as in start for mvtbifct(full=F) function
#simnestfact0= function(nn, grsize, param, cop=5, qcond=qcondfrk, ivect=T)
#{ d = sum(grsize);
#  mgrp = length(grsize);
#  if(cop==3) { qcond=qcondgum; } 
#  if(cop==1) { z0 = rnorm(nn); } else z0 = runif(nn);
#  z = matrix(0,nrow = nn,ncol = mgrp);
#  zdata = matrix(0,nrow = nn,ncol = d)
#  # change interface later with ivect and "frank" etc
#  if(cop == 1) 
#  { for (jg in 1:mgrp) { z[,jg] = z0*param[jg] + sqrt(1-param[jg]^2)*rnorm(nn); }}
#  else if(ivect) 
#  { for (jg in 1:mgrp) { z[,jg] = qcond(runif(nn),z0,param[jg]); } } 
#  else # ivect=F
#  { for (jg in 1:mgrp)  
#    { for (i in 1:nn) { z[i,jg] = qcond(runif(1),z0[i],param[jg]); } }
#  }
#
#  ind = 0;
#  for (jg in 1:mgrp)  
#  { ind1 = ind + 1;  
#    ind2 = ind + grsize[jg];
#    ind  = ind + grsize[jg];  
#    for(ij in ind1:ind2)  
#    { ijm=ij+mgrp
#      if(cop==1) { zdata[,ij]=z[,jg]*param[ijm]+sqrt(1-param[ijm]^2)*rnorm(nn); }
#      else if(ivect) { zdata[,ij]=qcond(runif(nn),z[,jg],param[ijm]); }
#      else 
#      { for(i in 1:nn) { zdata[i,ij]=qcond(runif(1),z[i,jg],param[ijm]); } } 
#      #else print("this copula family is not implemented");
#    }
#  }
#  zdata  
#}


# Simulate data from 2-nested copula model; second version
# nn = sample size 
# grsize = vector of group sizes 
# cop = number code: 1: Gaussian; 2: t; 3: Gumbel; 5: Frank; 10: gum+bb1
#  if cop = 1, data have standard normal marginals
#  if cop > 1, data have uniform marginals
# param = vector of parameters, length is d+mgrp(+1 for cop==2): 
#     starting with mgrp
#     length(param) = d+d+mgrp if cop == 10
# The order in param is the same as in start for mvtbifct(full=F) function.
# Output: d-dimensional random sample with U(0,1) margins
simnestfact=function(nn,grsize,cop,param)
{ d=sum(grsize);
  mgrp=length(grsize);
  if(cop==2) { df0=param[mgrp+d+1]; }
  if(cop==1 || cop==2) { z0=rnorm(nn); } else { z0=runif(nn); }
  z=matrix(0,nrow=nn,ncol=mgrp);
  zdata=matrix(0,nrow=nn,ncol=d)
  if(cop==1 || cop==2) 
  { for (jg in 1:mgrp) { z[,jg]= z0*param[jg]+ sqrt(1-param[jg]^2)*rnorm(nn); }}
  else if(cop==5) 
  { for (jg in 1:mgrp) { z[,jg]= qcondfrk(runif(nn),z0,param[jg]); } } 
  else if(cop==3 || cop==10)  
  { for (jg in 1:mgrp)  
    { for (i in 1:nn) { z[i,jg]= qcondgum(runif(1),z0[i],param[jg]); } }
  }
  
  ind = 0;
  for (jg in 1:mgrp)  
  { ind1=ind+1;  
    ind2=ind+grsize[jg];
    ind =ind+grsize[jg];  
    for(ij in ind1:ind2)  
    { ijm=ij+mgrp
      if(cop==1 || cop==2) 
      { zdata[,ij]=z[,jg]*param[ijm]+sqrt(1-param[ijm]^2)*rnorm(nn); }
      else if(cop==5) { zdata[,ij]=qcondfrk(runif(nn),z[,jg],param[ijm]); }
      else if(cop==3)
      { for(i in 1:nn) { zdata[i,ij]=qcondgum(runif(1),z[i,jg],param[ijm]); } } 
      else if(cop==10)
      { for(i in 1:nn) 
        { zdata[i,ij]=qcondbb1(runif(1),z[i,jg],param[c(ijm,ijm+d)]); } 
      } 
      else print("this copula family is not implemented");
    }
  }
  if(cop==2) 
  { for(i in 1:nn)
    { zdata[i,]=zdata[i,]/sqrt(rchisq(1,df=df0)/df0); 
    }
  }
  zdata  
}

# Simulate data from structured factor copula model
# nn = sample size 
# grsize = vector of group sizes 
# cop = number code: 1: Gaussian; 2: t; 3: Gumbel; 5: Frank; 9: bb1+frank
#  if cop = 1, data have standard normal marginals
#  if cop = 2, data have t marginals
#  if cop > 2, data have uniform marginals
# param = vector of parameters (those for the common factor go first)
# The order in param is the same as in start for mvtbifct(full=T) function.
# Output: d-dimensional random sample with U(0,1) or N(0,1) or t(df) margins
simbifact=function(nn,grsize,cop=5,param)
{ ivect=F;
  d=sum(grsize);
  mgrp=length(grsize);
  th1=param[1:d];
  th2=param[(d+1):(2*d)];
  if(cop==2) { df0=param[2*d+1];}
  if(cop==3) { qcond=qcondgum; ivect=T } 
  #  ivect=T if qcond can be vectorized
  if(cop==5) { qcond=qcondfrk; ivect=T } 
  if(cop==9) { qcond1=qcondbb1; qcond2=qcondfrk; th3=param[(2*d+1):(3*d)];} 
  #if(cop==1) th2 = th2*sqrt(1-th1^2);
  if(cop==1 || cop==2) lm2 = th2*sqrt(1-th1^2);
  zdata = matrix(0,nrow = nn,ncol = d);
  if(cop==1 || cop==2) # normal or Student t
  { z0=rnorm(nn);
    z=matrix(rnorm(nn*mgrp), ncol=mgrp);
  } 
  else 
  { z0=runif(nn);
    z=matrix(runif(nn*mgrp), ncol=mgrp);
  }
  ind=0;
  for (jg in 1:mgrp)  
  { ind1=ind+1;  
    ind2=ind+grsize[jg];
    ind =ind+grsize[jg];  
    for(ij in ind1:ind2) 
    { if(cop==1 || cop==2)  # normal or Student t
      { zdata[,ij]=th1[ij]*z0+z[,jg]*lm2[ij]+sqrt(1-th1[ij]^2-lm2[ij]^2)*rnorm(nn); }
      else if(ivect) 
      { q1=qcond(runif(nn),z[,jg],th2[ij]); 
        zdata[,ij]=qcond(q1,z0,th1[ij]); 
      }
      else if(cop==9)
      { for(i in 1:nn) 
        { q1=qcond2(runif(1),z[i,jg],th3[ij]); 
          zdata[i,ij]=qcond1(q1,z0[i],c(th1[ij],th2[ij])); 
        } 
      }
      else # ivect=F
      { for(i in 1:nn) 
        { q1=qcond(runif(1),z[i,jg],th2[ij]); 
          zdata[i,ij]=qcond(q1,z0[i],th1[ij]); 
        } 
      }
    } 
  }
  if(cop==2) 
  { for(i in 1:nn)
    { zdata[i,]=zdata[i,]/sqrt(rchisq(1,df=df0)/df0); 
    }
  }
  zdata
}


# Simulate data from 3-nested copula model (groups and subgroups)
# nn = sample size 
# grsize = vector of group sizes 
# sbgrsize = vector with sizes of subgroups
# cop = number code: 1: Gaussian; 3: Gumbel; 5: Frank; 
#   cop=2: t, for this copula the argument df is a vector of length 3 is used
#  if cop = 1, data have standard normal marginals
#  if cop > 1, data have uniform marginals (incl. t copula)
# param = vector of parameters, length is mgrp+msbgrp+d : starting with mgrp
# df = degree of freedom 3-vector if cop=2
# The order in param is the same as in start for mvttrifct(full=F) function
# Output: d-dimensional random sample with U(0,1) or N(0,1) or t(df) margins
simnest3fact= function(nn,grsize,sbgrsize,param,df,cop=5)
{ d= sum(grsize);
  mgrp= length(grsize);
  msbgrp = length(sbgrsize);
  if(cop==5) { qcond=qcondfrk; ivect=T} 
  if(cop==3) { qcond=qcondgum; ivect=T} 
  if(cop==2) { qcond=qcondt; ivect=T}
  if(cop==1) { z0=rnorm(nn); } else z0=runif(nn);
  z=matrix(0,nrow=nn,ncol=mgrp);
  zz=matrix(0,nrow=nn,ncol=msbgrp);

  zdata=matrix(0,nrow=nn,ncol=d)
  if(cop==1) 
  { for (jg in 1:mgrp) { z[,jg]= z0*param[jg]+ sqrt(1-param[jg]^2)*rnorm(nn); }}
  else if(cop == 2) 
  { for (jg in 1:mgrp) { z[,jg]= qcond(runif(nn),z0,param[jg],df[1]); } } 
  else if(ivect) 
  { for (jg in 1:mgrp) { z[,jg]= qcond(runif(nn),z0,param[jg]); } } 
  else # ivect=F
  { for (jg in 1:mgrp)  
    { for (i in 1:nn) { z[i,jg]= qcond(runif(1),z0[i],param[jg]); } }
  }
  sumv=0; ngr=rep(0,mgrp); sbind=1; temv=0
  for (jg in 1:mgrp)                     
  { temv=temv+grsize[jg]
    while(sumv<temv)
    { sumv=sumv+sbgrsize[sbind]
      ngr[jg]=ngr[jg]+1; sbind=sbind+1
    }  
  }       
  ind=0;
  for (jg in 1:mgrp)  
  { ind1=ind+1;  
    ind2=ind+ngr[jg];
    ind=ind+ngr[jg]; 
    for(ij in ind1:ind2)  
    { ijm=ij+mgrp
      if(cop==1) { zz[,ij]=z[,jg]*param[ijm]+sqrt(1-param[ijm]^2)*rnorm(nn); }
      else if(cop==2) { zz[,ij]=qcond(runif(nn),z[,jg],param[ijm],df[2]); }
      else if(ivect) { zz[,ij]=qcond(runif(nn),z[,jg],param[ijm]); }
      else 
      { for(i in 1:nn) { zz[i,ij]=qcond(runif(1),z[i,jg],param[ijm]); } } 
      #else print("this copula family is not implemented");
    }
  }
  indsb=0; jg=1; ind=0; 
  for (jsbg in 1:msbgrp)  
  { indsb1=indsb+1; 
    indsb2=indsb+sbgrsize[jsbg];
    indsb =indsb+sbgrsize[jsbg]; 
    for(ijk in indsb1:indsb2)  
    { ijkm=ijk+mgrp+msbgrp;  
      if(cop==1) 
      { zdata[,ijk]=zz[,jsbg]*param[ijkm]+sqrt(1-param[ijkm]^2)*rnorm(nn); }
      else if(cop==2) 
      { zdata[,ijk]=qcond(runif(nn),zz[,jsbg],param[ijkm],df[3]); }
      else if(ivect) { zdata[,ijk]=qcond(runif(nn),zz[,jsbg],param[ijkm]); }
      else 
      { for(i in 1:nn) { zdata[i,ijk]=qcond(runif(1),zz[i,jsbg],param[ijkm]); }}
      #else print("this copula family is not implemented");
    }
  }
  zdata  
}
