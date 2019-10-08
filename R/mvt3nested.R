# auxiliary functions for multivariate Student t with 3-nested structure
# used by mvttrifact() in f90mvtfactor.R
# Code written by Pavel Krupskii

# grsize = vector of group sizes 
# sbgrsize = vector of subgroup sizes 
# nestpar = vector of correlations in a nested model,
#   the order: global/common latent V0 and group latent Vg (size ng)
#     then Vg and subgroup latent Vkg (size nsbg) 
#     then Vkg and observed Uikg (size d)
# Output: thii (i=1,2,3) = auxiliary variables used in grad2nestgrad()
#         thic (i=1,2,3) = partial correlations in a tri-factor structure
#  Note: in bifct.R, ng is mgrp, nsbg is msbgrp
nest2condcor= function(nestpar,grsize,sbgrsize)
{ d = sum(grsize);
  ng = length(grsize);
  nsbg = length(sbgrsize);
  th1 = nestpar[1:ng];
  th2 = nestpar[(ng+1):(ng+nsbg)];
  th3 = nestpar[(ng+nsbg+1):(ng+nsbg+d)];
  th11 = rep(0,d);
  th22 = rep(0,d);
  th33 = th3;
  igr = 0;
  for(ig in 1:ng)
  { igr1 = igr+1; igr2 = igr+grsize[ig]; 
    th11[igr1:igr2]=th1[ig];
    igr = igr+grsize[ig];
  }
  isbgr = 0;
  for(isbg in 1:nsbg)
  { isbgr1 = isbgr+1; isbgr2 = isbgr+sbgrsize[isbg]; 
    th22[isbgr1:isbgr2]=th2[isbg];
    isbgr = isbgr+sbgrsize[isbg];
  }
  th23 = th22*th33;
  th1c = th11*th23;                          #condcor; common factor
  th2c = sqrt(th23^2-th1c^2)/sqrt(1-th1c^2); #condcor; group factor given common factor
  th3c = sqrt(th33^2-th23^2)/sqrt(1-th23^2); #condcor; subgroup factor given
                                             #  common and group factor
  list(th11=th11,th22=th22,th33=th33,th1c=th1c,th2c=th2c,th3c=th3c)
}

# th11, th22, th33, th1c, th2c, th3c = output from nest2condcor()
# lgrad = vector of derivatives (loglik wrt conditional correlations in a
#       tri-factor structure)
# grsize = vector of group sizes 
# sbgrsize = vector of subgroup sizes 
# Output: nestgrad = vector of derivatives 
#    (of loglik with respect to (correlation) parameters  in a 3-nested model) 
grad2nestgrad= function(th11,th22,th33,th1c,th2c,th3c,lgrad,grsize,sbgrsize)
{ d = sum(grsize);
  ng = length(grsize);
  nsbg = length(sbgrsize);
  th23 = th22*th33;
  tem1c = th1c/(1-th1c^2);
  tem11 = th11/(1-th11^2);
  tem22 = th22/(1-th22^2); 
  tem23 = th23/(1-th23^2); 
  lgrad1=lgrad[1:d];
  lgrad2=lgrad[(d+1):(2*d)];
  lgrad3=lgrad[(2*d+1):(3*d)];
  grad1th1 = lgrad1*th33*th22+lgrad2*th2c*(-tem11+th23*tem1c);
  grad1th2 = lgrad1*th33*th11+lgrad2*th2c*(1/th22+th11*th33*tem1c)+lgrad3*th3c*(-tem22+th33*tem23);
  grad1th3 = lgrad1*th22*th11+lgrad2*th2c*(1/th33+th11*th22*tem1c)+lgrad3*th3c*(1/th33+th22*tem23);
  nestgrad = rep(0,ng+nsbg+d);
  igr = 0;
  for(ig in 1:ng)
  { igr1 = igr+1; igr2 = igr+grsize[ig];
    nestgrad[ig] = sum(grad1th1[igr1:igr2]);
    igr = igr + grsize[ig]; 
  }
  isbgr = 0;
  for(isbg in 1:nsbg)
  { isbgr1 = isbgr+1; isbgr2 = isbgr+sbgrsize[isbg];
    nestgrad[isbg+ng] = sum(grad1th2[isbgr1:isbgr2]);
    isbgr = isbgr + sbgrsize[isbg];
  }
  nestgrad[(ng+nsbg+1):(ng+nsbg+d)] = grad1th3;
  nestgrad
}
