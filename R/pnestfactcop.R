# Functions for bivariate margins of nested factor copulas
# Code written by Pavel Krupskii

# 2-nested factor copula cdf
# u1, u2 = vectors of the same length; 
# dcop1 = function with copula pdf linking V0 and Vj
#  (common latent and group latent)
# pcondcop2 = function with conditional copula cdf linking Vj and Uij
# param1 = a vector of length 2 or a 2x2 matrix with parameters for dcop1
# param2 = a vector of length 2 or a 2x2 matrix with parameters for pcondcop2
# nq = number of quadrature points for Gauss-Legendre
# Output: bivariate cdf that is bivariate margin of a nested factor copula
pnest2cop=function(u1,u2,dcop1,pcondcop2,param1,param2,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;
  if(is.matrix(param1)) { par1.1=param1[1,]; par1.2=param1[2,]; }
  else { par1.1=param1[1]; par1.2=param1[2]; }         
  if(is.matrix(param2)) { par2.1=param2[1,]; par2.2=param2[2,]; }
  else { par2.1=param2[1]; par2.2=param2[2]; }         
  l1=length(u1);  l2=length(u2);
  l=max(l1,l2);   cdf=rep(0,l); 
  if(length(u1)==1)  u1=rep(u1,l);
  if(length(u2)==1)  u2=rep(u2,l);
  int0=rep(0,l);
  for(i in 1:l)
  { for (iq in 1:nq)
    { int1= sum(wl*dcop1(xl,xl[iq],par1.1)*pcondcop2(u1[i],xl,par2.1));
      int2= sum(wl*dcop1(xl,xl[iq],par1.2)*pcondcop2(u2[i],xl,par2.2));
      int0[i]= int0[i] + wl[iq]*int1*int2;
    }
  }
  int0
}

############################################################################

# u1 = vector of values in (0,1)
# u2 = vector of values in (0,1), same length as u1
# param = parameter matrix:
#  column 1 has parameters for global/common latent, column 2 
#  (and column 3 for pnest2tbb1,pnest2gumbb1) has parameters for group latent

pnest2frk=function(u1,u2,param)
{ f= pnest2cop(u1,u2,dfrk,pcondfrk,param[,1],param[,2],35);
  f= f*(f<=1)+(f>1)
  f
}

pnest2gum=function(u1,u2,param)
{ f= pnest2cop(u1,u2,dgum,pcondgum,param[,1],param[,2],35);
  f= f*(f<=1)+(f>1)
  f
}

pnest2t=function(u1,u2,param,df)
{ f= pnest2cop(u1,u2,dbvtcop,pcondbvtcop,cbind(param[,1],df[1]),
     cbind(param[,2],df[2]),35);
  f= f*(f<=1)+(f>1)
  f
}

pnest2tgum=function(u1,u2,param,df)
{ f= pnest2cop(u1,u2,dbvtcop,pcondgum,cbind(param[,1],df[1]),param[,2],35);
  f= f*(f<=1)+(f>1)
  f
}

pnest2tbb1=function(u1,u2,param,df)
{ f= pnest2cop(u1,u2,dbvtcop,pcondbb1,cbind(param[,1],df[1]),
     cbind(param[,2],param[,3]),35);
  f= f*(f<=1)+(f>1)
  f
}

pnest2gumbb1=function(u1,u2,param)
{ f= pnest2cop(u1,u2,dgum,pcondbb1,param[,1],cbind(param[,2],param[,3]),35);
  f= f*(f<=1)+(f>1)
  f
}

