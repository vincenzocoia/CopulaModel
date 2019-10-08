
# Code templates by Pavel Krupskii

# Functions for margins of factor copula models

# u1, u2 = vectors of the same length (values in (0,1))
# pcondcop = function for conditional cdf of copula family for factor 1
# param = vector of length 2 
#   (param[i] : dep. par. of the ith copula, i = 1,2) 
#  or 2x2 matrix 
#   (param[i,] : dep. par. of the ith copula, i = 1,2)
#    or one of these variables can be a constant 
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate cdf that is bivariate margin of a 1-factor copula
pfact1cop=function(u1,u2,pcondcop,param,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  if(is.matrix(param)) { par.1=param[1,]; par.2=param[2,]; } 
  else { par.1=param[1]; par.2=param[2]; }          
  m1=length(u1);  m2=length(u2);
  m=max(m1,m2);   cdf=rep(0,m); 
  if(length(u1)==1) u1=rep(u1, m);
  if(length(u2)==1) u2=rep(u2, m);
  for(i in 1:m) cdf[i]=sum(wl*pcondcop(u1[i],xl,par.1)*pcondcop(u2[i],xl,par.2)); 
  cdf
}


# u1, u2 = vectors of the same length (values in (0,1))
# pcondcop1 = function for conditional cdf of copula family for factor 1
# pcondcop2 = function for conditional cdf of copula family for factor 2
# param1 = vector of length 2 for factor 1 
#     (param1[i] : dep. par. of the ith copula, i = 1,2) 
#   or a 2x2 matrix 
#     (param1[i,] : dep. par. of the ith copula, i = 1,2)
#      or one of these variables can be a constant 
# param2 = vector of length 2 for factor 2 or a  2x2 matrix
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output; bivariate cdf that is bivariate margin of a 2-factor copula
pfact2cop=function(u1,u2,pcondcop1,pcondcop2,param1,param2,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes;  
  if(is.matrix(param1)) { par1.1=param1[1,]; par1.2=param1[2,]; }
  else { par1.1=param1[1]; par1.2=param1[2]; }         
  xl0=rep(xl,each=nq); 
  wl0=rep(wl,each=nq)*rep(wl,nq);
  if(is.matrix(param2)) { par2.1=param2[1,]; par2.2=param2[2,]; }
  else { par2.1=param2[1]; par2.2=param2[2]; }
  m1=length(u1);  m2=length(u2);
  m=max(m1,m2);   cdf=rep(0,m); 
  if(length(u1)==1) u1=rep(u1,m);
  if(length(u2)==1) u2=rep(u2,m);
  for(i in 1:m)
  { pc1= pcondcop1(u1[i],xl,par1.1); pc1=rep(pc1,nq);
    pc2= pcondcop1(u2[i],xl,par1.2); pc2=rep(pc2,nq);
    cdf[i]= sum(wl0*pcondcop2(pc1,xl0,par2.1)*pcondcop2(pc2,xl0,par2.2));
  }
  cdf
}

#======================================================================

# precision is up to three digits after the decimal point
# In 1-factor marginal bivariate cdf:
#  param is a vector of length 2, param[i] = dep. par. of the ith linking copula
# or a 2x2 matrix, param[i,] = dep. par of the ith copula
# In 2-factor marginal bivariate cdf:
#  pmatrix is a 2x2 matrix, pmatrix[,i] = dep. par. of ith linking copula
#  or a 2 by 3 matrix, pmatrix[,1:2] = dep. par. of the 1st factor 
#  (e.g. bb1 copula) and pmatrix[,3] = dep. par. of the 2nd factor 

#more cdfs can be added if necessary
# default nq changed from 55 and 45 to 35

pfact1gau=function(u1,u2,param)
{ pfact1cop(u1,u2,pcondbvncop,param,35); }

pfact1frk=function(u1,u2,param)
{ pfact1cop(u1,u2,pcondfrk,param,35); }

#pfact1pla= function(u1,u2,param)
#{ pfact1cop(u1,u2,pcondpla,param,35); }

pfact1gum=function(u1,u2,param)
{ pfact1cop(u1,u2,pcondgum,param,35); }

#pfact1mtcj= function(u1,u2,param)
#{ pfact1cop(u1,u2,pcondmtcj,param,35); }

pfact1bb1=function(u1,u2,param)
{ pfact1cop(u1,u2,pcondbb1,param,35); }

pfact2gau=function(u1,u2,pmatrix)
{ pfact2cop(u1,u2,pcondbvncop,pcondbvncop,pmatrix[,1],pmatrix[,2],35); }

pfact2frk=function(u1,u2,pmatrix)
{ pfact2cop(u1,u2,pcondfrk,pcondfrk,pmatrix[,1],pmatrix[,2],35); }

#pfact2pla= function(u1,u2,pmatrix)
#{ pfact2cop(u1,u2,pcondpla,pcondpla,pmatrix[,1],pmatrix[,2],35); }

pfact2gum=function(u1,u2,pmatrix)
{ pfact2cop(u1,u2,pcondgum,pcondgum,pmatrix[,1],pmatrix[,2],35); }

#pfact2mtcj= function(u1,u2,pmatrix)
#{ pfact2cop(u1,u2,pcondmtcj,pcondmtcj,pmatrix[,1],pmatrix[,2],35); }

pfact2bb1frk=function(u1,u2,pmatrix)
{ pfact2cop(u1,u2,pcondbb1,pcondfrk,pmatrix[,1:2],pmatrix[,3],35); }


#============================================================

# u1, u2 = vectors of the same length (values in (0,1))
#    or one of these variables can be a constant 
# param = vector of length 2 
#    (param[i] : dep. par. of the ith copula, i = 1,2) 
#  or 2x2 matrix 
#    (param[i,] : dep. par. of the ith copula, i = 1,2)
# pcondcop = function for conditional cdf of copula family for factor 1
# dcop = function for copula density
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: conditional cdf C_{2|1}(u2|u1) of bivariate margin of 1-factor copula
pcondfact1=function(u2,u1,pcondcop,dcop,param,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  if(is.matrix(param)) { par.1=param[1,]; par.2=param[2,]; } 
  else { par.1=param[1]; par.2=param[2]; }          
  m1=length(u1);  m2=length(u2);
  m=max(m1,m2);   ccdf=rep(0,m); 
  if(length(u1)==1) u1=rep(u1, m);
  if(length(u2)==1) u2=rep(u2, m);
  for(i in 1:m) ccdf[i]=sum(wl*dcop(u1[i],xl,par.1)*pcondcop(u2[i],xl,par.2))
  ccdf
}

# u1, u2 = vectors of the same length (values in (0,1))
#   or one of these variables can be a constant 
# dcop = function for copula density
# param = vector of length 2 
#    (param[i] : dep. par. of the ith copula, i = 1,2) 
#   or 2x2 matrix 
#    (param[i,] : dep. par. of the ith copula, i = 1,2)
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: copula density of bivariate margin of 1-factor copula
# could replace nq by gl=Gauss-Legendre object
dfact1cop=function(u1,u2,dcop,param,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  if(is.matrix(param)) { par.1=param[1,]; par.2=param[2,]; } 
  else { par.1=param[1]; par.2=param[2]; }          
  m1=length(u1);  m2=length(u2);
  m=max(m1,m2);   pdf=rep(0,m); 
  if(length(u1)==1) u1=rep(u1,m);
  if(length(u2)==1) u2=rep(u2,m);
  for(i in 1:m) pdf[i]=sum(wl*dcop(u1[i],xl,par.1)*dcop(u2[i],xl,par.2))
  pdf
}

# u1, u2 = vectors of the same length (values in (0,1))
#    or one of these variables can be a constant 
# pcondcop1 = function for conditional cdf of copula family for factor 1
# pcondcop2 = function for conditional cdf of copula family for factor 2
# dcop1 = function for copula density with factor1
# dcop2 = function for copula density with factor2
# param1 = vector of length 2 for factor 1 
#     (param1[i] : dep. par. of the ith copula, i = 1,2) 
#   or a 2x2 matrix 
#     (param1[i,] : dep. par. of the ith copula, i = 1,2)
#      or one of these variables can be a constant 
# param2 = vector of length 2 for factor 2 or a  2x2 matrix
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: conditional cdf that is bivariate margin of a 2-factor copula
pcondfact2=function(u2,u1,pcondcop1,pcondcop2,dcop1,dcop2,param1,param2,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes;  
  if(is.matrix(param1)) { par1.1=param1[1,]; par1.2=param1[2,]; }
  else { par1.1=param1[1]; par1.2=param1[2]; }         
  xl0=rep(xl,each=nq); 
  wl0=rep(wl,each=nq)*rep(wl,nq);
  if(is.matrix(param2)) { par2.1=param2[1,]; par2.2=param2[2,]; }
  else { par2.1=param2[1]; par2.2=param2[2]; }
  m1=length(u1);  m2=length(u2);
  m=max(m1, m2);  ccdf=rep(0,m); 
  if(length(u1)==1) u1=rep(u1,m);
  if(length(u2)==1) u2=rep(u2,m);
  for(i in 1:m)
  { pc1=pcondcop1(u1[i],xl,par1.1); pc1=rep(pc1,nq);
    pc2=pcondcop1(u2[i],xl,par1.2); pc2=rep(pc2,nq);
    pdf1=dcop1(u1[i],xl,par1.1);  pdf1=rep(pdf1,nq);
    tem=dcop2(pc1,xl0,par2.1)*pcondcop2(pc2,xl0,par2.2)
    tem=tem*pdf1;
    ccdf[i]=sum(wl0*tem);
  }
  ccdf
}

# u1, u2 = vectors of the same length (values in (0,1))
#    or one of these variables can be a constant 
# pcondcop1 = function for conditional cdf of copula family for factor 1
# dcop1 = function for copula density with factor1
# dcop2 = function for copula density with factor2
# param1 = vector of length 2 for factor 1 
#     (param1[i] : dep. par. of the ith copula, i = 1,2) 
#   or a 2x2 matrix 
#     (param1[i,] : dep. par. of the ith copula, i = 1,2)
#      or one of these variables can be a constant 
# param2 = vector of length 2 for factor 2 or a  2x2 matrix
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: copula density of bivariate margin of 2-factor copula
dfact2cop=function(u1,u2,pcondcop1,dcop1,dcop2,param1,param2,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes;  
  if(is.matrix(param1)) { par1.1=param1[1,]; par1.2=param1[2,]; }
  else { par1.1=param1[1]; par1.2=param1[2]; }         
  xl0=rep(xl,each=nq); 
  wl0=rep(wl,each=nq)*rep(wl,nq);
  if(is.matrix(param2)) { par2.1=param2[1,]; par2.2=param2[2,]; }
  else { par2.1=param2[1]; par2.2=param2[2]; }
  m1=length(u1);  m2=length(u2);
  m=max(m1, m2);  pdf=rep(0,m); 
  if(length(u1)==1) u1=rep(u1,m);
  if(length(u2)==1) u2=rep(u2,m);
  for(i in 1:m)
  { pc1=pcondcop1(u1[i],xl,par1.1); pc1=rep(pc1,nq);
    pc2=pcondcop1(u2[i],xl,par1.2); pc2=rep(pc2,nq);
    pdf1=dcop1(u1[i],xl,par1.1);  pdf1=rep(pdf1,nq);
    pdf2=dcop1(u2[i],xl,par1.2);  pdf2=rep(pdf2,nq);
    tem=dcop2(pc1,xl0,par2.1)*dcop2(pc2,xl0,par2.2)
    tem=tem*pdf1*pdf2;
    pdf[i]=sum(wl0*tem);
  }
  pdf
}

#============================================================


