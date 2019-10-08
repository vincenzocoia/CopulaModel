# Functions for Spearman rho for bivariate margins of factor copula
# Code written by Pavel Krupskii

# Spearman's rho for a given bivariate copula cdf
#   via Gauss-Legendre quadrature, for use with margins of factor copulas
# param = a vector of length npar or a npar/2 by 2 matrix of dep. parameters
# pcop = function for copula cdf
# nq = number of qudrature points for Gauss-Legendre quadrature
# Output: Spearman's rhos
cop2srho=function(param,pcop,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;
  xl01=rep(xl,each=nq);
  xl02=rep(xl,times=nq); 
  wl0= rep(wl,each=nq)*rep(wl,nq);
  rho= sum(wl0*pcop(xl01,xl02,param));
  rho=12*rho-3
  rho
}

# running time with nq = 25 and Gumbel copula:
#  ~0.02 sec
# running time with nq = 35 and Gumbel copula:
#  ~0.02 sec
# nq = 65 is recommended if t copula is used (time is ~0.06-0.09sec)
# for other copulas with nq = 35 the accuracy is OK (+-0.001)

# For bivariate margin of nested copula
# param1 = vector/matrix of copula parameters for linking copulas to global latent 
# param2 = vector of copula parameters for linking copulas to group latent 
# dcop1 = function for pdf of copula family for global latent variable
# pcondcop2 = function for conditional cdf of copula family for group latent variable 
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: Spearman's rho for biv margin (2 different groups) of nested copula
rho2nestfact=function(param1,param2,dcop1,pcondcop2,nq)
{ gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;
  if(is.matrix(param1)) { par1.1=param1[1,]; par1.2=param1[2,]; }
  else { par1.1=param1[1]; par1.2=param1[2]; }         
  if(is.matrix(param2)) { par2.1=param2[1,]; par2.2=param2[2,]; }
  else { par2.1=param2[1]; par2.2=param2[2]; }
  
  int0=0; 
  pcint1=rep(0,nq); pcint2=rep(0,nq);
  for (iq in 1:nq)
  { pcint1[iq]= sum(wl*pcondcop2(xl,xl[iq],par2.1));
    pcint2[iq]= sum(wl*pcondcop2(xl,xl[iq],par2.2));
  }
  for (iq in 1:nq)
  { int1= sum(wl*dcop1(xl,xl[iq],par1.1)*pcint1);
    int2= sum(wl*dcop1(xl,xl[iq],par1.2)*pcint2);
    int0= int0 + wl[iq]*int1*int2;
  }
  int0=12*int0-3
  int0
}

# much faster version of rho1fact (old version not exported)
# running time with nq = 25 and Gumbel copula:
# rho1fact: ~0.17 sec; srho1fact: ~0.02 sec
# running time with nq = 35 and Gumbel copula:
# rho1fact: ~0.33 sec; srho1fact: ~0.02 sec

# For bivariate margin of 1-factor copula
# param1 = vector/matrix of copula parameters for linking copulas to factor 1 
# pcondcop = function for conditional cdf of copula family for factor 1 
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: Spearman's rho for biv margin of 1-factor copula
srho1fact=function(param1,pcondcop,nq)
{ if(is.matrix(param1)) d=nrow(param1) else d=length(param1);
  rho=matrix(0,ncol=d,nrow=d);
  gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;
  #j1=1; j2=2;
  for(j1 in 1:d)
  { for(j2 in 1:(j1-1))
    { 
      if(is.matrix(param1)) { par.1=param1[j1,]; par.2=param1[j2,]; } 
      else { par.1=param1[j1]; par.2=param1[j2]; }
      pcint1=rep(0,nq);  pcint2=rep(0,nq);
      for(iq in 1:nq)
      { pcint1[iq]= sum(wl*pcondcop(xl,xl[iq],par.1));
        pcint2[iq]= sum(wl*pcondcop(xl,xl[iq],par.2));
      }
      rho[j1,j2]=sum(wl*pcint1*pcint2);
      rho[j2,j1]=rho[j1,j2];
    }
  }
  rho=12*rho-3;
  diag(rho)=1;
  rho
}

#much faster version of rho2fact
#running time with nq = 25 and Gumbel copula:
#rho2fact: ~1.74 sec; srho2fact: ~0.26 sec
#running time with nq = 35 and Gumbel copula:
#rho2fact: ~6.14 sec; srho2fact: ~0.58 se

# For bivariate margin of 1-factor copula
# param1 = vector/matrix of copula parameters for linking copulas to factor 1 
# param2 = vector of copula parameters for linking copulas to factor 2
# pcondcop1 = function for conditional cdf of copula family for factor 1 
# pcondcop2 = function for conditional cdf of copula family for factor 2
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: Spearman's rho for biv margin of 2-factor copula
srho2fact=function(param1,param2,pcondcop1,pcondcop2,nq)
{ if(is.matrix(param1)) d=nrow(param1) else d=length(param1);
  rho= matrix(0,ncol=d,nrow=d);
  gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;
  wl0=rep(wl,each=nq)*rep(wl,nq);
  xl1=rep(xl,nq);
  xl2=rep(xl,each=nq);
  #j1 = 1; j2 = 2;
  for(j1 in 1:d)
  { for(j2 in 1:(j1-1))
    { 
      if(is.matrix(param1)) { par1.1= param1[j1,]; par1.2= param1[j2,]; } 
      else { par1.1= param1[j1]; par1.2= param1[j2]; } 
      if(is.matrix(param2)) { par2.1= param2[j1,]; par2.2= param2[j2,]; }
      else { par2.1= param2[j1]; par2.2= param2[j2]; }
      pcint1=0;  pcint2=0;  
      for(iq1 in 1:nq)
      { for(iq2 in 1:nq)
        { pc1= pcondcop1(xl,xl[iq1],par1.1)
          pc2= pcondcop1(xl,xl[iq1],par1.2)
          pcint1= c(pcint1, sum(wl*pcondcop2(pc1,xl[iq2],par2.1)));
          pcint2= c(pcint2, sum(wl*pcondcop2(pc2,xl[iq2],par2.2)));
        }
      }
      rho[j1,j2]= sum(wl0*pcint1[-1]*pcint2[-1]);
      rho[j2,j1]=rho[j1,j2];
    }
  }
  rho=12*rho-3;
  diag(rho)=1;
  rho
}


# Matrix of bivariate dependence measures for factor copula models
# param = parameter vector of factor copulas
# param is dx1 for 1-factor copulas with 1-parameter bivariate linking copulas
# param is dx2 for 1-factor copulas with 2-parameter bivariate linking copulas
# pcondcop = function for conditional cdf of the linking copula family 
# nq = number of quadrature points for Gauss-Legendre quadrature
# output:
#  dxd matrix of Spearman's rho
rho1fact= function(param,pcondcop,nq)
{ if(is.matrix(param)) d=nrow(param) else d=length(param);
  rho=matrix(0,ncol=d,nrow=d);
  gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;
  for(j1 in 1:d)
  { for(j2 in 1:(j1-1))
    { for(i1 in 1:nq)
      { for(i2 in 1:nq)
        { if(is.matrix(param)) { par.1= param[j1,]; par.2= param[j2,]; } 
          else { par.1=param[j1]; par.2=param[j2]; } 
          rho[j1,j2]= rho[j1,j2] + wl[i1]*wl[i2]*
                 sum(wl*pcondcop(xl[i1],xl,par.1)*pcondcop(xl[i2],xl,par.2));
          rho[j2,j1]=rho[j1,j2];
        }
      }
    }
  }
  rho=12*rho-3;
  diag(rho)=1;
  rho
}

# param = parameter vector of factor copulas
# param1 = dx1 for 1-parameter bivariate linking copulas for factor1
# param1 = dx2 for 2-parameter bivariate linking copulas for factor1
# param2 = dx1 for 1-parameter bivariate linking copulas for factor2
# pcondcop1 = function for the conditional cdf of copula used for factor1
# pcondcop2 = function for the conditional cdf of copula used for factor2
# nq = number of quadrature points for Gauss-Legendre quadrature
# output: dxd matrix of Spearman's rho 
rho2fact= function(param1, param2, pcondcop1, pcondcop2, nq)
{ if(is.matrix(param1)) d=nrow(param1) else d=length(param1);
  rho=matrix(0,ncol=d,nrow=d);
  gl=gausslegendre(nq);
  wl=gl$weights;
  xl=gl$nodes;

  for(j1 in 1:d)
  { for(j2 in 1:(j1-1))
    { for(i1 in 1:nq)
      { for(i2 in 1:nq)
        { if(is.matrix(param1)) { par1.1=param1[j1,]; par1.2=param1[j2,]; } 
          else { par1.1=param1[j1]; par1.2=param1[j2]; }          
          pc1= pcondcop1(xl[i1],xl,par1.1); pc1=rep(pc1,nq);
          pc2= pcondcop1(xl[i2],xl,par1.2); pc2=rep(pc2,nq);
          xl0= rep(xl,each=nq); 
          wl0= rep(wl,each=nq)*rep(wl,nq);
          if(is.matrix(param2)) { par2.1=param2[j1,]; par2.2=param2[j2,]; }
          else { par2.1=param2[j1]; par2.2=param2[j2]; }
          rho[j1,j2]= rho[j1,j2] + wl[i1]*wl[i2]*
            sum(wl0*pcondcop2(pc1,xl0,par2.1)*pcondcop2(pc2,xl0,par2.2));
          rho[j2,j1]=rho[j1,j2];
        }
      }
    }
  }
  rho=12*rho-3;
  diag(rho)=1;
  rho
}

# fct = number of factors
# param = d x fct matrix of conditional or partial correlations
# Output: matrix of Spearman's rho for a MVN factor model
rhomvn=function(fct,param)
{ uncon=param;
  if(fct>1)
  { for(j in 2:fct)
    { for(k in 1:(j-1)) uncon[,j]=uncon[,j]*sqrt(1-param[,k]^2); }
  }
  sigma=uncon%*%t(uncon);
  diag(sigma)=1;
  srho=(6/pi)*asin(sigma/2);
  srho
}

