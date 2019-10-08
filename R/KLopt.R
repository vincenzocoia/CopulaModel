# Functions to find the parameter on bivariate copula family 2
# that is closest in Kullback-Leibler divergence to bivariate copula family 1
# with a specified parameter.

# KL divergence, true family is dcop1 with parameter par1;
# computations are based on Gauss-Legendre quadrature

# Given one family and parameter,
# use KL divergence to find closest member of another family
# report tau, beta, rhoS, rhoN for the two families 

# KL divergence family c2 relative to c1 (true) : 
# when parameter of c1 is par1  \int c1 * log c1/c2 = Delta
# par2 = copula parameter for family 2
# dcop2 = function name of econd bivariate copula density
# dcop1 = function name of first bivariate copula density
# par1 = copula parameter for family 1
# parlb = lower bound for copula parameter for family 2
# gl = Gauss-Legendre quadrature object with $nodes and $weights
# Output: KL divergence family c2 relative to c1  
KL12gl=function(par2,dcop2,par1,dcop1,parlb=0,gl)
{ if(any(par2<=parlb)) { return(1.e10); }
  Del=0.; 
  xl=gl$nodes
  wl=gl$weights
  nq=length(xl)
  for(iq in 1:nq)
  { x1=xl[iq]; w1=wl[iq];
    for(jq in 1:nq)
    { tem1=dcop1(x1,xl[jq],par1)
      tem2=dcop2(x1,xl[jq],par2)
      Del=Del+w1*wl[jq]*tem1*log(tem1/tem2);
    }
  }
  Del 
}


# Find parameter of cop2 closest to par1 of cop1 (true) in KL divergence
# dcop1 = function name of first bivariate copula density
# dcop2 = function name of econd bivariate copula density
# par1 = copula parameter for family 1
# par2 = copula parameter for family 2
# name1 = name of first bivariate copula density
# name2 = name of second bivariate copula density
# parlb = lower bound for copula parameter for family 2
# gl = Gauss-Legendre quadrature object with $nodes and $weights
# pcop1 = function name of first bivariate copula cdf
# pcop2 = function name of second bivariate copula cdf
# ccdf1a = C_{2|1} for pcop1
# ccdf1b = C_{1|2} for pcop1
# ccdf2a = C_{2|1} for pcop2
# ccdf2b = C_{1|2} for pcop2
# prlevel = print.level for nlm()
# Output: parameter from cop2 matching minimum Kl divergence from cop1
KLoptgl=function(dcop1,dcop2,par1,par2,name1,name2,par2lb=0,gl,pcop1,pcop2,
   ccdf1a,ccdf1b,ccdf2a,ccdf2b,prlevel=0)
{ 
  out=nlm(KL12gl,p=par2,dcop2=dcop2,par1=par1,dcop1=dcop1,parlb=par2lb,
    gl=gl,print.level=prlevel)
  nq=length(gl$weights)
  cat("nq=", nq, " KLdiv=", out$min,"\n")
  par2=out$estimate
  if(name1=="bvn")
  { be1=bvn.cpar2b(par1); tau1=bvn.cpar2tau(par1); rhoS1=bvn.cpar2rhoS(par1);
    rhoN1=par1;
  }
  else
  { be1=4*pcop1(.5,.5,par1)-1
    # allow asymmetric pcond21!=pcond12
    tau1=ktau(par1,pcond21=ccdf1a,pcond12=ccdf1b,zero=0)
    rhoS1=rhoS(par1,cop=pcop1,zero=0)
    # is next correct if copula is asymmetric? looks OK from code, check theory
    rhoN1=rhoN(par1,pcond=ccdf1a) 
  }
  if(name2=="bvn")
  { be2=bvn.cpar2b(par2); tau2=bvn.cpar2tau(par2); rhoS2=bvn.cpar2rhoS(par2);
    rhoN2=par2;
  }
  else
  { be2=4*pcop2(.5,.5,par2)-1
    # allow asymmetric pcond21!=pcond12
    tau2=ktau(par2,pcond21=ccdf2a,pcond12=ccdf2b,zero=0)
    rhoS2=rhoS(par2,cop=pcop2,zero=0)
    rhoN2=rhoN(par2,pcond=ccdf2a) 
  }
  if(prlevel>0)
  { cat(name1, "par1=", par1, "(be,tau,rS,rN)=", be1,tau1,rhoS1,rhoN1, "\n")
    cat(name2, "par2=", par2, "(be,tau,rS,rN)=", be2,tau2,rhoS2,rhoN2, "\n")
    cat("diff :", be1-be2,tau1-tau2,rhoS1-rhoS2,rhoN1-rhoN2,"\n")
  }
  if(prlevel>1) cat("============================================================\n\n")
  # more output
  list(cpar2=par2,depm1=c(be1,tau1,rhoS1,rhoN1),depm2=c(be2,tau2,rhoS2,rhoN2))
}


