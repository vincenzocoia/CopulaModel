# Functions for tail-weighted dependence measures,
# Model-based (factor copulas) and empirical.
# Code written by Pavel Krupskii

# numerical computation of the tail-weighted measure of dependence
# version 2, truncation after ranking
# very stable
# running time: < 1sec for 1-factor model and 3-4 sec for 2-factor model

# pcop = pfact1cop, where cop = normal/Gaussian, frk, pla etc.
# pcop = pfact2cop, where cop = normal/Gaussian, frk,...,bb1frk.
# see factorcopcdf.R for details

#============================================================

# pcop = function for copula cdf
# param = the dependence parameter of the copula pcop
# power = power to use for tail-weighted dependence measure, good choice is 6
# nq = number of quadrature points for Gauss-Legendre quadrature
# Output: lower and upper tail-weighted dependence measures
twdm=function(pcop,param,power,nq,tscore=F)
{ # reflected or survival copula
  pcopr= function(u1,u2,param)
  { cdf= -1+u1+u2+pcop(1-u1,1-u2,param);
    cdf
  }
  ql=0.5;
  tql=ql;
  gl=gausslegendre(nq);
  wl=ql*gl$weights; xl=ql*gl$nodes;
  tl=xl;
  if(tscore==T && param[2]<200)  
  { tl=qt(xl,param[2]); tql=qt(ql,param[2]); }
  if(tscore==T && param[2]>=200) 
  { tl=qnorm(xl); tql=qnorm(ql); param=param[1]; }
  
  den=ql^2*pcop(tql,tql,param);   
  denr=ql^2*pcopr(tql,tql,param); 
  if(tscore==T) { denr=den; }  
  xl1=rep(xl,each=nq); xl2=rep(xl,nq);
  tl1=rep(tl,each=nq); tl2=rep(tl,nq);
  xl0=((1-xl1/ql)*(1-xl2/ql))^(power-1);
  wl1=rep(wl,each=nq); wl2=rep(wl,nq);
  wl0=wl1*wl2;
  e12= power^2*sum(wl0*xl0*pcop(tl1,tl2,param))/den;
  e1.1= ql*power*sum(wl*(1-xl/ql)^(power-1)*pcop(tl,tql,param))/den;
  e2.1= 2*ql*power*sum(wl*(1-xl/ql)^(2*power-1)*pcop(tl,tql,param))/den;
  e1.2= ql*power*sum(wl*(1-xl/ql)^(power-1)*pcop(tql,tl,param))/den; 
  e2.2= 2*ql*power*sum(wl*(1-xl/ql)^(2*power-1)*pcop(tql,tl,param))/den;
  v1=e2.1-e1.1^2;
  v2=e2.2-e1.2^2;
  if(tscore==F)
  { e12r= power^2*(-ql^2/power/(power+1)+sum(wl0*xl0*pcop(1-tl1,1-tl2,param)))/denr;
    e1.1r= ql*power*(-ql^2/(power+1)+sum(wl*(1-xl/ql)^(power-1)*pcop(1-tl,1-tql,param)))/denr;
    e2.1r = 2*ql*power*(-ql^2/(2*power+1)+sum(wl*(1-xl/ql)^(2*power-1)*pcop(1-tl,1-tql,param)))/denr;
    e1.2r= ql*power*(-ql^2/(power+1)+sum(wl*(1-xl/ql)^(power-1)*pcop(1-tql,1-tl,param)))/denr;
    e2.2r= 2*ql*power*(-ql^2/(2*power+1)+sum(wl*(1-xl/ql)^(2*power-1)*pcop(1-tql,1-tl,param)))/denr;
    v1r=e2.1r-e1.1r^2;
    v2r=e2.2r-e1.2r^2;
    outr=(e12r-e1.1r*e1.2r)/sqrt(v1r*v2r);
  }  
  
  out=(e12-e1.1*e1.2)/sqrt(v1*v2);
  if(tscore==T) outr=out;  
  out0=c(out,outr);
  out0
}

# empirical version of twdm
# data = nxd data set
# power = power to use for tail-weighted dependence measure, good choice is 6
# Output: twdm in the lower tail for each pair of columns in data
twdm.emp=function(data,power)
{ d = ncol(data);  n = nrow(data);
  ltwdm = matrix(0,ncol=d,nrow=d);
  for(j1 in 2:d)
  { for(j2 in 1:(j1-1))
    { tem1 = data[,j1];
      tem2 = data[,j2];
      # ranking or uniform scores
      rnk1 = (rank(tem1)-0.5)/n;
      rnk2 = (rank(tem2)-0.5)/n;
      # truncation
      ind = (rnk1 < 0.5 & rnk2 < 0.5);
      ltwdm[j1,j2] = cor((1-2*rnk1[ind])^power,(1-2*rnk2[ind])^power);
      ltwdm[j2,j1] = ltwdm[j1,j2];
    }
  }
  diag(ltwdm) = 1;
  ltwdm
}

# empirical version of twdm: vectorized
# data = nxd data set
# power = power to use for tail-weighted dependence measure, good choice is 6
# Output: twdm in the lower tail for each pair of columns in data
twdm.emp.vec=function(data,power)
{ n = nrow(data);
  rk= apply(data,2,rank);
  rk= (rk-0.5)/n;
  l.rk = rk*(rk<0.5);
  l.rk[l.rk==0]=NA;
  out = cor((1-2*l.rk)^power,use="pairwise.complete.obs")
  diag(out)=1
  out
}

#############################################################################

# Much faster version of twdm() for nested factor copulas
# and even nq = 70 if Student copula is used (bad approximation of t-quantiles)
# running times with nq = 55 and diff. copulas for twd2:
# Frank: ~14.75sec; Gumbel: ~31.95sec; Gumbel+BB1: ~29.85sec
# running times with nq = 55 and diff. copulas for twdm.nestcop:
# Frank: ~0.05sec; Gumbel: ~0.1sec; Gumbel+BB1: ~0.09sec; Student: ~0.2sec.
# slight difference between two methods comes from pnest2cop functions, used in the old version
# these functions use numerical integration with nq = 35 (fixed)
# in the new function nq is set equal to 55 everywhere
# with not very large dep. pars, the diff. between two methods is very small (0.001 or less)

# dcop = function for density of copula linking latent variables
# pcondcop = function for conditional cdf linking latent and observed variables
# param1 = dependence parameter of dcop for two group variables
# param2 = dependence parameter of pcondcop for two observed variables
# power = power to use for tail-weighted dependence measure, good choice is 6
# nq = number of quadrature points for Gauss-Legendre
#     55 is recommended for tail dependent copulas for higher accuracy
# Output: vector of length 2 with
#   lower and upper tail-weighted dependence measure values
twdmnestcop=function(dcop,pcondcop,param1,param2,power,nq)
{ 
  if(is.matrix(param1)) { par1.1= param1[1,]; par1.2= param1[2,]; }
  else { par1.1= param1[1]; par1.2= param1[2]; }         
  if(is.matrix(param2)) { par2.1= param2[1,]; par2.2= param2[2,]; }
  else { par2.1= param2[1]; par2.2= param2[2]; }

  dcopr= function(u1,u2,param)
  { return(dcop(1-u1,1-u2,param)) }

  pcondcopr= function(u1,u2,param)
  { return(1-pcondcop(1-u1,1-u2,param)) }

  symm=F;
  tst=seq(0.01,0.99,0.01);
  dif1= max(abs(dcop(tst,tst,par1.1) - dcopr(tst,tst,par1.1)));
  dif2= max(abs(pcondcop(tst,tst,par2.1) - pcondcopr(tst,tst,par2.1)));
  dif=max(dif1,dif2);
  if(dif<1e-5) symm = T; #numerical check for reflection symmetric dependence 
  ql=0.5;
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes;
  wl2=ql*gl$weights; xl2=ql*gl$nodes;
  den= ql^2*pnest2cop(ql,ql,dcop,pcondcop,param1,param2,35);   
  denr=den; #C(0.5,0.5) = C_R(0.5,0.5)
  int0=0; int1ql=0; int21ql=0; int2ql=0; int22ql=0;
  pcint1=rep(0,nq); pcint2=rep(0,nq);
  pcint21=rep(0,nq); pcint22=rep(0,nq);
  int0r=0; int1qlr=0; int21qlr=0; int2qlr=0; int22qlr=0;
  pcint1r=rep(0,nq); pcint2r=rep(0,nq);
  pcint21r=rep(0,nq); pcint22r=rep(0,nq);
  for (iq in 1:nq)
  { pcint1[iq]= sum(wl2*(1-xl2/ql)^(power-1)*pcondcop(xl2,xl[iq],par2.1));
    pcint2[iq]= sum(wl2*(1-xl2/ql)^(power-1)*pcondcop(xl2,xl[iq],par2.2));
    pcint21[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop(xl2,xl[iq],par2.1));
    pcint22[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop(xl2,xl[iq],par2.2));
    #other tail, by reflection
    if(symm==F)
    { pcint1r[iq]= sum(wl2*(1-xl2/ql)^(power-1)*pcondcop(1-xl2,xl[iq],par2.1));
      pcint2r[iq]= sum(wl2*(1-xl2/ql)^(power-1)*pcondcop(1-xl2,xl[iq],par2.2));
      pcint21r[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop(1-xl2,xl[iq],par2.1));
      pcint22r[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop(1-xl2,xl[iq],par2.2));
    }
  }

  #. - integrate pcondcop wrt the 1st argument
  #1 - integrate pcondcop*(1-xl/ql)^(power-1) -...-
  #2 - integrate pcondcop*(1-xl/ql)^(2*power-1) -...-
  #! - integrate pcondcop(ql,xl,par)
  for (iq in 1:nq)
  { int1  = sum(wl*dcop(xl,xl[iq],par1.1)*pcint1);                 #1 1 0 0
    qlint1= sum(wl*dcop(xl,xl[iq],par1.1)*pcondcop(ql,xl,par2.1)); #! ! 0 0 
    int21 = sum(wl*dcop(xl,xl[iq],par1.1)*pcint21);                #2 2 0 0  
    int2  = sum(wl*dcop(xl,xl[iq],par1.2)*pcint2);                 #0 0 1 1
    qlint2= sum(wl*dcop(xl,xl[iq],par1.2)*pcondcop(ql,xl,par2.2)); #0 0 ! !
    int22 = sum(wl*dcop(xl,xl[iq],par1.2)*pcint22);                #0 0 2 2
    int0  = int0 + wl[iq]*int1*int2;
    int1ql= int1ql + wl[iq]*int1*qlint2;
    int21ql= int21ql + wl[iq]*int21*qlint2;
    int2ql = int2ql + wl[iq]*int2*qlint1;
    int22ql= int22ql + wl[iq]*int22*qlint1;

    #other tail, by reflection
    if(symm==F)
    { int1r  = sum(wl*dcop(xl,xl[iq],par1.1)*pcint1r);                
      qlint1r= sum(wl*dcop(xl,xl[iq],par1.1)*pcondcop(1-ql,xl,par2.1)); 
      int21r = sum(wl*dcop(xl,xl[iq],par1.1)*pcint21r);            
      int2r  = sum(wl*dcop(xl,xl[iq],par1.2)*pcint2r);                 
      qlint2r= sum(wl*dcop(xl,xl[iq],par1.2)*pcondcop(1-ql,xl,par2.2));
      int22r = sum(wl*dcop(xl,xl[iq],par1.2)*pcint22r);            
      int0r  = int0r + wl[iq]*(int1r)*(int2r);
      int1qlr= int1qlr + wl[iq]*int1r*qlint2r;
      int21qlr= int21qlr + wl[iq]*int21r*qlint2r;
      int2qlr = int2qlr + wl[iq]*int2r*qlint1r;
      int22qlr= int22qlr + wl[iq]*int22r*qlint1r;
    }
  } 

  #for the other tail straightforward computation is used, without using dcopr and pcondcopr; 
  #to increase stability of the integration if there is a strong dependence in the tails
  if(symm==F)
  { int0r= int0r - 0.25/power/(power+1)
    int1qlr= int1qlr - 0.25/(power+1)
    int21qlr= int21qlr - 0.25/(2*power+1)
    int2qlr = int2qlr - 0.25/(power+1)
    int22qlr= int22qlr - 0.25/(2*power+1)
    e12r=power^2*int0r/denr; 
    e1.1r=ql*power*int1qlr/denr;
    e2.1r=2*ql*power*int21qlr/denr;
    e1.2r=ql*power*int2qlr/denr;
    e2.2r=2*ql*power*int22qlr/denr; 
    v1r=e2.1r-e1.1r^2;
    v2r=e2.2r-e1.2r^2;
  }

  e12=power^2*int0/den; 
  e1.1=ql*power*int1ql/den;
  e2.1=2*ql*power*int21ql/den;
  e1.2=ql*power*int2ql/den;
  e2.2=2*ql*power*int22ql/den; 
  v1=e2.1-e1.1^2;
  v2=e2.2-e1.2^2;
  out =(e12-e1.1*e1.2)/sqrt(v1*v2);
  outr=out;
  if(symm==F) outr=(e12r-e1.1r*e1.2r)/sqrt(v1r*v2r);
  out0=c(out,outr);
  out0
}

# Much faster version of twdm() for 1 factor copulas
# running times with nq = 35 and diff. copulas for twdm2:
# BB1: ~0.22sec.; Gumbel: ~0.24sec.; Frank: ~0.11sec.
# running times with nq = 35 and diff. copulas for twdm1factcop:
# BB1: ~0.03sec.; Gumbel: ~0.04sec.; Frank: ~0.02sec.; t: ~0.07sec.
# difference between two methods is < 0.0005 for all copulas

# pcondcop = conditional copula cdf linking V1 and 2 observed variables
# param = dependence parameter in pcondcop for the 2 variables
# power = power to use for tail-weighted dependence measure, good choice is 6
# nq = number of quadrature points for Gauss-Legendre
#    nq = 35 gives very good accuracy for all tail dependent copulas 
# Output: vector of length 2 with
#   lower and upper tail-weighted dependence measure values
twdm1factcop=function(pcondcop,param,power,nq)
{ 
  if(is.matrix(param)) { par1=param[1,]; par2=param[2,]; }
  else { par1=param[1]; par2=param[2]; }   

  ql=0.5;
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes;
  wl2=ql*gl$weights; xl2=ql*gl$nodes;
   
  den=ql^2*pfact1cop(ql,ql,pcondcop,param,35);   
  denr=den; #C(0.5,0.5) = C_R(0.5,0.5)      

  pcondcopr= function(u1,u2,param)
  { return(1-pcondcop(1-u1,1-u2,param)) }

  symm=F;
  tst=seq(0.01,0.99,0.01);
  dif1= max(abs(pcondcop(tst,tst,par1) - pcondcopr(tst,tst,par1)));
  dif2= max(abs(pcondcop(tst,tst,par2) - pcondcopr(tst,tst,par2)));
  dif=max(dif1,dif2);

  if(dif<1e-5) symm=T;  #numerical check for reflection symmetric dependence 
  int0=0; int1ql=0; int21ql=0; int2ql=0; int22ql=0;
  pcint1=rep(0,nq); pcint2=rep(0,nq);
  pcint21=rep(0,nq); pcint22=rep(0,nq);
  
  int0r=0; int1qlr=0; int21qlr=0; int2qlr=0; int22qlr=0;
  pcint1r=rep(0,nq); pcint2r=rep(0,nq);
  pcint21r=rep(0,nq); pcint22r=rep(0,nq);

  for (iq in 1:nq)
  {
    pcint1[iq] = sum(wl2*(1-xl2/ql)^(power-1)*pcondcop(xl2,xl[iq],par1));
    pcint2[iq] = sum(wl2*(1-xl2/ql)^(power-1)*pcondcop(xl2,xl[iq],par2));
    pcint21[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop(xl2,xl[iq],par1));
    pcint22[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop(xl2,xl[iq],par2));

    #other tail, by reflection
    if(symm==F)
    { pcint1r[iq] = sum(wl2*(1-xl2/ql)^(power-1)*pcondcopr(xl2,xl[iq],par1));
      pcint2r[iq] = sum(wl2*(1-xl2/ql)^(power-1)*pcondcopr(xl2,xl[iq],par2));
      pcint21r[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcopr(xl2,xl[iq],par1));
      pcint22r[iq]= sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcopr(xl2,xl[iq],par2));
    }
  } 
  int0 = sum(wl*pcint1*pcint2);
  int1ql = sum(wl*pcint1*pcondcop(ql,xl,par2));
  int21ql= sum(wl*pcint21*pcondcop(ql,xl,par2));
  int2ql = sum(wl*pcint2*pcondcop(ql,xl,par1));
  int22ql= sum(wl*pcint22*pcondcop(ql,xl,par1));
    
  #other tail, by reflection
  if(symm==F)
  { int0r = sum(wl*pcint1r*pcint2r);
    int1qlr = sum(wl*pcint1r*pcondcopr(ql,xl,par2));
    int21qlr= sum(wl*pcint21r*pcondcopr(ql,xl,par2));
    int2qlr = sum(wl*pcint2r*pcondcopr(ql,xl,par1));
    int22qlr= sum(wl*pcint22r*pcondcopr(ql,xl,par1));
    e12r=power^2*int0r/denr; 
    e1.1r=ql*power*int1qlr/denr;
    e2.1r=2*ql*power*int21qlr/denr;
    e1.2r=ql*power*int2qlr/denr;
    e2.2r=2*ql*power*int22qlr/denr; 
    v1r=e2.1r-e1.1r^2;
    v2r=e2.2r-e1.2r^2;
  }
  
  e12=power^2*int0/den; 
  e1.1=ql*power*int1ql/den;
  e2.1=2*ql*power*int21ql/den;
  e1.2=ql*power*int2ql/den;
  e2.2=2*ql*power*int22ql/den; 
  v1=e2.1-e1.1^2;
  v2=e2.2-e1.2^2;
  out=(e12-e1.1*e1.2)/sqrt(v1*v2);
  outr=out;
  if(symm==F) outr=(e12r-e1.1r*e1.2r)/sqrt(v1r*v2r);
  out0=c(out,outr);
  out0
}

# Much faster version of twdm() for 2 factor copulas
# running times with nq = 35 and diff. copulas for twdm2:
# BB1+FRK: ~1.51sec.; Gumbel: ~4.92sec.; Frank: ~1.42sec.
# running times with nq = 35 and diff. copulas for twdm1factcop:
# BB1+FRK: ~0.47sec.; Gumbel: ~0.80sec.; Frank: ~0.17sec.; t: ~1.33sec.
# difference between two methods is < 0.001 for all copulas

# pcondcop1 = conditional copula cdf linking V1 and observed variables
# pcondcop2 = conditional copula cdf linking V2 and observed variables
# param1 = dependence parameter of pcondcop1 for two observed variables
# param2 = dependence parameter of pcondcop2 for two observed variables
# power = power to use for tail-weighted dependence measure, good choice is 6
# nq = number of quadrature points for Gauss-Legendre
#   nq = 35 gives very good accuracy for all tail dependent copulas 
#   nq = 25 twice as faster and accuracy is reasonable (+-0.001) 
# Output: vector of length 2 with
#   lower and upper tail-weighted dependence measure values
twdm2factcop=function(pcondcop1,pcondcop2,param1,param2,power,nq)
{ 
  if(is.matrix(param1)) { par1.1=param1[1,]; par1.2=param1[2,]; }
  else { par1.1=param1[1]; par1.2=param1[2]; }         
  if(is.matrix(param2)) { par2.1=param2[1,]; par2.2=param2[2,]; }
  else { par2.1=param2[1]; par2.2=param2[2]; }

  ql=0.5;
  gl=gausslegendre(nq);
  wl=gl$weights;
  wl0=rep(wl,each=nq)*rep(wl,nq);
  xl=gl$nodes;
  wl2=ql*gl$weights;
  xl2=ql*gl$nodes;
  den=  ql^2*pfact2cop(ql,ql,pcondcop1,pcondcop2,param1,param2,35);   
  denr=den; #C(0.5,0.5) = C_R(0.5,0.5)

  pcondcop1r= function(u1,u2,param)
  { return(1-pcondcop1(1-u1,1-u2,param)) }

  pcondcop2r= function(u1,u2,param)
  { return(1-pcondcop2(1-u1,1-u2,param)) }

  symm=F;
  tst=seq(0.01,0.99,0.01);
  dif1= max(abs(pcondcop1(tst,tst,par1.1) - pcondcop1r(tst,tst,par1.1)));
  dif2= max(abs(pcondcop2(tst,tst,par2.1) - pcondcop2r(tst,tst,par2.1)));
  dif=max(dif1,dif2);

  if(dif<1e-5) symm=T; #numerical check for reflection symmetric dependence 
  int0=0; int1ql=0; int21ql=0; int2ql=0; int22ql=0;
  pcint1=0; pcint2=0;
  pcint21=0; pcint22=0;
  int0r=0; int1qlr=0; int21qlr=0; int2qlr=0; int22qlr=0;
  pcint1r=0; pcint2r=0;
  pcint21r=0; pcint22r=0;

  for(iq1 in 1:nq)
  { for(iq2 in 1:nq)
    { pc1= pcondcop1(xl2,xl[iq1],par1.1)
      pc2= pcondcop1(xl2,xl[iq1],par1.2)
      pcint1= c(pcint1, sum(wl2*(1-xl2/ql)^(power-1)*pcondcop2(pc1,xl[iq2],par2.1)));
      pcint2= c(pcint2, sum(wl2*(1-xl2/ql)^(power-1)*pcondcop2(pc2,xl[iq2],par2.2)));
      pcint21= c(pcint21, sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop2(pc1,xl[iq2],par2.1)));
      pcint22= c(pcint22, sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop2(pc2,xl[iq2],par2.2)));

      if(symm==F)
      { pc1r= pcondcop1r(xl2,xl[iq1],par1.1)
        pc2r= pcondcop1r(xl2,xl[iq1],par1.2)
        pcint1r= c(pcint1r, sum(wl2*(1-xl2/ql)^(power-1)*pcondcop2r(pc1r,xl[iq2],par2.1)));
        pcint2r= c(pcint2r, sum(wl2*(1-xl2/ql)^(power-1)*pcondcop2r(pc2r,xl[iq2],par2.2)));
        pcint21r= c(pcint21r, sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop2r(pc1r,xl[iq2],par2.1)));
        pcint22r= c(pcint22r, sum(wl2*(1-xl2/ql)^(2*power-1)*pcondcop2r(pc2r,xl[iq2],par2.2)));
      } 
    }
  }
  pcint1=pcint1[-1]; pcint2=pcint2[-1]; 
  pcint21=pcint21[-1]; pcint22=pcint22[-1];
  pc1ql= rep(pcondcop1(ql,xl,par1.1),each=nq);  
  pc2ql= rep(pcondcop1(ql,xl,par1.2),each=nq);  

  int0 = sum(wl0*pcint1*pcint2);
  int1ql = sum(wl0*pcint1*pcondcop2(pc2ql,rep(xl,nq),par2.2));
  int21ql = sum(wl0*pcint21*pcondcop2(pc2ql,rep(xl,nq),par2.2));
  int2ql = sum(wl0*pcint2*pcondcop2(pc1ql,rep(xl,nq),par2.1));
  int22ql = sum(wl0*pcint22*pcondcop2(pc1ql,rep(xl,nq),par2.1));
 
  if(symm==F)
  { pcint1r=pcint1r[-1]; pcint2r=pcint2r[-1]; 
    pcint21r=pcint21r[-1]; pcint22r=pcint22r[-1];
    pc1qlr= rep(pcondcop1r(ql,xl,par1.1),each=nq);  
    pc2qlr= rep(pcondcop1r(ql,xl,par1.2),each=nq);  
    int0r = sum(wl0*pcint1r*pcint2r);
    int1qlr = sum(wl0*pcint1r*pcondcop2r(pc2qlr,rep(xl,nq),par2.2));
    int21qlr = sum(wl0*pcint21r*pcondcop2r(pc2qlr,rep(xl,nq),par2.2));
    int2qlr = sum(wl0*pcint2r*pcondcop2r(pc1qlr,rep(xl,nq),par2.1));
    int22qlr = sum(wl0*pcint22r*pcondcop2r(pc1qlr,rep(xl,nq),par2.1));
    e12r=power^2*int0r/denr; 
    e1.1r=ql*power*int1qlr/denr;
    e2.1r=2*ql*power*int21qlr/denr;
    e1.2r=ql*power*int2qlr/denr;
    e2.2r=2*ql*power*int22qlr/denr; 
    v1r=e2.1r-e1.1r^2;
    v2r=e2.2r-e1.2r^2;
  }
  
  e12=power^2*int0/den; 
  e1.1=ql*power*int1ql/den;
  e2.1=2*ql*power*int21ql/den;
  e1.2=ql*power*int2ql/den;
  e2.2=2*ql*power*int22ql/den; 
  v1=e2.1-e1.1^2;
  v2=e2.2-e1.2^2;
  out=(e12-e1.1*e1.2)/sqrt(v1*v2);
  outr=out;
  if(symm==F) outr=(e12r-e1.1r*e1.2r)/sqrt(v1r*v2r);
  out0=c(out,outr);
  out0 
}
