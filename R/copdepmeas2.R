# More functions for transforms of copula parameter to tau and rhoS
# BB2-BB9 copulas are included here if appropriate

# cpar is usually a 2-vector in these functions
#   some of these functions can be converted to take a 2-column matrix 

# cpar = copula parameter with th>1, de>0
bb5.cpar2tau=function(cpar)
{ # B(w) = [ w^th+w1^th - (w^(-th*de) + w1^(-th*de))^(-1/de) ]^(1/th)
  # B'(w) = [ ]^(1/th-1) { w^(th-1) - w1^(th-1)
  #  +  (w^(-th*de) + w1^(-th*de))^(-1/de-1) * (-w^(-th*de-1) + w1^(-th*de-1))
  th=cpar[1]; de=cpar[2]; td=th*de
  taubb5= function(w)
  { w1=1-w; 
    wt=w^th; w1t=w1^th
    wtd=w^(-td); w1td=w1^(-td); 
    tem=wtd+w1td; temd=(wtd+w1td)^(-1/de)
    sm=wt+w1t - temd
    B=sm^(1/th)
    Bp= (B/sm)* ( wt/w - w1t/w1 + (temd/tem) * (-wtd/w + w1td/w1) )   
    ((2*w-1)*Bp*B + w*w1*Bp*Bp )/(B*B)
  }
  tem=integrate(taubb5,0,1,rel.tol=1.e-6)
  tau=tem$value
  tau
}

# cpar = copula parameter with th>1, de>0
bb5.cpar2rhoS=function(cpar)
{ th=cpar[1]; de=cpar[2]; td=th*de
  spbb5= function(w,cpar)
  { w1=1-w; 
    wt=w^th; w1t=w1^th
    wtd=w^(-td); w1td=w1^(-td); 
    tem=wtd+w1td; temd=(wtd+w1td)^(-1/de)
    sm=wt+w1t - temd
    B=sm^(1/th)
    tem=1/(B+1)
    tem^2
  }
  tem=integrate(spbb5,0,1,rel.tol=1.e-6)
  rho=12*tem$value-3
  rho
}

# cpar = copula parameter with th>0, de>0
bb2.cpar2tau=function(cpar)
{ th=cpar[1]; de=cpar[2]; 
  psider= function(s)
  { s1=1+s; ls1=log(s1)
    psid=((1+ls1/de)^(-1/th-1))/s1/th/de  # negative deriv
    s*psid^2
  }
  tem=integrate(psider,0,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau

}

# cpar = copula parameter with th>1, de>0
bb3.cpar2tau=function(cpar)
{ th=cpar[1]; de=cpar[2]; deth=de^(-1/th) 
  psider= function(s)
  { s1=1+s; ls1=log(s1); tem=ls1^(1/th)
    psid= exp(-deth*tem) * deth* tem/ls1/s1/th # negative deriv
    s*psid^2  # infinite at 0 
  }
  tem=integrate(psider,0.,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau
}

# cpar = copula parameter with th>1, de>1
bb6.cpar2tau=function(cpar)
{ th=cpar[1]; de=cpar[2]; 
  psider= function(s)
  { spd=s^(1/de);  espd=exp(-spd)
    #psid=((1-espd)^(1/th-1)) * espd*spd  /s/th/de  # negative deriv
    #s*psid^2  # infinite at 0 
    psid=((1-espd)^(1/th-1)) *espd*spd /th/de  # alternative
    psid^2/s  # infinite at 0 
    
  }
  tem=integrate(psider,0.,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau
}

# cpar = copula parameter with th>1, de>0
# split into 2 cases 1<th<2 and th>=2
# cpar is a 2-vector not a matrix
bb7.cpar2tau=function(cpar)
{ th=cpar[1]; de=cpar[2]; 
  if(th<2) 
  { tau=1-2/de/(2-th)+4*beta(de+2,2/th-1)/de/th^2
    return(tau)
  }
  # else th>=2
  psider= function(s)
  { s1=1+s; tem=s1^(-1/de)
    psid=((1-tem)^(1/th-1))*tem/s1/th/de   # negative deriv
    s*psid^2
  }
  tem=integrate(psider,0,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau
}

# cpar = copula parameter with th>0, 0<p<1
bb8.cpar2tau=function(cpar)
{ vth=cpar[1]; de=cpar[2]; eta=1-(1-de)^vth
  psider= function(s)
  { s1=exp(-s); 
    psid= (eta/vth/de) * ((1-eta*s1)^(1/vth-1)) * s1 # negative deriv
    s*psid^2  
  }
  tem=integrate(psider,0.,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau
}

# cpar = copula parameter with th>1, ga>0
bb9.cpar2tau=function(cpar)
{ th1=1/cpar[1]; ga=cpar[2]; gath=ga^(-cpar[1]) 
  psider= function(s)
  { s1=gath+s; tem=s1^th1
    psid=exp(-tem+1/ga) * th1* tem/s1  # negative deriv
    s*psid^2
  }
  tem=integrate(psider,0,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau

}

# cpar = copula parameter with th>0, 0<p<1
bb10.cpar2tau=function(cpar)
{ th=cpar[1]; p=cpar[2]; qq=(1-p)^(1/th)
  psider= function(s)
  { s1=exp(-s); sth=exp(-s/th)
    psid= (qq/th) * ((1-p*s1)^(-1/th-1)) * sth # negative deriv
    s*psid^2  
  }
  tem=integrate(psider,0.,Inf,rel.tol=1.e-6)
  tau=1-4*tem$value
  tau
}


