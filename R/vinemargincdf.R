
# Functions for bivariate marginal cdfs, pdfs, conditionals of vine copulas
# Can be used for Spearman's rho and Kendall's tau 
#  (it is better to compute gl function calls to cdf etc)
# Possible later change so that last argument is gl and not nq

# Inputs in form of C-vine with nodes a,b,c:
# param = vector or 3-row matrix
#  par,ab,par.ac,par.bc : parameters of pcondba, pcondca, pcopbc
#  ub, uc = vectors of the same length
#    or one of these variables can be a constant 
#  pcondba = function for conditional cdf of copula for first tree1 edge
#  pcondca = function for conditional cdf of copula for second tree1 edge
#  pcopbc = copula cdf for tree2 edge joining the above two
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate cdf that is bivariate margin of edge of tree 2
ptree2cop=function(ub,uc,param,pcondba,pcondca,pcopbc,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.ac=param[2,]; par.bc=param[3,] }
  else
  { par.ab=param[1]; par.ac=param[2]; par.bc=param[3] }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  m1=length(ub); m2=length(uc);
  m=max(m1,m2);  cdf=rep(0,m); 
  if(length(ub)==1) ub=rep(ub, m);
  if(length(uc)==1) uc=rep(uc, m);
  for(i in 1:m) 
  { tem=pcopbc(pcondba(ub[i],xl,par.ab),pcondca(uc[i],xl,par.ac),par.bc)
    cdf[i]=sum(wl*tem)
  }
  cdf
}

# Inputs in form of C-vine with nodes a,b,c:
# param = vector or 3-row matrix
#  par,ab,par.ac,par.bc : parameters of pcondba, pcondca, pcopbc
#  ub, uc = vectors of the same length
#    or one of these variables can be a constant 
#  pcondba = function for conditional cdf of copula for first tree1 edge
#  pcondca = function for conditional cdf of copula for second tree1 edge
#  pcondcb= function for conditional cdf of copula  for tree2 edge 
#  dcopab = function for copula density for first tree1 edge
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: conditional cdf of bivariate margin of edge of tree 2
pcondtree2=function(uc,ub,param,pcondba,pcondca,pcondcb,dcopab,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.ac=param[2,]; par.bc=param[3,] }
  else
  { par.ab=param[1]; par.ac=param[2]; par.bc=param[3] }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  m1=length(ub); m2=length(uc);
  m=max(m1,m2);  ccdf=rep(0,m); 
  if(length(ub)==1) ub=rep(ub, m);
  if(length(uc)==1) uc=rep(uc, m);
  for(i in 1:m) 
  { tem=pcondcb(pcondca(uc[i],xl,par.ac),pcondba(ub[i],xl,par.ab),par.bc)
    tem=tem*dcopab(xl,ub[i],par.ab)
    ccdf[i]=sum(wl*tem)
  }
  ccdf
}

# Inputs in form of C-vine with nodes a,b,c:
# param = vector or 3-row matrix
#  par,ab,par.ac,par.bc : parameters of pcondba, pcondca, pcopbc
#  ub, uc = vectors of the same length
#    or one of these variables can be a constant 
#  pcondba = function for conditional cdf of copula for first tree1 edge
#  pcondca = function for conditional cdf of copula for second tree1 edge
#  dcopbc = function for copula density for tree2 edge 
#  dcopab = function for copula density for first tree1 edge
#  dcopac = function for copula density for second tree1 edge
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate pdf that is bivariate margin of edge of tree 2
dtree2cop=function(ub,uc,param,pcondba,pcondca,dcopbc,dcopab,dcopac,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.ac=param[2,]; par.bc=param[3,] }
  else
  { par.ab=param[1]; par.ac=param[2]; par.bc=param[3] }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  m1=length(ub); m2=length(uc);
  m=max(m1, m2); pdf=rep(0,m); 
  if(length(ub)==1) ub=rep(ub, m);
  if(length(uc)==1) uc=rep(uc, m);
  for(i in 1:m) 
  { tem=dcopbc(pcondba(ub[i],xl,par.ab),pcondca(uc[i],xl,par.ac),par.bc)
    tem=tem*dcopab(xl,ub[i],par.ab)*dcopac(xl,uc[i],par.ac)
    pdf[i]=sum(wl*tem)
  }
  pdf
}

#============================================================

# Inputs in form of C-vine with nodes a,b,c,d
# param = vector or 6-row matrix
#  uc and ud are vectors of the same length
#    or one of these variables can be a constant 
#  pcondba = function for conditional cdf of copula for first tree1 edge
#  pcondca = function for conditional cdf of copula for second tree1 edge
#  pcondda = function for conditional cdf of copula for third tree1 edge
#  pcondcb = function for conditional cdf of copula for first tree2 edge
#  pconddb = function for conditional cdf of copula for second tree2 edge
#  pcopcd = copula cdf for tree3 edge joining the above two
#  dcopab = copula pdf for tree1 edge 
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate cdf that is bivariate margin of edge of tree 3
ptree3cop.cvine=function(uc,ud,param,pcondba,pcondca,pcondda,pcondcb,pconddb,
  pcopcd,dcopab,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.ac=param[2,]; par.ad=param[3,];
    par.bc=param[4,]; par.bd=param[5,]; par.cd=param[6,];
  }
  else
  { par.ab=param[1]; par.ac=param[2]; par.ad=param[3] 
    par.bc=param[4]; par.bd=param[5]; par.cd=param[6]
  }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  xl2=rep(xl,each=nq); 
  xl1=rep(xl,nq)
  wl0= rep(wl,each=nq)*rep(wl,nq);
  m1=length(uc); m2=length(ud);
  m=max(m1, m2); cdf=rep(0,m); 
  if(length(uc)==1) uc=rep(uc, m);
  if(length(ud)==1) ud=rep(ud, m);
  for(i in 1:m) 
  { ccdfca= pcondca(uc[i],xl,par.ac); ccdfca= rep(ccdfca, nq);
    ccdfda= pcondda(ud[i],xl,par.ad); ccdfda= rep(ccdfda, nq);
    ccdfba= pcondba(xl2,xl1,par.ab); 
    pdfab= dcopab(xl1,xl2,par.ab); 
    ccdfcb= pcondcb(ccdfca,ccdfba,par.bc); 
    ccdfdb= pconddb(ccdfda,ccdfba,par.bd);
    cdfcd= pcopcd(ccdfcb,ccdfdb,par.cd)
    tem= cdfcd* pdfab
    cdf[i]=sum(wl0*tem)
  }
  cdf
}

# Inputs in form of C-vine with nodes a,b,c,d
# param = vector or 6-row matrix
#  uc and ud are vectors of the same length
#    or one of these variables can be a constant 
#  pcondba = function for conditional cdf of copula for first tree1 edge
#  pcondca = function for conditional cdf of copula for second tree1 edge
#  pcondda = function for conditional cdf of copula for third tree1 edge
#  pcondcb = function for conditional cdf of copula for first tree2 edge
#  pconddb = function for conditional cdf of copula for second tree2 edge
#  pconddc = function for conditional cdf for tree3 edge joining the above two
#  dcopab = function for copula density for first tree1 edge
#  dcopac = function for copula density for second tree1 edge
#  dcopbc = function for copula density for tree2 edge
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: conditional cdf of bivariate margin of edge of tree 3
pcondtree3.cvine=function(ud,uc,param,pcondba,pcondca,pcondda,pcondcb,pconddb,
  pconddc,dcopab,dcopac,dcopbc,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.ac=param[2,]; par.ad=param[3,];
    par.bc=param[4,]; par.bd=param[5,]; par.cd=param[6,];
  }
  else
  { par.ab=param[1]; par.ac=param[2]; par.ad=param[3] 
    par.bc=param[4]; par.bd=param[5]; par.cd=param[6]
  }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  xl2=rep(xl,each=nq); 
  xl1=rep(xl,nq)
  wl0=rep(wl,each=nq)*rep(wl,nq);
  m1=length(uc);  m2=length(ud);
  m=max(m1,m2);   ccdf=rep(0,m); 
  if(length(uc)==1) uc=rep(uc, m);
  if(length(ud)==1) ud=rep(ud, m);
  for(i in 1:m) 
  { ccdfca= pcondca(uc[i],xl,par.ac); ccdfca= rep(ccdfca, nq);
    pdfac= dcopac(xl,uc[i],par.ac); pdfac= rep(pdfac, nq);
    ccdfda= pcondda(ud[i],xl,par.ad); ccdfda= rep(ccdfda, nq);
    ccdfba= pcondba(xl2,xl1,par.ab); 
    pdfab= dcopab(xl1,xl2,par.ab); 
    pdfbc= dcopbc(ccdfba,ccdfca,par.bc)
    ccdfcb= pcondcb(ccdfca,ccdfba,par.bc); 
    ccdfdb= pconddb(ccdfda,ccdfba,par.bd);
    ccdfdc= pconddc(ccdfdb,ccdfcb,par.cd)
    tem= ccdfdc* pdfab* pdfac* pdfbc
    ccdf[i]=sum(wl0*tem)
  }
  ccdf
}


# Inputs in form of C-vine with nodes a,b,c,d
# param = vector or 6-row matrix
#  uc, ud = vectors of the same length
#    or one of these variables can be a constant 
#  pcondba = function for conditional cdf of copula for first tree1 edge
#  pcondca = function for conditional cdf of copula for second tree1 edge
#  pcondda = function for conditional cdf of copula for third tree1 edge
#  pcondcb = function for conditional cdf of copula for first tree2 edge
#  pconddb = function for conditional cdf of copula for second tree2 edge
#  dcopcd = function for conditional cdf for tree3 edge joining the above two
#  dcopab = function for copula density for first tree1 edge
#  dcopac = function for copula density for second tree1 edge
#  dcopad = function for copula density for third tree1 edge
#  dcopbc = function for copula density for tree2 edge
#  dcopbd = function for copula density for tree2 edge
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate pdf that is bivariate margin of edge of tree 3
dtree3cop.cvine=function(uc,ud,param,pcondba,pcondca,pcondda,pcondcb,pconddb,
  dcopcd,dcopab,dcopac,dcopad,dcopbc,dcopbd,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.ac=param[2,]; par.ad=param[3,];
    par.bc=param[4,]; par.bd=param[5,]; par.cd=param[6,];
  }
  else
  { par.ab=param[1]; par.ac=param[2]; par.ad=param[3] 
    par.bc=param[4]; par.bd=param[5]; par.cd=param[6]
  }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  xl2=rep(xl,each=nq); 
  xl1=rep(xl,nq)
  wl0=rep(wl,each=nq)*rep(wl,nq);
  m1=length(uc); m2=length(ud);
  m=max(m1,m2);  pdf=rep(0,m); 
  if(length(uc)==1) uc=rep(uc, m);
  if(length(ud)==1) ud=rep(ud, m);
  for(i in 1:m) 
  { ccdfca= pcondca(uc[i],xl,par.ac); ccdfca= rep(ccdfca, nq);
    pdfac= dcopac(xl,uc[i],par.ac); pdfac= rep(pdfac, nq);
    pdfad= dcopad(xl,ud[i],par.ad); pdfad= rep(pdfad, nq);
    ccdfda= pcondda(ud[i],xl,par.ad); ccdfda= rep(ccdfda, nq);
    ccdfba= pcondba(xl2,xl1,par.ab); 
    pdfab= dcopab(xl1,xl2,par.ab); 
    pdfbc= dcopbc(ccdfba,ccdfca,par.bc)
    ccdfcb= pcondcb(ccdfca,ccdfba,par.bc); 
    pdfbd= dcopbd(ccdfba,ccdfda,par.bd)
    ccdfdb= pconddb(ccdfda,ccdfba,par.bd);
    pdfcd= dcopcd(ccdfcb,ccdfdb,par.cd)
    tem= pdfcd* pdfab* pdfac* pdfbc* pdfad* pdfbd
    pdf[i]=sum(wl0*tem)
  }
  pdf
}

#============================================================

# Inputs in form of D-vine with nodes a,b,c,d
# param = vector or 6-row matrix
#  ua, ud = vectors of the same length
#    or one of these variables can be a constant 
#  pcondab = function for conditional cdf of copula for first tree1 edge
#  pcondbc = function for conditional cdf of copula for second tree1 edge
#  pcondcb = function for conditional cdf of copula for second tree1 edge
#  pconddc = function for conditional cdf of copula for third tree1 edge
#  pcondac = function for conditional cdf of copula for first tree2 edge
#  pconddb = function for conditional cdf of copula for second tree2 edge
#  pcopad = copula cdf for tree3 edge joining the above two
#  dcopbc = copula pdf for tree1 edge 
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate cdf that is bivariate margin of edge of tree 3
ptree3cop.dvine=function(ua,ud,param,pcondab,pcondbc,pcondcb,pconddc,
  pcondac,pconddb,pcopad,dcopbc,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.bc=param[2,]; par.cd=param[3,];
    par.ac=param[4,]; par.bd=param[5,]; par.ad=param[6,];
  }
  else
  { par.ab=param[1]; par.bc=param[2]; par.cd=param[3] 
    par.ac=param[4]; par.bd=param[5]; par.ad=param[6]
  }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  xl2=rep(xl,each=nq); 
  xl1=rep(xl,nq)
  wl0=rep(wl,each=nq)*rep(wl,nq);
  m1=length(ua); m2=length(ud);
  m=max(m1,m2);  cdf=rep(0,m); 
  if(length(ua)==1) ua=rep(ua,m);
  if(length(ud)==1) ud=rep(ud,m);
  for(i in 1:m) 
  { ccdfab= pcondab(ua[i],xl,par.ab); ccdfab= rep(ccdfab, nq);
    ccdfdc= pconddc(ud[i],xl,par.cd); ccdfdc= rep(ccdfdc, each=nq);
    ccdfcb= pcondcb(xl2,xl1,par.bc);  
    ccdfbc= pcondbc(xl1,xl2,par.bc); 
    pdfbc= dcopbc(xl1,xl2,par.bc); 
    ccdfac= pcondac(ccdfab,ccdfcb,par.ac); 
    ccdfdb= pconddb(ccdfdc,ccdfbc,par.bd);
    cdfad= pcopad(ccdfac,ccdfdb,par.ad)
    tem= cdfad* pdfbc
    cdf[i]=sum(wl0*tem)
  }
  cdf
}


# Inputs in form of D-vine with nodes a,b,c,d
# param = vector or 6-row matrix
#  ua, ud = vectors of the same length
#    or one of these variables can be a constant 
#  pcondab = function for conditional cdf of copula for first tree1 edge
#  pcondbc = function for conditional cdf of copula for second tree1 edge
#  pcondcb = function for conditional cdf of copula for second tree1 edge
#  pconddc = function for conditional cdf of copula for third tree1 edge
#  pcondac = function for conditional cdf of copula for first tree2 edge
#  pconddb = function for conditional cdf of copula for second tree2 edge
#  pcondda = function for conditional cdf for tree3 edge joining the above two
#  dcopbc = copula pdf for tree1 edge 
#  dcopab = copula pdf for tree1 edge 
#  dcopac = copula pdf for tree2 edge 
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: conditional cdf of bivariate margin of edge of tree 3
pcondtree3.dvine=function(ud,ua,param,pcondab,pcondbc,pcondcb,pconddc,
  pcondac,pconddb,pcondda,dcopab,dcopbc,dcopac,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.bc=param[2,]; par.cd=param[3,];
    par.ac=param[4,]; par.bd=param[5,]; par.ad=param[6,];
  }
  else
  { par.ab=param[1]; par.bc=param[2]; par.cd=param[3] 
    par.ac=param[4]; par.bd=param[5]; par.ad=param[6]
  }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  xl2=rep(xl,each=nq); 
  xl1=rep(xl,nq)
  wl0=rep(wl,each=nq)*rep(wl,nq);
  m1=length(ua);  m2=length(ud);
  m=max(m1,m2);   ccdf=rep(0,m); 
  if(length(ua)==1) ua=rep(ua, m);
  if(length(ud)==1) ud=rep(ud, m);
  for(i in 1:m) 
  { ccdfab= pcondab(ua[i],xl,par.ab); ccdfab= rep(ccdfab, nq);
    pdfab= dcopab(ua[i],xl,par.ab); pdfab= rep(pdfab, nq);
    ccdfdc= pconddc(ud[i],xl,par.cd); ccdfdc= rep(ccdfdc, each=nq);
    ccdfcb= pcondcb(xl2,xl1,par.bc);  
    ccdfbc= pcondbc(xl1,xl2,par.bc); 
    pdfbc= dcopbc(xl1,xl2,par.bc); 
    pdfac= dcopac(ccdfab,ccdfcb,par.ac)
    ccdfac= pcondac(ccdfab,ccdfcb,par.ac); 
    ccdfdb= pconddb(ccdfdc,ccdfbc,par.bd);
    ccdfda= pcondda(ccdfdb,ccdfac,par.ad)
    tem= ccdfda* pdfbc* pdfab* pdfac
    ccdf[i]=sum(wl0*tem)
  }
  ccdf
}


# Inputs in form of D-vine with nodes a,b,c,d
# param = vector or 6-row matrix
#  ua, ud = vectors of the same length
#    or one of these variables can be a constant 
#  pcondab = function for conditional cdf of copula for first tree1 edge
#  pcondbc = function for conditional cdf of copula for second tree1 edge
#  pcondcb = function for conditional cdf of copula for second tree1 edge
#  pconddc = function for conditional cdf of copula for third tree1 edge
#  pcondac = function for conditional cdf of copula for first tree2 edge
#  pconddb = function for conditional cdf of copula for second tree2 edge
#  dcopad = copula pdf for tree3 edge joining the above two
#  dcopbc = copula pdf for tree1 edge 
#  dcopab = copula pdf for tree1 edge 
#  dcopcd = copula pdf for tree1 edge 
#  dcopac = copula pdf for tree2 edge 
#  dcopbd = copula pdf for tree2 edge 
#  nq = number of quadrature points for Gauss-Legendre quadrature
# Output: bivariate pdf that is bivariate margin of edge of tree 3
dtree3cop.dvine=function(ua,ud,param,pcondab,pcondbc,pcondcb,pconddc,
  pcondac,pconddb,dcopad,dcopab,dcopbc,dcopcd,dcopac,dcopbd,nq) 
{ if(is.matrix(param))
  { par.ab=param[1,]; par.bc=param[2,]; par.cd=param[3,];
    par.ac=param[4,]; par.bd=param[5,]; par.ad=param[6,];
  }
  else
  { par.ab=param[1]; par.bc=param[2]; par.cd=param[3] 
    par.ac=param[4]; par.bd=param[5]; par.ad=param[6]
  }
  gl=gausslegendre(nq);
  wl=gl$weights; xl=gl$nodes; 
  xl2=rep(xl,each=nq); 
  xl1=rep(xl,nq)
  wl0=rep(wl,each = nq)*rep(wl,nq);
  m1=length(ua);  m2=length(ud);
  m=max(m1,m2);   pdf=rep(0,m); 
  if(length(ua)==1) ua=rep(ua, m);
  if(length(ud)==1) ud=rep(ud, m);
  for(i in 1:m) 
  { ccdfab= pcondab(ua[i],xl,par.ab); ccdfab= rep(ccdfab, nq);
    pdfab= dcopab(ua[i],xl,par.ab); pdfab= rep(pdfab, nq);
    ccdfdc= pconddc(ud[i],xl,par.cd); ccdfdc= rep(ccdfdc,each=nq);
    pdfcd= dcopcd(xl,ud[i],par.cd); pdfcd= rep(pdfcd, each=nq)
    ccdfcb= pcondcb(xl2,xl1,par.bc);  
    ccdfbc= pcondbc(xl1,xl2,par.bc); 
    pdfbc= dcopbc(xl1,xl2,par.bc); 
    pdfac= dcopac(ccdfab,ccdfcb,par.ac)
    ccdfac= pcondac(ccdfab,ccdfcb,par.ac); 
    pdfbd= dcopbd(ccdfbc,ccdfdc,par.bd);
    ccdfdb= pconddb(ccdfdc,ccdfbc,par.bd);
    pdfad= dcopad(ccdfac,ccdfdb,par.ad)
    tem= pdfad* pdfbc* pdfab* pdfac * pdfcd* pdfbd
    pdf[i]=sum(wl0*tem)
  }
  pdf
}

