# check for problem with Frank copula for cpar>35

library(CopulaModel)

chkfrk=function(cpar)
{ cat("\ncpar=", cpar,"\n")
  u=seq(.1,.6,.1)
  v=seq(.4,.9,.1)
  print(pfrk(u,v,cpar))
  print(pcondfrk(v,u,cpar))
  print(dfrk(u,v,cpar))
  print(dfrk(u,u,cpar))
  invisible(0)
}

chkfrk(10)
chkfrk(20)
chkfrk(30)
chkfrk(35)
chkfrk(37)
chkfrk(38)
chkfrk(40)
chkfrk(50)

cat("\nchecking pcop/pcond/dcop\n")
u=.3
v=seq(.4,.9,.1)
for(be in seq(.1,.9,.2))
{ cpar=frk.b2cpar(be)
  chkcopderiv(u,v,cpar,bcdf=pfrk,pcond=pcondfrk,bpdf=dfrk,str="frk",eps=1.e-5)
}


cat("\nchecking pcond/qcond\n")
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
be=.5
cpar=frk.b2cpar(be)
chkcopcond(u,v,cpar,pcondfrk,qcondfrk,"frk")
be=.9
cpar=frk.b2cpar(be)
chkcopcond(u,v,cpar,pcondfrk,qcondfrk,"frk")
