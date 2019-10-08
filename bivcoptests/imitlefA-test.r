# Archimedean copula with integrated Mittag-Leffler LT:
# LT involves incomplete beta
library(CopulaModel)

# param=(vth,de), vth>0, de>1
param=c(0.5,1.3)
param=c(0.5,1.4)
#param=c(1.6,1.4)
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
print(pimitlefA(u,v,param))
print(pcondimitlefA(v,u,param))
print(dimitlefA(u,v,param))

cat("\ncheck pcop, pcond, dcop\n")
param=c(2.2,1.4)
param=c(1/2.2,1.4)
u=.3
#u=.8
v=seq(.4,.9,.1)
chkcopderiv(u,v,param,bcdf=pimitlefA,pcond=pcondimitlefA,bpdf=dimitlefA,str="imitlefA",eps=1.e-5)

cat("\ncheck pcond, qcond\n")
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
chkcopcond(u,v,param,pcondimitlefA,qcondimitlefA,"imitlefA")

# This function in imitlefA.R but not exported
#   set con=beta(de,ze-de) in calling routine
# s = positive value
# ga = 1/vth = ze-de >0, 
# de = second parameter >1
# Output: first derivative of LT for the imitlefA copula
ibetalt1=function(s,ga,de,con)
{ s1=s^(1/de)
  ze=de+ga
  -(1+s1)^(-ze) /(de*con)
}

# integrand for Kendall tau formula
#   set con=beta(de,ze-de) in calling routine
# s = positive value
# vth = first parameter = 1/ga = 1/(ze-de) >0, 
# de = second parameter >1
psiderfn=function(s,vth,de)
{ ga=1/vth
  ze=ga+de
  con=beta(de,ga)
  der=ibetalt1(s,ga,de,con)
  s*der^2
}

# vth = first parameter >0 
# de = second parameter >1
# Output: compare the integral version of tau with the closed form formula
chktau=function(vth,de)
{ out=integrate(psiderfn,0,Inf,vth=vth,de=de)
  tau=imitlefA.cpar2tau(c(vth,de))
  print(c(vth,de,1-4*out$value,tau))  
  tau
}

cat("\ncheck on Kendall's tau\n")
tem=chktau(2,1.5)
tem=chktau(2,4.4)
tem=chktau(1/1.5,4.4)
tem=chktau(1/6.5,4.4)
tem=chktau(10,4.4)
# OK matches


