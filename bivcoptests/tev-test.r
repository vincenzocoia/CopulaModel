# checks for tEV cdf, pcond, pdf

library(CopulaModel)
#source("../R/tev.R")

# first parameter rho, second nu
#cpar=c(.5,5)
#cpar=c(.5,15)
#cpar=c(-.5,5)
#cpar=c(-.9,3)
#cpar=c(.8,2)

u=.3
#u=.8
v=seq(.4,.9,.1)
cat("check of ptev, pcondtev, dtev\n")
chkcopderiv(u,v,cpar=c(.5,5),bcdf=ptev,pcond=pcondtev,bpdf=dtev,str="tev",eps=1.e-4)
chkcopderiv(u,v,cpar=c(.5,15),bcdf=ptev,pcond=pcondtev,bpdf=dtev,str="tev",eps=1.e-4)
chkcopderiv(u,v,cpar=c(-.5,5),bcdf=ptev,pcond=pcondtev,bpdf=dtev,str="tev",eps=1.e-4)
chkcopderiv(u,v,cpar=c(-.9,3),bcdf=ptev,pcond=pcondtev,bpdf=dtev,str="tev",eps=1.e-4)
chkcopderiv(u,v,cpar=c(.8,2),bcdf=ptev,pcond=pcondtev,bpdf=dtev,str="tev",eps=1.e-4)
# looks OK

# nu = degree of freedom parameter >0
# Output within function: tau and rhoS given rho parameter in (-1,1)
wraptaurho=function(nu)
{ cat("\nnu=",nu,"\n")
  for(rho in seq(-.9,.9,.2))
  { ktau=tev.cpar2tau(c(rho,nu))
    srho=tev.cpar2rhoS(c(rho,nu))
    cat(rho,ktau,srho,"\n")
  }
  invisible(0)
}

cat("\nKendall's tau and Spearman's rho given nu,rho\n")
wraptaurho(.5)
wraptaurho(2)
wraptaurho(5)
wraptaurho(10)
wraptaurho(20)
wraptaurho(30)

cat("\nCompare Huesler-Reiss with de=1/sqrt(logn*(1-rho)) for large n\n")
n=1e8; logn=log(n)
for(rho in seq(-.9,.9,.2))
{ de=1/sqrt(logn*(1-rho))
  ktau=hr.cpar2tau(de)
  srho=hr.cpar2rhoS(de)
  cat(rho,de,ktau,srho,"\n")
}

cat("\ncompare lambda and beta\n")
for(rho in c(0,.3,.5,.7))
{ for(nu in c(1:5,10,30))
  { lm=bvt.cpar2lm(c(rho,nu))
    be=tev.cpar2b(c(rho,nu))
    print(c(rho,nu,lm,be))
  }
}

