# First step for creating the function depmeas2cpar
# Dependence tables (beta, tau, rhoS, rhoN, lm) for bivariate copula familes

library(CopulaModel)

bvec=seq(0,.9,.02)  
bvec=c(bvec,.95)
bvec=c(bvec,1)
np=length(bvec)

# output cpar, be, tau, rhoS. rhoN, lm 

cat("\nGalambos\n")
gal.deptab=makedeptable(bvec,bfn=gal.b2cpar,pcop=pgal,
  pcond12=pcondgal,pcond21=pcondgal,LBcpar=0,
  UBcpar=Inf,itaildep=T,lmfn=gal.cpar2lm,zero=1.e-6,zbd=7,iprint=T)
# replace by code with 1-dimensional numerical integrals (instead of 2-dim)
for(i in 2:(np-1))
{ gal.deptab[i,3]=gal.cpar2tau(gal.deptab[i,1])
  gal.deptab[i,4]=gal.cpar2rhoS(gal.deptab[i,1])
}
print(gal.deptab)

cat("\nGumbel\n")
gum.deptab=makedeptable(bvec,bfn=gum.b2cpar,pcop=pgum,
  pcond12=pcondgum,pcond21=pcondgum,LBcpar=1,
  UBcpar=Inf,itaildep=T,lmfn=gum.cpar2lm,zero=1.e-6,zbd=7,iprint=T)
tauinteg=gum.deptab[,3]
# replace by code with exact or 1-dim numerical integrals (instead of 2-dim)
for(i in 2:(np-1))
{ gum.deptab[i,3]=gum.cpar2tau(gum.deptab[i,1])
  gum.deptab[i,4]=gum.cpar2rhoS(gum.deptab[i,1]) 
}
print(gum.deptab)
cat("accuracy of tau from 2-dim integrals\n")
print(abs(tauinteg-gum.deptab[,3]))

cat("\nHuesler-Reiss\n")
hr.deptab=makedeptable(bvec,bfn=hr.b2cpar,pcop=phr,
  pcond12=pcondhr,pcond21=pcondhr,LBcpar=0,
  UBcpar=Inf,itaildep=T,lmfn=hr.cpar2lm,zero=1.e-6,zbd=7,iprint=T)
# replace by code with 1-dimensional numerical integrals (instead of 2-dim)
for(i in 2:(np-1))
{ hr.deptab[i,3]=hr.cpar2tau(hr.deptab[i,1])
  hr.deptab[i,4]=hr.cpar2rhoS(hr.deptab[i,1])
}
print(hr.deptab)

cat("\nPlackett\n")
pla.deptab=makedeptable(bvec,bfn=pla.b2cpar,pcop=ppla,
  pcond12=pcondpla,pcond21=pcondpla,LBcpar=0,
  UBcpar=Inf,itaildep=F,zbd=7,iprint=T)
rhointeg=pla.deptab[,4]
# replace by code with exact (instead of 2-dim integrals)
pla.deptab[,4]=pla.cpar2rhoS(pla.deptab[,1])
pla.deptab[1,4]=0
pla.deptab[np,4]=1
print(pla.deptab)
cat("accuracy of rhoS from 2-dim integrals\n")
print(abs(rhointeg-pla.deptab[,4]))

cat("\nJoeB5\n")
joe.deptab=makedeptable(bvec,bfn=joe.b2cpar,pcop=pjoe,
  pcond12=pcondjoe,pcond21=pcondjoe,LBcpar=1,
  UBcpar=Inf,itaildep=T,lmfn=joe.cpar2lm,zero=0,zbd=7,iprint=T)
tauinteg=joe.deptab[,3]
# replace by code with exact (instead of 2-dim integrals)
for(i in 2:(np-1))
{ joe.deptab[i,3]=joe.cpar2tau(joe.deptab[i,1]) }
print(joe.deptab)
cat("accuracy of tau from 2-dim integrals\n")
print(abs(tauinteg-joe.deptab[,3]))

cat("\nMTCJ\n") # fails for beta=.95
bvec=bvec[-(np-1)]
np=length(bvec)
mtcj.deptab=makedeptable(bvec,bfn=mtcj.b2cpar,pcop=pmtcj,
  pcond12=pcondmtcj,pcond21=pcondmtcj,LBcpar=0,
  UBcpar=Inf,itaildep=T,lmfn=mtcj.cpar2lm,zero=0,zbd=7,iprint=T)
tauinteg=mtcj.deptab[,3]
# replace by code with exact (instead of 2-dim integrals)
for(i in 2:(np-1))
{ mtcj.deptab[i,3]=mtcj.cpar2tau(mtcj.deptab[i,1]) }
print(mtcj.deptab)
cat("accuracy of tau from 2-dim integrals\n")
print(abs(tauinteg-mtcj.deptab[,3]))

cat("\nFrank\n") # fails for beta=.95
frk.deptab=makedeptable(bvec,bfn=frk.b2cpar,pcop=pfrk,
  pcond12=pcondfrk,pcond21=pcondfrk,LBcpar=0,
  UBcpar=Inf,itaildep=F,zero=0,zbd=7,iprint=T)
tauinteg=frk.deptab[,3]
rhointeg=frk.deptab[,4]
# replace by code with tau and rhoS as 1-dimensional integrals
for(i in 2:(np-1))
{ frk.deptab[i,3]=frk.cpar2tau(frk.deptab[i,1]) 
  frk.deptab[i,4]=frk.cpar2rhoS(frk.deptab[i,1])
}
print(frk.deptab)
cat("accuracy of tau and rhoS from 2-dim integrals\n")
print(abs(tauinteg-frk.deptab[,3]))
print(abs(rhointeg-frk.deptab[,4]))

#save(file="deptab.RData", pla.deptab,frk.deptab,mtcj.deptab,joe.deptab,gum.deptab,gal.deptab,hr.deptab)

