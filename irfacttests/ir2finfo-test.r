# check of Fisher info for item response factor models
# compare standalone R versus f90 code

library(CopulaModel)
library(abind)
source("fisherinfo-IR.r") # for non-f90 approach
source("cop-deriv-IR.r")
################################################################################
K=4
d=5
ucut=seq(1:(K-1))/K
ucuts=matrix(rep(ucut,d),ncol=d)
nn=1000
nq=15
################################################################################
# conditional independence model with t
#pcondcop1=pcondbvtcop
#pcondcop2=pcondbvtcop
#pconddotcop1=pconddott
#pconddotcop2=pconddott
#dcop1=dtcop
#dcop2=dtcop
################################################################################

# gumbel
cat("Gumbel\n")
theta=c(1.5,1.6,1.5,1.6,2.4)
delta=c(1.3,1.4,1.3,1.4,1.2)
cat("\ntheta\n"); print(theta); cat("delta\n"); print(delta)
# SEs with 2d copulas
Finfgum=irfisherinfo2(nq,ucuts,theta,delta,pcondcop1=pcondgum,pcondcop2=pcondgum,
  pconddotcop1=pconddotgum,pconddotcop2=pconddotgum,dcop1=dgum,dcop2=dgum,
  iindep=F)
cat("SEs with 2d linking copulas\n")
se1=sqrt(diag(solve(Finfgum))/nn)
print(se1)
ifixed=rep(F,2*d)
finfgum=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="gumbel",nq=nq,
  ifixed,nn=nn)
print(finfgum$SE)
print(max(abs(finfgum$finfo - Finfgum)))

cat("\nGumbel: fix one parameter for factor 2 at cond independence\n")
delta[d]=1.001
Finfgumb=irfisherinfo2(nq,ucuts,theta,delta,pcondcop1=pcondgum,pcondcop2=pcondgum,
  pconddotcop1=pconddotgum,pconddotcop2=pconddotgum,dcop1=dgum,dcop2=dgum,
  iindep=T)
cat("SEs with 2d-1 linking copulas\n")
se1=sqrt(diag(solve(Finfgumb))/nn)
print(se1)
ifixed=rep(F,2*d)
ifixed[2*d]=TRUE 
finfgumb=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="gumbel",nq=15,
  ifixed,nn=nn)
print(finfgumb$SE)
print(max(abs(finfgumb$finfo - Finfgumb)))

#============================================================

dfdefault=5

cat("\nt5/t5\n")
theta=c(0.5,0.6,0.5,0.6,0.4)
delta=c(0.3,0.4,0.3,0.4,0.2)
cat("\ntheta\n"); print(theta); cat("delta\n"); print(delta)
ifixed=rep(F,2*d)
Finft=irfisherinfo2(nq,ucuts,theta,delta,pcondcop1=pcondt,pcondcop2=pcondt,
  pconddotcop1=pconddott,pconddotcop2=pconddott,dcop1=dtcop,dcop2=dtcop,
  iindep=F)
cat("SEs with 2d linking copulas\n")
se1=sqrt(diag(solve(Finft))/nn)
print(se1)
finft=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="t",nq=nq,ifixed,nu1=5,nu2=5,nn=nn)
print(finft$SE)
print(max(abs(finft$finfo - Finft)))
# Ok matches

#============================================================

cat("\nt5/Gumbel\n")
delta=c(1.3,1.4,1.3,1.4,1.2)
theta=c(0.5,0.6,0.5,0.6,0.4)
cat("\ntheta\n"); print(theta); cat("delta\n"); print(delta)
ifixed=rep(F,2*d)
Finftgum=irfisherinfo2(nq,ucuts,theta,delta,pcondcop1=pcondt,pcondcop2=pcondgum,
  pconddotcop1=pconddott,pconddotcop2=pconddotgum,dcop1=dtcop,dcop2=dgum,
  iindep=F)
cat("SEs with 2d linking copulas\n")
se1=sqrt(diag(solve(Finftgum))/nn)
print(se1)
finftgum=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="tgumbel",nq=nq,
  ifixed,nu1=5,nn=nn)
print(finftgum$SE)
print(max(abs(finftgum$finfo - Finftgum)))
# Ok matches

cat("\nGumbel/t5\n")
theta=c(1.5,1.6,1.5,1.6,2.4)
delta=c(0.3,0.4,0.3,0.4,0.2)
cat("\ntheta\n"); print(theta); cat("delta\n"); print(delta)
ifixed=rep(F,2*d)
Finfgumt=irfisherinfo2(nq,ucuts,theta,delta,pcondcop1=pcondgum,pcondcop2=pcondt,
  pconddotcop1=pconddotgum,pconddotcop2=pconddott,dcop1=dgum,dcop2=dtcop,
  iindep=F)
cat("SEs with 2d linking copulas\n")
se1=sqrt(diag(solve(Finfgumt))/nn)
print(se1)
finfgumt=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="gumbelt",nq=nq,
  ifixed,nu2=5,nn=nn)
print(finfgumt$SE)
print(max(abs(finfgumt$finfo - Finfgumt)))
# Ok matches

#============================================================

cat("\nGaussian\n")
theta=c(0.5,0.6,0.5,0.6,0.4)
delta=c(0.3,0.4,0.3,0.4,0.2)
cat("\ntheta\n"); print(theta); cat("delta\n"); print(delta)
ifixed=rep(F,2*d)
Finfn=irfisherinfo2(nq,ucuts,theta,delta,pcondcop1=pcondbvncop,pcondcop2=pcondbvncop,
  pconddotcop1=pconddotbvncop,pconddotcop2=pconddotbvncop,dcop1=dbvncop,dcop2=dbvncop,
  iindep=F)
cat("SEs with 2d linking copulas\n")
se1=sqrt(diag(solve(Finfn))/nn)
print(se1)
finfn=f90irfisherinfo2(c(theta,delta),ucutp=ucuts,copname="gaussian",nq=nq,
  ifixed,nu1=5,nu2=5,nn=nn)
print(finfn$SE)
print(max(abs(finfn$finfo - Finfn)))
# Ok matches
