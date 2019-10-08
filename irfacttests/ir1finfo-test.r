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
#pcondcop=pcondbvtcop
#pconddotcop=pconddott
#dcop=dtcop
################################################################################

# 1-factor copula models

# gumbel
cat("Gumbel\n")
theta=c(1.5,1.6,1.5,1.6,2.4)
cat("\ntheta\n"); print(theta)
Finfgum=irfisherinfo1(nq,ucuts,theta,pcondcop=pcondgum,
  pconddotcop=pconddotgum,dcop=dgum)
cat("SEs with d copulas\n")
se1=sqrt(diag(solve(Finfgum))/nn)
print(se1)
finfgum=f90irfisherinfo1(theta,ucutp=ucuts,copname="gumbel",nq=nq,nn=nn)
cat(finfgum$SE,"\n")
print(max(abs(finfgum$finfo - Finfgum)))

#============================================================

dfdefault=5
cat("\nt5\n")
theta=c(0.5,0.6,0.5,0.6,0.4)
cat("\ntheta\n"); print(theta)
Finft=irfisherinfo1(nq,ucuts,theta,pcondcop=pcondbvtcop,
  pconddotcop=pconddott,dcop=dtcop)
cat("SEs with d copulas\n")
se1=sqrt(diag(solve(Finft))/nn)
print(se1)
finft=f90irfisherinfo1(theta,ucutp=ucuts,copname="t",nq=nq,nu1=5,nn=nn)
print(finft$SE)
print(max(abs(finft$finfo - Finft)))


#============================================================

cat("\nGaussian\n")
theta=c(0.5,0.6,0.5,0.6,0.4)
cat("\ntheta\n"); print(theta)
Finfg=irfisherinfo1(nq,ucuts,theta,pcondcop=pcondbvncop,
  pconddotcop=pconddotbvncop,dcop=dbvncop)
cat("SEs with d copulas\n")
se1=sqrt(diag(solve(Finfg))/nn)
print(se1)
finfg=f90irfisherinfo1(theta,ucutp=ucuts,copname="gaussian",nq=nq,nn=nn)
print(finfg$SE)
print(max(abs(finfg$finfo - Finfg)))

