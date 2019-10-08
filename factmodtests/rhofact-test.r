# spearman rho matrix for 1-factor and 2-factor copulas
# compare integration and simulation

library(CopulaModel)

bevec=c(.8,.7,.6,.5,.5)
thfrk=frk.b2cpar(bevec)
bevec2=c(.0001,.6,.6,.6,.7)
thfrk2=frk.b2cpar(bevec2)

n=10000  # within .016
n=40000  # within .003
set.seed(123)
frk1dat=sim1fact(n,thfrk,qcondfrk,"frank")

set.seed(123)
frk2dat=sim2fact(n,thfrk,thfrk2,qcond1=qcondfrk,qcond2=qcondfrk,"frank","frank")

options(digits=5)

# 1-factor
cat("1-factor with Frank\n")
r1=cor(frk1dat)
r1s=cor(frk1dat,method="spearman")
cat("empirical correlation matrix\n")
print(r1)
cat("empirical Spearman correlation matrix\n")
print(r1s) # should be close to above
# theoretical
cat("theoretical Spearman correlation matrix\n")
r1t=srho1fact(thfrk,pcondfrk,nq=15)
print(r1t)
print(max(abs(r1-r1t)))

# 2-factor
cat("\n2-factor with Frank\n")
r2=cor(frk2dat)
r2s=cor(frk2dat,method="spearman")
cat("empirical correlation matrix\n")
print(r2)
cat("empirical Spearman correlation matrix\n")
print(r2s)
cat("theoretical Spearman correlation matrix\n")
r2t=srho2fact(thfrk,thfrk2,pcondfrk,pcondfrk,nq=15)
print(r2t)
print(max(abs(r2-r2t)))

