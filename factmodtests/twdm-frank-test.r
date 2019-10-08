# tail-weighted dependence measure for factor copulas with Frank 
library(CopulaModel)

# Frank
th1=frk.b2cpar(.7)
th2=frk.b2cpar(.6)
frk0.tw=twdm(pfrk,th1,power=6,nq=15)
# default in pfact1frk and pfact2frk changed to 35
frk1.tw=twdm(pfact1frk,c(th1,th1),power=6,nq=25)
frk2.tw=twdm(pfact2frk,matrix(c(th1,th1,th2,th2),2,2),power=6,nq=25)

# theoretical
cat("twdm with power=6\n")
cat("theoretical bivariate Frank\n")
print(frk0.tw)
set.seed(12345)
n=4000
frkdat0=rfrk(n,th1)
cat("empirical bivariate Frank\n")
frk0.ltw=twdm.emp(frkdat0,pow=6)
frk0.utw=twdm.emp(1-frkdat0,pow=6)
print(c(frk0.ltw[1,2],frk0.utw[1,2]))

cat("\ntheoretical 1-factor, 2-factor\n")
print(frk1.tw)
print(frk2.tw)

set.seed(123)
frkdat1=sim1fact(n,c(th1,th1),qcondfrk,"frank")
set.seed(124)
frkdat2=sim2fact(n,c(th1,th1),c(th2,th2),qcondfrk,qcondfrk,"frank","frank")

cat("empirical 1-factor, 2-factor\n")
frk1.ltw=twdm.emp(frkdat1,pow=6)
frk1.utw=twdm.emp(1-frkdat1,pow=6)
print(c(frk1.ltw[1,2],frk1.utw[1,2]))

frk2.ltw=twdm.emp(frkdat2,pow=6)
frk2.utw=twdm.emp(1-frkdat2,pow=6)
print(c(frk2.ltw[1,2],frk2.utw[1,2]))

#============================================================

# asymmetric
th1=frk.b2cpar(.7)
th2=frk.b2cpar(.3)
# default in pfact1frk and pfact2frk changed to 35
frk1.tw=twdm(pfact1frk,c(th1,th2),power=6,nq=25)
frk2.tw=twdm(pfact2frk,matrix(c(th1,th2,th1,th2),2,2),power=6,nq=25)

cat("\nperm. asymmetric\n")
# theoretical
cat("theoretical 1-factor, 2-factor\n")
print(frk1.tw)
print(frk2.tw)

set.seed(123)
frkdat1=sim1fact(n,c(th1,th2),qcondfrk,"frank")
set.seed(124)
frkdat2=sim2fact(n,c(th1,th2),c(th1,th2),qcondfrk,qcondfrk,"frank","frank")

# empirical
cat("empirical 1-factor, 2-factor\n")
frk1.ltw=twdm.emp(frkdat1,pow=6)
frk1.utw=twdm.emp(1-frkdat1,pow=6)
print(c(frk1.ltw[1,2],frk1.utw[1,2]))

frk2.ltw=twdm.emp(frkdat2,pow=6)
frk2.utw=twdm.emp(1-frkdat2,pow=6)
print(c(frk2.ltw[1,2],frk2.utw[1,2]))
