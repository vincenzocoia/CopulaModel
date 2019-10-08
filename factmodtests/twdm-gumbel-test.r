# tail-weighted dependence measure for factor copulas with Gumbel
library(CopulaModel)

# Gumbel
th1=gum.b2cpar(.7)
th2=gum.b2cpar(.6)
gum0.tw=twdm(pgum,th1,power=6,nq=35)
# default in pfact1gum and pfact2gum changed to 35
gum1.tw=twdm(pfact1gum,c(th1,th1),power=6,nq=35)
gum2.tw=twdm(pfact2gum,matrix(c(th1,th1,th2,th2),2,2),power=6,nq=35)

# theoretical
cat("twdm with power=6\n")
cat("theoretical bivariate Gumbel\n")
print(gum0.tw)
set.seed(12345)
n=4000
gumdat0=rgum(n,th1)
cat("empirical bivariate Gumbel\n")
gum0.ltw=twdm.emp(gumdat0,pow=6)
gum0.utw=twdm.emp(1-gumdat0,pow=6)
print(c(gum0.ltw[1,2],gum0.utw[1,2]))

cat("\ntheoretical 1-factor, 2-factor\n")
print(gum1.tw)
print(gum2.tw)

set.seed(123)
gumdat1=sim1fact(n,c(th1,th1),qcondgum,"gumbel")
set.seed(124)
gumdat2=sim2fact(n,c(th1,th1),c(th2,th2),qcondgum,qcondgum,"gumbel","gumbel")

cat("empirical 1-factor, 2-factor\n")
gum1.ltw=twdm.emp(gumdat1,pow=6)
gum1.utw=twdm.emp(1-gumdat1,pow=6)
print(c(gum1.ltw[1,2],gum1.utw[1,2]))

gum2.ltw=twdm.emp(gumdat2,pow=6)
gum2.utw=twdm.emp(1-gumdat2,pow=6)
print(c(gum2.ltw[1,2],gum2.utw[1,2]))

#============================================================

# asymmetric
th1=gum.b2cpar(.7)
th2=gum.b2cpar(.3)
# default in pfact1gum and pfact2gum changed to 35
gum1.tw=twdm(pfact1gum,c(th1,th2),power=6,nq=35)
gum2.tw=twdm(pfact2gum,matrix(c(th1,th2,th1,th2),2,2),power=6,nq=35)

cat("\nperm. asymmetric\n")
# theoretical
cat("theoretical 1-factor, 2-factor\n")
print(gum1.tw)
print(gum2.tw)

set.seed(123)
gumdat1=sim1fact(n,c(th1,th2),qcondgum,"gumbel")
set.seed(124)
gumdat2=sim2fact(n,c(th1,th2),c(th1,th2),qcondgum,qcondgum,"gumbel","gumbel")

# empirical
cat("empirical 1-factor, 2-factor\n")
gum1.ltw=twdm.emp(gumdat1,pow=6)
gum1.utw=twdm.emp(1-gumdat1,pow=6)
print(c(gum1.ltw[1,2],gum1.utw[1,2]))

gum2.ltw=twdm.emp(gumdat2,pow=6)
gum2.utw=twdm.emp(1-gumdat2,pow=6)
print(c(gum2.ltw[1,2],gum2.utw[1,2]))

