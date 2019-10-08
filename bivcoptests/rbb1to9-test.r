# random pairs, 2-parameter copulas: difficult cases with strong dependence
library(CopulaModel)

n=1000

cat("\nBB1\n")
set.seed(123)
param=c(2.5,3.1)
udat=rbb1(n,param,icheck=T)
tau=bb1.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB2\n")
set.seed(123)
param=c(2,1.5)
udat=rbb2(n,param,icheck=T)
tau=bb2.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB3\n")
set.seed(123)
param=c(4,2)
udat=rbb3(n,param,icheck=T)
tau=bb3.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB7\n")
set.seed(123)
param=c(1.8,2.1)
udat=rbb7(n,param,icheck=T)
tau=bb7.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB9\n")
set.seed(123)
param=c(4.,3)
udat=rbb9(n,param,icheck=T)
tau=bb9.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB4\n")
set.seed(123)
param=c(1.5,2.3)
udat=rbb4(n,param,icheck=T)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tauo))

cat("\nBB6\n")
set.seed(123)
param=c(4.5,4.5)
udat=rbb6(n,param,icheck=T)
tau=bb6.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB8\n")
set.seed(123)
param=c(5.5,.9) # fails for bb8.cpar2tau
udat=rbb8(n,param,icheck=T)
tau=bb8.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

cat("\nBB5\n")
set.seed(123)
#param=c(1.5,3.5) # problems C code hasn't been fixed yet like qcondbb5
#param=c(1.5,2.5) # problems
param=c(1.2,2.5) # problems
udat=rbb5(n,param,icheck=T)
tau=bb5.cpar2tau(param)
tauo=taucor(udat[,1],udat[,2])
print(c(param,tau,tauo))

