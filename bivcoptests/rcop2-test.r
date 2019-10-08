# random pairs, 2-parameter copulas
library(CopulaModel)

n=100000
n=1000

# check kendall tau for BB1
set.seed(123)
param=c(.3,1.1)
out=rbb1(n,param,icheck=T)
tau=bb1.cpar2tau(param)
tauo=taucor(out[,1],out[,2])
print(c(param,tau,tauo))

set.seed(123)
param=c(.8,1.1)
out=rbb1(n,param,icheck=T)
tau=bb1.cpar2tau(param)
tauo=taucor(out[,1],out[,2])
print(c(param,tau,tauo))

set.seed(123)
param=c(1.5,1.1)
out=rbb1(n,param,icheck=T)
tau=bb1.cpar2tau(param)
tauo=taucor(out[,1],out[,2])
print(c(param,tau,tauo))

set.seed(123)
param=c(1.5,2.1)
out=rbb1(n,param,icheck=T)
tau=bb1.cpar2tau(param)
tauo=taucor(out[,1],out[,2])
print(c(param,tau,tauo))


set.seed(123)
param=c(1.3,.1)
out=rbb7(n,param,icheck=T)

set.seed(123)
param=c(1.8,.1)
out=rbb7(n,param,icheck=T)

set.seed(123)
param=c(2.5,.1)
out=rbb7(n,param,icheck=T)

set.seed(123)
param=c(1.3,1.1)
out=rbb7(n,param,icheck=T)

