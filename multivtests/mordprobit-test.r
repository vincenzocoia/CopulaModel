# check that probabilities from pmfmordprobit add to 1
# when pmvnorm in library(mvtnorm) with quasi-random method is used.

library(CopulaModel)
library(mvtnorm)  
#source("../R/mprobit.R")

ncateg=3
ucuts3=matrix(c(.4,.5,.4, .7,.8,.6),ncateg-1,3,byrow=T)
ucuts4=matrix(c(.4,.5,.4,.3, .7,.8,.6,.6),ncateg-1,4,byrow=T)
ucuts5=matrix(c(.4,.5,.4,.3,.2, .7,.8,.6,.6,.6),ncateg-1,5,byrow=T)
zcuts3=qnorm(ucuts3)
zcuts3=rbind(rep(-6,3),zcuts3,rep(6,3))
zcuts4=qnorm(ucuts4)
zcuts4=rbind(rep(-6,4),zcuts4,rep(6,4))
zcuts5=qnorm(ucuts5)
zcuts5=rbind(rep(-6,5),zcuts5,rep(6,5))

# results below: examples without fixed seed for pmvnorm within pmfmordprobit
cat("\nfirst run of AR1, rh=.6\n")
rh=.6  # first run of AR1
rmat3=toeplitz(rh^(0:2))
rmat4=toeplitz(rh^(0:3))
rmat5=toeplitz(rh^(0:4))
pmf3=pmfmordprobit(zcuts3,rmat3,iprint=T,ifixseed=T)  # 0.999996 = sum(pmf3) if seed=12345
pmf4=pmfmordprobit(zcuts4,rmat4,iprint=T,ifixseed=T)  # 0.9999504
pmf5=pmfmordprobit(zcuts5,rmat5,iprint=T,ifixseed=T)  # 1.000018
print(c(sum(pmf3),sum(pmf4),sum(pmf5))) # 0.9998280 0.9999901 0.9999627
# Ok to 4 decimal places

cat("\nsecond run of AR1, rh=.7\n")
rh=.7  # second run
rmat3=toeplitz(rh^(0:2))
rmat4=toeplitz(rh^(0:3))
rmat5=toeplitz(rh^(0:4))
pmf3=pmfmordprobit(zcuts3,rmat3,iprint=F,ifixseed=T) # 0.9999956
pmf4=pmfmordprobit(zcuts4,rmat4,iprint=F,ifixseed=T) # 0.999957
pmf5=pmfmordprobit(zcuts5,rmat5,iprint=F,ifixseed=T) # 1.000014
print(c(sum(pmf3),sum(pmf4),sum(pmf5))) # 0.9999405 1.0000088 1.0000022
# Ok to 4 decimal places

cat("\nthird run of AR1, rh=.8\n")
rh=.8  # third run
rmat3=toeplitz(rh^(0:2))
rmat4=toeplitz(rh^(0:3))
rmat5=toeplitz(rh^(0:4))
pmf3=pmfmordprobit(zcuts3,rmat3,iprint=F,ifixseed=T) # 1.000019
pmf4=pmfmordprobit(zcuts4,rmat4,iprint=F,ifixseed=T) # 0.9999326
pmf5=pmfmordprobit(zcuts5,rmat5,iprint=F,ifixseed=T) # 0.9999966
print(c(sum(pmf3),sum(pmf4),sum(pmf5))) # 0.9999470 0.9999428 1.0000392

cat("\nfirst run of AR2, rh1=.5, rh2=.6\n")
rh1=.5; rh2=.6  # first run of AR2
acf10=ar2acf(rh1,rh2,10)
r2mat3=toeplitz(acf10[1:3])
r2mat4=toeplitz(acf10[1:4])
r2mat5=toeplitz(acf10[1:5])
pmf3=pmfmordprobit(zcuts3,r2mat3,iprint=F,ifixseed=T) # 1.000022
pmf4=pmfmordprobit(zcuts4,r2mat4,iprint=F,ifixseed=T) # 0.9999544
pmf5=pmfmordprobit(zcuts5,r2mat5,iprint=F,ifixseed=T) # 0.9999747
print(c(sum(pmf3),sum(pmf4),sum(pmf5))) # 0.9999462 0.9999589 0.9999777

cat("\nsecond run of AR2, rh1=.6, rh2=.5\n")
rh1=.6; rh2=.5  # second run
acf10=ar2acf(rh1,rh2,10)
r2mat3=toeplitz(acf10[1:3])
r2mat4=toeplitz(acf10[1:4])
r2mat5=toeplitz(acf10[1:5])
pmf3=pmfmordprobit(zcuts3,r2mat3,iprint=F,ifixseed=T) # 1.000009
pmf4=pmfmordprobit(zcuts4,r2mat4,iprint=F,ifixseed=T) # 0.9999107
pmf5=pmfmordprobit(zcuts5,r2mat5,iprint=F,ifixseed=T) # 0.9999716
print(c(sum(pmf3),sum(pmf4),sum(pmf5))) # 1.0000354 0.9999981 1.0000152

cat("\nthird run of AR2, rh1=.7, rh2=.6\n")
rh1=.7; rh2=.6  # third run
acf10=ar2acf(rh1,rh2,10)
r2mat3=toeplitz(acf10[1:3])
r2mat4=toeplitz(acf10[1:4])
r2mat5=toeplitz(acf10[1:5])
pmf3=pmfmordprobit(zcuts3,r2mat3,iprint=F,ifixseed=T) # 1.000011
pmf4=pmfmordprobit(zcuts4,r2mat4,iprint=F,ifixseed=T) # 0.999906
pmf5=pmfmordprobit(zcuts5,r2mat5,iprint=F,ifixseed=T) # 0.9999607
print(c(sum(pmf3),sum(pmf4),sum(pmf5))) # 1.000011 1.000024 1.000068

