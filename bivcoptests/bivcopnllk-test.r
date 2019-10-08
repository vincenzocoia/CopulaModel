# test of 
# bivcopnllk(param,udat,logdcop,ivect=T,LB=0,UB=1000)
# bivmodnllk(param,xdat,logdcop,logpdf1,cdf1,np1,logpdf2,cdf2,np2,ivect=T,LB,UB)
# with ivect=T and ivect=F

library(CopulaModel)
cpar=depmeas2cpar(.5,type="rhoS","frk")
set.seed(123)
udat=rfrk(200,cpar)
nllk1=bivcopnllk(cpar,udat,logdfrk,ivect=F,LB=0,UB=35)
nllk2=bivcopnllk(cpar,udat,logdfrk,LB=-20,UB=35)
cat(nllk1,nllk2,"\n")

# log of density of exponential distribution
logdexp=function(x,scale) { -log(scale)-x/scale } 
xdat=-log(1-udat) # exponential margin
fnllk1=bivmodnllk(c(1,1,cpar),xdat,logdfrk,logdexp,pexp,1,logdexp,pexp,1,
  ivect=F,LB=c(0,0,-20),UB=c(20,20,35))
fnllk2=bivmodnllk(c(1,1,cpar),xdat,logdfrk,logdexp,pexp,1,logdexp,pexp,1,
  LB=c(0,0,-20),UB=c(20,20,35))
cat(fnllk1,fnllk2,"\n")
