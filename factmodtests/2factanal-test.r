# simulation from 2-factor model given a parametric bivariate copula family
# MLE for 2-factor copula log-likelihood with different starting points
#  based on factanal()

library(CopulaModel)

d=7
be1=c(.7,.6,.7,.6,.7,.6,.7)
be2=c(.4,.4,.4,.4,.3,.3,.3)

cpar1.gum=gum.b2cpar(be1)
cpar2.gum=gum.b2cpar(be2)
cat("copula parameters\n")
cat(cpar1.gum,cpar1.gum,"\n")

d=length(cpar1.gum)
n=300
set.seed(123)
gumdat=sim2fact(n,cpar1.gum,cpar2.gum,qcondgum,qcondgum,"gum","gum")

zdat=nscore(gumdat,iopt=T)
cat("corr of normal scores\n")
rr=cor(zdat)
print(rr)
mvn2f=factanal(covmat=rr,factors=2)
amat=c(mvn2f$loadings)
amat=matrix(amat,d,2)
cat("\nloadings\n")
print(amat)

gl21 = gausslegendre(21)

# MLE
pcmat=load2pcor(amat)
pcmat=c(pcmat)
pcmat[pcmat<=0]=0.01
cparstart=depmeas2cpar(pcmat,"rhoN","gumbel")
dstrgum=list(copname="gumbel",data=gumdat,quad=gl21,repar=0)
cat("start from factanal with transform to Gumbel parameters\n")

out=pdhessminb(cparstart,f90cop2nllk,ifixed=rep(F,2*d),
  dstruct=dstrgum, LB=rep(1,2*d),UB=rep(10,2*d),iprint=T,eps=1.e-5)

np=2*d
for(j in 1:d)
{ rhmat=grotate2(amat,row=j)
  rhmat=c(load2pcor(rhmat))  # second column -> partial corr
  cat("\nj=", j, "loadings rotated and partial\n")
  print(rhmat)
  rhmat[rhmat<=0]=0.01
  stgum=depmeas2cpar(rhmat,"rhoN","gumbel")
  cat("start=", stgum,"\n")
  ifixed=rep(F,np)
  cat("\nQN: copula for variable ", j, " and factor 2 = Cindep\n")
  stgum[d+j]=1.01
  ifixed[d+j]=TRUE
  out=ml2factb(21,stgum,ifixed,gumdat,"gumbel",LB=1,UB=20,prlevel=1,mxiter=20)

  cat("\npdhessmin: copula for variable ", j, " and factor 2 = Cindep\n")
  outj=pdhessminb(stgum,f90cop2nllk,ifixed=ifixed,
    dstruct=dstrgum, LB=rep(1,2*d),UB=rep(10,2*d),iprint=T,eps=1.e-5)
}

