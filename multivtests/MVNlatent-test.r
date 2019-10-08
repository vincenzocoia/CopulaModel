# estimation of univariate and bivariate parameters for
# latent multivariate normal/Gaussian model with discrete response

library(CopulaModel)
#source("../R/discreteresponse.R")

options(digits=5)

# count regression model with common betas for all margins: longitudinal count
cat("NB1 regression model with common betas\n")
# data(rwmsubset)
#rwm=read.table("../data/rwmsubset.tab",header=T)
data(rwmsubset)
rwm=rwmsubset
rwm$agec=(rwm$age-50)/10
rwm$ageq=(rwm$agec)^2
rwm$handfra=rwm$handper/100
xdat=cbind(rwm$sex,rwm$agec,rwm$ageq,rwm$hsat,rwm$handfra,rwm$univ)
xdat=as.matrix(xdat)
ydat=rwm$docvis
nc=ncol(xdat)
nc1=nc+1
out=MVNlatent1(ydat,xdat,nrep=5,upmf=nb1pmf,ucdf=nb1cdf,upmfcdf=nb1pmfcdf,
  mx=7,ustart=c(1.7,.3,.2,.1,-.2,.7,-.5,1.5),
  LB=c(rep(-20,nc1),0),UB=rep(10,nc1+1),prlevel=1)
print(out$uparam)

cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("\nGP1 regression model with uncommon betas for different margins\n")
# GP1 count regression models for each margin: uncommon regression coefficients 
#data(kzrepmeas)
kz=read.table("../data/kzrepmeas.tab",header=T)
kz$agehun=kz$age/100
xdat=cbind(kz$agehun,kz$sex,kz$msmok)
ydat=kz[,6]
nrep=4
outgp1=MVNlatent2(ydat,xdat,nrep,unllks=rep("gp1nllk",4),ucdfs=rep("gp1cdf",4),
  upmfcdfs=rep("gp1pmfcdf",4),
  mx=3,ustart=c(0,0,0,0,1),LB=c(-20,-20,-20,-20,0),UB=rep(10,5),prlevel=0)
print(outgp1)
