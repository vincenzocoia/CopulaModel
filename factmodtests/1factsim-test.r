# simulation from 1-factor model given a parametric bivariate copula family
# MLE for 1-factor copula log-likelihood

library(CopulaModel)

# create data sets to be use for comparing maximum likelihood estimation
# takes subsets of columns and rows for faster computations

# dependence parameters chosen so that Spearman rho with latent variable are:
#  .5, .5, .6, .6, .7, .7, .8, .8, .9, .9
#cpar.frk=c(3.45,3.45,4.47,4.47,5.82,5.82,7.90,7.90,12.2,12.2)
#cpar.mtcj=c(1.06,1.06,1.51,1.51,2.14,2.14,3.19,3.19,5.56,5.56)
#cpar.gum=c(1.54,1.54,1.75,1.75,2.07,2.07,2.58,2.58,3.73,3.73)
#cpar.gumr=c(1.54,1.54,1.75,1.75,2.07,2.07,2.58,2.58,3.73,3.73)
#  .9, .5, .6, .6, .7
cpar.frk=c(12.2,3.45,4.47,4.47,5.82)
cpar.mtcj=c(5.56,1.06,1.51,1.51,2.14)
cpar.gum=c(3.73,1.54,1.75,1.75,2.07)
cpar.gumr=c(3.73,1.54,1.75,1.75,2.07)
lmbb1=matrix(c(.3,.4,.5,.3,.5, .6,.6,.6,.7,.7),5,2)
cpar.bb1=lmbb1
for(i in 1:nrow(lmbb1))
{ cpar.bb1[i,]=bb1.lm2cpar(lmbb1[i,]) }

d=length(cpar.gum)
n=300
set.seed(123)
frkdat=sim1fact(n,cpar.frk,qcondfrk,"frk")
set.seed(123)
mtcjdat=sim1fact(n,cpar.mtcj,qcondmtcj,"mtcj")
set.seed(123)
gumdat=sim1fact(n,cpar.gum,qcondgum,"gum")
set.seed(123)
gumrdat=sim1fact(n,cpar.gumr,qcondgumr,"gumr")
set.seed(123)
bb1dat=sim1fact(n,cpar.bb1,qcondbb1,"bb1")

cat("\nfrkdat correlations\n")
print(cor(frkdat))
cat("\nmtcjdat correlations\n")
print(cor(mtcjdat))
cat("\ngumdat correlations\n")
print(cor(gumdat))
cat("\ngumrdat correlations\n")
print(cor(gumrdat))
cat("\nbb1dat correlations\n")
print(cor(bb1dat))

#pairs(qnorm(gumdat))
#pairs(qnorm(gumrdat))
#pairs(qnorm(frkdat))
#pairs(qnorm(mtcjdat))
#pairs(qnorm(bb1dat))

gl21 = gausslegendre(21)
# MLE
# standalone R interface

cat("\nGumbel 1-factor, pure R and then f90\n")
out.gum=ml1fact(nq=21,cpar.gum,gumdat,dgum,LB=1,UB=20,prlevel=1,mxiter=100)

#iteration = 12
#Parameter:
#[1] 3.829970 1.454046 1.936581 1.841941 1.948722
#Function Value
#[1] -264.3501
#Gradient:
#[1]  1.893806e-05  1.465998e-05 -1.943133e-05  7.545420e-05  8.756711e-05
#
#Relative gradient close to zero.
#Current iterate is probably solution.
#
#MLE: 
#[1] 3.829970 1.454046 1.936581 1.841941 1.948722
#SEs: 
#[1] 0.52472201 0.07345783 0.11245366 0.10150970 0.11457481
#nllk: 
#[1] -264.3501

# R interface via f90
dstrgum=list(copname="gumbel",data=gumdat,quad=gl21,repar=0)
out=pdhessminb(cpar.gum,f90cop1nllk,ifixed=rep(F,d),dstruct=dstrgum,
  LB=rep(1,d),UB=rep(20,d),iprint=T,eps=1.e-5)
cat(cpar.gum,"\n")
cat("\n============================================================\n")

#        param              
#[1,] 3.829965 -6.689663e-08
#[2,] 1.454047  1.115680e-08
#[3,] 1.936582 -2.230966e-10
#[4,] 1.841941  8.614938e-09
#[5,] 1.948722 -1.460217e-08
#5 -264.3501 1.823201e-08 

cat("\nreflected Gumbel 1-factor, pure R and then f90\n")
out.gumr=ml1fact(nq=21,cpar.gumr,gumrdat,dgumr,LB=1,UB=20,prlevel=1,mxiter=100)
dstrgumr=list(copname="gumbel",data=1-gumrdat,quad=gl21,repar=0)
out=pdhessminb(cpar.gumr,f90cop1nllk,ifixed=rep(F,d),dstruct=dstrgumr,
  LB=rep(1,d),UB=rep(20,d),iprint=T,eps=1.e-5)
cat(cpar.gumr,"\n")
cat("\n============================================================\n")

cat("\nMTCJ 1-factor\n")
out.mtcj=ml1fact(nq=21,cpar.mtcj,mtcjdat,dmtcj,LB=0,UB=20,prlevel=1,mxiter=100)
cat(cpar.mtcj,"\n")
cat("\n============================================================\n")

cat("\nFrank 1-factor, pure R and then f90\n")
out.frk=ml1fact(nq=21,cpar.frk,frkdat,dfrk,LB=-30,UB=30,prlevel=1,mxiter=100)
dstrfrk=list(copname="frk",data=frkdat,quad=gl21,repar=0)
out=pdhessminb(cpar.frk,f90cop1nllk,ifixed=rep(F,d),dstruct=dstrfrk,
  LB=rep(-30,d),UB=rep(30,d),iprint=T,eps=1.e-5)
cat(cpar.frk,"\n")
cat("\n============================================================\n")

cat("\nBB1 1-factor, pure R and then f90\n")
cpar.bb1=c(t(cpar.bb1))
print(cpar.bb1)
out.bb1=ml1fact(nq=21,cpar.bb1,bb1dat,dbb1,LB=rep(c(0,1),d),UB=10,prlevel=1,mxiter=100)
dstrbb1=list(copname="bb1",data=bb1dat,quad=gl21,repar=0)
out=pdhessminb(cpar.bb1,f90cop1nllk,ifixed=rep(F,2*d),dstruct=dstrbb1,
  LB=rep(c(0,1),d),UB=rep(10,2*d),iprint=T,eps=1.e-5)
cat(cpar.bb1,"\n")
cat("\n============================================================\n")


