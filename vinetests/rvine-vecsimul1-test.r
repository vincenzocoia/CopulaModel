# simulation for R-vine and fit with maximum likelihood 
# all pair-copulas with 1-dimensional parameter 

library(CopulaModel)

# Example 1
# Gumbel/Frank
d=7
qcondnames1=c("qcondgum","qcondfrk","qcondfrk")
pcondnames1=c("pcondgum","pcondfrk","pcondfrk")
parvec1=c(1.4,1.4,1.6,1.6,2.,2., 2.6,3.4,4.4,4.4,5.8, 2,2,2,2)
np1=matrix(0,d,d)
np1[1,2:d]=1
np1[2,3:d]=1
np1[3,4:d]=1
parmat=matrix(0,d,d)
parmat[1,2:d]=parvec1[1:6]
parmat[2,3:d]=parvec1[7:11]
parmat[3,4:d]=parvec1[12:15]

A1=vnum2array(d,300)
print(A1)
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    2    1    2    1
#[2,]    0    2    2    1    3    1    3
#[3,]    0    0    3    3    2    4    5
#[4,]    0    0    0    4    4    3    2
#[5,]    0    0    0    0    5    5    4
#[6,]    0    0    0    0    0    6    6
#[7,]    0    0    0    0    0    0    7
set.seed(123)
nsim=200
nsim=2000
udat1=rvinesimvec(nsim,A1,parvec1,np1,qcondnames1,pcondnames1,iprint=F)
#print(udat1)

cat("\nempirical correlation and tau for tree 1\n")
for(j in 2:d)
{ erho=cor(udat1[,A1[1,j]],udat1[,j])
  etau=taucor(udat1[,A1[1,j]],udat1[,j])
  cat(j,erho,etau, "\n")
}
cat("\n") 
cat("theoretical rhoS and tau for tree 1\n")
for(j in 2:d) 
{ rho=rhoS(parvec1[j-1],pgum)
  tau=ktau(parvec1[j-1],icond=T,pcondgum,pcondgum)
  cat(j,rho,tau, "\n")
}
cat("\n") 

cat("empirical correlation and tau for tree 2\n")
for(j in 3:d)
{ erho=cor(udat1[,A1[2,j]],udat1[,j])
  etau=taucor(udat1[,A1[2,j]],udat1[,j])
  cat(j,erho,etau, "\n")
}
cat("\n") 

# copula cdf, conditional cdf C_{2|1}, C_{1|2}
# for bivariate margin from tree 2 of a vine when edges are
# Gumbel and Frank pair-copulas
ptree2gumfrk=function(u,v,param)
{ ptree2cop(u,v,param,pcondgum,pcondgum,pfrk,nq=35) }
pcondtree2gumfrk21=function(v,u,param)
{ pcondtree2(v,u,param,pcondgum,pcondgum,pcondfrk,dgum,nq=35) }
pcondtree2gumfrk12=function(v,u,param)
{ pcondtree2(v,u,param[c(2,1,3)],pcondgum,pcondgum,pcondfrk,dgum,nq=35) }

cat("theoretical rhoS and tau for tree 2\n")
for(j in 3:d)
{ j0=A1[1,j]
  j1=A1[2,j]; j2=A1[j,j]
  jmin=min(j0,j1); jmax=max(j0,j1)
  param3=c(parmat[1,jmax],parmat[1,j],parmat[2,j])
  rho=rhoS(param3,ptree2gumfrk)
  tau=ktau(param3,pcond12=pcondtree2gumfrk12,pcond21=pcondtree2gumfrk21)
  cat(j,rho,tau,"\n")
}

cat("\nmaximum likelihood\n");
ifixed1=rep(F,3*d-6)
parfixed=NA
lb1=c(rep(1,d-1),rep(-20,2*d-5))
ub1=c(rep(10,d-1),rep(30,2*d-5))
lognvec1=c("logdgum","logdfrk","logdfrk")
pcnvec1=c("pcondgum","pcondfrk","pcondfrk")
start1=c(rep(1.5,d-1),rep(2,2*d-5))

vineml1=nlm(rvinenllk.trunc,start1,udat=udat1,A=A1,logdcopnames=lognvec1,
      pcondnames=pcnvec1,np1,ifixed1,parfixed,
      LB=lb1,UB=ub1,hessian=T,print.level=1,iterlim=70)
acov=solve(vineml1$hessian)
cat("SEs\n")
print(sqrt(diag(acov)))

