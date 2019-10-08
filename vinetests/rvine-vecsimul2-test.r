# simulation for R-vine and fit with maximum likelihood 
# pair-copulas with 2-dimensional parameter in tree 1 and
# 1-dimensional copula in trees 2,3

library(CopulaModel)

# Example 2  
# BB1/Frank/Frank
d=7
qcondnames2=c("qcondbb1","qcondfrk","qcondfrk")
pcondnames2=c("pcondbb1","pcondfrk","pcondfrk")
parvec2=c(.8,2,.7,1.6,.3,2.1,.6,1.5,.9,2.3,.9,2.4, 3,4,3,3,3, 2,2,2,2)
np2=matrix(0,d,d)
np2[1,2:d]=2
np2[2,3:d]=1
np2[3,4:d]=1
par1=matrix(0,d,d)
par2=matrix(0,d,d)
par1[1,2:d]=parvec2[seq(1,11,2)]
par2[1,2:d]=parvec2[seq(2,12,2)]
par1[2,3:d]=parvec2[13:17]
par1[3,4:d]=parvec2[18:21]


#A1=vnum2array(d,300)
A2=vnum2array(d,100)
print(A2)
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    1    1    2    3    1    3
#[2,]    0    2    2    1    1    2    1
#[3,]    0    0    3    3    2    4    5
#[4,]    0    0    0    4    4    3    2
#[5,]    0    0    0    0    5    5    4
#[6,]    0    0    0    0    0    6    6
#[7,]    0    0    0    0    0    0    7
set.seed(123)
nsim=200
nsim=2000
udat2=rvinesimvec(nsim,A2,parvec2,np2,qcondnames2,pcondnames2,iprint=F)
#print(udat2)

cat("\nmaximum likelihood\n");
ifixed2=rep(F,3*d-6+d-1)
parfixed=NA
lb2=rep(c(0,1),d-1)
lb2=c(lb2,rep(-20,d-2))
lb2=c(lb2,rep(-20,d-3))
ub2=rep(10,2*d-2)
ub2=c(ub2,rep(30,d-2))
ub2=c(ub2,rep(30,d-3))

lognvec2=c("logdbb1","logdfrk","logdfrk")
pcnvec2=c("pcondbb1","pcondfrk","pcondfrk")
start2=rep(c(.7,2.4),d-1)
start2=c(start2,rep(2,d-2))
start2=c(start2,rep(1.2,d-3))

vineml2=nlm(rvinenllk.trunc,start2,udat=udat2,A=A2,logdcopnames=lognvec2,
      pcondnames=pcnvec2,np2,ifixed2,parfixed,
      LB=lb2,UB=ub2,hessian=T,print.level=1,iterlim=60)
# looks OK compared with parvec2
acov=solve(vineml2$hessian)
cat("SEs\n")
print(sqrt(diag(acov)))

#============================================================

cat("\nempirical correlation and tau for tree 1\n")
for(j in 2:d)
{ erho=cor(udat2[,A2[1,j]],udat2[,j])
  etau=taucor(udat2[,A2[1,j]],udat2[,j])
  cat(j,erho,etau, "\n")
}
cat("\n") 
cat("theoretical rhoS and tau for tree 1\n")  # Ok
for(j in 2:d) 
{ j1=j-1; jj=c(2*j1-1,2*j1)
  rho=rhoS(parvec2[jj],pbb1)
  tau=ktau(parvec2[jj],icond=T,pcondbb1,pcondbb1)
  cat(j,rho,tau, "\n")
}
cat("\n") 

cat("empirical correlations for tree 2\n")
for(j in 3:d)
{ erho=cor(udat2[,A2[2,j]],udat2[,j])
  etau=taucor(udat2[,A2[2,j]],udat2[,j])
  cat(j,erho,etau, "\n")
}
cat("\n") 

# copula cdf, conditional cdf C_{2|1}, C_{1|2}
# for bivariate margin from tree 2 of a vine when edges are
# BB1 and Frank pair-copulas
pcondfrk2=function(v,u,param2) { pcondfrk(v,u,param2[1]) }
pfrk2=function(u,v,param2) { pfrk(u,v,param2[1]) }
ptree2bb1frk=function(u,v,param)
{ ptree2cop(u,v,param,pcondbb1,pcondbb1,pfrk2,nq=35) }
pcondtree2bb1frk21=function(v,u,param)
{ pcondtree2(v,u,param,pcondbb1,pcondbb1,pcondfrk2,dbb1,nq=35) }
pcondtree2bb1frk12=function(v,u,param)
{ pcondtree2(v,u,param[c(2,1,3),],pcondbb1,pcondbb1,pcondfrk2,dbb1,nq=35) } 

cat("theoretical rhoS and tau for tree 2\n") # OK
for(j in 3:d)
{ j0=A2[1,j]
  j1=A2[2,j]; j2=A2[j,j]
  jmin=min(j0,j1); jmax=max(j0,j1)
  param3=matrix(c(par1[1,jmax],par1[1,j],par1[2,j],
    par2[1,jmax],par2[1,j],par2[2,j]),3,2)
  rho=rhoS(param3,ptree2bb1frk)
  tau=ktau(param3,pcond12=pcondtree2bb1frk12,pcond21=pcondtree2bb1frk21)
  cat(j,rho,tau,"\n")
}

