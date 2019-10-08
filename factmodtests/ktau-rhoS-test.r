# bivariate marginal cdf and its conditionals of 1-factor and 2-factor copulas
# check Spearman's rho and Kendall's tau
# then compare with vine simulation

library(CopulaModel)
#source("../R/fact1cop.R")
#source("../R/fact2cop.R")

pcondfact1gum21=function(v,u,param)
{ pcondfact1(v,u,pcondgum,dgum,param,nq=35) }

pcondfact1gum12=function(v,u,param)
{ pcondfact1(v,u,pcondgum,dgum,param[c(2,1)],nq=35) }

dfact1gum=function(u,v,param)
{ dfact1cop(u,v,dgum,param,nq=35) }

pcondfact2gum21=function(v,u,param)
{ pcondfact2(v,u,pcondgum,pcondgum,dgum,dgum,param[,1],param[,2],nq=35) }

pcondfact2gum12=function(v,u,param)
{ pcondfact2(v,u,pcondgum,pcondgum,dgum,dgum,param[c(2,1),1],param[c(2,1),2],nq=35) }

dfact2gum=function(u,v,param)
{ dfact2cop(u,v,pcondgum,dgum,dgum,param[,1],param[,2],nq=35) }

# ============================================================
# 1-factor and 2-factor checks


par1=c(2,1.5)
par2=matrix(c(2,1.5,1.6,1.3),2,2)
u=.8
v=seq(.4,.9,.1)
cat("copderiv for pfact1gum\n")
chkcopderiv(u,v,par1,bcdf=pfact1gum,pcond=pcondfact1gum21,bpdf=dfact1gum,
  str="fact1gum",eps=1.e-5)

# factor k parameters in column k, k=1,2
cat("copderiv for pfact2gum\n")
chkcopderiv(u,v,par2,bcdf=pfact2gum,pcond=pcondfact2gum21,bpdf=dfact2gum,
  str="fact2gum",eps=1.e-5)

# Spearman's rho and Kendall's tau
cat("\nSpearman's rho and Kendall's tau for 1-factor\n")
spear1=rhoS(par1,cop=pfact1gum,zero=0.0001)
tau1=ktau(par1,pcond12=pcondfact1gum12,pcond21=pcondfact1gum21,zero=0.0001)
print(c(spear1,tau1))
# 0.3388759 0.2185302
cat("\nSpearman's rho and Kendall's tau for 2-factor\n")
spear2=rhoS(par2,cop=pfact2gum,zero=0.0001)
tau2=ktau(par2,pcond12=pcondfact2gum12,pcond21=pcondfact2gum21,zero=0.0001)
print(c(spear2,tau2))
# 0.4610530 0.3160495

# compare simulations
n=10000
set.seed(123)
cat("\nSimulations for comparisons for 1-factor\n")
dat1=sim1fact(n,par1,qcondgum,"gumbel",ivect=T)
sspear1=cor(dat1,method="spearman")
stau1=taucor(dat1)
print(c(sspear1[1,2],stau1[1,2])) 

set.seed(123)
cat("\nSimulations for comparisons for 2-factor\n")
dat2=sim2fact(n,par2[,1],par2[,2],qcondgum,qcondgum,ivect=T)
sspear2=cor(dat2,method="spearman")
stau2=taucor(dat2)
print(c(sspear2[1,2],stau2[1,2]))  

