# Check bivariate marginal cdf and its conditionals for tree 3 of
# C-vine and D-vine,
# then check Spearman's rho and Kendall's tau
# then compare with vine simulation

library(CopulaModel)

# copula cdf, conditional cdf C_{2|1}, C_{1|2}, copula pdf
# for bivariate margin from trees 2 and 3 of a C or D vine when edges are
# all Gumbel pair-copulas

ptree2gum=function(u,v,param)
{ ptree2cop(u,v,param,pcondgum,pcondgum,pgum,nq=35) }

pcondtree2gum21=function(v,u,param)
{ pcondtree2(v,u,param,pcondgum,pcondgum,pcondgum,dgum,nq=35) }

pcondtree2gum12=function(v,u,param)
{ pcondtree2(v,u,param[c(2,1,3)],pcondgum,pcondgum,pcondgum,dgum,nq=35) }

dtree2gum=function(u,v,param)
{ dtree2cop(u,v,param,pcondgum,pcondgum,dgum,dgum,dgum,nq=35) }

ptree3gum=function(u,v,param)
{ ptree3cop.cvine(u,v,param,pcondgum,pcondgum,pcondgum,pcondgum,pcondgum,
    pgum,dgum,nq=35) 
}

pcondtree3gum21=function(v,u,param)
{ pcondtree3.cvine(v,u,param,pcondgum,pcondgum,pcondgum,pcondgum,pcondgum,
    pcondgum, dgum,dgum,dgum,nq=35)
}

pcondtree3gum12=function(v,u,param)
{ pcondtree3.cvine(v,u,param[c(1,3,2,5,4,6)],pcondgum,pcondgum,pcondgum,
    pcondgum,pcondgum,pcondgum, dgum,dgum,dgum,nq=35) }

dtree3gum=function(u,v,param)
{ dtree3cop.cvine(u,v,param,pcondgum,pcondgum,pcondgum,pcondgum,pcondgum,
    dgum,dgum,dgum,dgum,dgum,dgum,nq=35) 
}

ptrde3gum=function(u,v,param)
{ ptree3cop.dvine(u,v,param,pcondgum,pcondgum,pcondgum,pcondgum,pcondgum,
    pcondgum,pgum,dgum,nq=35) 
}

pcondtrde3gum21=function(v,u,param)
{ pcondtree3.dvine(v,u,param,pcondgum,pcondgum,pcondgum,pcondgum,pcondgum,
    pcondgum,pcondgum, dgum,dgum,dgum,nq=35)
}

pcondtrde3gum12=function(v,u,param)
{ pcondtree3.dvine(v,u,param[c(3,2,1,5,4,6)],pcondgum,pcondgum,pcondgum,
    pcondgum,pcondgum,pcondgum,pcondgum, dgum,dgum,dgum,nq=35) }

dtrde3gum=function(u,v,param)
{ dtree3cop.dvine(u,v,param,pcondgum,pcondgum,pcondgum,pcondgum,pcondgum,
    pcondgum, dgum,dgum,dgum,dgum,dgum,dgum,nq=35) 
}

# to replace nq with gl object ??

# ============================================================

par2=c(2,1.5,1.3)
par3=c(2,1.7,1.6,1.4,1.3,1.1)
u=.3
v=seq(.4,.9,.1)
cat("copderiv for ptree2gum\n")
chkcopderiv(u,v,par2,bcdf=ptree2gum,pcond=pcondtree2gum21,bpdf=dtree2gum,
   str="tree2gum",eps=1.e-5)
cat("C-vine: copderiv for tree 3\n")
chkcopderiv(u,v,par3,bcdf=ptree3gum,pcond=pcondtree3gum21,bpdf=dtree3gum,
  str="tree3gum.cvine",eps=1.e-5)

cat("D-vine: copderiv for tree 3\n")
#print(ptrde3gum(.1,.1,par3))
#print(ptrde3gum(.1,.2,par3))
#print(pcondtrde3gum21(.1,.2,par3))
#print(pcondtrde3gum21(.3,.2,par3))
#print(pcondtrde3gum21(.3,.4,par3))
chkcopderiv(u,v,par3,bcdf=ptrde3gum,pcond=pcondtrde3gum21,bpdf=dtrde3gum,
  str="tree3gum.dvine",eps=1.e-5)

# Spearman's rho and Kendall's tau
cat("\nSpearman's rho and Kendall's tau for tree 2\n")
spear2=rhoS(par2,cop=ptree2gum,zero=0.0001)
tau2=ktau(par2,pcond12=pcondtree2gum12,pcond21=pcondtree2gum21,zero=0.0001)
print(c(spear2,tau2))

cat("\nSpearman/Kendall for tree 3 C-vine\n")
spear3=rhoS(par3,cop=ptree3gum,zero=0.0001)
tau3=ktau(par3,pcond12=pcondtree3gum12,pcond21=pcondtree3gum21,zero=0.0001)
print(c(spear3,tau3))

cat("\nSpearman/Kendall for tree 3 D-vine\n")
spear3d=rhoS(par3,cop=ptrde3gum,zero=0.0001)
tau3d=ktau(par3,pcond12=pcondtrde3gum12,pcond21=pcondtrde3gum21,zero=0.0001)
print(c(spear3d,tau3d))

# compare simulations
n=10000
cat("\nTree 2 simulations for comparisons\n")
set.seed(123)
udat=matrix(0,n,3)
th3=matrix(0,3,3)
th3[1,2]=par2[1]; th3[1,3]=par2[2]; th3[2,3]=par2[3]
for(i in 1:n)
{ udat[i,]=cvinesim(runif(3),th3,qcondgum,pcondgum) }
sspear=cor(udat,method="spearman")
stau=taucor(udat)
print(c(sspear[2,3],stau[2,3]))

cat("\nC-vine tree 3 simulations for comparisons\n")
set.seed(123)
udat=matrix(0,n,4)
th4=matrix(0,4,4)
th4[1,2]=par3[1]; th4[1,3]=par3[2]; th4[1,4]=par3[3]
th4[2,3]=par3[4]; th4[2,4]=par3[5]; th4[3,4]=par3[6]
for(i in 1:n)
{ udat[i,]=cvinesim(runif(4),th4,qcondgum,pcondgum) }
sspear=cor(udat,method="spearman")
stau=taucor(udat)
print(c(sspear[3,4],stau[3,4]))

cat("\nD-vine tree 3 simulations for comparisons\n")
set.seed(123)
for(i in 1:n)
{ udat[i,]=dvinesim(runif(4),th4,qcondgum,pcondgum) }
sspear=cor(udat,method="spearman")
stau=taucor(udat)
print(c(sspear[1,4],stau[1,4]))

