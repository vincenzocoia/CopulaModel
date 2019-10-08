# bivariate margins for vine for trees 2 and 3 
# twdm = tail-weighted dependence measure
# compare with vine simulation

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

# ============================================================

par2=c(2,1.5,1.3)
par3=c(2,1.7,1.6,1.4,1.3,1.1)

cat("twdm (lower and upper) for tree 2 C-vine\n")
tem=twdm(ptree2gum,par2,power=6,nq=25)
print(tem)
# 0.1978552 0.5829328

cat("\ntwdm (lower and upper) for tree 3 C-vine\n")
tem=twdm(ptree3gum,par3,power=6,nq=25)
print(tem)
# 0.1482317 0.5531348

cat("\ntwdm (lower and upper) for tree 3 D-vine\n")
tem=twdm(ptrde3gum,par3,power=6,nq=25)
print(tem)
# 0.1712516 0.5810486

cat("\nsimulations\n")
n=10000
set.seed(123)
udat=matrix(0,n,3)
th3=matrix(0,3,3)
th3[1,2]=par2[1]; th3[1,3]=par2[2]; th3[2,3]=par2[3]
for(i in 1:n)
{ udat[i,]=cvinesim(runif(3),th3,qcondgum,pcondgum) }
ltail=twdm.emp(udat,power=6)
utail=twdm.emp(1-udat,power=6)
cat("\nempirical twdm (lower and upper) for tree 2 C-vine\n")
print(ltail)
print(utail)

set.seed(123)
udat=matrix(0,n,4)
th4=matrix(0,4,4)
th4[1,2]=par3[1]; th4[1,3]=par3[2]; th4[1,4]=par3[3]
th4[2,3]=par3[4]; th4[2,4]=par3[5]; th4[3,4]=par3[6]
for(i in 1:n)
{ udat[i,]=cvinesim(runif(4),th4,qcondgum,pcondgum) }
ltail=twdm.emp(udat,power=6)
utail=twdm.emp(1-udat,power=6)
cat("\nempirical twdm (lower and upper) for tree 3 C-vine\n")
print(ltail)
print(utail)

set.seed(123)
for(i in 1:n)
{ udat[i,]=dvinesim(runif(4),th4,qcondgum,pcondgum) }
ltail=twdm.emp(udat,power=6)
utail=twdm.emp(1-udat,power=6)
cat("\nempirical twdm (lower and upper) for tree 3 D-vine\n")
print(ltail)
print(utail)


