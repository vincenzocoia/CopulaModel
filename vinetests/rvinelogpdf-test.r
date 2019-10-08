# check that R-vine becomes specialized code for C-vine and D-vine
# with appropriate vine arrays

# rvinelogpdf,cvinelogpdf,dvinelogpdf
# rvinelogpdf=function(uu,A,th,logdcop,pcond,iprint=F)
# cvinelogpdf=function(uu,th,logdcop,pcond,iprint=F)
# dvinelogpdf=function(uu,th,logdcop,pcond,iprint=F)
# these code only for 1-parameter copulas; see 
# rvinenllk.trunc for general interface

library(CopulaModel)
#source("../R/logcopden.R")

C=Cvinearray(5)
D=Dvinearray(5)

uu=c(.1,.3,.4,.5,.7)
th= matrix(c(0,0,0,0,0, 1.5,0,0,0,0, 1.5,1.2,0,0,0, 1.5,1.2,1.3,0,0, 1.5,1.2,1.3,1.4,0), 5,5)

pcondcop=pcondgum; logdcop=logdgum
pcondcop=pcondfrk; logdcop=logdfrk
pcondcop=pcondpla; logdcop=logdpla
pcondcop=pcondmtcj; logdcop=logdmtcj
pcondcop=pcondgal; logdcop=logdgal
pcondcop=pcondhr; logdcop=logdhr
pcondcop=pcondipsA; logdcop=logdipsA
rh= matrix(c(0,0,0,0,0, .5,0,0,0,0, .5,.2,0,0,0, .5,.2,.3,0,0, .5,.2,.3,.4,0), 5,5)
#pcondcop=pcondbvncop; logdcop=logdbvncop; th=rh
pcondcop=pcondjoe; logdcop=logdjoe

logdgumr=function(u,v,cpar) { logdgum(1-u,1-v,cpar) }
pcondcop=pcondgumr; logdcop=logdgumr

# c-vine
cat("\nC-vine\n")
out1=rvinelogpdf(uu,C,th,logdcop=logdcop,pcond=pcondcop,iprint=T)
print(out1)
out2=cvinelogpdf(uu,th,logdcop=logdcop,pcond=pcondcop,iprint=F)
print(out2)

# d-vine
cat("\nD-vine\n")
out1=rvinelogpdf(uu,D,th,logdcop=logdcop,pcond=pcondcop,iprint=T)
print(out1)
out2=dvinelogpdf(uu,th,logdcop=logdcop,pcond=pcondcop,iprint=F)
print(out2)

