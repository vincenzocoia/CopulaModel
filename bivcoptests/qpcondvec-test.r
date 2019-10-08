# check that qcond and pcond for bvn,t,bb1,frk work as vectorized

library(CopulaModel)

set.seed(123)
pp=runif(10)
print(pp)

cat("check for qcondfrk\n")
vv=seq(.1,.9,.1)
vv=c(vv,.45)
uu=qcondfrk(pp,vv,3.2)
print(uu)
p2=pcondfrk(uu,vv,3.2)
print(max(abs(p2-pp)))

cat("\ncheck for qcondbb1\n")
uu=qcondbb1(pp,vv,c(.5,2.2))
print(uu)
p2=pcondbb1(uu,vv,c(.5,2.2))
print(max(abs(p2-pp)))

cat("\ncheck for qcondbvncop\n")
uu=qcondbvncop(pp,vv,c(.52))
print(uu)
p2=pcondbvncop(uu,vv,c(.52))
print(max(abs(p2-pp)))

cat("\ncheck for qcondbvtcop\n")
dfdefault=5
uu=qcondbvtcop(pp,vv,c(.52))
print(uu)
p2=pcondbvtcop(uu,vv,c(.52))
print(max(abs(p2-pp)))

uu=qcondbvtcop(pp,vv,c(.52,5))
print(uu)
p2=pcondbvtcop(uu,vv,c(.52,5))
print(max(abs(p2-pp)))

cat("\ncheck for qcondgum\n")
uu=qcondgum(pp,vv,2.2)
print(uu)
p2=pcondgum(uu,vv,2.2)
print(max(abs(p2-pp)))


