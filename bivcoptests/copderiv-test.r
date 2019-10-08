
# derivatives of log copula pdf and conditional cdf C_{2|1}
#  wrt u,v and copula parameters for Frank, Gumbel, t (fixed df), BB1.
# check numerically  for correctness

library(CopulaModel)
#source("../R/copderiv.R")

u=.4; v=.7; cpar=2
eps=1.e-4

cat("numerical derivatives and then analytic\n")

cat("\nFrank logdfrk.deriv\n")
lpdf=logdfrk(u,v,cpar)
lpdfu=logdfrk(u+eps,v,cpar)
lpdfv=logdfrk(u,v+eps,cpar)
lpdfcpar=logdfrk(u,v,cpar+eps)
cat(lpdf,(lpdfu-lpdf)/eps,(lpdfv-lpdf)/eps,(lpdfcpar-lpdf)/eps,"\n")
tem=logdfrk.deriv(u,v,cpar) # OK
print(tem)

cat("\nFrank pcondfrk.deriv\n")
ccdf=pcondfrk(v,u,cpar)
ccdfx=pcondfrk(v,u+eps,cpar)
ccdfy=pcondfrk(v+eps,u,cpar)
ccdfcpar=pcondfrk(v,u,cpar+eps)
cat(ccdf,(ccdfx-ccdf)/eps,(ccdfy-ccdf)/eps,(ccdfcpar-ccdf)/eps,"\n")
temc=pcondfrk.deriv(v,u,cpar) # OK
print(temc)

#========================

cat("\nGumbel logdgum.deriv\n")
lpdf=logdgum(u,v,cpar)
lpdfu=logdgum(u+eps,v,cpar)
lpdfv=logdgum(u,v+eps,cpar)
lpdfcpar=logdgum(u,v,cpar+eps)
cat(lpdf,(lpdfu-lpdf)/eps,(lpdfv-lpdf)/eps,(lpdfcpar-lpdf)/eps,"\n")
tem=logdgum.deriv(u,v,cpar) # OK
print(tem)

cat("\nGumbel pcondgum.deriv\n")
ccdf=pcondgum(v,u,cpar)
ccdfx=pcondgum(v,u+eps,cpar)
ccdfy=pcondgum(v+eps,u,cpar)
ccdfcpar=pcondgum(v,u,cpar+eps)
cat(ccdf,(ccdfx-ccdf)/eps,(ccdfy-ccdf)/eps,(ccdfcpar-ccdf)/eps,"\n")
temc=pcondgum.deriv(v,u,cpar) # OK
print(temc)

#========================

# t with fixed nu=df
dfdefault=5
cpar=.6
cat("\nt logdbvtcop.deriv\n")
lpdf=logdbvtcop(u,v,cpar)
lpdfu=logdbvtcop(u+eps,v,cpar)
lpdfv=logdbvtcop(u,v+eps,cpar)
lpdfcpar=logdbvtcop(u,v,cpar+eps)
cat(lpdf,(lpdfu-lpdf)/eps,(lpdfv-lpdf)/eps,(lpdfcpar-lpdf)/eps,"\n")
tem=logdbvtcop.deriv(u,v,cpar) # OK
print(tem)

cat("\nt pcondbvtcop.deriv\n")
ccdf=pcondbvtcop(v,u,cpar)
ccdfx=pcondbvtcop(v,u+eps,cpar)
ccdfy=pcondbvtcop(v+eps,u,cpar)
ccdfcpar=pcondbvtcop(v,u,cpar+eps)
cat(ccdf,(ccdfx-ccdf)/eps,(ccdfy-ccdf)/eps,(ccdfcpar-ccdf)/eps,"\n")
temc=pcondbvtcop.deriv(v,u,cpar) # 
print(temc)

#========================

# BB1
cpar=c(.5,1.3)
cat("\nBB1 logdbb1.deriv\n")
lpdf=logdbb1(u,v,cpar)
lpdfu=logdbb1(u+eps,v,cpar)
lpdfv=logdbb1(u,v+eps,cpar)
lpdfth=logdbb1(u,v,cpar+c(eps,0))
lpdfdl=logdbb1(u,v,cpar+c(0,eps))
cat(lpdf,(lpdfu-lpdf)/eps,(lpdfv-lpdf)/eps,(lpdfth-lpdf)/eps,
    (lpdfdl-lpdf)/eps, "\n")
tem=logdbb1.deriv(u,v,cpar) # OK
print(tem)


cat("\nBB1 pcondbb1.deriv\n")
ccdf=pcondbb1(v,u,cpar)
ccdfx=pcondbb1(v,u+eps,cpar)
ccdfy=pcondbb1(v+eps,u,cpar)
ccdfth=pcondbb1(v,u,cpar+c(eps,0))
ccdfdl=pcondbb1(v,u,cpar+c(0,eps))
cat(ccdf,(ccdfx-ccdf)/eps,(ccdfy-ccdf)/eps,(ccdfth-ccdf)/eps,
   (ccdfdl-ccdf)/eps,"\n")
temc=pcondbb1.deriv(v,u,cpar) # OK
print(temc)

