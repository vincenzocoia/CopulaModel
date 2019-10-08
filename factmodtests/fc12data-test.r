# tests for 1-factor and 2-factor MLE

library(CopulaModel)
data(euro07gf)
n=nrow(euro07gf$uscore) # 239
logn=log(n)
d=ncol(euro07gf$uscore) # 7

# 1-factor

gl=gausslegendre(25)

cat("\nFrank 1-factor\n")
cpar.frk=depmeas2cpar(.8,"rhoS","frank")
st.frk=rep(cpar.frk,d)
dstrfrk = list(data=euro07gf$uscore, copname="frank", quad=gl,repar=0);
out.frk=pdhessminb(st.frk,f90cop1nllk, ifixed = rep(F,d), dstruct=dstrfrk,
   LB=rep(-20,d),UB=rep(35,d),mxiter=30,eps=1.e-4,iprint=T)
cat("iposdef=", out.frk$iposdef,"\n") 
nllk2=2*out.frk$fnval
cat("nllk*2, aic, bic: ", c(nllk2,nllk2+2*(d),nllk2+logn*d),"\n")
cat("\n============================================================\n")

cat("\nt 1-factor\n")
st.t=rep(.85,d)
for(nu in c(5,10,15,20,200))
{ cat("nu=", nu,"\n")
  dstrt=list(data=euro07gf$uscore, copname="t", quad=gl, repar=0, nu=nu)
  out.t=pdhessminb(st.t,f90cop1nllk, ifixed = rep(F,d), dstruct=dstrt,
     LB=rep(-1,d),UB=rep(1,d),mxiter=30,eps=1.e-4,iprint=T)
  cat("iposdef=", out.t$iposdef,"\n") 
  nllk2=2*out.t$fnval
  cat("nllk*2, aic, bic: ", c(nllk2,nllk2+2*(d+1),nllk2+logn*(d+1)),"\n")
  cat("nu=", nu,"\n")
  cat("\n============================================================\n")
}

cat("\nBB1 1-factor\n")
cpar.bb1=bb1.lm2cpar(c(.6,.6))
st.bb1=rep(cpar.bb1,d)
dstrbb1=list(data=euro07gf$uscore, copname="bb1", quad=gl, repar=0)
#f90cop1nllk(st.bb1, dstrbb1, iprfn=F)
out.bb1=pdhessminb(st.bb1,f90cop1nllk, ifixed = rep(F,2*d), dstruct=dstrbb1,
   LB=rep(c(0,1),d),UB=rep(10,2*d),mxiter=30,eps=1.e-4,iprint=T)
cat("iposdef=", out.bb1$iposdef,"\n") 
nllk2=2*out.bb1$fnval
#print(c(nllk2,nllk2+2*(d*2),nllk2+logn*(d*2)))
cat("nllk*2, aic, bic: ", c(nllk2,nllk2+2*(d*2),nllk2+logn*(d*2)),"\n")
cat("\n============================================================\n")

#============================================================

# 2-factor copulas

cat("\nFrank 2-factor\n")
st.frk2=c(rep(15,d),rep(3,d))
dstrfrk = list(data=euro07gf$uscore, copname="frank", quad=gl,repar=0);
#f90cop1nllk(st.frk2, dstrfrk, iprfn=T)
out.frk2=pdhessminb(st.frk2,f90cop2nllk, ifixed = rep(F,2*d), dstruct=dstrfrk,
   LB=rep(-20,2*d),UB=rep(35,2*d),mxiter=30,eps=1.e-4,iprint=T)
cat("iposdef=", out.frk2$iposdef,"\n") 
nllk2=2*out.frk2$fnval
#print(c(nllk2,nllk2+2*(2*d),nllk2+logn*(2*d)))
cat("nllk*2, aic, bic: ", c(nllk2,nllk2+2*(d*2),nllk2+logn*(d*2)),"\n")
cat("\n============================================================\n")

cat("\nBB1/Frank 2-factor\n")
cpar.bb1=bb1.lm2cpar(c(.6,.6))
st.bb1frk=c(rep(cpar.bb1,d),rep(3,d))
st.bb1frk=c(.2,1.3,.6,1.7,.5,1.7,.8,2.1,.8,1.7,.8,1.8,.4,1.4, 7,12,11,17,6,10,8)
dstrbb1frk=list(data=euro07gf$uscore, copname="bb1frank", quad=gl,repar=0)
#f90cop2nllk(st.bb1frk, dstrbb1frk, iprfn=F)
out.bb1frk=pdhessminb(st.bb1frk,f90cop2nllk, ifixed = rep(F,3*d), 
  dstruct=dstrbb1frk, LB=c(rep(c(0,1),d),rep(-20,d)),
  UB=c(rep(10,2*d),rep(30,d)),mxiter=50,eps=1.e-4,iprint=T)
cat("iposdef=", out.bb1frk$iposdef,"\n") 
nllk2=2*out.bb1frk$fnval
#print(c(nllk2,nllk2+2*(3*d),nllk2+logn*(3*d)))
cat("nllk*2, aic, bic: ", c(nllk2,nllk2+2*(d*3),nllk2+logn*(d*3)),"\n")


cat("\n============================================================\n")

#cat("\nt 2-factor\n")
# use tapprox in separate file
