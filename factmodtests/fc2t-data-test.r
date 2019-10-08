# ML for 2-factor using bivariate t copulas and BB1

library(CopulaModel)
data(euro07gf)
n=nrow(euro07gf$uscore) # 239
logn=log(n)
d=ncol(euro07gf$uscore) # 7

#============================================================

# 2-factor copulas t with approx via monotone interpolation for pt()

gl=gausslegendre(25)

# wrapper function for different nq, nu1,nu2
# nq = number of quadrature points for Gauss-Legendre
# nu1 = degree of freedom parameter for t linking copula for factor 1
# nu2 = degree of freedom parameter for t linking copula for factor 2
# start = (2xd)-vector starting point of correlation parameters 
myrun= function(nq,nu1,nu2,start)
{ gl=gausslegendre(nq)
  cat("\nnq=", nq, " nu1=", nu1, " nu2=", nu2,"\n")
  dstr=list(data=euro07gf$uscore, copname="tapprox", quad=gl,repar=0,nu1=nu1,nu2=nu2)
  out=pdhessminb(start,f90cop2nllk, ifixed = rep(F,2*d), dstruct=dstr,
     LB=rep(-1,2*d),UB=rep(1,2*d),mxiter=50,eps=1.e-4,iprint=T)
  print(out$iposdef) 
  nllk2=2*out$fnval
  nnu=(nu1<100)+(nu2<100)
  print(c(nllk2,nllk2+2*(2*d+nnu),nllk2+logn*(2*d+nnu)))
  out
}

#st.t2=c(rep(.85,d),rep(.3,d))
cpar.frk= c(
  4.705524,10.953168,10.112137,17.585899,8.624347,11.631441,5.290361,
  5.423608,6.255184,7.104717,10.197778,4.443289,7.704374,6.586640)
tem=frk.cpar2tau(cpar.frk)
st.t2=bvn.tau2cpar(tem)
# [1] 0.6346786 0.8822358 0.8660946 0.9476766 0.8287296 0.8932819 0.6748619
# [8] 0.6835379 0.7326087 0.7737457 0.8678790 0.6160535 0.7979571 0.7496907

out=myrun(25,5,5,st.t2)
out=myrun(25,10,10,st.t2)
out=myrun(25,5,10,st.t2)
out=myrun(25,10,5,st.t2)
out=myrun(25,10,15,st.t2)
out=myrun(25,10,20,st.t2)
out=myrun(25,10,25,st.t2)
out=myrun(25,200,200,st.t2)

# also try BB1Frank with nq=31 (it failed to converge at nq=25)
cat("\nBB1/Frank 2-factor\n")
gl=gausslegendre(31)
st.bb1frk=c(.3,1.3,1.1,1.2,.9,1.3,1.1,1.4,.9,1.4,.8,1.6,.6,1.3, 7,12,12,17,6,12,6)
dstrbb1frk=list(data=euro07gf$uscore, copname="bb1frank", quad=gl,repar=0)
#f90cop2nllk(st.bb1frk, dstrbb1frk, iprfn=F)
out.bb1frk=pdhessminb(st.bb1frk,f90cop2nllk, ifixed = rep(F,3*d), 
  dstruct=dstrbb1frk, LB=c(rep(c(0,1),d),rep(-20,d)),
  UB=c(rep(10,2*d),rep(30,d)),mxiter=50,eps=1.e-4,iprint=T)
cat("iposdef=", out.bb1frk$iposdef,"\n") 
nllk2=2*out.bb1frk$fnval
#print(c(nllk2,nllk2+2*(3*d),nllk2+logn*(3*d)))
cat("nllk*2, aic, bic: ", c(nllk2,nllk2+2*(d*3),nllk2+logn*(d*3)),"\n")

