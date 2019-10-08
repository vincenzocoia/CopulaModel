# compare truncated vine/MST algorithms 
library(CopulaModel)
library(igraph)

rrvec=c(
-0.04,-0.08, 0.47,-0.07, 0.67, 0.50,-0.02, 0.52, 0.61, 0.46, 0.00,-0.06,-0.08,
-0.07,-0.10, 0.15,-0.03,-0.14,-0.15,-0.06, 0.19, 0.03,-0.01,-0.04,-0.03,-0.03,
 0.10, 0.22, 0.19, 0.10,-0.09,-0.11, 0.02,-0.02, 0.62, 0.18, 0.10, 0.04, 0.02,
-0.03, 0.05, 0.07, 0.20, 0.27, 0.16, 0.10, 0.18,-0.09,-0.10, 0.01, 0.00, 0.32,
 0.04, 0.53, 0.05,-0.16, 0.03, 0.12, 0.15, 0.06,-0.25,-0.83,-0.24,-0.52,-0.23,
-0.26,-0.11,-0.49, 0.05, 0.01,-0.14,-0.03,-0.38,-0.11,-0.57,-0.18,-0.43, 0.35,
 0.15,-0.02,-0.14,-0.13,-0.05, 0.17, 0.95, 0.21, 0.63, 0.18, 0.33,-0.81,-0.38,
-0.05, 0.07, 0.13, 0.10, 0.08,-0.33,-0.32,-0.15,-0.13,-0.21,-0.06, 0.42, 0.13,
-0.31, 0.40,-0.03,-0.11,-0.08,-0.05,-0.07, 0.13, 0.01, 0.25, 0.06, 0.22,-0.12,
-0.17, 0.15, 0.02, 0.66, 0.00,-0.11,-0.06,-0.02,-0.02, 0.18, 0.04, 0.26, 0.11,
 0.16,-0.18,-0.18, 0.19,-0.03, 0.52, 0.29,-0.01,-0.02,-0.04,-0.01,-0.02, 0.05,
 0.01, 0.09, 0.07, 0.04,-0.05,-0.05, 0.05, 0.02, 0.15, 0.30, 0.35,-0.22,-0.12,
-0.17,-0.09,-0.07, 0.10,-0.02, 0.15, 0.03, 0.14,-0.10, 0.01, 0.12, 0.02, 0.80,
 0.50, 0.19)

rmat=corvec2mat(rrvec)

nn=1187
mxtrunc=10
d=ncol(rmat)

# truncated vine
bestmst=gausstrvine.mst(rmat,ntrunc=mxtrunc,iprint=F)
print(bestmst)


perm=diag(bestmst$RVM$VineA)
Aori=bestmst$RVM$VineA
iperm=order(perm)
AA=varrayperm(Aori,iperm)
pcmat=bestmst$RVM$pc
rperm=rmat[perm,perm]
npar=0
for(ell in 1:mxtrunc)
{ Rtrun=pcor2cor.truncvine(pcmat,AA,ntrunc=ell)
  npar=npar+(d-ell)
  if(ell>1) Rtrun=Rtrun$rmat
  cat("\nMST",ell,"-truncated R-vine with npar=", npar,"\n")
  cat("max abs difference =", max(abs(rperm-Rtrun)),"\n")
  cat("ave abs difference =", sum(abs(rperm-Rtrun))/d/(d-1),"\n")
  outtr=corDis(Rtrun,rperm,nn,npar=npar)
  cat("corDis:Dfit, 2*nllk, AIC, BIC=", outtr,"\n")
}

cat("\n============================================================\n")
cat("\nnon-greedy, seed=123\n")
set.seed(123)
kbmst=gausstrvine.galg(rmat,selectcrit="prop",method="nfi",n=nn,iprint=T)
print(names(kbmst))
print(names(kbmst$RVM))
print(kbmst)

perm=diag(kbmst$RVM$VineA)
Aori=kbmst$RVM$VineA
iperm=order(perm)
AA=varrayperm(Aori,iperm)
pcmat=kbmst$RVM$pc
rperm=rmat[perm,perm]
npar=0
for(ell in 1:mxtrunc)
{ Rtrun=pcor2cor.truncvine(pcmat,AA,ntrunc=ell)
  npar=npar+(d-ell)
  if(ell>1) Rtrun=Rtrun$rmat
  cat("\nMST",ell,"-truncated R-vine with npar=", npar,"\n")
  cat("max abs difference =", max(abs(rperm-Rtrun)),"\n")
  cat("ave abs difference =", sum(abs(rperm-Rtrun))/d/(d-1),"\n")
  outtr=corDis(Rtrun,rperm,nn,npar=npar)
  cat("corDis:Dfit, 2*nllk, AIC, BIC=", outtr,"\n")
}

cat("\n============================================================\n")
cat("\nnon-greedy, seed=12345\n")
set.seed(12345)
kbmst2=gausstrvine.galg(rmat,selectcrit="prop",method="nfi",n=nn,iprint=T)
print(kbmst2)
cat("\n============================================================\n")

perm=diag(kbmst2$RVM$VineA)
Aori=kbmst2$RVM$VineA
iperm=order(perm)
AA=varrayperm(Aori,iperm)
pcmat=kbmst2$RVM$pc
rperm=rmat[perm,perm]
npar=0
for(ell in 1:mxtrunc)
{ Rtrun=pcor2cor.truncvine(pcmat,AA,ntrunc=ell)
  npar=npar+(d-ell)
  if(ell>1) Rtrun=Rtrun$rmat
  cat("\nMST",ell,"-truncated R-vine with npar=", npar,"\n")
  cat("max abs difference =", max(abs(rperm-Rtrun)),"\n")
  cat("ave abs difference =", sum(abs(rperm-Rtrun))/d/(d-1),"\n")
  outtr=corDis(Rtrun,rperm,nn,npar=npar)
  cat("corDis:Dfit, 2*nllk, AIC, BIC=", outtr,"\n")
}
