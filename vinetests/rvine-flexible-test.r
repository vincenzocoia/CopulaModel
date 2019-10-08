# example with different bivariate copula family for each edge
# sequential estimate by tree with check of pseudo-observations
#  for tail dependence etc

library(CopulaModel)

parmat1=matrix(c(0,1.5,2,2.5,2.2,0,0,1.5,1.5,3,0,0,0,1.2,1.2,0,0,0,1.2,1.2,0,0,0,0,1.2),5,5,byrow=T)
parvec1=c(parmat1[1,2:5],parmat1[2,3:5],parmat1[3,4:5],parmat1[4,5])
parmat2=matrix(c(0,1.5,2,2.5,2.2,0,0,.5,.4,.6,0,0,0,.2,.2,0,0,0,.2,.2,0,0,0,0,.2),5,5,byrow=T)
parvec2=c(parmat2[1,2:5],parmat2[2,3:5],parmat2[3,4:5],parmat2[4,5])

C5=Cvinearray(5)

pcondnames1=rep("pcondgum",4)
pcondmat=matrix(c("",rep("pcondgum",4),"","",rep("pcondbvncop",3),
  "","","",rep("pcondbvncop",2),"","","","","pcondbvncop",rep("",5)),
  5,5,byrow=T)
qcondmat=matrix(c("",rep("qcondgum",4),"","",rep("qcondbvncop",3),
  "","","",rep("qcondbvncop",2),"","","","","qcondbvncop",rep("",5)),
  5,5,byrow=T)
logdcopmat=matrix(c("",rep("logdgum",4),"","",rep("logdbvncop",3),
  "","","",rep("logdbvncop",2),"","","","","logdbvncop",rep("",5)),5,5,byrow=T)
qcondnames1=rep("qcondgum",4)
pcondnames2=c("pcondgum","pcondbvncop","pcondbvncop","pcondbvncop")
qcondnames2=c("qcondgum","qcondbvncop","qcondbvncop","qcondbvncop")

np=matrix(1,5,5)

nsim=300
set.seed(123)
udat1c=rvinesimvec(nsim,C5,parvec1,np,qcondnames1,pcondnames1)
set.seed(123)
udat2c=rvinesimvec(nsim,C5,parvec2,np,qcondnames2,pcondnames2)
set.seed(123)
udat2cflex=rvinesimvec2(nsim,C5,ntrunc=4,parvec2,np,qcondmat,pcondmat)
# OK same as udat2c
print(max(abs(udat2c-udat2cflex)))
cat("\n============================================================\n")

# put a BVN edge in tree 1 to mix with Gumbel
#  and a Gumbel edge in tree 2 to mix with BVN
pcondmat3=pcondmat
qcondmat3=qcondmat
logdcopmat3=logdcopmat
pcondmat3[1,3]="pcondbvncop"
qcondmat3[1,3]="qcondbvncop"
logdcopmat3[1,3]="logdbvncop"
pcondmat3[2,3]="pcondgum"
qcondmat3[2,3]="qcondgum"
logdcopmat3[2,3]="logdgum"
parvec3=parvec2
parvec3[2]=0.5;  parvec3[5]=2.0
set.seed(123)
udat3c=rvinesimvec2(nsim,C5,ntrunc=4,parvec3,np,qcondmat3,pcondmat3)
zdat3c=nscore(udat3c) 
# pairs(zdat3c) # (1,3) margin more elliptical in shape than others

#============================================================
# 1-truncation
cat("\n1-truncation\n")

cat("Gumbel/BVN C-vine, mixed at levels 1 and 2\n")

cat("fit all Gumbel at tree 1\n")
mle3cmis=nlm(rvinenllk.trunc2,p=parvec1[1:4],
    np=np,ifixed=rep(F,4),parfixed=NULL,
    udat=udat3c,A=C5,ntrunc=1,logdcopmat=logdcopmat,pcondmat=pcondmat,
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
# nllk=-322.1337
cat("fit Gumbel/BVN/Gumbel/Gumbel at tree 1\n")
mle3c=nlm(rvinenllk.trunc2,p=parvec3[1:4],
    np=np,ifixed=rep(F,4),parfixed=NULL,
    udat=udat3c,A=C5,ntrunc=1,logdcopmat=logdcopmat3,pcondmat=pcondmat3,
    hessian=T,iterlim=30,print.level=1,LB=c(1,-1,1,1),UB=c(20,1,20,20))
# nllk=-328.9861
pseud3c=rvinenllkpseud2(mle3c$estimate,udat3c,C5,ntrunc=1,logdcopmat3,
  pcondmat3,np)

# should be 2|1 3|1 4|1 5|1 in $condforw
zdat3c=nscore(pseud3c$condforw)
# pairs(zdat3c) # 2,3|1 shows  upper tail dependence
out=semicortable(zdat3c,inscore=T)
out[,1]=out[,1]+1
out[,2]=out[,2]+1
print(out)
#  j1 j2     ncorr     lcorr     ucorr  bvnsemic
#1  2  3 0.6717155 0.2387019 0.7090820 0.4301044
#2  2  4 0.4691834 0.1131556 0.2155522 0.2450662
#3  3  4 0.4803906 0.1929006 0.2913079 0.2535240
#4  2  5 0.6357750 0.3629089 0.4981312 0.3915010
#5  3  5 0.5250522 0.1946253 0.5013440 0.2890220
#6  4  5 0.4487680 0.0772523 0.2735893 0.2300923
cat("\n============================================================\n")

cat("compare different 2-truncated vines")

cat("\npair-copulas for first model\n")
pcondtmp=pcondmat3
logdcoptmp=logdcopmat3
print(logdcoptmp)
try1=nlm(rvinenllk.trunc2,p=c(1.5,.5,2.5,2.2,2,.4,.6),
    np=np,ifixed=rep(F,7),parfixed=NULL,
    udat=udat3c,A=C5,ntrunc=2,logdcopmat=logdcoptmp,pcondmat=pcondtmp,
    hessian=T,iterlim=30,print.level=1,
    LB=c(1,-1,1,1,1,-1,-1),UB=c(20,1,20,20,20,1,1))
# nllk=-559.6373

cat("\npair-copulas for second model\n")
pcondtmp[2,5]="pcondgum"
logdcoptmp[2,5]="logdgum"
print(logdcoptmp)
try2=nlm(rvinenllk.trunc2,p=c(1.5,.5,2.5,2.2,2,.4,1.6),
    np=np,ifixed=rep(F,7),parfixed=NULL,
    udat=udat3c,A=C5,ntrunc=2,logdcopmat=logdcoptmp,pcondmat=pcondtmp,
    hessian=T,iterlim=30,print.level=1,
    LB=c(1,-1,1,1,1,-1,1),UB=c(20,1,20,20,20,1,20))
# nllk=-553.6895
# first one is better for nllk

cat("\npseudo-observations for first model\n")
pcondtmp[2,5]="pcondbvncop"
logdcoptmp[2,5]="logdbvncop"
pseud3c2=rvinenllkpseud2(try1$estimate,udat3c,C5,ntrunc=2,logdcoptmp,
  pcondtmp,np)

# should be 3|12 4|12 5|12 in $condforw
zdat3c2=nscore(pseud3c2$condforw)
# pairs(zdat3c2) 
out=semicortable(zdat3c2,inscore=T)
out[,1]=out[,1]+2
out[,2]=out[,2]+2
print(out)
#  j1 j2     ncorr     lcorr     ucorr  bvnsemic
#1  3  4 0.2826826  0.07308482 0.081108575 0.12596357
#2  3  5 0.1589957  0.15278148 0.051072115 0.06449663
#3  4  5 0.2138694 -0.09535638 0.006865775 0.09035861

