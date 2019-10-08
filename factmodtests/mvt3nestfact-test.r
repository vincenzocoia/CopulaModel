# mvt with 3-nested-factor structure, simulation example
# full=F

# if full = F, reduced model is used
# if full = F, in param: values for (Vg and V0) go first (size ng)
#  then V_{kg} and Vg: 1st group, then second group etc (size nsbg)
#  then Z_{ikg} and V_{kg}: 1st subgroup, then second subgroup etc (size d) 

library(CopulaModel)

grsize=c(6,6)
sbgrsize=c(3,3,3,3)
d=sum(grsize)
nestpar=c(.6,.7, .5,.6,.5,.6, .8,.8,.8,.9,.9,.9,.7,.7,.7,.7,.5,.5)
mgrp=length(grsize);
msbgrp=length(sbgrsize);

# code from mvttrifct in  file bifct.R
# convert to parameters for tri-factor structure
ind=0; indsb=0;
a3=nestpar[(mgrp+msbgrp+1):(mgrp+msbgrp+d)]; a2=rep(0,d); a1=rep(0,d);
for (jsbg in 1:msbgrp)
{ indsb1=indsb+1;  indsb2=indsb+sbgrsize[jsbg];
  a2[indsb1:indsb2]=nestpar[jsbg+mgrp]; 
  indsb=indsb+sbgrsize[jsbg]; 
}
for (jg in 1:mgrp)
{ ind1=ind+1; ind2=ind+grsize[jg];
  a1[ind1:ind2]=nestpar[jg];
  ind=ind+grsize[jg]; 
}
rh11= a1*a2*a3;
rh21= a2*a3*sqrt(1-a1^2)/sqrt(1-rh11^2);
rh31= a3*sqrt(1-a2^2)/sqrt(1-rh11^2)/sqrt(1-rh21^2);
rh1=rh11; rh2=rh21; rh3=rh31;
cat("parameters from 3-nested factor converted to parameters for tri-factor\n")
print(cbind(rh1,rh2,rh3))

#n=100
n=600
df=10
# correlation matrix and inverse, determinant
robj=trifct(grsize,sbgrsize,rh1,rh2,rh3)
achol=chol(robj$fctmat)

# simulate data
set.seed(123)
z=matrix(rnorm(n*d),n,d)
z=z%*%achol
udata=uscore(z)
tdata=qt(udata,df)
tripar=c(rh1,rh2,rh3)
cat("\nnllk at ", tripar,"\n")
out=mvttrifactnllk(tripar,grsize,sbgrsize,tdata,df)
print(out)
# check with ML
ml=mvttrifact(tdata,start=nestpar,grsize,sbgrsize,df,prlevel=1,full=F,mxiter=150)
#MLE: 
# [1] 0.4514109 0.5754071 0.4703391 0.5776824 0.4847008 0.4575413 0.7516366
# [8] 0.8135473 0.7731043 0.8756071 0.8921518 0.8879098 0.6719751 0.6110680
#[15] 0.6846689 0.8155632 0.3886726 0.4400177
#nllk:  
#[1] -1213.062

cat("\nparameters for simulation\n")
print(nestpar)
# [1] 0.6 0.7 0.5 0.6 0.5 0.6 0.8 0.8 0.8 0.9 0.9 0.9 0.7 0.7 0.7 0.7 0.5 0.5


