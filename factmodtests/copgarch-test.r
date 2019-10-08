# test of rgarch1fact, rgarch2fact, rgarchbifact, rgarchnestfact

library(CopulaModel)

garchpar=matrix(c(0.09475955, 0.12552862, 0.06885825, 0.21206611, 0.09284262,
 -0.06241426,-0.04088966,-0.08266326,-0.02045153,-0.05177665,
  0.01906628, 0.01485872, 0.01698609, 0.08236219, 0.01437761,
  0.07118317, 0.07564840, 0.08851696, 0.09204580, 0.08004666,
  0.90983966, 0.91433214, 0.88253580, 0.84222487, 0.90219596,
 10.41578901, 8.30321969,15.54704352, 7.91622415,11.69870286), 6,5,byrow=T)
sigma0=c(1.113026,1.286220,0.8577638,1.176423,0.9684971)
garchpar6=cbind(garchpar,c(0.1, -0.05,0.01,0.08,0.9,12.))
sigma6=c(sigma0,1)
grsize=c(2,2,2)

cpar=seq(1.1,1.5,.1)
cparbvn=.9*(1:5)/5
cparbvt=c(cparbvn,4)
cparbb1=c(.2,1.1,.2,1.2,.2,1.3,.2,1.4,.2,1.5)

cpar2=c(seq(1.1,1.5,.1),rep(1.1,5))
cpar2bb1=c(cparbb1,seq(1.1,1.5,.1))
cpar2bvn=c(.8,.7,.6,.5,.5,.4,.4,.4,.4,.4)
cpar2bvt=c(cpar2bvn,4)

cparbi=c(seq(1.1,1.6,.1),rep(1.1,6))
cparbibb1=c(cparbb1,.2,1.6,seq(1.1,1.6,.1))
cparbibvn=c(.8,.7,.6,.5,.5,.6,.4,.4,.4,.4,.4,.3)
cparbibvt=c(cparbibvn,4)

cparne=c(rep(1.1,3),seq(1.1,1.6,.1))
cparnebb1=c(seq(0.6,1.1,.1),cparbb1)
cparnebvn=c(.7,.7,.7,.8,.7,.6,.5,.5,.4)
cparnebvt=c(cparnebvn,4)

#============================================================

cat("\ntesting r1fact\n")
set.seed(123)
tem=r1fact(3,5,cpar,3); print(tem)
tem=r1fact(3,5,cpar,3); print(tem)
set.seed(123)
tem=r1fact(3,5,cpar,3); print(tem)
tem=r1fact(3,5,cpar,-3); print(tem)
tem=r1fact(3,5,cpar,5); print(tem)
tem=r1fact(3,5,cparbb1,9); print(tem)
tem=r1fact(3,5,cparbvn,1); print(tem)
tem=r1fact(3,5,cparbvt,2); print(tem)

cat("\ntesting r2fact\n")
set.seed(123)
tem=r2fact(3,5,cpar2,3); print(tem)
tem=r2fact(3,5,cpar2,3); print(tem)
set.seed(123)
tem=r2fact(3,5,cpar2,3); print(tem)
tem=r2fact(3,5,cpar2,-3); print(tem)
tem=r2fact(3,5,cpar2,5); print(tem)
tem=r2fact(3,5,cpar2bb1,9); print(tem)
tem=r2fact(3,5,cpar2bvn,1); print(tem)
tem=r2fact(3,5,cpar2bvt,2); print(tem)

cat("\ntesting rbifact\n")
set.seed(123)
tem=rbifact(3,grsize,cparbi,3); print(tem)
tem=rbifact(3,grsize,cparbi,3); print(tem)
set.seed(123)
tem=rbifact(3,grsize,cparbi,3); print(tem)
tem=rbifact(3,grsize,cparbi,-3); print(tem)
tem=rbifact(3,grsize,cparbi,5); print(tem)
tem=rbifact(3,grsize,cparbibb1,9); print(tem)
tem=rbifact(3,grsize,cparbibvn,1); print(tem)
tem=rbifact(3,grsize,cparbibvt,2); print(tem)

cat("\ntesting rnestfact\n")
set.seed(123)
tem=rnestfact(3,grsize,cparne,3); print(tem)
tem=rnestfact(3,grsize,cparne,3); print(tem)
set.seed(123)
tem=rnestfact(3,grsize,cparne,3); print(tem)
tem=rnestfact(3,grsize,cparne,-3); print(tem)
tem=rnestfact(3,grsize,cparne,5); print(tem)
tem=rnestfact(3,grsize,cparnebb1,9); print(tem)
tem=rnestfact(3,grsize,cparnebvn,1); print(tem)
tem=rnestfact(3,grsize,cparnebvt,2); print(tem)

#============================================================

# more tests of r1fact
cat("\nmore tests of r1fact with mvn, mvt\n")
rhvec=c(.8,.7,.6,.5,.5)
out=r1fact(3,5,rhvec,1)
out=r1fact(3,5,c(rhvec,4),2)
set.seed(123)
n=3000
dat=r1fact(n,5,rhvec,1)
print(cov(dat))
print(cor(dat))
print(outer(rhvec,rhvec))
dat=r1fact(n,5,c(rhvec,4),2)
print(cov(dat))
print(cor(dat))

#============================================================

cat("\ntesting rgarch1fact\n")
set.seed(123)
out=rgarch1fact(3,garchpar,cpar,sigma0,copcode=3)
print(out)
out=rgarch1fact(3,garchpar,cpar,sigma0,copcode=3)
print(out)
set.seed(123)
out=rgarch1fact(3,garchpar,cpar,sigma0,copcode=3)
set.seed(123)
out=rgarch1fact(3,garchpar,cpar,sigma0,copcode=-3)
print(out)
set.seed(123)
out=rgarch1fact(3,garchpar,cpar,sigma0,copcode=5)
print(out)
set.seed(123)
out=rgarch1fact(3,garchpar,cparbvn,sigma0,copcode=1)
print(out)
set.seed(123)
out=rgarch1fact(3,garchpar,cparbvt,sigma0,copcode=2)
print(out)
set.seed(123)
out=rgarch1fact(3,garchpar,cparbb1,sigma0,copcode=9)
print(out)

cat("\ntesting rgarch2fact\n")
set.seed(123)
out=rgarch2fact(3,garchpar,cpar2,sigma0,copcode=3)
print(out)

cat("\ntesting rgarchbifact\n")
set.seed(123)
out=rgarchbifact(3,grsize,garchpar6,cparbi,sigma6,copcode=3)
print(out)

cat("\ntesting rgarchnestfact\n")
set.seed(123)
out=rgarchnestfact(3,grsize,garchpar6,cparne,sigma6,copcode=3)
print(out)

