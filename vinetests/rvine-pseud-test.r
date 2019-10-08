# sequential estimate by tree with check of pseudo-observations
#  for tail dependence etc

library(CopulaModel)

parmat1=matrix(c(0,1.5,2,2.5,2.2,0,0,1.5,1.5,3,0,0,0,1.2,1.2,0,0,0,1.2,1.2,0,0,0,0,1.2),5,5,byrow=T)
parvec1=c(parmat1[1,2:5],parmat1[2,3:5],parmat1[3,4:5],parmat1[4,5])
parmat2=matrix(c(0,1.5,2,2.5,2.2,0,0,.5,.4,.6,0,0,0,.2,.2,0,0,0,.2,.2,0,0,0,0,.2),5,5,byrow=T)
parvec2=c(parmat2[1,2:5],parmat2[2,3:5],parmat2[3,4:5],parmat2[4,5])

C5=Cvinearray(5)
B1=vnum2array(5,5)
# B1b=vnum2array(5,3)

pcondnames1=rep("pcondgum",4)
qcondnames1=rep("qcondgum",4)
pcondnames2=c("pcondgum","pcondbvncop","pcondbvncop","pcondbvncop")
qcondnames2=c("qcondgum","qcondbvncop","qcondbvncop","qcondbvncop")
np=matrix(1,5,5)

nsim=300
set.seed(123)
udat1c=rvinesimvec(nsim,C5,parvec1,np,qcondnames1,pcondnames1)
set.seed(123)
udat1b=rvinesimvec(nsim,B1,parvec1,np,qcondnames1,pcondnames1)
set.seed(123)
udat2c=rvinesimvec(nsim,C5,parvec2,np,qcondnames2,pcondnames2)
set.seed(123)
udat2b=rvinesimvec(nsim,B1,parvec2,np,qcondnames2,pcondnames2)

#============================================================
# 1-truncation
cat("\n1-truncation\n")

cat("Gumbel C-vine\n")
logdcopnames1="logdgum"
mle1c=nlm(rvinenllk1.trunc,p=parvec1[1:4],
    udat=udat1c,A=C5,logdcopnames=logdcopnames1,pcondnames=pcondnames1[1],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud1c=rvinenllkpseud(mle1c$estimate,udat1c,C5,logdcopnames1,
  pcondnames1[1],np)
# should be 2|1 3|1 4|1 5|1 in $condforw
zdat1c=nscore(pseud1c$condforw)
# pairs(zdat1c)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
for(i2 in 2:4)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat1c[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.4805246 0.1208167 0.6378915 0.2536262 
#1 3 0.5731644 0.1196589 0.5118941 0.330791 
#2 3 0.5328685 0.1228473 0.533721 0.2955467 
#1 4 0.8552503 0.5300056 0.8644127 0.6859402 
#2 4 0.5187759 0.1797216 0.6059007 0.2838526 
#3 4 0.6138828 0.1093684 0.6042844 0.3693834 
cat("\n------------------------------------------------------------\n")

cat("Gumbel/BVN C-vine\n")
logdcopnames1="logdgum"
mle2c=nlm(rvinenllk1.trunc,p=parvec2[1:4],
    udat=udat2c,A=C5,logdcopnames=logdcopnames1,pcondnames=pcondnames2[1],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud2c=rvinenllkpseud(mle2c$estimate,udat2c,C5,logdcopnames1,
  pcondnames2[1],np)
# should be 2|1 3|1 4|1 5|1 in $condforw
zdat2c=nscore(pseud2c$condforw)
# pairs(zdat2c)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
for(i2 in 2:4)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat2c[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.480421 0.1527189 0.3382402 0.2535473 
#1 3 0.4691834 0.1131556 0.2155522 0.2450662 
#2 3 0.4450296 0.1968744 0.2181895 0.2274091 
#1 4 0.635775 0.3629089 0.4981312 0.3915010 
#2 4 0.4185118 0.1635124 0.3340247 0.2088761 
#3 4 0.448768 0.0772523 0.2735893 0.2300923
cat("\n------------------------------------------------------------\n")

cat("Gumbel B1-vine\n")
logdcopnames1="logdgum"
mle1b=nlm(rvinenllk1.trunc,p=parvec1[1:4],
    udat=udat1b,A=B1,logdcopnames=logdcopnames1,pcondnames=pcondnames1[1],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud1b=rvinenllkpseud(mle1b$estimate,udat1b,B1,logdcopnames1,
  pcondnames1[1],np)
# should be 2|1 3|1 4|1 5|2 in $condforw
# 1|2 in $condbackw
zdat1b=nscore(pseud1b$condforw[,1:3])
zzdat1b=nscore(cbind(pseud1b$condforw[,4],pseud1b$condbackw[,1]))
# pairs(zdat1b)
cat("i1 i2 ncor lcor ucor bvnsemic: forward\n")
for(i2 in 2:3)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat1b[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.4805246 0.1208167 0.6378915 0.2536262 
#1 3 0.5731644 0.1196589 0.5118941 0.330791 
#2 3 0.5328685 0.1228473 0.533721 0.2955467
# pairs(zzdat1b)
cat("i1 i2 ncor lcor ucor bvnsemic: backward\n")
semic=semicor(zzdat1b)
ze=bvnsemic(semic[1])
cat(i1,i2,semic,ze,"\n")
cat("\n------------------------------------------------------------\n")

cat("Gumbel/BVN B1-vine\n")
logdcopnames1="logdgum"
mle2b=nlm(rvinenllk1.trunc,p=parvec2[1:4],
    udat=udat2b,A=B1,logdcopnames=logdcopnames1,pcondnames=pcondnames2[1],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud2b=rvinenllkpseud(mle2b$estimate,udat2b,B1,logdcopnames1,
  pcondnames2[1],np)
# should be 2|1 3|1 4|1 5|2 in $condforw
# 1|2 in $condbackw
zdat2b=nscore(pseud2b$condforw[,1:3])
zzdat2b=nscore(cbind(pseud2b$condforw[,4],pseud2b$condbackw[,1]))
# pairs(zdat2b)
cat("i1 i2 ncor lcor ucor bvnsemic: forward\n")
for(i2 in 2:3)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat2b[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.480421 0.1527189 0.3382402 0.2535473 
#1 3 0.4691834 0.1131556 0.2155522 0.2450662 
#2 3 0.4450296 0.1968744 0.2181895 0.2274091
# pairs(zzdat2b)
cat("i1 i2 ncor lcor ucor bvnsemic: backward\n")
semic=semicor(zzdat2b)
ze=bvnsemic(semic[1])
cat(i1,i2,semic,ze,"\n")
# 2 3 0.5761853 0.3409091 0.3453862 0.333547

cat("\n============================================================\n")

#============================================================

# 2-truncation

cat("\n2-truncation\n")

cat("Gumbel C-vine\n")
logdcopnames2=c("logdgum","logdgum")
mle1c=nlm(rvinenllk1.trunc,p=parvec1[1:7],
    udat=udat1c,A=C5,logdcopnames=logdcopnames2,pcondnames=pcondnames1[1:2],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud1c=rvinenllkpseud(mle1c$estimate,udat1c,C5,logdcopnames2,
  pcondnames1[1:2],np)
# should be 3|12 4|12 5|12 in $condforw
zdat1c=nscore(pseud1c$condforw)
# pairs(zdat1c)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
for(i2 in 2:3)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat1c[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.3509610 0.1069968 0.349834 0.1653167 
#1 3 0.2160991 0.1159900 0.3391459 0.09145468 
#2 3 0.2995734 -0.09751323 0.4049379 0.1353011 
cat("\n------------------------------------------------------------\n")

cat("Gumbel/BVN C-vine\n")
logdcopnames2=c("logdgum","logdbvncop")
mle2c=nlm(rvinenllk1.trunc,p=parvec2[1:7],
    udat=udat2c,A=C5,logdcopnames=logdcopnames2,pcondnames=pcondnames2[1:2],
    hessian=T,iterlim=30,print.level=1,LB=c(rep(1,4),rep(-1,3)),
    UB=c(rep(20,4),rep(1,3)))
pseud2c=rvinenllkpseud(mle2c$estimate,udat2c,C5,logdcopnames2,
  pcondnames2[1:2],np)
# should be 3|12 4|12 5|12 in $condforw
zdat2c=nscore(pseud2c$condforw)
# pairs(zdat2c)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
for(i2 in 2:3)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat2c[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.2834874 0.07062834 0.09127613 0.1264028 
#1 3 0.1611091 0.1530898 0.04949042 0.0654546 
#2 3 0.2143484 -0.0944196 0.006465896 0.09059374 
cat("\n------------------------------------------------------------\n")

cat("Gumbel B1-vine\n")
logdcopnames2=c("logdgum","logdgum")
mle1b=nlm(rvinenllk1.trunc,p=parvec1[1:7],
    udat=udat1c,A=B1,logdcopnames=logdcopnames2,pcondnames=pcondnames1[1:2],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud1b=rvinenllkpseud(mle1b$estimate,udat1c,B1,logdcopnames2,
  pcondnames1[1:2],np)
# should be 3|12 4|12 5|12 in $condforw
zdat1b=nscore(pseud1b$condforw)
# pairs(zdat1b)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
for(i2 in 2:3)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat1b[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.3513768 0.1072861 0.3486435 0.1655699 
#1 3 0.1971193 0.03413146 0.3866341 0.08224073 
#2 3 0.2801743 -0.09804698 0.3096864 0.1245980
cat("\n------------------------------------------------------------\n")

cat("Gumbel/BVN B1-vine\n")
logdcopnames2=c("logdgum","logdbvncop")
mle2b=nlm(rvinenllk1.trunc,p=parvec2[1:7],
    udat=udat2b,A=B1,logdcopnames=logdcopnames2,pcondnames=pcondnames2[1:2],
    hessian=T,iterlim=30,print.level=1,LB=c(rep(1,4),rep(-1,3)),
    UB=c(rep(20,4),rep(1,3)))
pseud2b=rvinenllkpseud(mle2b$estimate,udat2b,B1,logdcopnames2,
  pcondnames2[1:2],np)
# should be 3|12 4|12 5|12 in $condforw
zdat2b=nscore(pseud2b$condforw)
# pairs(zdat2b)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
for(i2 in 2:3)
{ for(i1 in 1:(i2-1))
  { semic=semicor(zdat2b[,c(i1,i2)])
    ze=bvnsemic(semic[1])
    cat(i1,i2,semic,ze,"\n")
  }
}
#1 2 0.2843866 0.07154063 0.09199461 0.1268943 
#1 3 0.1606774 0.1294251 0.08974877 0.0652587 
#2 3 0.2153023 -0.05918468 0.02401360 0.09106255

cat("\n============================================================\n")

#============================================================

# 3-truncation

cat("\n3-truncation\n")

cat("Gumbel C-vine\n")
logdcopnames3=c("logdgum","logdgum","logdgum")
mle1c=nlm(rvinenllk1.trunc,p=parvec1[1:9],
    udat=udat1c,A=C5,logdcopnames=logdcopnames3,pcondnames=pcondnames1[1:3],
    hessian=T,iterlim=30,print.level=1,LB=1,UB=20)
pseud1c=rvinenllkpseud(mle1c$estimate,udat1c,C5,logdcopnames3,
  pcondnames1[1:3],np)
# should be 4|123 5|123 in $condforw
zdat1c=nscore(pseud1c$condforw)
# plot(zdat1c)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
semic=semicor(zdat1c)
ze=bvnsemic(semic[1])
cat(i1,i2,semic,ze,"\n")
# 2 3 0.2300463 0.04999838 0.2922819 0.09839489
cat("\n------------------------------------------------------------\n")

cat("Gumbel/BVN C-vine\n")
logdcopnames3=c("logdgum","logdbvncop","logdbvncop")
mle2c=nlm(rvinenllk1.trunc,p=parvec2[1:9],
    udat=udat2c,A=C5,logdcopnames=logdcopnames3,pcondnames=pcondnames2[1:3],
    hessian=T,iterlim=30,print.level=1,LB=c(rep(1,4),rep(-1,5)),
    UB=c(rep(20,4),rep(1,5)))
pseud2c=rvinenllkpseud(mle2c$estimate,udat2c,C5,logdcopnames3,
  pcondnames2[1:3],np)
# should be 4|123 5|123 in $condforw
zdat2c=nscore(pseud2c$condforw)
# plot(zdat2c)
cat("i1 i2 ncor lcor ucor bvnsemic\n")
semic=semicor(zdat2c)
ze=bvnsemic(semic[1])
cat(i1,i2,semic,ze,"\n")
#2 3 0.1770966 0.01828018 0.02235871 0.07279782

cat("\n============================================================\n")
