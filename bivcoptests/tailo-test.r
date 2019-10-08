# validate code via tail orders and
# check if computational problems for near (0,0) and (1,1)

library(CopulaModel)

cpar=2
cpar=3
uvec=c(.1,.05,.01,.005,.001,.0001)
uvec1=1-uvec

# limit should be de
cat("Plackett, cpar=",cpar, "\n")
print(ppla(uvec,uvec,cpar)/uvec^2)
print((2*uvec-1+ppla(uvec1,uvec1,cpar))/uvec^2)
print(dpla(uvec,uvec,cpar))
print(dpla(uvec1,uvec1,cpar))
cat("limits should be ", cpar,"\n")

# limit should be cpar/(1-exp(-cpar))
cat("\nFrank, cpar=",cpar, "\n")
print(pfrk(uvec,uvec,cpar)/uvec^2)
print((2*uvec-1+pfrk(uvec1,uvec1,cpar))/uvec^2)
print(dfrk(uvec,uvec,cpar))
print(dfrk(uvec1,uvec1,cpar))
cat("limits should be ", cpar/(1-exp(-cpar)),"\n")

kap=2^(1/cpar)
cat("\nGumbel, cpar=",cpar,"\n")
#print(pgum(uvec,uvec,cpar)/uvec^kap)
print(log(pgum(uvec,uvec,cpar))/log(uvec))
cat("kappaL= ", kap, "\n")
print((1-2*uvec1+pgum(uvec1,uvec1,cpar))/uvec)
cat("lambdaU=", gum.cpar2lm(cpar),"\n")
print(dgum(uvec,uvec,cpar)*uvec^(2-kap))
cat("density at (0,0) should go to ", (kap/2)^2,"\n")
print(dgum(uvec1,uvec1,cpar)*uvec)
cat("(1-u)*density at (1,1) should go to ", (cpar-1)*2^(1/cpar-2),"\n")

cat("\nGalambos, cpar=",cpar, "\n")
kap=2-2^(-1/cpar)
#print(pgal(uvec,uvec,cpar)/uvec^kap)
print(log(pgal(uvec,uvec,cpar))/log(uvec))
cat("kappaL= ", kap, "\n")
print((1-2*uvec1+pgal(uvec1,uvec1,cpar))/uvec)
cat("lambdaU=", gal.cpar2lm(cpar),"\n")
print(dgal(uvec,uvec,cpar)*uvec^(2-kap))
cat("density at (0,0) should go to ", (kap/2)^2,"\n")
print(dgal(uvec1,uvec1,cpar)*uvec)
cat("(1-u)*density at (1,1) should go to ", (cpar+1)*2^(-1/cpar-2),"\n")

cat("\nMTCJ, cpar=",cpar, "\n")
print(pmtcj(uvec,uvec,cpar)/uvec)
cat("lambdaL=", mtcj.cpar2lm(cpar),"\n")
print((1-2*uvec1+pmtcj(uvec1,uvec1,cpar))/uvec^2) # 1+cpar
cat("limit should be ", 1+cpar,"\n")
print(dmtcj(uvec,uvec,cpar)*uvec) # lower tail
cat("limit should be ", (cpar+1)*2^(-1/cpar-2),"\n")
print(dmtcj(1-uvec,1-uvec,cpar))  # 1+cpar
cat("limit should be ", 1+cpar,"\n")

cat("\nJoe/B5, cpar=",cpar, "\n")
print(pjoe(uvec,uvec,cpar)/uvec^2)   # cpar
cat("limit should be ", cpar,"\n")
print((1-2*uvec1+pjoe(uvec1,uvec1,cpar))/uvec)
cat("lambdaU=",joe.cpar2lm(cpar),"\n")
print(djoe(uvec,uvec,cpar)) # cpar
cat("limit should be ", cpar,"\n")
print(djoe(uvec1,uvec1,cpar)*uvec) # upper tail
cat("(1-u)*density at (1,1) should go to ", (cpar-1)*2^(1/cpar-2),"\n")

