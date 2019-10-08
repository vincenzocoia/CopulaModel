# check of integrals for tau and rhoS for bivariate Frank copula

library(CopulaModel)
bvec=seq(0,.9,.02)  
bvec=c(bvec,.95)
bvec=c(bvec,1)
np=length(bvec)
bvec=bvec[-(np-1)]
cpar=frk.b2cpar(bvec)
cat("\npositive dependence\n")
for(i in 2:(np-2))
{ rho=frk.cpar2rhoS(cpar[i]); 
  tau=frk.cpar2tau(cpar[i])
  cat(i,bvec[i],cpar[i],tau,rho,"\n") 
}

cat("\nnegative dependence\n")
cpar=frk.b2cpar(-bvec)
for(i in 2:(np-2))
{ rho=frk.cpar2rhoS(cpar[i]); 
  tau=frk.cpar2tau(cpar[i])
  cat(i,bvec[i],cpar[i],tau,rho,"\n") 
}
