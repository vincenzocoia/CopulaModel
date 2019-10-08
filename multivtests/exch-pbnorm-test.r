# compare exchmvn and pbnorm in some extreme cases
# code in pbnorm now used exchmvn when the original returns 2.e-9 or less

library(CopulaModel)

# wrapper function to compare exchmvn and pbnorm for bivariate normal cdf
# z1,z2 = two real values
wrapcomp=function(z1,z2)
{ rhvec=seq(0,.9,.1)
  rhvec=c(-.9,-.6,-.3,rhvec)
  cat("\nz1,z2 = ", z1,z2,"\n")
  out= pbnorm(z1,z2,rhvec)
  nn=length(rhvec)
  eout=rep(0,nn)
  for (i in 1:nn) 
  { rh=rhvec[i]
    if(rh>=0)  
    { lb=c(min(z1-5,-6),min(z2-5,-6))
       ub=c(z1,z2)
       eout[i]=exchmvn(lb,ub,rh) 
    }
    else 
    { lb=c(min(z1-5,-6),min(-z2-5,-6))
      ub=c(z1,-z2)
      eout[i]=pnorm(z1)-exchmvn(lb,ub,-rh)
    }
  }
  diff=out-eout
  print(cbind(rhvec,out,eout,diff))
  invisible(0)
}

wrapcomp(1,1)
wrapcomp(1,2)
wrapcomp(-1,2)
wrapcomp(2,3)
wrapcomp(2,5)
wrapcomp(2,6)
wrapcomp(-1,7)

wrapcomp(2,-4)  
wrapcomp(-1,-4)
wrapcomp(-2,-4)
wrapcomp(-3,-4) 
wrapcomp(-1,-5) 
wrapcomp(-2,-5)
wrapcomp(-3,-5)

wrapcomp(-2,-2)
wrapcomp(-2,-3)
wrapcomp(-3,-3)
wrapcomp(-3,-4)
wrapcomp(-3.5,-4)
wrapcomp(-3.9,-4)
wrapcomp(-4,-4)
wrapcomp(-4,-5) 
wrapcomp(-5,-5) 
wrapcomp(-6,-6)

