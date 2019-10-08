# checks for 2-parameter copula families
library(CopulaModel)

# check monotonicity of bivariate cdf
# u = value in (0,1) and v = increasing vector of values in (0,1), or
# v = value in (0,1) and u = increasing vector of values in (0,1)
# cpar = copula parameter 
# bcdf = function for bivariate copula cdf
# str = string for name of copula
# Output: cdf values at (u,v)
chkmono=function(u,v,cpar,bcdf,str=" ")
{ tem=bcdf(u,v,cpar)
  cat("\n",str, cpar,"\n")
  print(tem)
  invisible(0)
}

cpar=c(1.4,1.2)
cpar8=c(1.4,.2)
cpar10=c(1.4,.4)
vvec=seq(.1,.9,.1)
u=.8
u=.4

cat("check monotonicity of cdf\n")
chkmono(u,vvec,cpar,pbb1,"BB1")
chkmono(u,vvec,cpar,pbb2,"BB2")
chkmono(u,vvec,cpar,pbb3,"BB3")
chkmono(u,vvec,cpar,pbb4,"BB4")
chkmono(u,vvec,cpar,pbb5,"BB5")
chkmono(u,vvec,cpar,pbb6,"BB6")
chkmono(u,vvec,cpar,pbb7,"BB7")
chkmono(u,vvec,cpar8,pbb8,"BB8")
chkmono(u,vvec,cpar,pbb9,"BB9")
chkmono(u,vvec,cpar10,pbb10,"BB10")

cat("\ncheck pcop, pcond, dcop\n")
chkcopderiv(u,vvec,cpar,pbb1,pcondbb1,dbb1,"BB1")
chkcopderiv(u,vvec,cpar,pbb2,pcondbb2,dbb2,"BB2")
chkcopderiv(u,vvec,cpar,pbb3,pcondbb3,dbb3,"BB3")
chkcopderiv(u,vvec,cpar,pbb4,pcondbb4,dbb4,"BB4")
chkcopderiv(u,vvec,cpar,pbb5,pcondbb5,dbb5,"BB5")
chkcopderiv(u,vvec,cpar,pbb6,pcondbb6,dbb6,"BB6")
chkcopderiv(u,vvec,cpar,pbb7,pcondbb7,dbb7,"BB7")
chkcopderiv(u,vvec,cpar8,pbb8,pcondbb8,dbb8,"BB8")
chkcopderiv(u,vvec,cpar,pbb9,pcondbb9,dbb9,"BB9")
chkcopderiv(u,vvec,cpar10,pbb10,pcondbb10,dbb10,"BB10")

cat("\ncheck pcond, qcond\n")
u=seq(.1,.6,.1)
v=seq(.4,.9,.1)
chkcopcond(u,v,cpar,pcondbb1,qcondbb1,"BB1")
chkcopcond(u,v,cpar,pcondbb2,qcondbb2,"BB2")
chkcopcond(u,v,cpar,pcondbb3,qcondbb3,"BB3")
chkcopcond(u,v,cpar,pcondbb4,qcondbb4,"BB4")
chkcopcond(u,v,cpar,pcondbb5,qcondbb5,"BB5")
chkcopcond(u,v,cpar,pcondbb6,qcondbb6,"BB6")
chkcopcond(u,v,cpar,pcondbb7,qcondbb7,"BB7")
chkcopcond(u,v,cpar8,pcondbb8,qcondbb8,"BB8")
chkcopcond(u,v,cpar,pcondbb9,qcondbb9,"BB9")
chkcopcond(u,v,cpar10,pcondbb10,qcondbb10,"BB10")


