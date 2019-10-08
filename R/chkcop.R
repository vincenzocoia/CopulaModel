# Functions used to check copula derivatives

# check copcdf, coppdf, copcond for correctness via numerical derivatives
# u = scalar in interval (0,1), 
# vvec = vector of numbers between 0 and 1
# cpar = parameter for bivariate copula
# bcdf = copula cdf, 
# pcond = C_{2|1} conditional cdf,
# bpdf = copula density, 
# str = string for copula name
# eps = increment for numerical derivatives
# Output is null, but diagnostic output is printed
chkcopderiv=function(u,vvec,cpar,bcdf,pcond,bpdf,str=" ",eps=1.e-4)
{ nn=length(vvec)
  cat("\n",str, " with parameter ", cpar,"\n")
  cat("u   v    cdf      numccdf     ccdf\n")
  cat("u   v    ccdf     numpdf      pdf\n")
  for(i in 1:nn)
  { v=vvec[i]
    cdf=bcdf(u,v,cpar)
    cdf2=bcdf(u+eps,v,cpar)
    ccdf=pcond(v,u,cpar)
    #cat(u,v,cdf,(cdf2-cdf)/eps,ccdf,"\n")
    tem=(cdf2-cdf)/eps
    cat(u,v,cdf,tem,ccdf,"\n")
    if(abs(tem-ccdf)>eps*30) cat("*** ccdf roundoff???\n")
    ccdf2=pcond(v+eps,u,cpar)
    pdf=bpdf(u,v,cpar)
    #cat(u,v,ccdf,(ccdf2-ccdf)/eps,pdf,"\n")
    tem=(ccdf2-ccdf)/eps
    cat(u,v,ccdf,tem,pdf,"\n")
    if(abs(tem-pdf)>eps*30) cat("*** dcop roundoff???\n")
  }
  invisible(0)
}

# uvec,vvec = vectors of length nn with values in (0,1) 
# cpar = copula parameter, 
# pcond = function for conditional cdf 
# qcond = function for inverse conditional cdf
# tol = tolerance pcond/qcond composition versus identity
# Output is null, but diagnostic output is printed
chkcopcond=function(uvec,vvec,cpar,pcond,qcond,str=" ",tol=1.e-5)
{ nn=length(uvec)
  cat("\n", str, " with parameter ", cpar, "\n")
  cat("u   v    pcond  vv\n")
  cat("u   p    qcond  pp\n")
  for(i in 1:nn)
  { pp=pcond(vvec[i],uvec[i],cpar)
    vv=qcond(pp,uvec[i],cpar)
    cat(uvec[i],vvec[i],pp,vv,"\n")
    if(abs(vvec[i]-vv)>tol) cat("*** v roundoff???\n")
    p=vvec[i]
    vv=qcond(p,uvec[i],cpar)
    pp=pcond(vv,uvec[i],cpar)
    cat(uvec[i],p,vv,pp,"\n")
    if(abs(p-pp)>tol) cat("*** p roundoff???\n")
  }
  invisible(0)
}


# check if a vector of numbers is in increasing order 
# for example, dependence measures computed via numerical integration
# as a function of increasing copula parameter 
# mat = matrix or dataframe
# Output: 1 if OK, 0 if not
chkincrease=function(mat)
{ nn=ncol(mat)
  iok=1
  for(j in 1:nn)
  { tem=diff(mat[,j])
    nneg=sum(tem<0)
    if(nneg>0) { cat("out of order in column ", j,"\n"); iok=0; }
  }
  iok
}
 
