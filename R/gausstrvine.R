# R interface to gausstrvine for 
# best Gaussian truncated vines with 1,2,3,...d-2 trees

# rmat = positive definite dxd correlation matrix (4<=d<=8)
# iprint = print flag for intermediate results in the f90 code
# Output: 
# best ell-truncated vine for ell=1,...,d-2 based on min logdet
#    bnum[ell] : number of the vine array for best 
#    permmat[,ell] : permutation for best
#    pcarr[,.ell]   : partial correlations for best stored as d*d matrix
#    logdet[ell] : log determinant of best
gausstrvine=function(rmat,iprint=F)
{ d=nrow(rmat)
  if(d>8 | d!=ncol(rmat)) { cat("dimension<=8\n"); return(NA) }
  # no checking if rmat is symmetric
  if(!isposdef(rmat)) { cat("not positive definite\n"); return(NA) }
  d2=d-2; d1=d-1; dd=(d*d1)/2; 
  
  out= .Fortran("gausstrvine",
      as.integer(d), as.integer(dd), as.double(rmat), as.integer(iprint),
      bnumv=as.integer(rep(0,d)),
      logdetmin=as.double(rep(0,d)),
      permmat=as.integer(rep(0,d*d)),
      pcmat=as.double(rep(0,dd*d))  )
  pcmat=out$pcmat[1:(dd*d2)]; pcmat=matrix(pcmat,dd,d2)
  pcarr=array(0,c(d,d,d2))
  for(ell in 1:d2)
  { i1=1; i2=d1
    for(i in 1:d1)
    { pcarr[i,(i+1):d, ell]=pcmat[i1:i2,ell]
      i1=i2+1; i2=i2+(d-i-1)
    }
  }
  if(iprint)
  { cat("\nbest truncated 1,...d-2: bnumv, logdetmin, permmat, pcarr[,,1:(d-2)]\n")
    print(out$bnumv[1:d2])
    print(out$logdetmin[1:d1])
    print(matrix(out$permmat,d,d))
    print(pcarr)
  }
  list(bnum=out$bnumv[1:d2], logdetmin=out$logdetmin[1:d1], 
       permmat=matrix(out$permmat[1:(d2*d)],d,d2),
       pcarr=pcarr) 
}

# rmat = positive definite dxd correlation matrix (4<=d<=8)
# jtrunc = truncation level to check on non-unique best truncated vines
# eps = tolerance to check on degree of non-uniqueness
# Output:
# The returned numbers are not the distinct numbers of non-equivalent
# truncated vines. The printed output can be used to retrieve some of the
# non-unique best truncated vines via (i) the single binary number for the
# vine array and (ii) the permutation of the variable indices.
gausstrvine.nonuniq=function(rmat, jtrunc=3,eps=1.e-7,iprint=F)
{ d=nrow(rmat)
  if(d>8 | d!=ncol(rmat)) return(NA)
  # no checking if rmat is symmetric
  if(!isposdef(rmat)) return(NA)
  d2=d-2; d1=d-1; dd=(d*d1)/2; 
  
  out= .Fortran("gausstrvinenonuniq",
      as.integer(d), as.integer(dd), as.double(rmat), 
      as.double(eps), as.integer(jtrunc), as.integer(iprint),
      bnumv=as.integer(rep(0,d)),
      logdetmin=as.double(rep(0,d)),
      permmat=as.integer(rep(0,d*d)),
      pcmat=as.double(rep(0,dd*d))  )
  pcmat=out$pcmat[1:(dd*d2)]; pcmat=matrix(pcmat,dd,d2)
  pcarr=array(0,c(d,d,d2))
  for(ell in 1:d2)
  { i1=1; i2=d1
    for(i in 1:d1)
    { pcarr[i,(i+1):d, ell]=pcmat[i1:i2,ell]
      i1=i2+1; i2=i2+(d-i-1)
    }
  }
  if(iprint)
  { cat("\nbest truncated 1,...d-2: bnumv, logdetmin, permmat, pcarr[,,1:(d-2)]\n")
    print(out$bnumv[1:d2])
    print(out$logdetmin[1:d1])
    print(matrix(out$permmat,d,d))
    print(pcarr)
  }
  list(bnum=out$bnumv[1:d2], logdetmx=out$logdetmx[1:d1], 
       permmat=matrix(out$permmat[1:(d2*d)],d,d2),
       pcarr=pcarr) 
}

