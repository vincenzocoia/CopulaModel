# check if a square symmetric matrix is positive definite

# amat = symmetrix matrix
# Output: T if amat is positive definite, F otherwise
isposdef=function(amat)
{ tt=try(chol(amat), silent=T)
  ifelse(class(tt)=="matrix",T,F)
}

#a1=matrix(c(1,.5,.5,1),2,2)
#a2=matrix(c(1,1.5,1.5,1),2,2)
#t1=try(chol(a1))
#t2=try(chol(a2))
#print(isposdef(a1))
#print(isposdef(a2))
