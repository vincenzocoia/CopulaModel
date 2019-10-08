# functions for partial correlations to correlations and 
# regression coefficients for vine array A;
# check inverse operations

library(CopulaModel)

d=6
bnum=2^((d-2)*(d-3)/2)-1  # largest bnum is for C-vine
set.seed(123)
b=sample(0:bnum,1)
A=vnum2array(d,b)
print(A)

pp=matrix(0,d,d)
for(ell in 1:(d-1))
{ pp[ell,(ell+1):d]=runif(d-ell) }
print(pp)

rmat=pcor2cor.rvine(pp,A)
print(rmat)
regmat=cor2reg(rmat,A,iprint=F) 
print(regmat)
rmat2=reg2cor(regmat,A)
print(max(abs(rmat2-rmat)))
pc=cor2pcor.rvine(rmat,A)
print(max(abs(pc$pctree-pp)))

# now choose more random inputs
set.seed(1234)
for(iarray in 1:10)
{ b=sample(0:bnum,1)
  A=vnum2array(d,b)
  cat("checking array ", iarray,"\n")
  for(ii in 1:10)
  { pp=matrix(0,d,d)
    for(ell in 1:(d-1))
    { pp[ell,(ell+1):d]=runif(d-ell) }
    rmat=pcor2cor.rvine(pp,A)
    regmat=cor2reg(rmat,A,iprint=F) 
    rmat2=reg2cor(regmat,A)
    tem1=max(abs(rmat2-rmat))
    pc=cor2pcor.rvine(rmat,A)
    tem2=max(abs(pc$pctree-pp))
    cat("  checking partial correlation array ", ii, " : ", tem1,tem2,"\n")
  }
}

