# check code for discrete R-vine and D-vine

# n = dimension 
# ii = integer in [0,2^(n-1))
# return binary vector jj of dimension n corresponding to integer ii
d2b=function(n,ii)
{ dd=ii
  jj=rep(0,n)
  for(i in n:1)
  { jj[i]=dd%%2; dd=floor(dd/2) }
  jj
}

library(CopulaModel)
#source("../R/rvinediscrete.R")

upar=c(.3,.4,.5,.6,.7,.8,.9)
parmat=matrix(c(0,1,1.5,1,1.5,1.5,1.5, 0,0,2,1.5,2,2,2, 0,0,0,2,1.5,1.5,1.5,
   0,0,0,0,1.6,1.6,1.6, 0,0,0,0,0,1.4,1.5, 0,0,0,0,0,0,1.4, 0,0,0,0,0,0,0), 
   7,7, byrow=T)
# select a value of d
d=3
#d=4
#d=5
#d=6
#d=7
if(d==3) iprint=T else iprint=F
# select one of the tow below
pcopnames=rep("pfrk",d)
#pcopnames=rep("pgum",d)

# wrapper function
# d = dimension
# strcop = d-vector of names of pair-copulas
wraprun=function(d,strcop)
{ cat("d=", d, "pcop=",strcop,"\n")
  if(d==3) iprint=T else iprint=F
  pcopnames=rep(strcop,d)
  cat("\nrunning dvinepmf.discrete\n")
  sm=0
  pmf=rep(0,2^d)
  for(ii in 0:(2^d-1))
  { yy=d2b(d,ii)
    u2vec=pbinom(yy,1,upar[1:d])
    u1vec=pbinom(yy-1,1,upar[1:d])
    out=dvinepmf.discrete(parmat[1:d,1:d],u1vec,u2vec,pcopnames,iprint)
    sm=sm+out
    pmf[ii+1]=out
  }
  print(sm)
  
  DD=Dvinearray(d)
  out=varray2M(DD); M=out$mxarray
  cat("\nrunning rvinepmf.discrete with A=D\n")
  smr=0
  pmfr=rep(0,2^d)
  for(ii in 0:(2^d-1))
  { yy=d2b(d,ii)
    u2vec=pbinom(yy,1,upar[1:d])
    u1vec=pbinom(yy-1,1,upar[1:d])
    out=rvinepmf.discrete(parmat[1:d,1:d],u1vec,u2vec,DD,M,pcopnames,iprint)
    smr=smr+out
    pmfr[ii+1]=out
  }
  print(smr)

  print(max(abs(pmf-pmfr)))
  cat("============================================================\n")
  invisible(0)
}

wraprun(3,"pfrk")
wraprun(4,"pfrk")
wraprun(5,"pfrk")
wraprun(6,"pfrk")
wraprun(7,"pfrk")
# sm and smr
# Frank
# d=3: 1
# d=4: 1
# d=5: 1
# d=6: 1
# d=7: 1

wraprun(3,"pgum")
wraprun(4,"pgum")
wraprun(5,"pgum")
wraprun(6,"pgum")
wraprun(7,"pgum")
# Gumbel
# d=3: 0.9999997
# d=4: 0.9999995
# d=5: 0.9999994
# d=6: 0.9999992
# d=7: 0.999999

