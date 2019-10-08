# Functions for simulation from factor copula models or multivariate t
# with factor structure. Included are bi-factor and nested-factor.
# Link to C code for speed.

# 1-factor copula: random generation
# n = sample size 
# d = #variables
# cpar = copula parameter
# copcode = copula code, currently a few options but more can be added later
#   1 for Gaussian, 2 for t, 3 for Gumbel,
#   -3 for reflected Gumbel, 5 for Frank, 9 for BB1
# Output: nxd matrix of d-variate random vectors in (0,1) 
r1fact=function(n,d,cpar,copcode)
{ dat=rep(0,n*d)
  # cpar should have length 2*d for BB1
  # cpar should have length d+1 for BVt
  # cpar should have length d for Gumbel, Frank and other 1-parameter biv 
  if(copcode==1 | copcode==2)
  { out= .C("r1factmvt", as.integer(n), as.integer(d), as.double(cpar),
           as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)  # random d-variate normal or t vector
  }
  else
  { out = .C("r1fact", as.integer(n), as.integer(d), as.double(cpar),
           as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)  # random uniform
  }
  dat
}

# 2-factor copula: random generation
# n = sample size 
# d = #variables
# cpar = copula parameter
# copcode = copula code, currently a few options but more can be added later
#   1 for Gaussian, 2 for t, 3 for Gumbel,
#   -3 for reflected Gumbel, 5 for Frank, 9 for BB1/Frank
# Output: nxd matrix of d-variate random vectors in (0,1) 
r2fact=function(n,d,cpar,copcode)
{ dat=rep(0,n*d)
  # cpar should have length 3*d for BB1+Frank
  # cpar should have length 2*d+1 for BVt
  # cpar should have length 2*d for Gumbel, Frank and other 1-parameter biv 
  if(copcode==1 | copcode==2)
  { out= .C("r2factmvt", as.integer(n), as.integer(d), as.double(cpar),
           as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)  # random d-variate normal or t vector
  }
  else
  { out= .C("r2fact", as.integer(n), as.integer(d), as.double(cpar),
           as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)  # random uniform
  }
  dat
}

# bifactor copula: random generation
# n = sample size 
# grsize = vector of group sizes for m groups with sum(grsize)=d
#    d = #variables
# cpar = copula parameter
#  order is factor1 (var[1],...,var[d]), then factor2 (var[1],...,var[d])
#  the order of parameters in cpar is different from simbifact for BB1
# copcode = copula code, currently a few options but more can be added later
#   1 for Gaussian, 2 for t, 3 for Gumbel,
#   -3 for reflected Gumbel, 5 for Frank, 9 for BB1/Frank
# Output: nxd matrix of d-variate random vectors in (0,1) 
rbifact=function(n,grsize,cpar,copcode)
{ d=sum(grsize); m=length(grsize)
  dat=rep(0,n*d)
  # cpar should have length 3*d for BB1+Frank
  # cpar should have length 2*d+1 for BVt
  # cpar should have length 2*d for Gumbel, Frank and other 1-parameter biv 
  if(copcode==1 | copcode==2)
  { out= .C("rbifactmvt", as.integer(n), as.integer(m), as.integer(grsize),
           as.double(cpar), as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)  # random d-variate normal or t vector
  }
  else
  { out= .C("rbifact", as.integer(n), as.integer(m), as.integer(grsize),
           as.double(cpar), as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)
  }
  dat
}

# nested factor copula: random generation
# n = sample size 
# grsize = vector of group sizes for m groups with sum(grsize)=d
#    d = #variables
# cpar = copula parameter
#  order is (grp[1],...,grp[m]), then (var[1],...,var[d])
# copcode = copula code, currently a few options but more can be added later
#   1 for Gaussian, 2 for t, 3 for Gumbel,
#   -3 for reflected Gumbel, 5 for Frank, 11 for Gumbel/BB1
# Output: nxd matrix of d-variate random vectors in (0,1) 
rnestfact=function(n,grsize,cpar,copcode)
{ d=sum(grsize); m=length(grsize)
  dat=rep(0,n*d)
  # cpar should have length m+d for Gumbel, Frank and other 1-parameter biv 
  # cpar should have length m+2*d for BB1+Frank
  if(copcode==1 | copcode==2)
  { out= .C("rnestfactmvt", as.integer(n), as.integer(m), as.integer(grsize),
           as.double(cpar), as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)  # random d-variate normal or t vector
  }
  else
  { out= .C("rnestfact", as.integer(n), as.integer(m), as.integer(grsize),
           as.double(cpar), as.integer(copcode), dat=as.double(dat)) 
    dat=matrix(out$dat,n,d)
  }
  dat
}

#============================================================
