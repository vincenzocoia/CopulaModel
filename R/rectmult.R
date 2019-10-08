# Functions for rectangle probabilities from a multivariate copula cdf  

# u1vec = lower vector of hyperrectangle in [0,1]^d
# u2vec = upper vector of hyperrectangle in [0,1]^d 
# cpar = copula parameter 
# pmcop = function for multivariate copula
# Output: multivariate rectangle probability
rectmult=function(u1vec,u2vec,cpar,pmcop)
{ d=length(u1vec)
  rect=0
  for(ii in 0:(2^d-1))
  { jj=d2b(d,ii)
    tem=(sum(jj==0))
    sg=ifelse(tem%%2==0,1,-1) 
    uu=u2vec
    kk=(jj==0)
    uu[kk]=u1vec[kk]
    cdf=pmcop(uu,cpar)
    rect=rect+sg*cdf
  }
  rect
}

# u1vec = lower vector of hyperrectangle in [0,1]^d
# u2vec = upper vector of hyperrectangle in [0,1]^d 
# cpar = copula parameter >0 or <0 
#   d=length(u1vec)=length(u2vec)
# Output: multivariate Frank copula rectangle probability 
rectmfrk=function(u1vec,u2vec,cpar)
{ rectmult(u1vec,u2vec,cpar,pmfrk) }

# u1vec = lower vector of hyperrectangle in [0,1]^d
# u2vec = upper vector of hyperrectangle in [0,1]^d 
# cpar = copula parameter>1 
#   d=length(u1vec)=length(u2vec)
# Output: multivariate Gumbel copula rectangle probability 
rectmgum=function(u1vec,u2vec,cpar)
{ rectmult(u1vec,u2vec,cpar,pmgum) }

# u1vec = lower vector of hyperrectangle in [0,1]^d
# u2vec = upper vector of hyperrectangle in [0,1]^d 
# cpar = copula parameter >0
#   d=length(u1vec)=length(u2vec)
# Output: multivariate Galambos copula rectangle probability 
rectmgal=function(u1vec,u2vec,cpar)
{ rectmult(u1vec,u2vec,cpar,pmgal) }

# u1vec = lower vector of hyperrectangle in [0,1]^d
# u2vec = upper vector of hyperrectangle in [0,1]^d 
# cpar = copula parameter in (-1,1)
#   d=length(u1vec)=length(u2vec)
# Output: exchangeable multivariate normal rectangle probability
rectemvn=function(u1vec,u2vec,cpar)
{ z1=qnorm(u1vec); z1[z1< -6]= -6
  z2=qnorm(u2vec); z2[z2>6]=6
  pr=exchmvn(z1,z2,cpar) # extracted from library(mprobit)
  pr
}


#d=4; cpar=3
#out=rectmgum(rep(.1,4),rep(.9,4),cpar)
#print(out)
#rectmgum(rep(.1,8),rep(.9,8),cpar) 
#rectmgum(rep(.1,10),rep(.9,10),cpar) 
#rectmgum(rep(.1,14),rep(.9,14),cpar) 

#d=4; cpar=3
#out=rectmfrk(rep(.1,4),rep(.9,4),cpar)
#print(out)
#rectmfrk(rep(.1,8),rep(.9,8),cpar)
#rectmfrk(rep(.1,10),rep(.9,10),cpar) 
#rectmfrk(rep(.1,14),rep(.9,14),cpar) 

###############################################################################
